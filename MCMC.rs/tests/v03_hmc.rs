use mcmc_rs::{
    DifferentiableLogDensity, EuclideanState, LogDensity, SamplingPhase, StaticHmc,
    TransitionKernel,
};
use rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

#[derive(Clone, Copy)]
struct StandardNormal;

impl LogDensity<[f64]> for StandardNormal {
    fn log_density(&mut self, position: &[f64]) -> f64 {
        -0.5 * position.iter().map(|value| value * value).sum::<f64>()
    }
}

impl DifferentiableLogDensity for StandardNormal {
    fn log_density_and_gradient(&mut self, position: &[f64], gradient: &mut [f64]) -> f64 {
        for (gradient, position) in gradient.iter_mut().zip(position.iter().copied()) {
            *gradient = -position;
        }
        self.log_density(position)
    }
}

#[test]
fn static_hmc_recovers_standard_normal_moments() {
    let mut target = StandardNormal;
    let mut state = EuclideanState::initialize(&mut target, vec![2.0]).unwrap();
    let mut kernel = StaticHmc::unit(1, 0.23, 7).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(91);

    for _ in 0..200 {
        kernel
            .transition(&mut target, &mut state, &mut rng, SamplingPhase::Warmup)
            .unwrap();
    }
    assert!(kernel.adaptation_is_frozen());

    let mut sum = 0.0;
    let mut square_sum = 0.0;
    let draws = 4_000;
    for _ in 0..draws {
        let report = kernel
            .transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
            .unwrap();
        assert_eq!(report.leapfrog_steps, 7);
        sum += state.position()[0];
        square_sum += state.position()[0] * state.position()[0];
    }
    let mean = sum / f64::from(draws);
    let variance = square_sum / f64::from(draws) - mean * mean;
    assert!(mean.abs() < 0.08, "mean={mean}");
    assert!((variance - 1.0).abs() < 0.12, "variance={variance}");
}

#[test]
fn divergent_trajectory_rejects_without_corrupting_state() {
    let mut target = StandardNormal;
    let mut state = EuclideanState::initialize(&mut target, vec![0.5]).unwrap();
    let before_position = state.position().clone();
    let before_log_density = state.log_density();
    let before_iteration = state.iteration();
    let mut kernel = StaticHmc::unit(1, 1.0e200, 2)
        .unwrap()
        .with_max_energy_error(10.0)
        .unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(7);

    let report = kernel
        .transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
        .unwrap();
    assert!(report.divergent);
    assert_eq!(report.accepted, Some(false));
    assert_eq!(state.position(), &before_position);
    assert_eq!(state.log_density(), before_log_density);
    assert_eq!(state.iteration(), before_iteration + 1);
    state.validate().unwrap();

    let encoded = serde_json::to_string(&kernel).unwrap();
    let _: StaticHmc<mcmc_rs::UnitMetric> = serde_json::from_str(&encoded).unwrap();
}

#[test]
fn gradient_cache_avoids_recomputing_the_accepted_gradient_after_rejection() {
    #[derive(Default)]
    struct CountingTarget {
        calls: usize,
    }
    impl LogDensity<[f64]> for CountingTarget {
        fn log_density(&mut self, position: &[f64]) -> f64 {
            -0.5 * position[0] * position[0]
        }
    }
    impl DifferentiableLogDensity for CountingTarget {
        fn log_density_and_gradient(&mut self, position: &[f64], gradient: &mut [f64]) -> f64 {
            self.calls += 1;
            gradient[0] = -position[0];
            self.log_density(position)
        }
    }

    let mut target = CountingTarget::default();
    let mut state = EuclideanState::initialize(&mut target, vec![0.25]).unwrap();
    let mut kernel = StaticHmc::unit(1, 100.0, 4)
        .unwrap()
        .with_max_energy_error(1.0)
        .unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(11);
    let first = kernel
        .transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
        .unwrap();
    let second = kernel
        .transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
        .unwrap();
    assert!(first.divergent && second.divergent);
    assert_eq!(first.gradient_evaluations, first.target_evaluations);
    assert_eq!(second.gradient_evaluations, second.target_evaluations);
    assert_eq!(second.gradient_evaluations, second.leapfrog_steps);
    assert_eq!(first.gradient_evaluations, second.gradient_evaluations + 1);
}

#[test]
fn hmc_warmup_checkpoint_roundtrip_preserves_future_trajectory() {
    use mcmc_rs::DiagonalMetric;

    let mut target = StandardNormal;
    let mut state = EuclideanState::initialize(&mut target, vec![1.25]).unwrap();
    let mut kernel = StaticHmc::new(DiagonalMetric::unit(1).unwrap(), 0.18, 5)
        .unwrap()
        .with_diagonal_adaptation(50, 0.8, 1.0e-3)
        .unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0xA11CE);

    for _ in 0..20 {
        kernel
            .transition(&mut target, &mut state, &mut rng, SamplingPhase::Warmup)
            .unwrap();
    }

    let encoded = serde_json::to_string(&(kernel, state, rng)).unwrap();
    let (mut restored_kernel, mut restored_state, mut restored_rng): (
        StaticHmc<DiagonalMetric>,
        EuclideanState,
        Xoshiro256PlusPlus,
    ) = serde_json::from_str(&encoded).unwrap();
    let (mut reference_kernel, mut reference_state, mut reference_rng): (
        StaticHmc<DiagonalMetric>,
        EuclideanState,
        Xoshiro256PlusPlus,
    ) = serde_json::from_str(&encoded).unwrap();

    let mut restored_target = StandardNormal;
    let mut reference_target = StandardNormal;
    for index in 0..55 {
        let phase = if index < 30 {
            SamplingPhase::Warmup
        } else {
            SamplingPhase::Sampling
        };
        let restored_report = restored_kernel
            .transition(
                &mut restored_target,
                &mut restored_state,
                &mut restored_rng,
                phase,
            )
            .unwrap();
        let reference_report = reference_kernel
            .transition(
                &mut reference_target,
                &mut reference_state,
                &mut reference_rng,
                phase,
            )
            .unwrap();
        assert_eq!(restored_report, reference_report);
        assert_eq!(restored_state.position(), reference_state.position());
        assert_eq!(restored_state.log_density(), reference_state.log_density());
        assert_eq!(restored_state.iteration(), reference_state.iteration());
    }
}
