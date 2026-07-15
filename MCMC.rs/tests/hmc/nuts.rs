use mcmc_rs::{
    DiagonalMetric, DifferentiableLogDensity, EuclideanState, LogDensity, Nuts, SamplingPhase,
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
fn nuts_recovers_standard_normal_moments() {
    let mut target = StandardNormal;
    let mut state = EuclideanState::initialize(&mut target, vec![2.0]).unwrap();
    let mut kernel = Nuts::new(DiagonalMetric::unit(1).unwrap(), 0.25, 8)
        .unwrap()
        .with_diagonal_adaptation(300, 0.8, 1.0e-3)
        .unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x4E55);

    for _ in 0..300 {
        kernel
            .transition(&mut target, &mut state, &mut rng, SamplingPhase::Warmup)
            .unwrap();
    }
    assert!(kernel.adaptation_is_frozen());

    let draws = 5_000_u32;
    let mut sum = 0.0;
    let mut square_sum = 0.0;
    let mut divergences = 0_u32;
    for _ in 0..draws {
        let report = kernel
            .transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
            .unwrap();
        assert!(report.tree_depth.is_some());
        assert!(report.tree_depth.unwrap() <= 8);
        assert!(report.leapfrog_steps > 0);
        assert!(report.energy.is_some());
        divergences += u32::from(report.divergent);
        let draw = state.position()[0];
        sum += draw;
        square_sum += draw * draw;
    }

    let mean = sum / f64::from(draws);
    let variance = square_sum / f64::from(draws) - mean * mean;
    assert!(mean.abs() < 0.08, "mean={mean}");
    assert!((variance - 1.0).abs() < 0.12, "variance={variance}");
    assert!(divergences < 10, "divergences={divergences}");
}

#[test]
fn nuts_reports_tree_depth_exhaustion() {
    let mut target = StandardNormal;
    let mut state = EuclideanState::initialize(&mut target, vec![0.0]).unwrap();
    let mut kernel = Nuts::unit(1, 0.01, 2).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(5);

    let report = kernel
        .transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
        .unwrap();

    assert_eq!(report.tree_depth, Some(2));
    assert_eq!(report.leapfrog_steps, 3);
    assert!(report.max_tree_depth_reached);
}

#[test]
fn divergent_nuts_trajectory_rejects_atomically() {
    let mut target = StandardNormal;
    let mut state = EuclideanState::initialize(&mut target, vec![0.5]).unwrap();
    let before_position = state.position().clone();
    let before_log_density = state.log_density();
    let before_iteration = state.iteration();
    let mut kernel = Nuts::unit(1, 1.0e200, 4)
        .unwrap()
        .with_max_energy_error(10.0)
        .unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(7);

    let report = kernel
        .transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
        .unwrap();

    assert!(report.divergent);
    assert_eq!(report.accepted, None);
    assert_eq!(state.position(), &before_position);
    assert_eq!(state.log_density(), before_log_density);
    assert_eq!(state.iteration(), before_iteration + 1);
    state.validate().unwrap();

    let encoded = serde_json::to_string(&kernel).unwrap();
    let _: Nuts<mcmc_rs::UnitMetric> = serde_json::from_str(&encoded).unwrap();
}

#[test]
fn nuts_warmup_checkpoint_roundtrip_preserves_future_trajectory() {
    let mut target = StandardNormal;
    let mut state = EuclideanState::initialize(&mut target, vec![1.25]).unwrap();
    let mut kernel = Nuts::new(DiagonalMetric::unit(1).unwrap(), 0.18, 7)
        .unwrap()
        .with_diagonal_adaptation(50, 0.8, 1.0e-3)
        .unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0xA11CE);

    for _ in 0..20 {
        kernel
            .transition(&mut target, &mut state, &mut rng, SamplingPhase::Warmup)
            .unwrap();
    }

    let state_json = serde_json::to_string(&state).unwrap();
    let kernel_json = serde_json::to_string(&kernel).unwrap();
    let rng_json = serde_json::to_string(&rng).unwrap();

    let mut restored_state: EuclideanState = serde_json::from_str(&state_json).unwrap();
    let mut restored_kernel: Nuts<DiagonalMetric> = serde_json::from_str(&kernel_json).unwrap();
    let mut restored_rng: Xoshiro256PlusPlus = serde_json::from_str(&rng_json).unwrap();
    let mut restored_target = StandardNormal;

    for index in 0..55 {
        let phase = if index < 30 {
            SamplingPhase::Warmup
        } else {
            SamplingPhase::Sampling
        };
        let report = kernel
            .transition(&mut target, &mut state, &mut rng, phase)
            .unwrap();
        let restored_report = restored_kernel
            .transition(
                &mut restored_target,
                &mut restored_state,
                &mut restored_rng,
                phase,
            )
            .unwrap();
        assert_eq!(report, restored_report);
        assert_eq!(state.position(), restored_state.position());
        assert_eq!(state.log_density(), restored_state.log_density());
        assert_eq!(state.iteration(), restored_state.iteration());
    }
}
