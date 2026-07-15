use mcmc_rs::{
    DiagonalMetric, DifferentiableLogDensity, DualAveraging, EuclideanState, LogDensity,
    SamplingPhase, StaticHmc, TransitionKernel, WarmupWindowConfig,
};
use rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

#[derive(Clone, Copy)]
struct AnisotropicGaussian;

impl LogDensity<[f64]> for AnisotropicGaussian {
    fn log_density(&mut self, position: &[f64]) -> f64 {
        -0.5 * (position[0] * position[0] / 0.04 + position[1] * position[1] / 9.0)
    }
}

impl DifferentiableLogDensity for AnisotropicGaussian {
    fn log_density_and_gradient(&mut self, position: &[f64], gradient: &mut [f64]) -> f64 {
        gradient[0] = -position[0] / 0.04;
        gradient[1] = -position[1] / 9.0;
        self.log_density(position)
    }
}

#[test]
fn dual_averaging_freezes_to_a_positive_finite_step_size() {
    let mut adaptation = DualAveraging::new(0.5, 0.8).unwrap();
    for index in 0..100 {
        let statistic = if index % 3 == 0 { 0.5 } else { 0.9 };
        adaptation.observe(statistic).unwrap();
    }
    let step_size = adaptation.freeze();
    assert!(step_size.is_finite() && step_size > 0.0);
    assert!(adaptation.is_frozen());
    assert!(adaptation.observe(0.8).is_err());
}

#[test]
fn windowed_diagonal_adaptation_updates_geometry_and_freezes() {
    let mut target = AnisotropicGaussian;
    let mut state = EuclideanState::initialize(&mut target, vec![0.1, 1.0]).unwrap();
    let metric = DiagonalMetric::unit(2).unwrap();
    let windows = WarmupWindowConfig {
        initial_buffer: 20,
        terminal_buffer: 20,
        initial_window: 20,
    };
    let mut kernel = StaticHmc::new(metric, 0.08, 6)
        .unwrap()
        .with_warmup_adaptation(
            120,
            0.8,
            mcmc_rs::MetricAdaptation::Diagonal {
                regularization: 1.0e-3,
            },
            windows,
        )
        .unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(123);

    for _ in 0..120 {
        kernel
            .transition(&mut target, &mut state, &mut rng, SamplingPhase::Warmup)
            .unwrap();
    }
    kernel
        .on_phase_end(&mut target, SamplingPhase::Warmup, &state)
        .unwrap();
    assert!(kernel.adaptation_is_frozen());
    assert!(kernel.step_size().is_finite() && kernel.step_size() > 0.0);
    let inverse_mass = kernel.metric().inverse_mass();
    assert!(inverse_mass
        .iter()
        .all(|value| value.is_finite() && *value > 0.0));
    assert_ne!(inverse_mass, &[1.0, 1.0]);
}

#[test]
fn entering_sampling_before_configured_warmup_is_rejected() {
    let mut target = AnisotropicGaussian;
    let state = EuclideanState::initialize(&mut target, vec![0.0, 0.0]).unwrap();
    let mut kernel = StaticHmc::new(DiagonalMetric::unit(2).unwrap(), 0.1, 4)
        .unwrap()
        .with_dual_averaging(50, 0.8)
        .unwrap();
    assert!(kernel
        .on_phase_start(&mut target, SamplingPhase::Sampling, &state)
        .is_err());
}
