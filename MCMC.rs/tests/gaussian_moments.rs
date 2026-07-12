use mcmc_rs::{
    EuclideanState, FnLogDensity, MemoryTrace, RandomWalkMetropolis, SamplingPhase, TraceStore,
    TransitionKernel,
};
use rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

#[test]
fn adaptive_random_walk_recovers_standard_normal_moments() {
    let mut target = FnLogDensity::new(|position: &[f64]| -0.5 * position[0] * position[0]);
    let mut state = EuclideanState::initialize(&mut target, vec![4.0]).unwrap();
    let mut kernel = RandomWalkMetropolis::isotropic(1, 0.7)
        .unwrap()
        .with_scale_adaptation(0.44)
        .unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(11);
    kernel
        .on_phase_start(SamplingPhase::Warmup, &state)
        .unwrap();
    for _ in 0..4_000 {
        kernel
            .transition(&mut target, &mut state, &mut rng, SamplingPhase::Warmup)
            .unwrap();
    }
    kernel.on_phase_end(SamplingPhase::Warmup, &state).unwrap();
    kernel
        .on_phase_start(SamplingPhase::Sampling, &state)
        .unwrap();
    assert!(kernel.adaptation_is_frozen());

    let mut trace = MemoryTrace::new(1, 1).unwrap();
    for _ in 0..30_000 {
        let report = kernel
            .transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
            .unwrap();
        trace.record(0, &state, &report).unwrap();
    }
    let values = trace.parameter(0).unwrap().collect::<Vec<_>>();
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    let variance = values
        .iter()
        .map(|value| (value - mean).powi(2))
        .sum::<f64>()
        / (values.len() - 1) as f64;
    assert!(mean.abs() < 0.08, "mean={mean}");
    assert!((variance - 1.0).abs() < 0.12, "variance={variance}");
}
