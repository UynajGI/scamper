use mcmc_rs::{
    EuclideanState, FnLogDensity, MemoryTrace, SamplingPhase, SliceSampler, TraceStore,
    TransitionKernel,
};
use rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

#[test]
fn slice_sampler_recovers_normal_mean() {
    let mut target = FnLogDensity::new(|position: &[f64]| -0.5 * position[0].powi(2));
    let mut state = EuclideanState::initialize(&mut target, vec![3.0]).unwrap();
    let mut kernel = SliceSampler::new(vec![1.0]).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(71);
    for _ in 0..500 {
        kernel
            .transition(&mut target, &mut state, &mut rng, SamplingPhase::Warmup)
            .unwrap();
    }
    let mut trace = MemoryTrace::new(1, 1).unwrap();
    for _ in 0..5_000 {
        let report = kernel
            .transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
            .unwrap();
        trace.record(0, &state, &report).unwrap();
    }
    let values = trace.parameter(0).unwrap().collect::<Vec<_>>();
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    assert!(mean.abs() < 0.08, "mean={mean}");
}
