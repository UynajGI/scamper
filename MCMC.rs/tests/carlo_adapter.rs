use carlo_rs::{Context, Run, RunConfig, RunId, TaskId};
use mcmc_rs::{FnLogDensity, McmcSampler, MemoryTrace, RandomWalkMetropolis, TraceStore};
use rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

#[test]
fn carlo_run_returns_sampler_trace() {
    let target = FnLogDensity::new(|position: &[f64]| -0.5 * position[0].powi(2));
    let kernel = RandomWalkMetropolis::isotropic(1, 1.0).unwrap();
    let trace = MemoryTrace::new(1, 1).unwrap();
    let sampler = McmcSampler::new(target, kernel, vec![0.0], trace, 0).unwrap();
    let config = RunConfig {
        thermalization_sweeps: 20,
        measurement_sweeps: 100,
        binsize: 10,
        base_seed: 44,
        progress_interval: 0,
        checkpoint_interval: 0,
    };
    let context = Context::new_with_binsize(
        Xoshiro256PlusPlus::seed_from_u64(config.base_seed),
        config.thermalization_sweeps,
        config.binsize,
    );
    let mut run = Run::from_parts(
        context,
        sampler,
        TaskId::new(0),
        RunId::new(0),
        config.clone(),
    );
    while !run.is_complete() {
        run.step();
    }
    let (results, sampler) = run.finalize_with_mc(config.base_seed);
    assert_eq!(sampler.trace().len(), 100);
    assert!(results.get("LogDensity").is_some());
}
