use mcmc_rs::{
    ChainCheckpoint, EuclideanState, FnLogDensity, MemoryTrace, RandomWalkMetropolis,
    SamplingPhase, TargetFingerprint, TraceStore, TransitionKernel, CHECKPOINT_FORMAT,
};
use rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

#[test]
fn json_checkpoint_preserves_exact_future_trajectory() {
    let mut target = FnLogDensity::new(|position: &[f64]| -0.5 * position[0].powi(2));
    let mut state = EuclideanState::initialize(&mut target, vec![0.4]).unwrap();
    let mut kernel = RandomWalkMetropolis::isotropic(1, 0.9).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(888);
    let mut trace = MemoryTrace::new(1, 1).unwrap();
    let fingerprint = TargetFingerprint {
        name: "standard-normal".to_string(),
        version: "1".to_string(),
        dimension: 1,
        parameter_names: vec!["x".to_string()],
    };

    for _ in 0..100 {
        let report = kernel
            .transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
            .unwrap();
        trace.record(0, &state, &report).unwrap();
    }
    let checkpoint = ChainCheckpoint {
        format: CHECKPOINT_FORMAT.to_string(),
        chain_id: 0,
        phase: SamplingPhase::Sampling,
        target: fingerprint.clone(),
        state: state.clone(),
        kernel: kernel.clone(),
        rng: rng.clone(),
        trace: trace.clone(),
        last_report: mcmc_rs::TransitionReport::default(),
    };
    let path = std::env::temp_dir().join(format!(
        "mcmc-rs-checkpoint-{}-{}.json",
        std::process::id(),
        state.iteration()
    ));
    checkpoint.save_json(&path).unwrap();
    let restored: ChainCheckpoint<RandomWalkMetropolis, Xoshiro256PlusPlus, MemoryTrace> =
        ChainCheckpoint::load_json(&path).unwrap();
    restored.validate_target(&fingerprint).unwrap();
    std::fs::remove_file(path).unwrap();

    let mut uninterrupted_state = state;
    let mut uninterrupted_kernel = kernel;
    let mut uninterrupted_rng = rng;
    let mut uninterrupted_trace = trace;
    for _ in 0..200 {
        let report = uninterrupted_kernel
            .transition(
                &mut target,
                &mut uninterrupted_state,
                &mut uninterrupted_rng,
                SamplingPhase::Sampling,
            )
            .unwrap();
        uninterrupted_trace
            .record(0, &uninterrupted_state, &report)
            .unwrap();
    }

    let mut resumed_target = FnLogDensity::new(|position: &[f64]| -0.5 * position[0].powi(2));
    let mut resumed_state = restored.state;
    let mut resumed_kernel = restored.kernel;
    let mut resumed_rng = restored.rng;
    let mut resumed_trace = restored.trace;
    for _ in 0..200 {
        let report = resumed_kernel
            .transition(
                &mut resumed_target,
                &mut resumed_state,
                &mut resumed_rng,
                SamplingPhase::Sampling,
            )
            .unwrap();
        resumed_trace.record(0, &resumed_state, &report).unwrap();
    }

    assert_eq!(uninterrupted_state.position(), resumed_state.position());
    assert_eq!(
        uninterrupted_state.log_density(),
        resumed_state.log_density()
    );
    assert_eq!(uninterrupted_trace, resumed_trace);
}
