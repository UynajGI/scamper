use mcmc_rs::proposal::standard_normal;
use mcmc_rs::{diagnose, EuclideanState, FnLogDensity, MemoryTrace, TraceStore, TransitionReport};
use rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

fn synthetic_trace(chain_id: usize, shift: f64, seed: u64) -> MemoryTrace {
    let mut target = FnLogDensity::new(|position: &[f64]| -0.5 * position[0].powi(2));
    let mut state = EuclideanState::initialize(&mut target, vec![0.0]).unwrap();
    let mut trace = MemoryTrace::new(1, 1).unwrap();
    let report = TransitionReport {
        accepted: Some(true),
        proposals: 1,
        acceptances: 1,
        ..TransitionReport::default()
    };
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    for _ in 0..2_000 {
        let value = standard_normal(&mut rng) + shift;
        state.replace(vec![value], -0.5 * value * value);
        trace.record(chain_id, &state, &report).unwrap();
    }
    trace
}

#[test]
fn iid_chains_have_good_rhat_and_large_ess() {
    let traces = (0..4)
        .map(|chain| synthetic_trace(chain, 0.0, chain as u64 + 5))
        .collect::<Vec<_>>();
    let result = diagnose(&traces, &["x".to_string()]).unwrap();
    let parameter = &result.parameters[0];
    assert!(parameter.rhat < 1.02, "rhat={}", parameter.rhat);
    assert!(parameter.ess_bulk > 3_000.0, "ess={}", parameter.ess_bulk);
}

#[test]
fn shifted_chain_is_detected() {
    let mut traces = (0..4)
        .map(|chain| synthetic_trace(chain, 0.0, chain as u64 + 19))
        .collect::<Vec<_>>();
    traces[3] = synthetic_trace(3, 1.5, 99);
    let result = diagnose(&traces, &["x".to_string()]).unwrap();
    assert!(result.parameters[0].rhat > 1.05);
}
