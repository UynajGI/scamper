use mcmc_rs::{
    diagnose, energy_bfmi, EuclideanState, FnLogDensity, MemoryTrace, TraceStore, TransitionReport,
};

#[test]
fn energy_bfmi_matches_a_hand_series_and_rejects_invalid_input() {
    let energies = [Some(0.0), Some(1.0), Some(0.0), Some(1.0)];
    assert_eq!(energy_bfmi(&energies), Some(3.0));
    assert_eq!(energy_bfmi(&[Some(1.0), Some(1.0), Some(1.0)]), None);
    assert_eq!(energy_bfmi(&[Some(0.0), None, Some(1.0)]), None);
    assert_eq!(energy_bfmi(&[Some(0.0), Some(1.0)]), None);
}

fn diagnostic_trace(chain: usize) -> MemoryTrace {
    let mut target = FnLogDensity::new(|position: &[f64]| -0.5 * position[0] * position[0]);
    let mut state = EuclideanState::initialize(&mut target, vec![chain as f64 * 0.01]).unwrap();
    let mut trace = MemoryTrace::new(1, 1).unwrap();

    for draw in 0..20 {
        if draw > 0 {
            state.replace(vec![(draw as f64 * 0.37 + chain as f64).sin()], -0.5);
        }
        let report = TransitionReport {
            accepted: Some(true),
            acceptance_statistic: Some(0.8),
            proposals: 1,
            acceptances: 1,
            energy: Some((draw as f64 * 0.31 + chain as f64).sin() + 2.0),
            tree_depth: Some(3),
            max_tree_depth_reached: draw % 7 == 0,
            subtransitions: 1,
            ..TransitionReport::default()
        };
        trace.record(chain, &state, &report).unwrap();
    }
    trace
}

#[test]
fn multichain_diagnostics_include_ebfmi_and_depth_hits() {
    let traces = vec![diagnostic_trace(0), diagnostic_trace(1)];
    let diagnostics = diagnose(&traces, &["x".to_string()]).unwrap();

    assert_eq!(diagnostics.chain_ebfmi.len(), 2);
    assert!(diagnostics.chain_ebfmi.iter().all(Option::is_some));
    assert_eq!(diagnostics.max_tree_depth_hits, 6);
    assert!((diagnostics.mean_acceptance - 0.8).abs() < 1.0e-12);
}

#[test]
fn v03_json_trace_can_be_extended_with_v04_columns() {
    let trace = diagnostic_trace(0);
    let mut value = serde_json::to_value(&trace).unwrap();
    let object = value.as_object_mut().unwrap();
    object.remove("energy");
    object.remove("tree_depth");
    object.remove("max_tree_depth_reached");

    let mut restored: MemoryTrace = serde_json::from_value(value).unwrap();
    let mut target = FnLogDensity::new(|position: &[f64]| -0.5 * position[0] * position[0]);
    let state = EuclideanState::initialize(&mut target, vec![0.25]).unwrap();
    let report = TransitionReport {
        energy: Some(1.0),
        tree_depth: Some(2),
        max_tree_depth_reached: true,
        ..TransitionReport::default()
    };
    restored.record(0, &state, &report).unwrap();
    restored.validate().unwrap();

    assert_eq!(restored.energies().len(), restored.len());
    assert!(restored.energies()[..restored.len() - 1]
        .iter()
        .all(Option::is_none));
    assert_eq!(restored.energies().last(), Some(&Some(1.0)));
    assert_eq!(restored.tree_depths().last(), Some(&Some(2)));
    assert_eq!(restored.max_tree_depth_reached().last(), Some(&1));
}
