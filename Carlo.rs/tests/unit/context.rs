use carlo_rs::{Context, ContextCheckpoint, RunPhase};
use rand::SeedableRng;
use rand_core::Rng;
use rand_xoshiro::Xoshiro256PlusPlus;

#[test]
fn test_context_thermalization() {
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut ctx = Context::new(rng, 10);

    // Before thermalization
    assert!(!ctx.is_thermalized());
    assert_eq!(ctx.sweep_count(), 0);

    // Advance to thermalization
    for _ in 0..10 {
        ctx.advance_sweep();
    }

    // The warmup boundary is reached exactly after the configured count.
    assert!(ctx.is_thermalized());
}

#[test]
fn test_context_measurements() {
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut ctx = Context::new_with_binsize(rng, 0, 10);

    // Add measurements
    for i in 1..=10 {
        ctx.measure("Energy", i as f64);
    }

    // Finalize
    let estimates = ctx.finalize_measurements();
    assert!(estimates.contains_key("Energy"));

    let est = estimates.get("Energy").unwrap();
    assert!((est.mean - 5.5).abs() < 1e-10);
}

#[test]
fn test_context_rng_reproducibility() {
    let rng1 = Xoshiro256PlusPlus::seed_from_u64(12345);
    let rng2 = Xoshiro256PlusPlus::seed_from_u64(12345);

    let mut ctx1 = Context::new(rng1, 0);
    let mut ctx2 = Context::new(rng2, 0);

    // Same seed should give same random numbers
    for _ in 0..100 {
        let v1 = ctx1.rng.next_u64();
        let v2 = ctx2.rng.next_u64();
        assert_eq!(v1, v2);
    }
}

#[test]
fn test_context_checkpoint_state() {
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let ctx = Context::new(rng, 100);
    let state = ctx.checkpoint_state();
    assert_eq!(state.sweep_count, 0);
    assert_eq!(state.thermalization_sweeps, 100);
}

#[test]
fn test_context_multiple_observables() {
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut ctx = Context::new_with_binsize(rng, 0, 5);

    for i in 1..=10 {
        ctx.measure("Energy", i as f64);
        ctx.measure("Magnetization", i as f64 * 0.1);
    }

    let estimates = ctx.finalize_measurements();
    assert!(estimates.contains_key("Energy"));
    assert!(estimates.contains_key("Magnetization"));
}

#[test]
fn test_context_register_observable() {
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut ctx = Context::new_with_binsize(rng, 0, 10);

    ctx.register_observable("CustomObs", 20);

    ctx.measure("CustomObs", 1.0);
    ctx.measure("CustomObs", 2.0);

    let estimates = ctx.finalize_measurements();
    assert!(estimates.contains_key("CustomObs"));
}

#[test]
fn test_context_sweep_count_increment() {
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut ctx = Context::new(rng, 0);

    for i in 1..=100 {
        ctx.advance_sweep();
        assert_eq!(ctx.sweep_count(), i);
    }
}

#[test]
fn test_context_thermalization_sweeps_accessor() {
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let ctx = Context::new(rng, 50);

    assert_eq!(ctx.thermalization_sweeps(), 50);
}

#[test]
fn test_context_binsize() {
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let ctx = Context::new_with_binsize(rng, 0, 100);

    // Just verify creation with custom binsize works
    assert_eq!(ctx.sweep_count(), 0);
}

#[test]
fn test_context_measure_auto_register() {
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut ctx = Context::new_with_binsize(rng, 0, 10);

    // measure() should auto-register observable
    ctx.measure("AutoObs", 1.0);

    let estimates = ctx.finalize_measurements();
    assert!(estimates.contains_key("AutoObs"));
}

#[test]
fn test_context_register_with_shape() {
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut ctx = Context::new_with_binsize(rng, 0, 10);

    // Shape hint is currently ignored but should not error
    ctx.register_observable_with_shape("ArrayObs", 10, &[3, 3]);

    ctx.measure("ArrayObs", 1.0);
}

#[test]
fn test_context_explicit_phase() {
    let rng = Xoshiro256PlusPlus::seed_from_u64(7);
    let mut ctx = Context::new(rng, 100);
    assert_eq!(ctx.phase(), RunPhase::Initialization);
    ctx.enter_phase(RunPhase::Thermalization);
    assert!(ctx.phase().allows_adaptation());
    assert!(!ctx.is_thermalized());
    ctx.enter_phase(RunPhase::Measurement);
    assert!(ctx.phase().collects_measurements());
    assert!(ctx.is_thermalized());
}

#[test]
fn legacy_checkpoint_phase_is_inferred_from_counters() {
    let checkpoint = ContextCheckpoint {
        sweep_count: 7,
        thermalization_sweeps: 10,
        thermalized: false,
        phase: RunPhase::Initialization,
        attempted_updates: 0,
        accepted_moves: 0,
        event_time: 0.0,
    };
    let rng = Xoshiro256PlusPlus::seed_from_u64(8);
    let context = Context::restore_from_checkpoint(checkpoint, rng, 10);
    assert_eq!(context.phase(), RunPhase::Thermalization);

    let checkpoint = ContextCheckpoint {
        sweep_count: 10,
        thermalization_sweeps: 10,
        thermalized: true,
        phase: RunPhase::Initialization,
        attempted_updates: 0,
        accepted_moves: 0,
        event_time: 0.0,
    };
    let rng = Xoshiro256PlusPlus::seed_from_u64(9);
    let context = Context::restore_from_checkpoint(checkpoint, rng, 10);
    assert_eq!(context.phase(), RunPhase::Measurement);
}

#[test]
fn explicit_thermalization_phase_overrides_fixed_counter_boundary() {
    let rng = Xoshiro256PlusPlus::seed_from_u64(10);
    let mut context = Context::new(rng, 1);
    context.enter_phase(RunPhase::Thermalization);
    for _ in 0..3 {
        context.advance_sweep();
    }
    assert!(!context.is_thermalized());
    context.enter_phase(RunPhase::Measurement);
    assert!(context.is_thermalized());
}

#[test]
fn checkpoint_preserves_explicit_adaptation_phase() {
    let rng = Xoshiro256PlusPlus::seed_from_u64(11);
    let mut context = Context::new(rng, 1);
    context.enter_phase(RunPhase::Thermalization);
    context.advance_sweep();
    context.advance_sweep();
    let checkpoint = context.checkpoint_state();

    let rng = Xoshiro256PlusPlus::seed_from_u64(12);
    let restored = Context::restore_from_checkpoint(checkpoint, rng, 10);
    assert_eq!(restored.phase(), RunPhase::Thermalization);
    assert!(!restored.is_thermalized());
}

#[test]
fn explicit_simulation_clocks_round_trip_checkpoint() {
    let rng = Xoshiro256PlusPlus::seed_from_u64(13);
    let mut context = Context::new(rng, 0);
    context.record_attempts(17);
    context.record_accepted_moves(9);
    context.advance_event_time(2.5);
    context.advance_sweep();
    let checkpoint = context.checkpoint_state();

    let rng = Xoshiro256PlusPlus::seed_from_u64(14);
    let restored = Context::restore_from_checkpoint(checkpoint, rng, 10);
    assert_eq!(restored.sweep_count(), 1);
    assert_eq!(restored.attempted_updates(), 17);
    assert_eq!(restored.accepted_moves(), 9);
    assert_eq!(restored.event_time(), 2.5);
    assert_eq!(restored.simulation_clocks()[3].value(), 2.5);
}

#[test]
fn legacy_json_checkpoint_defaults_new_clocks_to_zero() {
    let checkpoint: ContextCheckpoint = serde_json::from_value(serde_json::json!({
        "sweep_count": 3,
        "thermalization_sweeps": 5,
        "thermalized": false,
        "phase": "Thermalization"
    }))
    .unwrap();
    assert_eq!(checkpoint.attempted_updates, 0);
    assert_eq!(checkpoint.accepted_moves, 0);
    assert_eq!(checkpoint.event_time, 0.0);
}

// ── Array measurements via Context ────────────────────────────────────────

#[test]
fn test_context_measure_array_basic() {
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut ctx = Context::new_with_binsize(rng, 0, 2);

    ctx.measure_array("Correlation", &[1.0, 2.0, 3.0]);
    ctx.measure_array("Correlation", &[4.0, 5.0, 6.0]);

    let estimates = ctx.finalize_measurements();
    assert!(estimates.contains_key("Correlation"));
}

#[test]
fn test_context_measure_array_finalize() {
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut ctx = Context::new_with_binsize(rng, 0, 1);

    // binsize=1, so each sample becomes one bin
    ctx.measure_array("Vec", &[10.0, 20.0]);
    ctx.measure_array("Vec", &[30.0, 40.0]);

    let estimates = ctx.finalize_measurements();
    let est = estimates.get("Vec").expect("Vec observable");
    assert_eq!(est.n_bins, 2);
}

// ── Measurement namespace ─────────────────────────────────────────────────

#[test]
fn test_context_set_measurement_namespace() {
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut ctx = Context::new_with_binsize(rng, 0, 5);

    ctx.set_measurement_namespace(Some("chain_0".into()));
    assert_eq!(ctx.measurement_namespace(), Some("chain_0"));

    ctx.measure("Energy", 1.0);
    ctx.set_measurement_namespace(None);
    ctx.measure("Energy", 2.0);

    let estimates = ctx.finalize_measurements();
    assert!(estimates.contains_key("chain_0/Energy"));
    assert!(estimates.contains_key("Energy"));
}

#[test]
fn test_context_measurement_namespace_default_none() {
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let ctx = Context::new(rng, 10);
    assert!(ctx.measurement_namespace().is_none());
}

#[test]
fn test_context_namespace_isolation() {
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut ctx = Context::new_with_binsize(rng, 0, 1);

    ctx.set_measurement_namespace(Some("chainA".into()));
    ctx.measure("Mag", 0.5);
    ctx.measure_array("Corr", &[1.0, 2.0]);

    ctx.set_measurement_namespace(Some("chainB".into()));
    ctx.measure("Mag", 0.9);

    ctx.set_measurement_namespace(None);

    let estimates = ctx.finalize_measurements();
    assert!(estimates.contains_key("chainA/Mag"));
    assert!(estimates.contains_key("chainB/Mag"));
    assert!(estimates.contains_key("chainA/Corr"));
    assert!(!estimates.contains_key("chainB/Corr"));
}

#[test]
fn test_context_empty_namespace_ignored() {
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut ctx = Context::new_with_binsize(rng, 0, 1);

    ctx.set_measurement_namespace(Some(String::new()));
    ctx.measure("Energy", 1.0);

    let estimates = ctx.finalize_measurements();
    assert!(estimates.contains_key("Energy"));
}
