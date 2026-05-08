use carlo_rs::Context;
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

    // Still not thermalized (sweeps > thermalization_sweeps)
    assert!(!ctx.is_thermalized());

    // One more sweep
    ctx.advance_sweep();
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
