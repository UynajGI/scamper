// File: Carlo.rs/tests/checkpoint_test.rs
use carlo_rs::Context;
use rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

#[test]
fn test_context_checkpoint_restore() {
    let rng1 = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut ctx = Context::new(rng1, 100);

    // Advance sweeps
    for _ in 0..150 {
        ctx.advance_sweep();
    }

    let state = ctx.checkpoint_state();
    assert!(state.thermalized);
    assert_eq!(state.sweep_count, 150);

    // Restore in new context
    let rng2 = Xoshiro256PlusPlus::seed_from_u64(123);
    let restored = Context::restore_from_checkpoint(state, rng2, 100);
    assert_eq!(restored.sweep_count(), 150);
    assert!(restored.is_thermalized());
}

#[test]
fn test_register_observable() {
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut ctx = Context::new(rng, 100);

    ctx.register_observable("custom_obs", 50);

    // Measure after registration
    ctx.measure("custom_obs", 1.0);
    ctx.measure("custom_obs", 2.0);
}

#[test]
fn test_register_observable_with_shape() {
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut ctx = Context::new(rng, 100);

    // Shape hint is currently ignored but should work
    ctx.register_observable_with_shape("array_obs", 100, &[3, 3]);

    ctx.measure("array_obs", 1.0);
}

#[test]
fn test_thermalization_sweeps_accessor() {
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let ctx = Context::new(rng, 500);

    assert_eq!(ctx.thermalization_sweeps(), 500);
}
