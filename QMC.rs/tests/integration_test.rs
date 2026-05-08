//! Integration test with Carlo.rs framework

use qmc_rs::{
    Context, MonteCarlo, Params, RayonBackend, RunConfig, Scheduler,
    HeisenbergModel, SSECore,
};
use qmc_rs::lattice::builders::build_chain;
use rand_xoshiro::Xoshiro256PlusPlus;
use rand_xoshiro::rand_core::SeedableRng;

#[test]
fn test_heisenberg_with_carlo_scheduler() {
    // Build lattice
    let lattice = build_chain(8, true);

    // Create model
    let model = HeisenbergModel::new(lattice, 1.0, 1.0);

    // Create SSE core
    let core = SSECore::new(model);

    // Verify MonteCarlo trait is implemented
    fn assert_monte_carlo<M: MonteCarlo>(_: &M) {}
    assert_monte_carlo(&core);

    // Create scheduler
    let backend = RayonBackend::new(1);
    let config = RunConfig {
        thermalization_sweeps: 100,
        measurement_sweeps: 500,
        binsize: 50,
        base_seed: 42,
        progress_interval: 0,
        checkpoint_interval: 0,
    };
    let scheduler = Scheduler::new(backend, config);

    // Run simulation with required params
    let mut params = Params::new();
    params.set("L", "8");
    params.set("beta", "1.0");
    let results = scheduler.run_one::<SSECore<HeisenbergModel>>(&params);

    // Check that Energy was measured
    let energy = results.get("Energy");
    assert!(energy.is_some(), "Energy should be measured");

    // Check that Magnetization was measured
    let mag = results.get("Magnetization");
    assert!(mag.is_some(), "Magnetization should be measured");

    println!("Energy: {:?}", energy);
    println!("Magnetization: {:?}", mag);
}

#[test]
fn test_sse_core_sweep() {
    let lattice = build_chain(4, true);
    let model = HeisenbergModel::new(lattice, 1.0, 1.0);
    let mut core = SSECore::new(model);

    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut ctx = Context::new(rng, 100);

    // Run a few sweeps
    for _ in 0..10 {
        core.sweep(&mut ctx);
    }

    // Should have advanced sweep count
    assert!(ctx.sweep_count() >= 10);
}

#[test]
fn test_sse_core_measure() {
    let lattice = build_chain(4, true);
    let model = HeisenbergModel::new(lattice, 1.0, 1.0);
    let mut core = SSECore::new(model);

    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut ctx = Context::new(rng, 100);

    // Run sweep and measure
    core.sweep(&mut ctx);
    core.measure(&mut ctx);

    // Context should have measurements
    // Note: ctx.measurements() returns reference to internal accumulator
}