//! Physical validation tests for SSE algorithm.

use qmc_rs::{
    MonteCarlo, Params, RayonBackend, RunConfig, Scheduler,
    HeisenbergModel, SSECore,
};
use qmc_rs::lattice::builders::build_chain;
use rand_xoshiro::rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

/// Bethe ansatz ground state energy per site for 1D Heisenberg chain.
/// E_0/N = 1/4 - ln(2)
const BETHE_ANSATZ_ENERGY: f64 = 0.25 - std::f64::consts::LN_2; // ~ -0.443147

/// Verifies that the SSE algorithm produces the correct energy for the 1D
/// Heisenberg chain, matching the Bethe ansatz ground state energy.
#[test]
fn test_heisenberg_chain_ground_state() {
    let n_sites = 16;
    let beta = 10.0;

    let mut params = Params::new();
    params.set("L", n_sites);
    params.set("beta", beta);
    params.set("J", 1.0);
    params.set("pbc", true);

    let backend = RayonBackend::new(1);
    let config = RunConfig {
        thermalization_sweeps: 10000,
        measurement_sweeps: 50000,
        binsize: 1000,
        base_seed: 42,
        progress_interval: 0,
        checkpoint_interval: 0,
    };
    let scheduler = Scheduler::new(backend, config);

    let results = scheduler.run_one::<SSECore<HeisenbergModel>>(&params);

    if let Some(energy) = results.get("Energy") {
        let expected = BETHE_ANSATZ_ENERGY;
        let tolerance = 3.0 * energy.stderr;

        println!("Energy: {:.6} +/- {:.6}", energy.mean, energy.stderr);
        println!("Expected (Bethe ansatz): {:.6}", expected);
        println!("Difference: {:.6}", (energy.mean - expected).abs());
        println!("Tolerance (3sigma): {:.6}", tolerance);

        assert!(
            (energy.mean - expected).abs() < tolerance,
            "Energy {:.6} not within {:.6} (3sigma) of expected {:.6}",
            energy.mean, tolerance, expected
        );
    } else {
        panic!("Energy not found in results");
    }
}

#[test]
fn test_sse_monte_carlo_trait() {
    let lattice = build_chain(4, true);
    let model = HeisenbergModel::new(lattice, 1.0, 1.0);
    let core = SSECore::new(model);

    fn assert_monte_carlo<M: MonteCarlo>(_: &M) {}
    assert_monte_carlo(&core);
}

#[test]
fn test_diagnostic_operator_count() {
    let n_sites = 16;
    let beta = 10.0;

    let lattice = build_chain(n_sites, true);
    let model = HeisenbergModel::new(lattice, beta, 1.0);
    let mut core = SSECore::new(model);

    println!("Diagnostic: max_length = {}", core.engine.op_seq.max_length);
    println!("Diagnostic: bond_list len = {}", core.engine.bond_list.len());
    println!("Diagnostic: weights = {:?}", core.engine.weights);
    println!("Diagnostic: diagonal_shift = {}", core.engine.diagonal_shift);
    let n_bonds = core.engine.bond_list.len();

    use qmc_rs::Context;
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut ctx = Context::new(rng, 100);

    for i in 0..50 {
        core.sweep(&mut ctx);
        if (i + 1) <= 5 || (i + 1) % 10 == 0 {
            let energy = core.engine.compute_energy();
            let aligned: usize = core.engine.bond_list.iter()
                .filter(|(i, j, _)| core.engine.spins[*i] == core.engine.spins[*j])
                .count();
            println!("Diagnostic: sweep {} -> n_ops = {}, energy = {:.4}, aligned = {}/{}",
                     i + 1, core.engine.op_seq.n_operators, energy, aligned, n_bonds);
        }
    }
}

#[test]
fn test_debug_loop_update() {
    let n_sites = 8;
    let beta = 4.0;

    let lattice = build_chain(n_sites, true);
    let model = HeisenbergModel::new(lattice, beta, 1.0);
    let mut core = SSECore::new(model);

    use qmc_rs::Context;
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut ctx = Context::new(rng, 100);

    for i in 0..50 {
        core.sweep(&mut ctx);
        if i < 5 || (i + 1) % 10 == 0 {
            let energy = core.engine.compute_energy();
            let aligned: usize = core.engine.bond_list.iter()
                .filter(|(si, sj, _)| core.engine.spins[*si] == core.engine.spins[*sj])
                .count();
            let n_diag: usize = core.engine.op_seq.vertices.iter()
                .filter(|v| v.vertex_idx >= 1 && v.vertex_idx <= 4)
                .count();
            let n_offdiag: usize = core.engine.op_seq.vertices.iter()
                .filter(|v| v.vertex_idx == 5 || v.vertex_idx == 6)
                .count();
            println!("Sweep {}: n_ops={}, n_diag={}, n_offdiag={}, aligned={}/{}, E={:.4}",
                     i + 1, core.engine.op_seq.n_operators, n_diag, n_offdiag, 
                     aligned, core.engine.bond_list.len(), energy);
        }
    }
    
    // Check that we have off-diagonal operators
    let n_offdiag: usize = core.engine.op_seq.vertices.iter()
        .filter(|v| v.vertex_idx == 5 || v.vertex_idx == 6)
        .count();
    assert!(n_offdiag > 0, "Should have off-diagonal operators after sweeps, but got none");
    
    // Check that spins are not all the same
    let n_up: usize = core.engine.spins.iter().filter(|&&s| s == 0).count();
    assert!(n_up > 0 && n_up < n_sites, "Spins should be mixed, not all up or all down");
}
