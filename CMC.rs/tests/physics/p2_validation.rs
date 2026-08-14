//! P2 validation tests: additional correctness checks for solvers
//! that already have basic validation but need more coverage.

use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use cmc_rs::{ClassicalMC, Hamiltonian, IsingModel, MultiSpinIsing};

// ═══════════════════════════════════════════════════════════════════════
// P2.1: MultiSpinIsing exact energy for 8-spin system
// ═══════════════════════════════════════════════════════════════════════

/// Exact ⟨E⟩ for small Ising chain by enumeration.
fn exact_energy(n: usize, j: f64, beta: f64, pbc: bool) -> f64 {
    let lattice = if pbc {
        cmc_rs::build_chain(n, true)
    } else {
        cmc_rs::build_chain(n, false)
    };
    let model = IsingModel::new(j);
    let mut z = 0.0;
    let mut we = 0.0;
    for mask in 0..(1u32 << n) {
        let spins: Vec<f64> = (0..n)
            .map(|i| if (mask >> i) & 1 == 1 { 1.0 } else { -1.0 })
            .collect();
        let e = model.compute_total_energy(&spins, &lattice, 1.0);
        let w = (-beta * e).exp();
        z += w;
        we += e * w;
    }
    we / z
}

#[test]
fn exact_enumeration_helper_is_self_consistent() {
    // Sanity check for the exact_energy helper used by other tests.
    // AFM Ising chain at finite β should have negative energy.
    let exact = exact_energy(3, 1.0, 0.5, true);
    assert!(exact < 0.0, "AFM energy should be negative, got {exact}");
    // Higher β should give lower (more negative) energy
    let cold = exact_energy(3, 1.0, 2.0, true);
    assert!(
        cold < exact,
        "colder β=2 energy {cold} should be < warmer β=0.5 energy {exact}"
    );
}

// P2.2: HybridCore needs explicit construction (no Default impl) —
// tested via integration/usage.rs smoke test instead.

#[test]
fn metropolis_8site_energy_matches_exact_enumeration() {
    // 8-site Ising chain: Metropolis MC vs exact 256-state enumeration.
    // Note: MultiSpinIsing is NOT tested here — it uses a different adapter
    // pattern (ParallelTemperingCompatible) and needs its own test.
    let exact = exact_energy(8, 1.0, 0.5, true);

    // 8 spins, β=0.5, J=1, PBC: exact ⟨E⟩ via 256-state enumeration
    // The exact value should be negative (AFM at warm T).
    assert!(exact < 0.0, "AFM energy should be negative, got {exact}");

    // Cross-check: Metropolis on same system should match
    let mut params = Params::new();
    params.set("L", 8);
    params.set("J", 1.0);
    params.set("beta", 0.5);
    let config = RunConfig {
        thermalization_sweeps: 5000,
        measurement_sweeps: 20000,
        binsize: 500,
        base_seed: 42,
        ..Default::default()
    };
    let scheduler = Scheduler::new(RayonBackend::new(1), config);
    let results = scheduler.run_one::<ClassicalMC<IsingModel, cmc_rs::MetropolisCore>>(&params);
    let e = results.get("Energy").expect("Energy");
    assert!(
        (e.mean - exact).abs() < 0.15,
        "Metropolis 8-site E={:.4}, exact={:.4}",
        e.mean,
        exact
    );
}

#[test]
fn multispin_ising_8site_energy_matches_exact_enumeration() {
    // 8-site Ising chain: MultiSpinIsing (64 packed replicas) vs exact
    // 256-state enumeration.  MultiSpinIsing implements both `MonteCarlo`
    // and `FromParams`, so it runs through the standard Carlo scheduler.
    // The scalar "Energy" observable is replica 0's energy; with enough
    // thermalization and measurement sweeps it converges to ⟨E⟩.
    let exact = exact_energy(8, 1.0, 0.5, true);
    assert!(exact < 0.0, "AFM energy should be negative, got {exact}");

    let mut params = Params::new();
    params.set("L", 8);
    params.set("J", 1.0);
    params.set("beta", 0.5);
    let config = RunConfig {
        thermalization_sweeps: 5000,
        measurement_sweeps: 20000,
        binsize: 500,
        base_seed: 42,
        ..Default::default()
    };
    let scheduler = Scheduler::new(RayonBackend::new(1), config);
    let results = scheduler.run_one::<MultiSpinIsing>(&params);
    let e = results.get("Energy").expect("Energy");
    assert!(
        (e.mean - exact).abs() < 0.15,
        "MultiSpinIsing 8-site E={:.4}, exact={:.4}",
        e.mean,
        exact
    );
}
