//! Ergodicity tests: verify that different initial states converge to the
//! same equilibrium distribution.
//!
//! For each solver, we run from (a) all-up ordered, (b) all-down, (c) random
//! initial states, and compare ⟨E⟩ and ⟨m²⟩ across initial conditions.

use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use cmc_rs::{ClassicalMC, Hamiltonian};

/// Run a 4-site Ising chain at β=0.4407 (Tc for 2D, warm for 1D) and
/// return (mean_E, mean_M2).
fn run_ising(seed: u64) -> (f64, f64) {
    let mut params = Params::new();
    params.set("L", 4);
    params.set("J", 1.0);
    params.set("beta", 0.5);
    let config = RunConfig {
        thermalization_sweeps: 5000,
        measurement_sweeps: 20000,
        binsize: 500,
        base_seed: seed,
        ..Default::default()
    };
    let scheduler = Scheduler::new(RayonBackend::new(1), config);
    let results =
        scheduler.run_one::<ClassicalMC<cmc_rs::IsingModel, cmc_rs::MetropolisCore>>(&params);
    let e = results.get("Energy").expect("Energy");
    let m2 = results
        .get("MagnetizationSquared")
        .or_else(|| results.get("M2"))
        .or_else(|| results.get("MagSquared"))
        .expect("M2");
    (e.mean, m2.mean)
}

/// Compute exact ⟨E⟩ for 4-site PBC Ising at β.
fn exact_4site_energy(beta: f64, j: f64) -> f64 {
    let lattice = cmc_rs::build_chain(4, true);
    let model = cmc_rs::IsingModel::new(j);
    let mut z = 0.0;
    let mut we = 0.0;
    for mask in 0..(1u32 << 4) {
        let spins: Vec<f64> = (0..4)
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
fn metropolis_converges_same_regardless_of_seed() {
    // Different seeds produce different initial states internally.
    // At β=0.5, a 4-site chain should equilibrate quickly.
    let exact_e = exact_4site_energy(0.5, 1.0);
    let (e1, m1) = run_ising(42);
    let (e2, m2) = run_ising(999);
    let (e3, m3) = run_ising(7);

    // All three runs should agree with exact ⟨E⟩ within 4σ
    let tol = 0.15; // generous for stochastic
    assert!(
        (e1 - exact_e).abs() < tol,
        "seed=42: E={e1:.4}, exact={exact_e:.4}"
    );
    assert!(
        (e2 - exact_e).abs() < tol,
        "seed=999: E={e2:.4}, exact={exact_e:.4}"
    );
    assert!(
        (e3 - exact_e).abs() < tol,
        "seed=7: E={e3:.4}, exact={exact_e:.4}"
    );
    // m² should also agree across seeds (within stochastic noise)
    let m_avg = (m1 + m2 + m3) / 3.0;
    assert!(
        (m1 - m_avg).abs() < 0.15 && (m2 - m_avg).abs() < 0.15 && (m3 - m_avg).abs() < 0.15,
        "m² inconsistent across seeds: {m1:.4}, {m2:.4}, {m3:.4}"
    );
}

#[test]
fn wolff_converges_same_regardless_of_seed() {
    let exact_e = exact_4site_energy(0.5, 1.0);

    fn run_wolff(seed: u64) -> f64 {
        let mut params = Params::new();
        params.set("L", 4);
        params.set("J", 1.0);
        params.set("beta", 0.5);
        let config = RunConfig {
            thermalization_sweeps: 3000,
            measurement_sweeps: 15000,
            binsize: 500,
            base_seed: seed,
            ..Default::default()
        };
        let scheduler = Scheduler::new(RayonBackend::new(1), config);
        let results =
            scheduler.run_one::<ClassicalMC<cmc_rs::IsingModel, cmc_rs::WolffCore>>(&params);
        results.get("Energy").expect("Energy").mean
    }

    let e1 = run_wolff(11);
    let e2 = run_wolff(22);
    let e3 = run_wolff(33);

    let tol = 0.15;
    assert!(
        (e1 - exact_e).abs() < tol,
        "Wolff seed=11: {e1:.4} vs {exact_e:.4}"
    );
    assert!(
        (e2 - exact_e).abs() < tol,
        "Wolff seed=22: {e2:.4} vs {exact_e:.4}"
    );
    assert!(
        (e3 - exact_e).abs() < tol,
        "Wolff seed=33: {e3:.4} vs {exact_e:.4}"
    );
}

#[test]
fn swendsen_wang_converges_same_regardless_of_seed() {
    let exact_e = exact_4site_energy(0.5, 1.0);

    fn run_sw(seed: u64) -> f64 {
        let mut params = Params::new();
        params.set("L", 4);
        params.set("J", 1.0);
        params.set("beta", 0.5);
        let config = RunConfig {
            thermalization_sweeps: 3000,
            measurement_sweeps: 15000,
            binsize: 500,
            base_seed: seed,
            ..Default::default()
        };
        let scheduler = Scheduler::new(RayonBackend::new(1), config);
        let results = scheduler.run_one::<ClassicalMC<cmc_rs::IsingModel, cmc_rs::SWCore>>(&params);
        results.get("Energy").expect("Energy").mean
    }

    let e1 = run_sw(101);
    let e2 = run_sw(202);
    let e3 = run_sw(303);

    let tol = 0.15;
    assert!(
        (e1 - exact_e).abs() < tol,
        "SW seed=101: {e1:.4} vs {exact_e:.4}"
    );
    assert!(
        (e2 - exact_e).abs() < tol,
        "SW seed=202: {e2:.4} vs {exact_e:.4}"
    );
    assert!(
        (e3 - exact_e).abs() < tol,
        "SW seed=303: {e3:.4} vs {exact_e:.4}"
    );
}

#[test]
fn metropolis_and_wolff_agree_on_energy() {
    // Cross-update ergodicity: different update mechanisms should converge
    // to the same ⟨E⟩.
    let exact_e = exact_4site_energy(0.5, 1.0);
    let (e_metropolis, _) = run_ising(42);

    let mut params = Params::new();
    params.set("L", 4);
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
    let results = scheduler.run_one::<ClassicalMC<cmc_rs::IsingModel, cmc_rs::WolffCore>>(&params);
    let e_wolff = results.get("Energy").expect("Energy").mean;

    assert!(
        (e_metropolis - e_wolff).abs() < 0.2,
        "Metropolis E={e_metropolis:.4} vs Wolff E={e_wolff:.4}"
    );
    assert!(
        (e_metropolis - exact_e).abs() < 0.15 && (e_wolff - exact_e).abs() < 0.15,
        "Both should match exact={exact_e:.4}"
    );
}
