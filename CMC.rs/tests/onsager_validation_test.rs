//! Onsager exact solution validation for 2D Ising model.
//!
//! Tests the 2D Ising model against the Onsager exact solution for the
//! infinite square lattice.

use cmc_rs::*;

// ─── Helpers ────────────────────────────────────────────────────────────────

fn run_simulation<T: MonteCarlo + FromParams>(
    params: &Params,
    thermalization: u64,
    measurements: u64,
    binsize: usize,
) -> Results {
    let backend = RayonBackend::new(1);
    let config = RunConfig {
        thermalization_sweeps: thermalization,
        measurement_sweeps: measurements,
        binsize,
        base_seed: 42,
        progress_interval: 0,
        checkpoint_interval: 0,
    };
    let scheduler = Scheduler::new(backend, config);
    scheduler.run_one::<T>(params)
}

fn make_params_2d(lx: usize, ly: usize, beta: f64, j: f64, pbc: bool) -> Params {
    let mut p = Params::new();
    p.set("Lx", lx);
    p.set("Ly", ly);
    p.set("beta", beta);
    p.set("J", j);
    p.set("pbc", pbc);
    p
}

// Onsager 2D Ising critical temperature: βc = ln(1+√2)/(2J)
fn ising_2d_critical_beta(j: f64) -> f64 {
    (1.0 + 2.0_f64.sqrt()).ln() / (2.0 * j)
}

// E/N at Tc = -coth(2*βc*J) = -√2 ≈ -1.41421 * J.
// At Tc, sinh(2βcJ) = 1, cosh(2βcJ) = √2, so coth(2βcJ) = √2.
// The full Onsager formula has a term (2/π)(2 tanh²(2βJ)-1)K(k₁) which
// vanishes exactly at Tc since 2 tanh²(2βcJ) - 1 = 0, leaving -coth(2βcJ).
fn ising_2d_energy_per_site_at_tc(j: f64) -> f64 {
    -2.0_f64.sqrt() * j
}

// ─── Test: 2D Ising at Tc with Metropolis ────────────────────────────────────

#[test]
fn test_2d_ising_at_tc_onsager() {
    // Run Metropolis on 16x16 PBC at βc — energy should match Onsager
    let j = 1.0;
    let beta_c = ising_2d_critical_beta(j);
    let expected_e_per_site = ising_2d_energy_per_site_at_tc(j);

    let params = make_params_2d(16, 16, beta_c, j, true);

    let results = run_simulation::<MetropolisCore<IsingModel2D>>(
        &params, 10000, 20000, 100,
    );

    let energy = results.get("Energy").unwrap();
    let e_per_site = energy.mean / 256.0; // 16*16 sites

    assert!(
        (e_per_site - expected_e_per_site).abs() < 0.05 * expected_e_per_site.abs(),
        "Energy per site {:.6} differs from Onsager value {:.6} by more than 5%",
        e_per_site,
        expected_e_per_site
    );
}

// ─── Test: 2D Ising temperature sweep ────────────────────────────────────────

#[test]
fn test_2d_ising_temperature_sweep() {
    // Run at several β values and verify monotonic energy decrease
    let j = 1.0;
    let betas = [0.2, 0.3, 0.4, 0.4407, 0.5, 0.6, 0.8];
    let mut energies = Vec::new();

    for &beta in &betas {
        let params = make_params_2d(16, 16, beta, j, true);
        let results = run_simulation::<MetropolisCore<IsingModel2D>>(
            &params, 5000, 10000, 100,
        );
        let e_per_site = results.get("Energy").unwrap().mean / 256.0;
        energies.push(e_per_site);
    }

    // At β=0.2 (high T), |E/N| should be small (< 0.5*J)
    assert!(
        energies[0].abs() < 0.5 * j,
        "At β=0.2 (high T), |E/N| = {} should be < 0.5",
        energies[0].abs()
    );

    // At β=0.8 (low T), E/N should be close to ground state E/N = -2.0*J
    let gs_per_site = -2.0 * j;
    assert!(
        (energies[energies.len() - 1] - gs_per_site).abs() < 0.15 * gs_per_site.abs(),
        "At β=0.8 (low T), E/N = {:.4} should be near ground state {:.4}",
        energies[energies.len() - 1],
        gs_per_site
    );

    // Energy should be monotonically decreasing with β (within statistical tolerance)
    for i in 1..energies.len() {
        assert!(
            energies[i] < energies[i - 1] + 0.3,
            "Energy should generally decrease with β at index {}: {:.4} vs {:.4}",
            i, energies[i], energies[i - 1]
        );
    }

    // Near Tc (β=0.4407), energy should be in the expected range
    let tc_idx = 3; // β=0.4407
    let expected_e_per_site = ising_2d_energy_per_site_at_tc(j);
    assert!(
        (energies[tc_idx] - expected_e_per_site).abs() < 0.10 * expected_e_per_site.abs(),
        "At Tc, E/N = {:.4} differs from Onsager {:.4} by more than 10%",
        energies[tc_idx],
        expected_e_per_site
    );
}

// ─── Test: 2D Ising at Tc with Wolff algorithm ──────────────────────────────

#[test]
fn test_2d_ising_wolff_at_tc() {
    // Same setup but with Wolff — should converge faster at Tc
    // due to reduced critical slowing down
    let j = 1.0;
    let beta_c = ising_2d_critical_beta(j);
    let expected_e_per_site = ising_2d_energy_per_site_at_tc(j);

    let params = make_params_2d(16, 16, beta_c, j, true);

    let results = run_simulation::<WolffCore<IsingModel2D>>(
        &params, 10000, 20000, 100,
    );

    let energy = results.get("Energy").unwrap();
    let e_per_site = energy.mean / 256.0; // 16*16 sites

    assert!(
        (e_per_site - expected_e_per_site).abs() < 0.05 * expected_e_per_site.abs(),
        "Wolff: Energy per site {:.6} differs from Onsager value {:.6} by more than 5%",
        e_per_site,
        expected_e_per_site
    );
}
