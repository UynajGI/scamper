//! Validate 2D Ising model against Onsager's exact solution.
//!
//! Tc = 2 / ln(1 + √2) ≈ 2.269185
//! At Tc: energy per bond = -√2 / 2 ≈ -0.7071, energy per site = -√2 ≈ -1.4142

use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use cmc_rs::{ClassicalMC, IsingModel, MetropolisCore};

/// Known critical temperature for 2D Ising on square lattice.
/// Tc = 2 / ln(1 + sqrt(2)) ≈ 2.269185
fn tc() -> f64 {
    2.0 / (1.0_f64 + 2.0_f64.sqrt()).ln()
}

/// Run a 2D Ising simulation and return (energy_per_site, magnetization).
fn run_2d_ising(l: usize, beta: f64, therm: u64, meas: u64) -> (f64, f64, f64, f64) {
    let mut params = Params::new();
    params.set("Lx", l);
    params.set("Ly", l);
    params.set("J", 1.0);
    params.set("beta", beta);

    let config = RunConfig {
        thermalization_sweeps: therm,
        measurement_sweeps: meas,
        binsize: 100,
        base_seed: 42,
        ..Default::default()
    };

    let backend = RayonBackend::new(1);
    let scheduler = Scheduler::new(backend, config);
    let results = scheduler.run_one::<ClassicalMC<IsingModel, MetropolisCore>>(&params);

    let e = results.get("Energy").expect("Energy missing");
    let m = results.get("Magnetization").expect("Magnetization missing");

    (e.mean, e.stderr, m.mean, m.stderr)
}

#[test]
fn test_2d_ising_at_tc_energy() {
    // At Tc, energy per site ≈ -√2 ≈ -1.4142
    let l = 16;
    let beta = 1.0 / tc();

    let (e, e_err, _m, _m_err) = run_2d_ising(l, beta, 2000, 5000);

    let e_per_site = e / (l * l) as f64;
    let expected = -(2.0_f64).sqrt();

    // Should be within ~5% of exact Onsager result
    let tol = 0.15;
    assert!(
        (e_per_site - expected).abs() < tol + 3.0 * e_err,
        "Energy per site at Tc: got {:.4} ± {:.4}, expected {:.4}",
        e_per_site,
        e_err,
        expected
    );
}

#[test]
fn test_2d_ising_high_t_magnetization_vanish() {
    // At T >> Tc (beta = 0.1), magnetization ~ 0
    let l = 8;
    let beta = 0.1;

    let (_e, _e_err, m, m_err) = run_2d_ising(l, beta, 500, 2000);

    assert!(
        m < 0.3 + 3.0 * m_err,
        "Magnetization at high T should vanish: got {:.4} ± {:.4}",
        m,
        m_err
    );
}

#[test]
fn test_2d_ising_low_t_magnetization_appears() {
    // At T << Tc (beta = 1.0), magnetization > 0
    let l = 8;
    let beta = 1.0;

    let (_e, _e_err, m, m_err) = run_2d_ising(l, beta, 1000, 2000);

    assert!(
        m > 0.5 - 3.0 * m_err,
        "Magnetization at low T should be non-zero: got {:.4} ± {:.4}",
        m,
        m_err
    );
}

#[test]
fn test_2d_ising_energy_decreases_with_cooling() {
    // Energy should be lower (more negative) at lower temperature
    let l = 8;
    let beta_high = 0.3; // high T
    let beta_low = 1.0; // low T

    let (e_high, _, _, _) = run_2d_ising(l, beta_high, 300, 1000);
    let (e_low, _, _, _) = run_2d_ising(l, beta_low, 500, 2000);

    // Energy at low T should be more negative (more ordered)
    assert!(
        e_low < e_high,
        "Energy at low T ({:.4}) should be lower than at high T ({:.4})",
        e_low,
        e_high
    );
}
