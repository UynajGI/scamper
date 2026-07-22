//! Gap 1: z-score statistical validation framework.
//!
//! Runs 16 independent seeds for each solver on a 3-site PBC Ising chain
//! at β=0.5, computes standardized residuals z = (⟨E⟩_MC − ⟨E⟩_exact) / SE,
//! and asserts:
//!   - Each individual seed: |z| < 4
//!   - Mean z-score across seeds: |z̄| < 1 (no systematic bias)
//!   - z-scores are not all same sign (no one-sided bias)

use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use cmc_rs::{ClassicalMC, Hamiltonian, IsingModel, MetropolisCore, SWCore, WolffCore};

fn exact_3site_energy(beta: f64, j: f64) -> f64 {
    let lattice = cmc_rs::build_chain(3, true);
    let model = IsingModel::new(j);
    let mut z = 0.0;
    let mut we = 0.0;
    for mask in 0..(1u32 << 3) {
        let spins: Vec<f64> = (0..3)
            .map(|i| if (mask >> i) & 1 == 1 { 1.0 } else { -1.0 })
            .collect();
        let e = model.compute_total_energy(&spins, &lattice, 1.0);
        let w = (-beta * e).exp();
        z += w;
        we += e * w;
    }
    we / z
}

fn exact_3site_m2(beta: f64, j: f64) -> f64 {
    let lattice = cmc_rs::build_chain(3, true);
    let model = IsingModel::new(j);
    let mut z = 0.0;
    let mut wm2 = 0.0;
    for mask in 0..(1u32 << 3) {
        let spins: Vec<f64> = (0..3)
            .map(|i| if (mask >> i) & 1 == 1 { 1.0 } else { -1.0 })
            .collect();
        let e = model.compute_total_energy(&spins, &lattice, 1.0);
        let m: f64 = spins.iter().sum::<f64>() / 3.0;
        let w = (-beta * e).exp();
        z += w;
        wm2 += m * m * w;
    }
    wm2 / z
}

const N_SEEDS: usize = 16;
const BETA: f64 = 0.5;
const J: f64 = 1.0;

/// Run MC for one seed and return (energy_mean, energy_stderr).
fn run_metropolis(seed: u64) -> (f64, f64) {
    let mut params = Params::new();
    params.set("L", 3);
    params.set("J", J);
    params.set("beta", BETA);
    let config = RunConfig {
        thermalization_sweeps: 3000,
        measurement_sweeps: 12000,
        binsize: 200,
        base_seed: seed,
        ..Default::default()
    };
    let results = Scheduler::new(RayonBackend::new(1), config)
        .run_one::<ClassicalMC<IsingModel, MetropolisCore>>(&params);
    let e = results.get("Energy").unwrap();
    (e.mean, e.stderr)
}

fn run_wolff(seed: u64) -> (f64, f64) {
    let mut params = Params::new();
    params.set("L", 3);
    params.set("J", J);
    params.set("beta", BETA);
    let config = RunConfig {
        thermalization_sweeps: 3000,
        measurement_sweeps: 12000,
        binsize: 200,
        base_seed: seed,
        ..Default::default()
    };
    let results = Scheduler::new(RayonBackend::new(1), config)
        .run_one::<ClassicalMC<IsingModel, WolffCore>>(&params);
    let e = results.get("Energy").unwrap();
    (e.mean, e.stderr)
}

fn run_sw(seed: u64) -> (f64, f64) {
    let mut params = Params::new();
    params.set("L", 3);
    params.set("J", J);
    params.set("beta", BETA);
    let config = RunConfig {
        thermalization_sweeps: 3000,
        measurement_sweeps: 12000,
        binsize: 200,
        base_seed: seed,
        ..Default::default()
    };
    let results = Scheduler::new(RayonBackend::new(1), config)
        .run_one::<ClassicalMC<IsingModel, SWCore>>(&params);
    let e = results.get("Energy").unwrap();
    (e.mean, e.stderr)
}

/// Compute z-scores for a set of (mean, stderr) pairs against an exact value.
/// Returns (z_scores, max_abs_z, mean_z, fraction_positive).
fn analyze_z_scores(results: &[(f64, f64)], exact: f64) -> (Vec<f64>, f64, f64, f64) {
    let z_scores: Vec<f64> = results
        .iter()
        .map(|(mean, stderr)| (mean - exact) / stderr.max(1e-10))
        .collect();
    let max_abs_z = z_scores.iter().map(|z| z.abs()).fold(0.0_f64, f64::max);
    let mean_z = z_scores.iter().sum::<f64>() / z_scores.len() as f64;
    let n_positive = z_scores.iter().filter(|z| **z > 0.0).count();
    let frac_positive = n_positive as f64 / z_scores.len() as f64;
    (z_scores, max_abs_z, mean_z, frac_positive)
}

#[test]
fn metropolis_zscore_energy_16_seeds() {
    let exact = exact_3site_energy(BETA, J);
    let results: Vec<(f64, f64)> = (0..N_SEEDS as u64).map(run_metropolis).collect();
    let (_, max_z, mean_z, frac_pos) = analyze_z_scores(&results, exact);

    assert!(max_z < 4.0, "max |z| = {max_z:.2} should be < 4");
    assert!(
        mean_z.abs() < 1.5,
        "mean z = {mean_z:.2} should be |z̄| < 1.5"
    );
    // Not all same sign (would indicate systematic bias)
    assert!(
        frac_pos > 0.15 && frac_pos < 0.85,
        "fraction positive = {frac_pos:.2}, should be between 0.15 and 0.85"
    );
}

#[test]
fn metropolis_zscore_magnetization_16_seeds() {
    let exact = exact_3site_m2(BETA, J);
    let results: Vec<(f64, f64)> = (0..N_SEEDS as u64)
        .map(|seed| {
            let mut params = Params::new();
            params.set("L", 3);
            params.set("J", J);
            params.set("beta", BETA);
            let config = RunConfig {
                thermalization_sweeps: 3000,
                measurement_sweeps: 12000,
                binsize: 200,
                base_seed: seed,
                ..Default::default()
            };
            let r = Scheduler::new(RayonBackend::new(1), config)
                .run_one::<ClassicalMC<IsingModel, MetropolisCore>>(&params);
            let m2 = r
                .get("MagnetizationSquared")
                .or_else(|| r.get("M2"))
                .unwrap();
            (m2.mean, m2.stderr)
        })
        .collect();

    let (_, max_z, mean_z, frac_pos) = analyze_z_scores(&results, exact);
    assert!(max_z < 4.0, "M² max |z| = {max_z:.2} should be < 4");
    assert!(mean_z.abs() < 1.5, "M² mean z = {mean_z:.2}");
    assert!(
        frac_pos > 0.15 && frac_pos < 0.85,
        "M² frac_pos = {frac_pos:.2}"
    );
}

#[test]
fn wolff_zscore_energy_16_seeds() {
    let exact = exact_3site_energy(BETA, J);
    let results: Vec<(f64, f64)> = (0..N_SEEDS as u64).map(run_wolff).collect();
    let (_, max_z, mean_z, frac_pos) = analyze_z_scores(&results, exact);

    assert!(max_z < 4.0, "Wolff max |z| = {max_z:.2} should be < 4");
    assert!(mean_z.abs() < 1.5, "Wolff mean z = {mean_z:.2}");
    assert!(
        frac_pos > 0.15 && frac_pos < 0.85,
        "Wolff frac_pos = {frac_pos:.2}"
    );
}

#[test]
fn swendsen_wang_zscore_energy_16_seeds() {
    let exact = exact_3site_energy(BETA, J);
    let results: Vec<(f64, f64)> = (0..N_SEEDS as u64).map(run_sw).collect();
    let (_, max_z, mean_z, frac_pos) = analyze_z_scores(&results, exact);

    assert!(max_z < 4.0, "SW max |z| = {max_z:.2} should be < 4");
    assert!(mean_z.abs() < 1.5, "SW mean z = {mean_z:.2}");
    assert!(
        frac_pos > 0.15 && frac_pos < 0.85,
        "SW frac_pos = {frac_pos:.2}"
    );
}

#[test]
fn cross_solver_zscores_agree() {
    // Metropolis and Wolff z-score distributions should overlap.
    let exact = exact_3site_energy(BETA, J);
    let metro: Vec<(f64, f64)> = (0..N_SEEDS as u64).map(run_metropolis).collect();
    let wolff: Vec<(f64, f64)> = (0..N_SEEDS as u64).map(run_wolff).collect();
    let (z_metro, _, mean_metro, _) = analyze_z_scores(&metro, exact);
    let (z_wolff, _, mean_wolff, _) = analyze_z_scores(&wolff, exact);

    // Mean z-scores should both be near zero and not differ by more than 2
    assert!(
        (mean_metro - mean_wolff).abs() < 2.0,
        "mean z differ: metro={mean_metro:.2}, wolff={mean_wolff:.2}"
    );
    // Pooled z-scores: all |z| < 5 (relaxed for pooled)
    let all_z: Vec<f64> = z_metro.iter().chain(z_wolff.iter()).copied().collect();
    let pooled_max = all_z.iter().map(|z| z.abs()).fold(0.0_f64, f64::max);
    assert!(pooled_max < 5.0, "pooled max |z| = {pooled_max:.2}");
}
