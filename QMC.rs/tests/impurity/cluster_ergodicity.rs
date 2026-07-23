//! Ergodicity tests for the longitudinal spin-boson cluster solver.
//!
//! Verifies that independent runs from different RNG seeds converge to
//! the same thermal expectation values for ⟨Sz⟩ and kink count,
//! confirming that the cluster update is ergodic across spin-flip
//! sectors and bosonic configurations.

use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use qmc_rs::LongitudinalSpinBosonClusterQmc;

fn run_cluster(seed: u64) -> (f64, f64, f64, f64) {
    let mut params = Params::new();
    params.set("beta", 4.0);
    params.set("bath", "single");
    params.set("omega0", 1.0);
    params.set("g", 0.3);
    params.set("tunnelling", 0.5);
    params.set("epsilon", 0.0);

    let run = RunConfig {
        thermalization_sweeps: 2000,
        measurement_sweeps: 10000,
        binsize: 100,
        base_seed: seed,
        ..Default::default()
    };
    let results = Scheduler::new(RayonBackend::new(1), run)
        .run_one::<LongitudinalSpinBosonClusterQmc>(&params);
    let mag = results
        .get("ClusterMagnetizationSz")
        .expect("ClusterMagnetizationSz");
    let kinks = results.get("ClusterKinkCount").expect("ClusterKinkCount");
    (mag.mean, mag.stderr, kinks.mean, kinks.stderr)
}

/// z-score check against the pooled mean.
///
/// For each observable the per-seed means should scatter around the
/// pooled (inverse-variance weighted) mean with z-scores well within
/// ±4, and the mean |z| should be under 2.
fn assert_z_scores(name: &str, values: &[f64], stderrs: &[f64]) {
    let n = values.len() as f64;

    // Pooled mean (inverse-variance weighted)
    let mut sum_w = 0.0;
    let mut sum_wm = 0.0;
    for i in 0..values.len() {
        let w = 1.0 / stderrs[i].max(0.01).powi(2);
        sum_w += w;
        sum_wm += w * values[i];
    }
    let pooled_mean = sum_wm / sum_w;

    // Pooled variance (sample variance of the per-seed means)
    let pooled_var = values
        .iter()
        .map(|&v| (v - pooled_mean).powi(2))
        .sum::<f64>()
        / (n - 1.0);

    // z-scores
    let mut z_values: Vec<f64> = Vec::new();
    for i in 0..values.len() {
        let denom = (pooled_var / n + stderrs[i].max(0.01).powi(2)).sqrt();
        let z = (values[i] - pooled_mean) / denom.max(1e-10);
        z_values.push(z);
        assert!(
            z.abs() < 4.0,
            "{name}: z-score for seed {i} = {z:.4} exceeds 4σ (values: {values:?})",
        );
    }

    let mean_abs_z: f64 = z_values.iter().map(|z| z.abs()).sum::<f64>() / n;
    assert!(
        mean_abs_z < 2.0,
        "{name}: mean |z| = {mean_abs_z:.4} exceeds 2 (z-values: {z_values:?})",
    );
}

#[test]
fn cluster_ergodicity_multi_seed_convergence() {
    let seeds = [42u64, 123, 456, 789];
    let results: Vec<(f64, f64, f64, f64)> = seeds.iter().map(|&s| run_cluster(s)).collect();

    // ⟨Sz⟩ consistency: max−min < 4 × max(stderr)
    let sz_values: Vec<f64> = results.iter().map(|r| r.0).collect();
    let sz_stderrs: Vec<f64> = results.iter().map(|r| r.1).collect();
    let sz_spread = sz_values.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
        - sz_values.iter().cloned().fold(f64::INFINITY, f64::min);
    let sz_max_stderr = sz_stderrs.iter().cloned().fold(0.0_f64, f64::max);
    assert!(
        sz_spread < 4.0 * sz_max_stderr.max(0.01),
        "⟨Sz⟩ spread={sz_spread:.6} exceeds 4σ={:.6} (values: {sz_values:?})",
        4.0 * sz_max_stderr.max(0.01)
    );

    // Kink count consistency: max−min < 4 × max(stderr)
    let kink_values: Vec<f64> = results.iter().map(|r| r.2).collect();
    let kink_stderrs: Vec<f64> = results.iter().map(|r| r.3).collect();
    let kink_spread = kink_values
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max)
        - kink_values.iter().cloned().fold(f64::INFINITY, f64::min);
    let kink_max_stderr = kink_stderrs.iter().cloned().fold(0.0_f64, f64::max);
    assert!(
        kink_spread < 4.0 * kink_max_stderr.max(0.1),
        "kink count spread={kink_spread:.6} exceeds 4σ={:.6} (values: {kink_values:?})",
        4.0 * kink_max_stderr.max(0.1)
    );

    // z-score checks against pooled mean
    assert_z_scores("⟨Sz⟩", &sz_values, &sz_stderrs);
    assert_z_scores("kink count", &kink_values, &kink_stderrs);
}
