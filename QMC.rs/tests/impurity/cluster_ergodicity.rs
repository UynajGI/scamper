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
}
