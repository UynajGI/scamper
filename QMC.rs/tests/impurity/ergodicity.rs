//! Impurity ergodicity tests: verify convergence from different initial states.
//!
//! The wormhole solver starts from an empty worldline. This test verifies
//! that the solver produces consistent results regardless of seed,
//! which implicitly tests sector accessibility.
//!
//! The z-score test honours `SCUTTLE_ZSCORE_SEEDS=<n>` for nightly
//! high-power monitoring (unset → the default 4 seeds, unchanged for CI).

use crate::zscore_seeds::zscore_seeds;
use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use qmc_rs::impurity::ImpurityQmc;

fn run_rabi(seed: u64) -> (f64, f64, f64, f64) {
    let mut params = Params::new();
    params.set("beta", 8.0);
    params.set("model", "rabi");
    params.set("bath", "single");
    params.set("omega0", 1.0);
    params.set("g", 0.3);
    params.set("tunnelling", 1.0);
    params.set("h_z", 0.0);
    let run = RunConfig {
        thermalization_sweeps: 3000,
        measurement_sweeps: 12000,
        binsize: 100,
        base_seed: seed,
        ..Default::default()
    };
    let results = Scheduler::new(RayonBackend::new(1), run).run_one::<ImpurityQmc>(&params);
    let mag = results
        .get("MagnetizationSigmaZ")
        .expect("MagnetizationSigmaZ");
    let order = results.get("ExpansionOrder").expect("ExpansionOrder");
    (mag.mean, mag.stderr, order.mean, order.stderr)
}

#[test]
fn wormhole_ergodicity_multi_seed_convergence() {
    // Run from 4 different seeds. All should converge to the same ⟨σz⟩ and ⟨n⟩.
    let seeds = [42u64, 123, 777, 2026];
    let results: Vec<(f64, f64, f64, f64)> = seeds.iter().map(|&s| run_rabi(s)).collect();

    // Check magnetization consistency
    let mag_mean: f64 = results.iter().map(|r| r.0).sum::<f64>() / results.len() as f64;
    for (i, &(mag, stderr, _, _)) in results.iter().enumerate() {
        assert!(
            (mag - mag_mean).abs() < 4.0 * stderr.max(0.02),
            "Seed {}: ⟨σz⟩={mag:.4}±{stderr:.4}, mean={mag_mean:.4}",
            seeds[i]
        );
    }

    // Check expansion order consistency
    let order_mean: f64 = results.iter().map(|r| r.2).sum::<f64>() / results.len() as f64;
    for (i, &(_, _, order, stderr)) in results.iter().enumerate() {
        assert!(
            (order - order_mean).abs() < 4.0 * stderr.max(0.1),
            "Seed {}: ⟨n⟩={order:.4}±{stderr:.4}, mean={order_mean:.4}",
            seeds[i]
        );
    }
}

#[test]
fn wormhole_ergodicity_zscore_4_seeds() {
    // z-score framework: 4 seeds (default), check |z| < 4 for magnetization
    let seeds = zscore_seeds(&[42u64, 123, 777, 2026]);
    let results: Vec<(f64, f64)> = seeds
        .iter()
        .map(|&s| {
            let (mag, stderr, _, _) = run_rabi(s);
            (mag, stderr)
        })
        .collect();

    let mean: f64 = results.iter().map(|r| r.0).sum::<f64>() / results.len() as f64;
    let z_scores: Vec<f64> = results
        .iter()
        .map(|(mag, stderr)| (mag - mean) / stderr.max(1e-10))
        .collect();

    for (i, &z) in z_scores.iter().enumerate() {
        assert!(
            z.abs() < 4.0,
            "Seed {}: z-score = {z:.2}, should be |z| < 4",
            seeds[i]
        );
    }

    // Mean z-score should be near zero (no systematic bias)
    let mean_z: f64 = z_scores.iter().sum::<f64>() / z_scores.len() as f64;
    assert!(
        mean_z.abs() < 2.0,
        "Mean z-score = {mean_z:.2}, should be |z̄| < 2"
    );
}
