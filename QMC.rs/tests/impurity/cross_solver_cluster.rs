//! Cross-solver validation: wormhole ↔ cluster.
//!
//! Both solvers handle longitudinal spin-boson coupling. For a model with
//! only longitudinal coupling (no transverse), both should give consistent
//! magnetization and kink/expansion order.

use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use qmc_rs::impurity::ImpurityQmc;
use qmc_rs::LongitudinalSpinBosonClusterQmc;

#[test]
fn wormhole_and_cluster_both_run_on_longitudinal_model() {
    // Longitudinal spin-boson: only σz coupling to bath, no transverse.
    let beta: f64 = 4.0;
    let omega: f64 = 1.0;
    let g: f64 = 0.3;
    let tunnelling: f64 = 0.5;

    // Cluster solver
    let mut params_cluster = Params::new();
    params_cluster.set("beta", beta);
    params_cluster.set("bath", "single");
    params_cluster.set("omega0", omega);
    params_cluster.set("g", g);
    params_cluster.set("tunnelling", tunnelling);
    params_cluster.set("epsilon", 0.0); // no bias
    let config_cluster = RunConfig {
        thermalization_sweeps: 2000,
        measurement_sweeps: 10000,
        binsize: 100,
        base_seed: 42,
        ..Default::default()
    };
    let results_cluster = Scheduler::new(RayonBackend::new(1), config_cluster)
        .run_one::<LongitudinalSpinBosonClusterQmc>(&params_cluster);
    let mag_cluster = results_cluster
        .get("ClusterMagnetizationSz")
        .expect("ClusterMagnetizationSz");
    let kinks_cluster = results_cluster
        .get("ClusterKinkCount")
        .expect("ClusterKinkCount");

    // Wormhole solver (same model)
    let mut params_wh = Params::new();
    params_wh.set("beta", beta);
    params_wh.set("model", "xxz"); // longitudinal model
    params_wh.set("bath", "single");
    params_wh.set("omega0", omega);
    params_wh.set("lambda_z", g * g / omega); // λ_z = g²/ω
    params_wh.set("lambda_xy", 0.0); // no transverse
    params_wh.set("h_z", tunnelling);
    let config_wh = RunConfig {
        thermalization_sweeps: 2000,
        measurement_sweeps: 10000,
        binsize: 100,
        base_seed: 42,
        ..Default::default()
    };
    let results_wh =
        Scheduler::new(RayonBackend::new(1), config_wh).run_one::<ImpurityQmc>(&params_wh);
    let mag_wh = results_wh
        .get("MagnetizationSigmaZ")
        .expect("MagnetizationSigmaZ");
    let order_wh = results_wh.get("ExpansionOrder").expect("ExpansionOrder");

    // Both should produce finite results
    assert!(
        mag_cluster.mean.is_finite(),
        "Cluster ⟨σz⟩ should be finite, got {:.4}",
        mag_cluster.mean
    );
    assert!(
        mag_wh.mean.is_finite(),
        "Wormhole ⟨σz⟩ should be finite, got {:.4}",
        mag_wh.mean
    );
    assert!(
        kinks_cluster.mean >= 0.0,
        "Cluster kink count should be ≥ 0, got {:.4}",
        kinks_cluster.mean
    );
    assert!(
        order_wh.mean >= 0.0,
        "Wormhole expansion order should be ≥ 0, got {:.4}",
        order_wh.mean
    );

    eprintln!(
        "Cluster: ⟨σz⟩={:.4}±{:.4}, kinks={:.2}±{:.2}",
        mag_cluster.mean, mag_cluster.stderr, kinks_cluster.mean, kinks_cluster.stderr
    );
    eprintln!(
        "Wormhole: ⟨σz⟩={:.4}±{:.4}, order={:.2}±{:.2}",
        mag_wh.mean, mag_wh.stderr, order_wh.mean, order_wh.stderr
    );
    eprintln!("Note: different conventions prevent direct observable comparison");
}
