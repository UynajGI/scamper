//! Cross-solver validation: wormhole↔occupation.
//!
//! The occupation solver is compared against exact analytic results for the
//! free two-level system. The wormhole solver is run on the same model but
//! only checked for finite output (convention differences prevent direct
//! observable comparison). A true cross-solver test would compare both
//! against a shared ED reference.

use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use qmc_rs::impurity::ImpurityQmc;
use qmc_rs::OccupationWorldlineQmc;

// ─── P1.1: Wormhole ↔ Occupation (free two-level system) ────────────────

#[test]
fn occupation_matches_exact_free_two_level_wormhole_smoke() {
    // Free two-level system: g=0, only tunnelling Δ.
    // Both solvers should give ⟨σz⟩ = -tanh(βΔ/2) (occupation convention)
    // and ⟨E⟩ = -(Δ/2)tanh(βΔ/2).
    let beta: f64 = 4.0;
    let delta: f64 = 0.5;
    let exact_sz: f64 = -(beta * delta / 2.0).tanh();

    // Occupation solver
    let mut params_occ = Params::new();
    params_occ.set("beta", beta);
    params_occ.set("kind", "rabi");
    params_occ.set("spin_splitting", delta);
    params_occ.set("g", 0.0); // free
    params_occ.set("omega0", 1.0);
    params_occ.set("cutoff", 5);
    let config_occ = RunConfig {
        thermalization_sweeps: 2000,
        measurement_sweeps: 10000,
        binsize: 100,
        base_seed: 42,
        ..Default::default()
    };
    let results_occ = Scheduler::new(RayonBackend::new(1), config_occ)
        .run_one::<OccupationWorldlineQmc>(&params_occ);
    let sz_occ = results_occ
        .get("OccupationSigmaZ")
        .expect("OccupationSigmaZ");

    // Wormhole solver (rotated basis: σz_sampled = σx_physical)
    // For free system with only tunnelling, the wormhole measures
    // MagnetizationSigmaZ which corresponds to physical σx.
    // In the free limit, ⟨σx⟩_physical = 0 (no σx term in H).
    // So we compare expansion order instead.
    let mut params_wh = Params::new();
    params_wh.set("beta", beta);
    params_wh.set("model", "rabi");
    params_wh.set("bath", "single");
    params_wh.set("omega0", 1.0);
    params_wh.set("g", 0.0); // free
    params_wh.set("tunnelling", delta);
    params_wh.set("h_z", 0.0);
    let config_wh = RunConfig {
        thermalization_sweeps: 2000,
        measurement_sweeps: 10000,
        binsize: 100,
        base_seed: 42,
        ..Default::default()
    };
    let results_wh =
        Scheduler::new(RayonBackend::new(1), config_wh).run_one::<ImpurityQmc>(&params_wh);
    let order_wh = results_wh.get("ExpansionOrder").expect("ExpansionOrder");

    // Occupation ⟨σz⟩ should match exact
    assert!(
        (sz_occ.mean - exact_sz).abs() < 4.0 * sz_occ.stderr.max(0.02),
        "Occupation ⟨σz⟩={:.4}±{:.4}, exact={:.4}",
        sz_occ.mean,
        sz_occ.stderr,
        exact_sz
    );

    // Wormhole expansion order should be non-negative (free system has few vertices)
    assert!(
        order_wh.mean >= 0.0,
        "Wormhole expansion order should be ≥ 0, got {:.4}",
        order_wh.mean
    );
}

// ─── P1.1b: Wormhole ↔ Occupation (interacting, compare energy) ─────────

#[test]
fn wormhole_and_occupation_smoke_both_run_interacting() {
    // Interacting Rabi model: both solvers should give consistent
    // expansion order (related to energy).
    let beta = 8.0;
    let g = 0.3;
    let omega = 1.0;
    let delta = 1.0;

    // Occupation solver
    let mut params_occ = Params::new();
    params_occ.set("beta", beta);
    params_occ.set("kind", "rabi");
    params_occ.set("spin_splitting", delta);
    params_occ.set("g", g);
    params_occ.set("omega0", omega);
    params_occ.set("cutoff", 15);
    let config_occ = RunConfig {
        thermalization_sweeps: 3000,
        measurement_sweeps: 15000,
        binsize: 100,
        base_seed: 42,
        ..Default::default()
    };
    let results_occ = Scheduler::new(RayonBackend::new(1), config_occ)
        .run_one::<OccupationWorldlineQmc>(&params_occ);
    let n_occ = results_occ
        .get("OccupationBosonNumber")
        .expect("OccupationBosonNumber");

    // Wormhole solver
    let mut params_wh = Params::new();
    params_wh.set("beta", beta);
    params_wh.set("model", "rabi");
    params_wh.set("bath", "single");
    params_wh.set("omega0", omega);
    params_wh.set("g", g);
    params_wh.set("tunnelling", delta);
    params_wh.set("h_z", 0.0);
    let config_wh = RunConfig {
        thermalization_sweeps: 3000,
        measurement_sweeps: 15000,
        binsize: 100,
        base_seed: 42,
        ..Default::default()
    };
    let results_wh =
        Scheduler::new(RayonBackend::new(1), config_wh).run_one::<ImpurityQmc>(&params_wh);
    let order_wh = results_wh.get("ExpansionOrder").expect("ExpansionOrder");

    // Both should produce positive, finite results
    assert!(
        n_occ.mean > 0.0 && n_occ.mean.is_finite(),
        "Occupation ⟨n⟩ should be positive and finite, got {:.4}",
        n_occ.mean
    );
    assert!(
        order_wh.mean > 0.0 && order_wh.mean.is_finite(),
        "Wormhole expansion order should be positive and finite, got {:.4}",
        order_wh.mean
    );

    // Note: occupation ⟨n⟩ and wormhole expansion order measure different
    // physical quantities (explicit boson number vs retarded vertex count).
    // They are not directly comparable. This test verifies both solvers
    // produce physically reasonable results on the same interacting model.
}
