//! Cross-solver validation: wormhole↔occupation.
//!
//! Both solvers are compared against the exact analytic result for the
//! free two-level system (g=0): |⟨σ⟩| = tanh(βΔ/2). The wormhole and
//! occupation solvers use different basis conventions (see
//! `cross_solver.rs` for the full catalogue), so the sign of the
//! measured observable differs, but the magnitude must agree and each
//! must match the exact result individually.

use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use qmc_rs::impurity::ImpurityQmc;
use qmc_rs::OccupationWorldlineQmc;

// ─── P1.1: Wormhole ↔ Occupation (free two-level system) ────────────────

#[test]
fn occupation_and_wormhole_match_exact_free_two_level() {
    // Free two-level system: g=0, only splitting Δ.
    //
    // The two solvers use different basis conventions (see cross_solver.rs):
    //   • Wormhole (rotated basis, σz_sampled = σx_physical):
    //       H = -(Δ/2)σx  →  MagnetizationSigmaZ = +tanh(βΔ/2)
    //   • Occupation (occupation basis):
    //       H = +(Δ/2)σz  →  OccupationSigmaZ = -tanh(βΔ/2)
    //
    // The sign flip is convention difference #4. Both solvers measure the
    // SAME physical splitting Δ, so their magnitudes must agree and each
    // must match the exact result individually.
    let beta: f64 = 4.0;
    let delta: f64 = 1.0;
    let exact_tanh: f64 = (beta * delta / 2.0).tanh(); // tanh(2.0) ≈ 0.964

    // ── Occupation solver ──────────────────────────────────────────
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

    // ── Wormhole solver ───────────────────────────────────────────
    // NOTE: the wormhole uses `tunnelling` (not `spin_splitting`).
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
    let sz_wh = results_wh
        .get("MagnetizationSigmaZ")
        .expect("MagnetizationSigmaZ");

    // ── 1. Occupation ⟨σz⟩ vs exact ──────────────────────────────
    // Occupation convention: ⟨σz⟩ = -tanh(βΔ/2)
    let exact_occ = -exact_tanh;
    assert!(
        (sz_occ.mean - exact_occ).abs() < 4.0 * sz_occ.stderr.max(0.02),
        "Occupation ⟨σz⟩={:.4}±{:.4}, exact={:.4}",
        sz_occ.mean,
        sz_occ.stderr,
        exact_occ
    );

    // ── 2. Wormhole ⟨σz⟩ (= physical ⟨σx⟩) vs exact ──────────────
    // Wormhole convention: MagnetizationSigmaZ = +tanh(βΔ/2)
    let exact_wh = exact_tanh;
    assert!(
        (sz_wh.mean - exact_wh).abs() < 4.0 * sz_wh.stderr.max(0.02),
        "Wormhole ⟨σz⟩={:.4}±{:.4}, exact={:.4}",
        sz_wh.mean,
        sz_wh.stderr,
        exact_wh
    );

    // ── 3. Cross-solver: occupation vs wormhole ──────────────────
    // Sign-corrected comparison: -⟨σz⟩_occ should agree with ⟨σz⟩_wh
    // within combined 4σ (sqrt of sum of squared stderrs).
    let combined_sigma = (sz_occ.stderr.powi(2) + sz_wh.stderr.powi(2))
        .sqrt()
        .max(0.02);
    assert!(
        ((-sz_occ.mean) - sz_wh.mean).abs() < 4.0 * combined_sigma,
        "Cross-solver mismatch: -occ⟨σz⟩={:.4}±{:.4}, wh⟨σz⟩={:.4}±{:.4}, combined 4σ={:.4}",
        -sz_occ.mean,
        sz_occ.stderr,
        sz_wh.mean,
        sz_wh.stderr,
        4.0 * combined_sigma
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
