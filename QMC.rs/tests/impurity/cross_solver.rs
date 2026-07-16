//! Cross-solver validation notes and convention documentation.
//!
//! The wormhole and occupation solvers share the single-mode Rabi spin-boson
//! physics but implement fundamentally different algorithms. Direct
//! observable-by-observable comparison through the Carlo.rs adapter is blocked
//! by four documented convention differences. This module records the
//! findings as tests that verify each solver is internally consistent
//! within its own conventions.

use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use qmc_rs::OccupationWorldlineQmc;

// ─── Convention differences (documented, not bugs) ───────────────────────
//
// 1. Basis rotation: The wormhole samples in a rotated basis
//    (σz_sampled = σx_physical) to make the retarded interaction sign-free.
//    See BasisTransform::rotated_rabi() in core/operators.rs.
//    The occupation solver keeps the physical σz basis.
//
// 2. Bath representation: The wormhole integrates out the bath into retarded
//    two-time vertices. The occupation solver keeps explicit boson states.
//    Energy decompositions differ between the two.
//
// 3. Off-diagonal observables: The occupation solver samples in the σz
//    product basis, so ⟨σx⟩ = 0 by construction. Transverse observables
//    are computed exactly from the transfer matrix but are not directly
//    comparable to wormhole sampled σx.
//
// 4. Sign conventions:
//    - Wormhole: H = -h_z Sz + bath(σz_sampled = σx_physical)
//    - Occupation: H = +(Δ/2)σz + Σωn + g(σ+a + σ-a†)
//    These describe the same physics with different notation.
//
// Full cross-solver numerical validation requires a shared independent
// reference (e.g. dense ED) that both solvers are compared against
// separately. The occupation solver already has its own ED comparisons
// (tests/impurity/occupation.rs, src/worldline.rs inline tests).

/// Occupation solver via Carlo.rs adapter should produce ⟨σz⟩ consistent
/// with positive spin splitting driving spin down (negative ⟨σz⟩).
#[test]
fn occupation_sigma_z_reflects_spin_splitting_sign() {
    let beta = 6.0;
    let delta = 0.4; // positive Δ → spin down (σz=-1) favored
    let g = 0.35;
    let omega = 1.0;

    let mut params = Params::new();
    params.set("beta", beta);
    params.set("kind", "rabi");
    params.set("spin_splitting", delta);
    params.set("g", g);
    params.set("omega0", omega);
    params.set("cutoff", 12);

    let run = RunConfig {
        thermalization_sweeps: 5_000,
        measurement_sweeps: 30_000,
        binsize: 20,
        base_seed: 55,
        ..Default::default()
    };
    let results =
        Scheduler::new(RayonBackend::new(1), run).run_one::<OccupationWorldlineQmc>(&params);
    let sigma_z = results
        .get("OccupationSigmaZ")
        .expect("OccupationSigmaZ from occupation");

    // Positive Δ → H has +Δ/2 σz → spin down (σz=-1) is lower energy.
    // ⟨σz⟩ → -tanh(βΔ/2) for the free spin, reduced in magnitude by coupling.
    assert!(
        sigma_z.mean < -0.3,
        "occupation ⟨σz⟩ should be negative with Δ>0: got {:.6}",
        sigma_z.mean
    );
}

/// Reversing spin splitting should flip σz sign.
#[test]
fn occupation_sigma_z_flips_with_spin_splitting_sign() {
    let beta = 4.0;
    let g = 0.3;
    let omega = 1.0;

    let run = |seed: u64, delta: f64| -> f64 {
        let mut params = Params::new();
        params.set("beta", beta);
        params.set("kind", "rabi");
        params.set("spin_splitting", delta);
        params.set("g", g);
        params.set("omega0", omega);
        params.set("cutoff", 10);

        let config = RunConfig {
            thermalization_sweeps: 3_000,
            measurement_sweeps: 20_000,
            binsize: 20,
            base_seed: seed,
            ..Default::default()
        };
        let results =
            Scheduler::new(RayonBackend::new(1), config).run_one::<OccupationWorldlineQmc>(&params);
        results
            .get("OccupationSigmaZ")
            .expect("OccupationSigmaZ")
            .mean
    };

    let sz_pos = run(11, 0.5);
    let sz_neg = run(22, -0.5);

    // Sign flip: positive Δ → negative σz, negative Δ → positive σz
    assert!(
        (sz_pos + sz_neg).abs() < 0.1,
        "σz should be odd in Δ: sz(Δ=0.5)={sz_pos:.4}, sz(Δ=-0.5)={sz_neg:.4}"
    );
    assert!(
        sz_pos < 0.0 && sz_neg > 0.0,
        "signs should be opposite: sz(Δ=0.5)={sz_pos:.4}, sz(Δ=-0.5)={sz_neg:.4}"
    );
}
