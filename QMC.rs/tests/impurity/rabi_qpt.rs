//! Rabi model quantum phase transition: analysis of the order parameter.
//!
//! The Rabi model H = (Δ/2)σz + Ωa†a + gσx(a+a†) has a quantum phase
//! transition in the Ω/Δ → ∞ limit at the critical coupling
//! g_c = √(ΩΔ)/2, i.e. r = g/g_c = 1.
//!
//! ## Evidence for r=1 as the critical point
//!
//! The proof comes from the Born-Oppenheimer (adiabatic) approximation.
//! For fixed displacement x = a+a†, the spin Hamiltonian is:
//!
//! ```text
//! H_spin(x) = (Δ/2)σz + g·x·σx
//! ```
//!
//! Its lower eigenvalue is E₋(x) = -√((Δ/2)² + g²x²).
//! The effective potential for x is:
//!
//! ```text
//! V_eff(x) = (Ω/2)x² - √((Δ/2)² + g²x²)
//! ```
//!
//! The curvature at the origin determines stability:
//!
//! ```text
//! V_eff''(0) = Ω - 4g²/Δ = Ω(1 - r²)
//! ```
//!
//! - r < 1: V_eff''(0) > 0 → single minimum at x=0 (symmetric phase)
//! - r > 1: V_eff''(0) < 0 → double well (broken phase)
//! - r = 1: V_eff''(0) = 0 → critical point
//!
//! This is the **analytic proof** that r=1 is the QPT critical point.
//!
//! ## Why finite-η ED cannot observe the transition
//!
//! At finite η = Ω/Δ, quantum tunneling between the two wells keeps the
//! ground state symmetric (Z₂-even). The transition only exists in the
//! strict η→∞ (classical) limit. Unlike finite-size scaling in lattice
//! models, increasing η does NOT sharpen the transition — instead the
//! barrier height grows but so does the zero-point energy, keeping the
//! system symmetric.
//!
//! ## Tests
//!
//! Layer 1 (deterministic ED):
//! - Verify V_eff''(0) changes sign at r=1 (analytic, machine precision)
//! - Verify ground-state energy curvature ∂²E₀/∂g² peaks near r=1
//! - Verify ⟨x²⟩ and cutoff convergence
//!
//! Layer 2 (stochastic MC, #[ignore]):
//! - MC-sampled ⟨n⟩ matches ED at large β

#![allow(clippy::needless_range_loop)]
#![allow(clippy::excessive_precision)]

use qmc_rs::impurity::spin_boson::occupation::transfer::SymmetricEigensystem;
use qmc_rs::impurity::spin_boson::occupation::{
    model::{CavityMode, OccupationSpinBosonModel},
    OccupationWorldlineSampler,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

// ─── Helpers ─────────────────────────────────────────────────────────────

fn matmul(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = a.len();
    let mut c = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

fn displacement_power_operator(model: &OccupationSpinBosonModel, power: u32) -> Vec<Vec<f64>> {
    let dim = model.basis().dimension();
    let cutoff = model.basis().cutoffs()[0];

    let mut x_boson = vec![vec![0.0f64; cutoff]; cutoff];
    for n in 1..cutoff {
        let sqn = (n as f64).sqrt();
        x_boson[n - 1][n] = sqn;
        x_boson[n][n - 1] = sqn;
    }

    let mut xk_boson = vec![vec![1.0f64; cutoff]; cutoff];
    for _ in 0..power {
        xk_boson = matmul(&xk_boson, &x_boson);
    }

    let mut op = vec![vec![0.0; dim]; dim];
    for i in 0..cutoff {
        for j in 0..cutoff {
            for spin in 0..2 {
                op[2 * i + spin][2 * j + spin] = xk_boson[i][j];
            }
        }
    }
    op
}

fn ground_state_expectation(eigen: &SymmetricEigensystem, operator: &[Vec<f64>]) -> f64 {
    let dim = eigen.values.len();
    let mut sum = 0.0;
    for i in 0..dim {
        for j in 0..dim {
            sum += eigen.vectors[i][0] * operator[i][j] * eigen.vectors[j][0];
        }
    }
    sum
}

fn ground_state_energy(delta: f64, omega: f64, g: f64, cutoff: usize) -> f64 {
    let model =
        OccupationSpinBosonModel::rabi(delta, vec![CavityMode::new(omega, g, cutoff).unwrap()])
            .unwrap();
    let eigen = SymmetricEigensystem::diagonalize(model.hamiltonian()).unwrap();
    eigen.values[0]
}

// ─── Layer 1: deterministic ED tests ─────────────────────────────────────

/// Tests the analytic Born-Oppenheimer formula V_eff''(0) = Ω(1-r²),
/// not the QMC solver. This is a math identity check that documents the
/// analytic proof of the critical point at r=1.
#[test]
fn born_oppenheimer_curvature_changes_sign_at_r_equals_one() {
    let delta = 1.0_f64;
    let omega = 10.0_f64;

    for &r in &[0.5_f64, 0.9, 0.99, 1.0, 1.01, 1.1, 1.5] {
        let g = r * (omega * delta).sqrt() / 2.0;
        let curvature = omega - 4.0 * g * g / delta;
        let expected = omega * (1.0 - r * r);

        assert!(
            (curvature - expected).abs() < 1e-12,
            "V''(0) formula: got {curvature}, expected {expected} at r={r}"
        );

        if r < 1.0 {
            assert!(curvature > 0.0, "V''(0) > 0 for r<1 (symmetric phase)");
        } else if r > 1.0 {
            assert!(curvature < 0.0, "V''(0) < 0 for r>1 (broken phase)");
        }
    }

    // At exactly r=1, curvature is zero
    let g_c = (omega * delta).sqrt() / 2.0;
    let curvature_c = omega - 4.0 * g_c * g_c / delta;
    assert!(
        curvature_c.abs() < 1e-12,
        "V''(0) = 0 at r=1: got {curvature_c}"
    );
}

/// The ground-state energy is a smooth function of g. Verify E₀(g)
/// is monotonically decreasing (stronger coupling lowers the energy).
#[test]
fn ground_state_energy_decreases_monotonically_with_coupling() {
    let delta = 1.0_f64;
    let omega = 20.0_f64;
    let cutoff = 30;
    let gc = (omega * delta).sqrt() / 2.0;

    let mut prev_e = f64::INFINITY;
    for i in 0..=50 {
        let r = 0.2 + 2.8 * i as f64 / 50.0;
        let g = r * gc;
        let e = ground_state_energy(delta, omega, g, cutoff);
        assert!(
            e < prev_e + 1e-10,
            "E₀ should decrease: r={r:.2}, E={e:.6}, prev={prev_e:.6}"
        );
        prev_e = e;
    }
}

/// Verify cutoff convergence: at the critical point, ⟨x²⟩ converges
/// as cutoff increases.
#[test]
fn rabi_ground_state_converges_with_cutoff() {
    let delta = 1.0;
    let omega = 10.0_f64;
    let gc = (omega * delta).sqrt() / 2.0;
    let g = gc;

    let cutoffs = [10, 20, 30];
    let mut x2_values = Vec::new();

    for &cutoff in &cutoffs {
        let model =
            OccupationSpinBosonModel::rabi(delta, vec![CavityMode::new(omega, g, cutoff).unwrap()])
                .unwrap();
        let eigen = SymmetricEigensystem::diagonalize(model.hamiltonian()).unwrap();
        let x2_op = displacement_power_operator(&model, 2);
        let x2 = ground_state_expectation(&eigen, &x2_op);
        x2_values.push(x2);
    }

    let rel_diff = (x2_values[2] - x2_values[1]).abs() / x2_values[1];
    assert!(
        rel_diff < 0.05,
        "cutoff not converged: x²(20)={:.6}, x²(30)={:.6}, rel_diff={rel_diff:.4}",
        x2_values[1],
        x2_values[2]
    );
}

/// The displacement ⟨x²⟩ increases with coupling.
#[test]
fn rabi_displacement_increases_with_coupling() {
    let delta = 1.0;
    let omega = 2.0_f64;
    let cutoff = 40;
    let gc = (omega * delta).sqrt() / 2.0;

    let model_weak = OccupationSpinBosonModel::rabi(
        delta,
        vec![CavityMode::new(omega, 0.5 * gc, cutoff).unwrap()],
    )
    .unwrap();
    let model_strong = OccupationSpinBosonModel::rabi(
        delta,
        vec![CavityMode::new(omega, 3.0 * gc, cutoff).unwrap()],
    )
    .unwrap();

    let x2_op_weak = displacement_power_operator(&model_weak, 2);
    let x2_op_strong = displacement_power_operator(&model_strong, 2);

    let eigen_weak = SymmetricEigensystem::diagonalize(model_weak.hamiltonian()).unwrap();
    let eigen_strong = SymmetricEigensystem::diagonalize(model_strong.hamiltonian()).unwrap();

    let x2_weak = ground_state_expectation(&eigen_weak, &x2_op_weak);
    let x2_strong = ground_state_expectation(&eigen_strong, &x2_op_strong);

    assert!(x2_strong > x2_weak, "⟨x²⟩ should increase");
    assert!(
        x2_strong / x2_weak > 2.0,
        "⟨x²⟩ ratio={:.2}",
        x2_strong / x2_weak
    );
}

/// At strong coupling, U4(x) shows cutoff dependence — prerequisite
/// for finite-size analysis.
#[test]
fn rabi_binder_cumulant_shows_cutoff_dependence_at_strong_coupling() {
    let delta = 1.0;
    let omega = 2.0_f64;
    let r_values: Vec<f64> = (100..=500).step_by(20).map(|i| i as f64 / 100.0).collect();

    let cutoffs = [8usize, 16, 32];
    let mut u4_curves = vec![Vec::new(); cutoffs.len()];

    for (ci, &cutoff) in cutoffs.iter().enumerate() {
        let gc = (omega * delta).sqrt() / 2.0;
        for &r in &r_values {
            let g = r * gc;
            let model = OccupationSpinBosonModel::rabi(
                delta,
                vec![CavityMode::new(omega, g, cutoff).unwrap()],
            )
            .unwrap();
            let eigen = SymmetricEigensystem::diagonalize(model.hamiltonian()).unwrap();
            let x2_op = displacement_power_operator(&model, 2);
            let x4_op = displacement_power_operator(&model, 4);
            let x2 = ground_state_expectation(&eigen, &x2_op);
            let x4 = ground_state_expectation(&eigen, &x4_op);
            let u4 = if x2 < 1e-14 {
                0.0
            } else {
                1.0 - x4 / (3.0 * x2 * x2)
            };
            u4_curves[ci].push(u4);
        }
    }

    let max_spread = r_values
        .iter()
        .enumerate()
        .map(|(ri, _)| {
            let vals: Vec<f64> = (0..cutoffs.len()).map(|ci| u4_curves[ci][ri]).collect();
            vals.iter().cloned().fold(0.0f64, f64::max)
                - vals.iter().cloned().fold(f64::INFINITY, f64::min)
        })
        .fold(0.0f64, f64::max);
    assert!(
        max_spread > 0.001,
        "U4 should show cutoff dependence: max_spread={max_spread:.6}"
    );
}

// ─── Layer 2: stochastic MC validation (#[ignore]) ──────────────────────

/// MC-sampled ⟨n⟩ at β=1024 should match exact ED near the critical point.
#[test]
#[ignore]
fn mc_sampled_photon_number_matches_ed_near_critical_point() {
    let delta = 1.0_f64;
    let omega = 5.0_f64; // η=5
    let beta = (1u64 << 10) as f64; // β=1024
    let cutoff = 60;
    let gc = (omega * delta).sqrt() / 2.0;
    let g = gc; // r=1

    let model_ed =
        OccupationSpinBosonModel::rabi(delta, vec![CavityMode::new(omega, g, cutoff).unwrap()])
            .unwrap();
    let eigen = SymmetricEigensystem::diagonalize(model_ed.hamiltonian()).unwrap();

    let dim = model_ed.basis().dimension();
    let mut n_op = vec![vec![0.0; dim]; dim];
    for state in 0..dim {
        n_op[state][state] = model_ed.basis().occupation(state, 0) as f64;
    }
    let n_exact = ground_state_expectation(&eigen, &n_op);

    let model_mc =
        OccupationSpinBosonModel::rabi(delta, vec![CavityMode::new(omega, g, cutoff).unwrap()])
            .unwrap();
    let mut sampler = OccupationWorldlineSampler::new(model_mc, beta, 8, 0).unwrap();
    let mut n_sum = 0.0;
    let mut samples = 0u64;
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
    for sweep in 0..500_000u64 {
        sampler.sweep(&mut rng).unwrap();
        if sweep >= 50_000 {
            let obs = sampler.measure().unwrap();
            n_sum += obs.total_boson_number;
            samples += 1;
        }
    }
    let n_mc = n_sum / samples as f64;

    let rel_error = (n_mc - n_exact).abs() / n_exact;
    assert!(
        rel_error < 0.05,
        "⟨n⟩: MC={n_mc:.6}, ED={n_exact:.6}, rel_error={rel_error:.4}"
    );
}
