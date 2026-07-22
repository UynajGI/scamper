//! Wormhole solver validation against exact diagonalization for the
//! **interacting** single-mode Rabi model.
//!
//! The wormhole samples in a rotated basis where `σz_sampled = σx_physical`
//! (see `BasisTransform::rotated_rabi`).  Its `MagnetizationSigmaZ` observable
//! therefore measures the physical `⟨σx⟩`.
//!
//! The physical (spin-boson) Hamiltonian in the wormhole's convention is:
//!
//!   H = Ω a†a − (Δ/2) σx + (g/2) σz (a + a†)
//!
//! where Δ = `tunnelling` and g = `coupling`.  This differs from the
//! quantum-optics Rabi model by a π/2 spin rotation; the ED helper below
//! uses the wormhole's convention so that `⟨σx⟩_ED` can be compared
//! directly to `MagnetizationSigmaZ`.
//!
//! Stochastic tests are `#[ignore]` — they take ~30 s each.

#![allow(clippy::needless_range_loop)]

use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use qmc_rs::impurity::spin_boson::occupation::transfer::SymmetricEigensystem;
use qmc_rs::impurity::ImpurityQmc;

// ── Exact diagonalization helper ──────────────────────────────────────────

/// Basis index: `|n, spin⟩ → 2n + spin`, spin 0 = ↓ (σz = −1), 1 = ↑ (σz = +1).
#[inline]
fn idx(n: usize, spin: usize) -> usize {
    2 * n + spin
}

/// Build the single-mode Rabi Hamiltonian in the wormhole's **physical**
/// (spin-boson) basis:
///
///   H = Ω a†a − (Δ/2) σx + (g/2) σz (a + a†)
///
/// This matches the `rabi_matrix(…, rotated=false)` convention used in the
/// inline tests of `observables.rs`.
fn rabi_hamiltonian(
    boson_states: usize,
    omega: f64,
    tunnelling: f64,
    coupling: f64,
) -> Vec<Vec<f64>> {
    let dim = 2 * boson_states;
    let mut h = vec![vec![0.0; dim]; dim];
    for n in 0..boson_states {
        for spin in 0..2 {
            let i = idx(n, spin);
            let sz = if spin == 0 { -0.5 } else { 0.5 };
            // Boson energy: Ω n
            h[i][i] += omega * n as f64;
            // Tunnelling: −(Δ/2) σx → off-diagonal in spin, same n
            let j = idx(n, 1 - spin);
            h[i][j] += -0.5 * tunnelling;
            // Coupling: (g/2) σz (a + a†) → diagonal in spin, changes n by ±1
            if n + 1 < boson_states {
                let amplitude = (n as f64 + 1.0).sqrt();
                let k = idx(n + 1, spin);
                h[i][k] += coupling * sz * amplitude;
                h[k][i] += coupling * sz * amplitude;
            }
        }
    }
    h
}

/// Build the σx Pauli operator (full Pauli matrix, eigenvalues ±1).
fn sigma_x_operator(boson_states: usize) -> Vec<Vec<f64>> {
    let dim = 2 * boson_states;
    let mut op = vec![vec![0.0; dim]; dim];
    for n in 0..boson_states {
        op[idx(n, 0)][idx(n, 1)] = 1.0;
        op[idx(n, 1)][idx(n, 0)] = 1.0;
    }
    op
}

/// Thermal expectation `⟨O⟩ = Σ_k ⟨k|O|k⟩ e^{−βE_k} / Z`.
fn thermal_expectation(eigen: &SymmetricEigensystem, operator: &[Vec<f64>], beta: f64) -> f64 {
    let dim = eigen.values.len();
    let ground = eigen.values[0];
    let mut numerator = 0.0;
    let mut partition = 0.0;
    for k in 0..dim {
        let exponent = -beta * (eigen.values[k] - ground);
        if exponent < -700.0 {
            continue;
        }
        let weight = exponent.exp();
        let mut expectation = 0.0;
        for i in 0..dim {
            for j in 0..dim {
                expectation += eigen.vectors[i][k] * operator[i][j] * eigen.vectors[j][k];
            }
        }
        numerator += weight * expectation;
        partition += weight;
    }
    numerator / partition
}

/// Thermal energy `⟨E⟩ = Σ_k E_k e^{−βE_k} / Z`.
fn thermal_energy(eigen: &SymmetricEigensystem, beta: f64) -> f64 {
    let ground = eigen.values[0];
    let mut numerator = 0.0;
    let mut partition = 0.0;
    for &energy in &eigen.values {
        let exponent = -beta * (energy - ground);
        if exponent < -700.0 {
            continue;
        }
        let weight = exponent.exp();
        numerator += weight * energy;
        partition += weight;
    }
    numerator / partition
}

/// Imaginary-time correlation at τ = β/2 via the Lehmann representation:
///
///   C(β/2) = (1/Z) Σ_{k,l} |⟨k|O|l⟩|² e^{−β(E_k + E_l)/2}
fn correlation_half(eigen: &SymmetricEigensystem, operator: &[Vec<f64>], beta: f64) -> f64 {
    let dim = eigen.values.len();
    let ground = eigen.values[0];
    let half_beta = 0.5 * beta;

    // Boltzmann half-weights: w_k = e^{−(β/2)(E_k − E_0)}
    let weights: Vec<f64> = eigen
        .values
        .iter()
        .map(|&e| {
            let exponent = -half_beta * (e - ground);
            if exponent < -700.0 {
                0.0
            } else {
                exponent.exp()
            }
        })
        .collect();
    // Z_shifted = Σ_k e^{−β(E_k − E_0)} = Σ_k w_k²
    let z: f64 = weights.iter().map(|w| w * w).sum();

    let mut result = 0.0;
    for k in 0..dim {
        if weights[k] == 0.0 {
            continue;
        }
        for l in 0..dim {
            if weights[l] == 0.0 {
                continue;
            }
            // ⟨k|O|l⟩ in the energy eigenbasis
            let mut matrix_element = 0.0;
            for i in 0..dim {
                for j in 0..dim {
                    matrix_element += eigen.vectors[i][k] * operator[i][j] * eigen.vectors[j][l];
                }
            }
            result += matrix_element * matrix_element * weights[k] * weights[l];
        }
    }
    result / z
}

/// Compute ED reference values for the interacting Rabi model.
///
/// Returns `(⟨σx⟩, ⟨E⟩, C(β/2))` where σx is the physical Pauli operator
/// that the wormhole's `MagnetizationSigmaZ` measures.
fn ed_rabi_observables(
    omega: f64,
    coupling: f64,
    tunnelling: f64,
    beta: f64,
    cutoff: usize,
) -> (f64, f64, f64) {
    let h = rabi_hamiltonian(cutoff, omega, tunnelling, coupling);
    let eigen = SymmetricEigensystem::diagonalize(h).expect("ED diagonalization failed");
    let sx = sigma_x_operator(cutoff);

    let sigma_x = thermal_expectation(&eigen, &sx, beta);
    let energy = thermal_energy(&eigen, beta);
    let corr_half = correlation_half(&eigen, &sx, beta);

    (sigma_x, energy, corr_half)
}

// ── Wormhole runner ───────────────────────────────────────────────────────

fn run_wormhole(beta: f64, omega0: f64, g: f64, tunnelling: f64, seed: u64) -> carlo_rs::Results {
    let mut params = Params::new();
    params.set("beta", beta);
    params.set("model", "rabi");
    params.set("bath", "single");
    params.set("omega0", omega0);
    params.set("g", g);
    params.set("tunnelling", tunnelling);
    params.set("h_z", 0.0);
    params.set("validate_each_sweep", true);
    let run = RunConfig {
        thermalization_sweeps: 5_000,
        measurement_sweeps: 20_000,
        binsize: 100,
        base_seed: seed,
        ..Default::default()
    };
    Scheduler::new(RayonBackend::new(1), run).run_one::<ImpurityQmc>(&params)
}

// ── Stochastic tests ──────────────────────────────────────────────────────

/// Compare wormhole MC against ED for the interacting Rabi model at
/// moderate coupling (g = 0.3, λ = g²/ω = 0.09).
///
/// Checks across three independent seeds:
/// - `MagnetizationSigmaZ` (= physical ⟨σx⟩) matches ED within 4σ or 0.05
/// - `ExpansionOrder` is positive, finite, and in a reasonable range
/// - `CorrelationSigmaZHalf` (= physical C(β/2)) matches ED within 4σ or 0.05
/// - Cross-seed consistency of the magnetization
#[test]
#[ignore = "long: wormhole interacting ED comparison"]
fn wormhole_interacting_matches_ed() {
    let omega = 1.0;
    let g = 0.3;
    let tunnelling = 1.0;
    let beta = 10.0;
    let cutoff = 20;

    let (ed_sigma_x, ed_energy, ed_corr_half) =
        ed_rabi_observables(omega, g, tunnelling, beta, cutoff);

    eprintln!("ED reference (cutoff={cutoff}, β={beta}):");
    eprintln!("  ⟨σx⟩     = {ed_sigma_x:.6}");
    eprintln!("  ⟨E⟩      = {ed_energy:.6}");
    eprintln!("  C(β/2)   = {ed_corr_half:.6}");
    eprintln!("  −β·⟨E⟩   = {:.6}", -beta * ed_energy);

    let seeds = [42u64, 137, 2026];
    let mut mag_values = Vec::new();
    let mut mag_errors = Vec::new();

    for &seed in &seeds {
        let results = run_wormhole(beta, omega, g, tunnelling, seed);

        let mag = results
            .get("MagnetizationSigmaZ")
            .expect("MagnetizationSigmaZ missing");
        let order = results
            .get("ExpansionOrder")
            .expect("ExpansionOrder missing");
        let corr = results
            .get("CorrelationSigmaZHalf")
            .expect("CorrelationSigmaZHalf missing");

        eprintln!(
            "\nSeed {seed}: ⟨σx⟩_MC = {:.4} ± {:.4},  ⟨k⟩ = {:.2} ± {:.2},  C(β/2) = {:.4} ± {:.4}",
            mag.mean, mag.stderr, order.mean, order.stderr, corr.mean, corr.stderr,
        );

        // ── Magnetization: MC ⟨σz_sampled⟩ = ED ⟨σx_physical⟩ ──
        let mag_tol = (4.0 * mag.stderr).max(0.05);
        assert!(
            (mag.mean - ed_sigma_x).abs() < mag_tol,
            "seed {seed}: MagnetizationSigmaZ {:.4} ± {:.4} vs ED ⟨σx⟩ {:.4} \
             (deviation {:.4}, tolerance {:.4})",
            mag.mean,
            mag.stderr,
            ed_sigma_x,
            (mag.mean - ed_sigma_x).abs(),
            mag_tol,
        );

        // ── Expansion order sanity ──
        // The expansion order counts retarded-interaction vertices and is
        // related to the retarded interaction energy, NOT the total energy
        // (which includes boson kinetic + tunnelling contributions).
        // We verify it is positive, finite, and in a physically reasonable
        // range for these parameters (λ = g²/ω = 0.09, β = 10).
        assert!(
            order.mean > 0.0 && order.mean.is_finite(),
            "seed {seed}: ExpansionOrder {:.2} should be positive and finite",
            order.mean,
        );
        assert!(
            order.mean > 1.0 && order.mean < 100.0,
            "seed {seed}: ExpansionOrder {:.2} outside reasonable range [1, 100]",
            order.mean,
        );

        // ── Correlation at β/2 ──
        let corr_tol = (4.0 * corr.stderr).max(0.05);
        assert!(
            (corr.mean - ed_corr_half).abs() < corr_tol,
            "seed {seed}: CorrelationSigmaZHalf {:.4} ± {:.4} vs ED C(β/2) {:.4} \
             (deviation {:.4}, tolerance {:.4})",
            corr.mean,
            corr.stderr,
            ed_corr_half,
            (corr.mean - ed_corr_half).abs(),
            corr_tol,
        );

        mag_values.push(mag.mean);
        mag_errors.push(mag.stderr);
    }

    // ── Cross-seed consistency ──
    for i in 0..seeds.len() {
        for j in (i + 1)..seeds.len() {
            let combined_err =
                (mag_errors[i] * mag_errors[i] + mag_errors[j] * mag_errors[j]).sqrt();
            let deviation = (mag_values[i] - mag_values[j]).abs();
            assert!(
                deviation < 4.0 * combined_err,
                "seeds {} and {} disagree: {:.4} vs {:.4} (combined σ = {:.4})",
                seeds[i],
                seeds[j],
                mag_values[i],
                mag_values[j],
                combined_err,
            );
        }
    }
}

/// Run three independent seeds and verify that the z-score of
/// `MagnetizationSigmaZ` relative to the ED reference stays below 4.
#[test]
#[ignore = "long: wormhole interacting ED comparison"]
fn wormhole_interacting_zscore_3_seeds() {
    let omega = 1.0;
    let g = 0.3;
    let tunnelling = 1.0;
    let beta = 10.0;
    let cutoff = 20;

    let (ed_sigma_x, _, _) = ed_rabi_observables(omega, g, tunnelling, beta, cutoff);
    eprintln!("ED ⟨σx⟩ = {ed_sigma_x:.6}");

    let seeds = [7u64, 99, 314];
    for &seed in &seeds {
        let results = run_wormhole(beta, omega, g, tunnelling, seed);
        let mag = results
            .get("MagnetizationSigmaZ")
            .expect("MagnetizationSigmaZ missing");

        let z = (mag.mean - ed_sigma_x) / mag.stderr;
        eprintln!(
            "seed {seed}: ⟨σx⟩_MC = {:.4} ± {:.4},  z = {z:.2}",
            mag.mean, mag.stderr,
        );
        assert!(
            z.abs() < 4.0,
            "seed {seed}: |z| = {:.2} ≥ 4  (MC = {:.4} ± {:.4}, ED = {:.4})",
            z.abs(),
            mag.mean,
            mag.stderr,
            ed_sigma_x,
        );
    }
}

// ── Deterministic ED sanity checks (fast, not ignored) ────────────────────

/// At g = 0 the boson decouples and H = Ωa†a − (Δ/2)σx.
/// Then ⟨σx⟩ = tanh(βΔ/2) and C(τ) = 1 for all τ (σx commutes with H).
#[test]
fn ed_free_spin_limit() {
    let beta = 5.0;
    let tunnelling = 1.0;
    let cutoff = 20; // large enough that boson truncation error is negligible

    let (sigma_x, energy, corr_half) = ed_rabi_observables(1.0, 0.0, tunnelling, beta, cutoff);

    let expected_mag = (0.5 * beta * tunnelling).tanh();
    assert!(
        (sigma_x - expected_mag).abs() < 1e-10,
        "free spin ⟨σx⟩ = {sigma_x:.10}, expected {expected_mag:.10}",
    );

    // Total energy = boson thermal energy + spin energy.
    // ⟨Ωa†a⟩ = Ω / (e^{βΩ} − 1) for a free boson at temperature 1/β.
    let omega = 1.0;
    let boson_energy = omega / ((beta * omega).exp() - 1.0);
    let expected_energy = boson_energy - 0.5 * tunnelling * expected_mag;
    assert!(
        (energy - expected_energy).abs() < 1e-10,
        "free spin ⟨E⟩ = {energy:.10}, expected {expected_energy:.10}",
    );

    // σx commutes with H = −(Δ/2)σx, so σx(τ) = σx(0) and C(τ) = ⟨σx²⟩ = 1.
    assert!(
        (corr_half - 1.0).abs() < 1e-10,
        "free spin C(β/2) = {corr_half:.10}, expected 1.0",
    );
}

/// The physical and rotated Hamiltonians must have identical spectra
/// (unitary equivalence under the wormhole's basis rotation).
#[test]
fn ed_spectrum_invariant_under_rotation() {
    let cutoff = 10;
    let omega = 1.0;
    let tunnelling = 1.0;
    let coupling = 0.3;

    let h_physical = rabi_hamiltonian(cutoff, omega, tunnelling, coupling);

    // Rotated (sampled) Hamiltonian: tunnelling becomes diagonal (σz),
    // coupling becomes off-diagonal (σx, flips spin).
    let dim = 2 * cutoff;
    let mut h_rotated = vec![vec![0.0; dim]; dim];
    for n in 0..cutoff {
        for spin in 0..2 {
            let i = idx(n, spin);
            let sz = if spin == 0 { -0.5 } else { 0.5 };
            h_rotated[i][i] += omega * n as f64;
            // Tunnelling: −(Δ/2)σz → diagonal
            h_rotated[i][i] += -tunnelling * sz;
            // Coupling: −(g/2)σx(a+a†) → off-diagonal, flips spin
            if n + 1 < cutoff {
                let amplitude = (n as f64 + 1.0).sqrt();
                let j = idx(n + 1, 1 - spin);
                h_rotated[i][j] += -0.5 * coupling * amplitude;
                h_rotated[j][i] += -0.5 * coupling * amplitude;
            }
        }
    }

    let eigen_phys = SymmetricEigensystem::diagonalize(h_physical).unwrap();
    let eigen_rot = SymmetricEigensystem::diagonalize(h_rotated).unwrap();

    for (k, (&e_phys, &e_rot)) in eigen_phys.values.iter().zip(&eigen_rot.values).enumerate() {
        assert!(
            (e_phys - e_rot).abs() < 1e-10,
            "eigenvalue {k}: physical {e_phys:.12} vs rotated {e_rot:.12}",
        );
    }
}

/// ED ⟨σx⟩ in the physical basis must equal ED ⟨σz⟩ in the rotated basis
/// (the observable mapping used by the wormhole).
#[test]
fn ed_observable_mapping_under_rotation() {
    let cutoff = 8;
    let omega = 1.3;
    let tunnelling = 0.7;
    let coupling = 0.45;
    let beta = 2.1;

    // Physical basis: measure σx
    let h_phys = rabi_hamiltonian(cutoff, omega, tunnelling, coupling);
    let eigen_phys = SymmetricEigensystem::diagonalize(h_phys).unwrap();
    let sx = sigma_x_operator(cutoff);
    let phys_sx = thermal_expectation(&eigen_phys, &sx, beta);

    // Rotated basis: measure σz
    let dim = 2 * cutoff;
    let mut h_rot = vec![vec![0.0; dim]; dim];
    for n in 0..cutoff {
        for spin in 0..2 {
            let i = idx(n, spin);
            let sz = if spin == 0 { -0.5 } else { 0.5 };
            h_rot[i][i] += omega * n as f64 - tunnelling * sz;
            if n + 1 < cutoff {
                let amplitude = (n as f64 + 1.0).sqrt();
                let j = idx(n + 1, 1 - spin);
                h_rot[i][j] += -0.5 * coupling * amplitude;
                h_rot[j][i] += -0.5 * coupling * amplitude;
            }
        }
    }
    let eigen_rot = SymmetricEigensystem::diagonalize(h_rot).unwrap();
    // σz operator (full Pauli matrix)
    let mut sz_op = vec![vec![0.0; dim]; dim];
    for n in 0..cutoff {
        sz_op[idx(n, 0)][idx(n, 0)] = -1.0;
        sz_op[idx(n, 1)][idx(n, 1)] = 1.0;
    }
    let rot_sz = thermal_expectation(&eigen_rot, &sz_op, beta);

    assert!(
        (phys_sx - rot_sz).abs() < 1e-10,
        "physical ⟨σx⟩ = {phys_sx:.12} vs rotated ⟨σz⟩ = {rot_sz:.12}",
    );
}
