//! Production validation of microcanonical over-relaxation (`MicrocanonicalCore`).
//!
//! Over-relaxation alone is deterministic and energy-conserving, hence not
//! ergodic — that is physics, not a defect. Production status is earned for
//! the solver **as operated, in composition** with an ergodic kernel:
//!
//! - **A/D (machine precision):** the reflection map `R(s) = 2 (s·ĥ) ĥ − s` is
//!   an exact involution (`R∘R = id`, 1e-15), norm-preserving (1e-15),
//!   preserves the local-field projection (1e-15) — the pair energy with every
//!   neighbor — and is an exact isometry with |Jacobian| = 1: for O(2) the
//!   angle form θ' = 2φ − θ gives |dθ'/dθ| = 1 (a reflection, det −1); for
//!   O(3) the map matrix 2 ĥ ĥᵀ − I is orthogonal with det +1 — a π rotation
//!   about the field axis (the field direction is fixed, the orthogonal plane
//!   reversed). A deterministic involutive isometry satisfies detailed balance
//!   against any rotation-invariant measure because T(s→s') = T(s'→s) = 1 on
//!   each orbit. The kernel itself is verified to apply exactly this map
//!   (bit-identical to the manual sequential reflection) with per-update
//!   |ΔE| < 1e-12 through the transactional evaluator.
//! - **B/E:** Hybrid(Metropolis, Microcanonical) on the XY 4-ring vs an exact
//!   spectral quadrature reference — ⟨E⟩, ⟨m²⟩, ⟨cos(θ1−θ3)⟩ — from both hot
//!   and cold initializations (multi-init convergence), 8 seeds, |z| < 4.
//! - **F:** the same composition vs Wolff on 8×8 O(2) and O(3) lattices
//!   (pooled cross-solver z on ⟨E⟩ and ⟨m²⟩).
//! - **C:** analytic limits in composition: β → 0 gives exactly ⟨E⟩ = 0 and
//!   ⟨m²⟩ = 1/N; strong coupling β = 8 approaches the harmonic (spin-wave)
//!   result ⟨E⟩ → −J·N_edges + N_edges/(2β) with ⟨m²⟩ → 1.

use super::common::{assert_close, exact_xy_ring4_moments, zscore_seed_count};
use cmc_rs::{
    build_chain, Algorithm, EnergyPatch, HybridCore, Initializable, LocalFieldModel, Measurable,
    MetropolisCore, MicrocanonicalCore, ONModel, OPSSStrategy, SimulationPhase, SiteSpinMove,
    System, TrialEvaluator, VisitSchedule, WolffCore,
};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

type Rng = Xoshiro256PlusPlus;

/// Binned point estimate (mean, stderr).
type Estimate = (f64, f64);

// ── A/D: the reflection map as an exact measure-preserving involution ───────

/// The kernel's reflection of `old` through the local field `field`
/// (mirrors `MicrocanonicalCore::sweep_with_phase` line for line).
fn reflect(old: &[f64], field: &[f64]) -> Vec<f64> {
    let norm_squared: f64 = field.iter().map(|value| value * value).sum();
    assert!(
        norm_squared > 1e-28,
        "zero-field sites are skipped by the kernel"
    );
    let projection = old
        .iter()
        .zip(field)
        .map(|(spin, value)| spin * value)
        .sum::<f64>()
        / norm_squared;
    (0..old.len())
        .map(|component| 2.0 * projection * field[component] - old[component])
        .collect()
}

fn random_spin_system<const D: usize>(coupling: f64, beta: f64, seed: u64) -> (System, ONModel<D>) {
    let mut system = System::new(build_chain(4, true), D, 0.0, beta);
    let model = ONModel::<D>::new(coupling);
    let mut rng = Rng::seed_from_u64(seed);
    for site in 0..system.n_sites() {
        let spin = model.random_spin(&mut rng);
        system.spin_at_mut(site, D).copy_from_slice(&spin);
    }
    system.recompute_energy(&model);
    (system, model)
}

#[test]
fn over_relaxation_reflection_is_an_exact_measure_preserving_involution() {
    for (seed, coupling) in [(0x11C, 0.9_f64), (0x11D, 1.7)] {
        // O(2): angle form θ' = 2φ − θ is an orientation-reversing isometry of
        // the circle (|dθ'/dθ| = 1): measure-preserving, deterministic DB.
        let (system, model) = random_spin_system::<2>(coupling, 0.7, seed);
        let mut field = [0.0f64; 2];
        for site in 0..system.n_sites() {
            model.local_field(&system.spins, &system.lattice, site, &mut field);
            let old = system.spin_at(site, 2);
            let reflected = reflect(old, &field);
            // Involution R∘R = id ⇒ deterministic detailed balance
            // T(s→s') = T(s'→s) = 1.
            let twice = reflect(&reflected, &field);
            for component in 0..2 {
                assert_close(twice[component], old[component], 1e-15);
            }
            // Norm preservation (the spin stays on the unit circle).
            assert_close(reflected.iter().map(|v| v * v).sum::<f64>(), 1.0, 1e-15);
            // Field projection preserved ⇒ bond energy with every neighbor
            // unchanged: R(s)·h = s·h.
            let old_projection: f64 = old.iter().zip(&field).map(|(a, b)| a * b).sum();
            let new_projection: f64 = reflected.iter().zip(&field).map(|(a, b)| a * b).sum();
            assert_close(new_projection, old_projection, 1e-15);
            // Angle form: θ' + θ = 2φ (mod 2π) — an exact reflection.
            let phi = field[1].atan2(field[0]);
            let theta = old[1].atan2(old[0]);
            let theta_prime = reflected[1].atan2(reflected[0]);
            let residual = (theta_prime + theta - 2.0 * phi)
                .rem_euclid(std::f64::consts::TAU)
                .min((2.0 * phi - theta - theta_prime).rem_euclid(std::f64::consts::TAU));
            assert!(residual < 1e-14, "angle reflection residual {residual:.3e}");
        }

        // O(3): the map matrix H = 2 ĥ ĥᵀ − I is orthogonal (|Jacobian| = 1,
        // measure-preserving) with det H = +1: it fixes the field axis and
        // reverses the orthogonal plane — a π rotation about ĥ, the textbook
        // over-relaxation move. The kernel's reflection is exactly H·s.
        let (system, model) = random_spin_system::<3>(coupling, 0.7, seed);
        let mut field = [0.0f64; 3];
        for site in 0..system.n_sites() {
            model.local_field(&system.spins, &system.lattice, site, &mut field);
            let old = system.spin_at(site, 3);
            let reflected = reflect(old, &field);
            let twice = reflect(&reflected, &field);
            for component in 0..3 {
                assert_close(twice[component], old[component], 1e-15);
            }
            assert_close(reflected.iter().map(|v| v * v).sum::<f64>(), 1.0, 1e-15);

            let norm = field.iter().map(|v| v * v).sum::<f64>().sqrt();
            let unit: Vec<f64> = field.iter().map(|v| v / norm).collect();
            let mut hh = [[0.0f64; 3]; 3];
            for row in 0..3 {
                for column in 0..3 {
                    hh[row][column] = 2.0 * unit[row] * unit[column] - f64::from(row == column);
                }
            }
            // HᵀH = I (orthogonal ⇒ |Jacobian| = 1).
            for row in 0..3 {
                for column in 0..3 {
                    let dot: f64 = (0..3).map(|k| hh[k][row] * hh[k][column]).sum();
                    assert_close(dot, f64::from(row == column), 1e-15);
                }
            }
            // det H = +1 in three dimensions (π rotation about ĥ; the
            // perpendicular plane is reversed, the field axis fixed).
            let determinant = hh[0][0] * (hh[1][1] * hh[2][2] - hh[1][2] * hh[2][1])
                - hh[0][1] * (hh[1][0] * hh[2][2] - hh[1][2] * hh[2][0])
                + hh[0][2] * (hh[1][0] * hh[2][1] - hh[1][1] * hh[2][0]);
            assert_close(determinant, 1.0, 1e-15);
            // The kernel's vector reflection is exactly H·s.
            for component in 0..3 {
                let matrix_product: f64 = (0..3).map(|k| hh[component][k] * old[k]).sum();
                assert_close(reflected[component], matrix_product, 1e-15);
            }
        }
    }
}

/// The kernel's sequential sweep must equal the manual per-site reflection
/// bit-for-bit, and every transactional update must report |ΔE| < 1e-12
/// (machine-precision per-update energy conservation).
fn kernel_equals_manual_reflection<const D: usize>(label: &str, coupling: f64, beta: f64) {
    let (mut kernel_system, _) = random_spin_system::<D>(coupling, beta, 0x5EED);
    let (mut manual_system, model) = random_spin_system::<D>(coupling, beta, 0x5EED);
    assert_eq!(kernel_system.spins, manual_system.spins);

    // Kernel sweep (deterministic sequential schedule; no randomness consumed).
    let mut kernel = MicrocanonicalCore::new().with_visit_schedule(VisitSchedule::Sequential);
    let mut rng = Rng::seed_from_u64(0);
    kernel.sweep_with_phase(
        &mut kernel_system,
        &model,
        &mut rng,
        SimulationPhase::Measurement,
    );

    // Manual replication through the transactional evaluator.
    let mut field = vec![0.0f64; D];
    let mut patch = EnergyPatch::default();
    for site in 0..manual_system.n_sites() {
        model.local_field(
            &manual_system.spins,
            &manual_system.lattice,
            site,
            &mut field,
        );
        if field.iter().map(|v| v * v).sum::<f64>() < 1e-28 {
            continue;
        }
        let old = manual_system.spin_at(site, D).to_vec();
        let reflected = reflect(&old, &field);
        let movement = SiteSpinMove::new(site, cmc_rs::Spin::from_slice(&reflected));
        <System as TrialEvaluator<ONModel<D>, SiteSpinMove>>::evaluate_trial(
            &manual_system,
            &model,
            &movement,
            &mut patch,
        );
        assert!(
            patch.delta_energy.abs() < 1e-12,
            "{label}: per-update |ΔE| = {:.3e} at site {site}",
            patch.delta_energy
        );
        <System as TrialEvaluator<ONModel<D>, SiteSpinMove>>::commit_trial(
            &mut manual_system,
            &movement,
            &patch,
        );
    }

    // Bit-identical state: the kernel is exactly the reflection map, and the
    // energy is conserved to machine precision against a full recomputation.
    assert_eq!(
        kernel_system.spins, manual_system.spins,
        "{label}: kernel sweep must equal the manual sequential reflection"
    );
    assert_close(manual_system.energy_error(&model), 0.0, 1e-12);
    assert_close(kernel_system.energy_error(&model), 0.0, 1e-12);
    assert_close(kernel_system.energy, manual_system.energy, 1e-15);
}

#[test]
fn over_relaxation_kernel_applies_exactly_this_reflection_with_zero_energy_delta() {
    kernel_equals_manual_reflection::<2>("O(2)", 1.3, 0.8);
    kernel_equals_manual_reflection::<3>("O(3)", 1.3, 0.8);
    kernel_equals_manual_reflection::<3>("O(3) weak coupling", 0.4, 0.05);
}

// ── B/E: composition vs exact quadrature on the XY 4-ring ──────────────────

const COUPLING: f64 = 1.0;
const THERM: usize = 2_000;
const MEAS: usize = 60_000;
const BIN: usize = 300;

fn binned_stats(samples: &[f64], binsize: usize) -> (f64, f64) {
    let usable = samples.len() / binsize * binsize;
    assert!(usable >= 2 * binsize, "not enough samples for binning");
    let bins: Vec<f64> = samples[..usable]
        .chunks(binsize)
        .map(|chunk| chunk.iter().sum::<f64>() / chunk.len() as f64)
        .collect();
    let n = bins.len() as f64;
    let mean = bins.iter().sum::<f64>() / n;
    let variance = bins.iter().map(|b| (b - mean).powi(2)).sum::<f64>() / (n - 1.0);
    (mean, (variance / n).sqrt())
}

/// Hybrid(Metropolis, Microcanonical) run on the XY 4-ring.
/// Returns ((E, se), (m², se), (cos(θ1−θ3), se)).
fn run_composition(beta: f64, seed: u64, cold_start: bool) -> ((f64, f64), (f64, f64), (f64, f64)) {
    let model = ONModel::<2>::new(COUPLING);
    let mut system = System::new(build_chain(4, true), 2, 0.0, beta);
    let mut rng = Rng::seed_from_u64(seed);
    for site in 0..system.n_sites() {
        let spin = if cold_start {
            [1.0f64, 0.0]
        } else {
            let theta = rng.random::<f64>() * std::f64::consts::TAU;
            [theta.cos(), theta.sin()]
        };
        system.spin_at_mut(site, 2).copy_from_slice(&spin);
    }
    system.recompute_energy(&model);

    let mut kernel = HybridCore::new(
        MetropolisCore::with_strategy(OPSSStrategy::new().with_sigma(0.4)),
        MicrocanonicalCore::new(),
    );
    for _ in 0..THERM {
        kernel.sweep_with_phase(&mut system, &model, &mut rng, SimulationPhase::Measurement);
    }
    let mut energies = Vec::with_capacity(MEAS);
    let mut m2 = Vec::with_capacity(MEAS);
    let mut c13 = Vec::with_capacity(MEAS);
    for _ in 0..MEAS {
        kernel.sweep_with_phase(&mut system, &model, &mut rng, SimulationPhase::Measurement);
        energies.push(system.energy);
        let magnetization = model.magnetization(&system.spins);
        m2.push(magnetization * magnetization);
        c13.push(system.spins[0] * system.spins[4] + system.spins[1] * system.spins[5]);
    }
    (
        binned_stats(&energies, BIN),
        binned_stats(&m2, BIN),
        binned_stats(&c13, BIN),
    )
}

#[test]
fn over_relaxation_composition_matches_exact_quadrature_from_two_inits() {
    // Quadrature convergence anchor: doubling the grid changes nothing at the
    // 1e-8 level (spectral accuracy), so the reference is effectively exact.
    for &beta in &[0.6_f64, 1.2] {
        let coarse = exact_xy_ring4_moments(beta, COUPLING, 96);
        let fine = exact_xy_ring4_moments(beta, COUPLING, 128);
        for (a, b) in [(coarse.0, fine.0), (coarse.1, fine.1), (coarse.2, fine.2)] {
            assert!(
                (a - b).abs() < 1e-8,
                "quadrature not converged at beta={beta}"
            );
        }

        let n_seeds = zscore_seed_count(8);
        for (label, cold) in [("hot", false), ("cold", true)] {
            let mut z_max = 0.0f64;
            let mut z_sum = [0.0f64; 3];
            for seed in 0..n_seeds as u64 {
                let ((e, e_se), (m2, m2_se), (c13, c13_se)) =
                    run_composition(beta, 0x0E25 + 131 * seed, cold);
                let zs = [
                    (e - coarse.0) / e_se.max(1e-12),
                    (m2 - coarse.1) / m2_se.max(1e-12),
                    (c13 - coarse.2) / c13_se.max(1e-12),
                ];
                for (index, z) in zs.iter().enumerate() {
                    z_max = z_max.max(z.abs());
                    z_sum[index] += z;
                }
            }
            let mean_z: Vec<f64> = z_sum.iter().map(|z| z / n_seeds as f64).collect();
            eprintln!(
                "[over-relaxation β={beta} {label}] max|z| = {z_max:.2}, mean z = {:.2}/{:.2}/{:.2}",
                mean_z[0], mean_z[1], mean_z[2]
            );
            assert!(
                z_max < 4.0,
                "β={beta} {label}: max |z| = {z_max:.2} vs exact ({:.5}, {:.5}, {:.5})",
                coarse.0,
                coarse.1,
                coarse.2
            );
            assert!(
                mean_z.iter().all(|z| z.abs() < 2.0),
                "β={beta} {label}: mean z = {mean_z:?}"
            );
        }
    }
}

// ── C: analytic limits in composition ───────────────────────────────────────

#[test]
fn over_relaxation_composition_analytic_limits() {
    // High temperature (β → 0): exactly ⟨E⟩ = 0 and ⟨m²⟩ = 1/N for independent
    // uniform spins (closed form, independent of any quadrature).
    let ((e, e_se), _, _) = run_composition(0.001, 0x11A7, false);
    let ((_, _), (m2, m2_se), _) = run_composition(0.001, 0x11A8, true);
    assert!(
        e.abs() < 4.0 * e_se.max(0.02) + 0.02,
        "β→0 ⟨E⟩ = {e:.4} ± {e_se:.4}, expected 0"
    );
    assert!(
        (m2 - 0.25).abs() < 4.0 * m2_se.max(0.005) + 0.01,
        "β→0 ⟨m²⟩ = {m2:.4} ± {m2_se:.4}, expected 1/N = 0.25"
    );

    // Strong coupling (β = 8): harmonic (spin-wave) expansion of the XY ring.
    // The relative angles δ are Gaussian with the single global zero mode
    // removed: N−1 independent quadratic modes × kT/2 give
    // ⟨E⟩ = −J·N + (N−1)/(2β); m² → 1. Anharmonic remainder O(1/β²) ≈ 0.016.
    let beta = 8.0;
    let n_sites = 4.0;
    let harmonic_energy = -COUPLING * n_sites + (n_sites - 1.0) / (2.0 * beta);
    let ((e, e_se), (m2, _), _) = run_composition(beta, 0x11B9, true);
    assert!(
        (e - harmonic_energy).abs() < 0.05,
        "β=8 ⟨E⟩ = {e:.4} ± {e_se:.4}, harmonic {harmonic_energy:.4}"
    );
    assert!(m2 > 0.9, "β=8 ⟨m²⟩ = {m2:.4}, expected → 1");
}

// ── F: cross-solver — composition vs Wolff on 8×8 O(N) ─────────────────────

fn run_8x8<const D: usize, A>(beta: f64, seed: u64) -> (Estimate, Estimate)
where
    cmc_rs::ClassicalMC<ONModel<D>, A>: carlo_rs::FromParams,
    A: Algorithm<ONModel<D>> + Default,
{
    let mut params = carlo_rs::Params::new();
    params.set("Lx", 8_usize);
    params.set("Ly", 8_usize);
    params.set("J", 1.0_f64);
    params.set("beta", beta);
    let config = carlo_rs::RunConfig {
        thermalization_sweeps: 3_000,
        measurement_sweeps: 2_000,
        binsize: 100,
        base_seed: seed,
        ..Default::default()
    };
    let results = carlo_rs::Scheduler::new(carlo_rs::RayonBackend::new(1), config)
        .run_one::<cmc_rs::ClassicalMC<ONModel<D>, A>>(&params);
    let e = results.get("Energy").expect("Energy observable");
    let m2 = results.get("M2").expect("M2 observable");
    ((e.mean, e.stderr), (m2.mean, m2.stderr))
}

fn pooled(results: &[(Estimate, Estimate)]) -> (Estimate, Estimate) {
    let pool_scalar = |values: &[(f64, f64)]| {
        let n = values.len() as f64;
        let mean = values.iter().map(|(m, _)| m).sum::<f64>() / n;
        let sem = values.iter().map(|(_, s)| s * s).sum::<f64>().sqrt() / n;
        (mean, sem)
    };
    let energy: Vec<_> = results.iter().map(|(e, _)| *e).collect();
    let m2: Vec<_> = results.iter().map(|(_, m)| *m).collect();
    (pool_scalar(&energy), pool_scalar(&m2))
}

fn assert_cross_solver<const D: usize>(beta: f64, label: &str) {
    let n_seeds = zscore_seed_count(8);
    type Composition = HybridCore<MetropolisCore<OPSSStrategy>, MicrocanonicalCore>;
    let hybrid: Vec<_> = (0..n_seeds as u64)
        .map(|seed| run_8x8::<D, Composition>(beta, seed))
        .collect();
    let wolff: Vec<_> = (0..n_seeds as u64)
        .map(|seed| run_8x8::<D, WolffCore>(beta, 0xA000 + seed))
        .collect();
    let ((hybrid_e, hybrid_e_se), (hybrid_m2, hybrid_m2_se)) = pooled(&hybrid);
    let ((wolff_e, wolff_e_se), (wolff_m2, wolff_m2_se)) = pooled(&wolff);
    let z_e = (hybrid_e - wolff_e) / (hybrid_e_se * hybrid_e_se + wolff_e_se * wolff_e_se).sqrt();
    let z_m2 =
        (hybrid_m2 - wolff_m2) / (hybrid_m2_se * hybrid_m2_se + wolff_m2_se * wolff_m2_se).sqrt();
    eprintln!(
        "[over-relaxation cross {label}] E: hybrid {hybrid_e:.4}±{hybrid_e_se:.4} wolff \
         {wolff_e:.4}±{wolff_e_se:.4} z={z_e:.2} | m²: hybrid {hybrid_m2:.4}±{hybrid_m2_se:.4} \
         wolff {wolff_m2:.4}±{wolff_m2_se:.4} z={z_m2:.2}"
    );
    assert!(z_e.abs() < 4.0, "{label}: ⟨E⟩ pooled-z = {z_e:.2}");
    assert!(z_m2.abs() < 4.0, "{label}: ⟨m²⟩ pooled-z = {z_m2:.2}");
}

#[test]
fn over_relaxation_composition_matches_wolff_8x8_xy() {
    assert_cross_solver::<2>(0.90, "XY β=0.90");
}

#[test]
fn over_relaxation_composition_matches_wolff_8x8_heisenberg() {
    assert_cross_solver::<3>(0.90, "O(3) β=0.90");
}
