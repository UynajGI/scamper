//! Validation of Swendsen-Wang cluster updates for continuous O(N) spins.
//!
//! `SWCore` requires `H: ClusterModel`, and `ONModel<D>` implements the SW
//! contract for continuous spins (global reflection auxiliary, per-cluster
//! independent flip-or-identity decisions, Wolff embedded-Ising bond
//! activation `1 − exp(−2βJ (s·r)(s'·r))`). This file validates that offered
//! path, closing the "SW continuous-spin cluster updates remain unvalidated"
//! residue with real evidence:
//!
//! - **Exact quadrature** (XY 4-ring): ⟨E⟩, ⟨m²⟩, ⟨cos(θ1−θ3)⟩ vs the
//!   spectral quadrature reference shared with `over_relaxation.rs`, 8 seeds,
//!   |z| < 4 per seed, |z̄| < 2 per observable.
//! - **Analytic limits** (O(3) 4-ring): β → 0 gives exactly ⟨E⟩ = 0 and
//!   ⟨m²⟩ = 1/N for independent uniform spins; strong coupling (β = 8)
//!   approaches the harmonic result ⟨E⟩ → −J·N_edges + N_edges/(2β).
//! - **Cross-solver** (8×8): SW vs Wolff on O(2) and O(3), pooled cross-solver
//!   z on ⟨E⟩ and ⟨m²⟩ — Wolff for O(N) is independently validated against
//!   the exact Langevin (O(3)) and Bessel-ratio (O(2)) results.

use super::common::{exact_xy_ring4_moments, zscore_seed_count};
use cmc_rs::{
    build_chain, Algorithm, ClassicalMC, Initializable, Measurable, ONModel, SWCore,
    SimulationPhase, System, WolffCore,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

type Rng = Xoshiro256PlusPlus;

/// Binned point estimate (mean, stderr).
type Estimate = (f64, f64);

const COUPLING: f64 = 1.0;
const THERM: usize = 2_000;
const MEAS: usize = 40_000;
const BIN: usize = 400;

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

/// SW run on the XY 4-ring: ((E, se), (m², se), (cos(θ1−θ3), se)).
fn run_sw_xy_ring(beta: f64, seed: u64) -> ((f64, f64), (f64, f64), (f64, f64)) {
    let model = ONModel::<2>::new(COUPLING);
    let mut system = System::new(build_chain(4, true), 2, 0.0, beta);
    let mut rng = Rng::seed_from_u64(seed);
    for site in 0..system.n_sites() {
        let spin = model.random_spin(&mut rng);
        system.spin_at_mut(site, 2).copy_from_slice(&spin);
    }
    system.recompute_energy(&model);

    let mut kernel = SWCore::new();
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
fn sw_continuous_matches_exact_quadrature_xy_ring() {
    for &beta in &[0.6_f64, 1.2] {
        let exact = exact_xy_ring4_moments(beta, COUPLING, 96);
        let n_seeds = zscore_seed_count(8);
        let mut z_max = 0.0f64;
        let mut z_sum = [0.0f64; 3];
        for seed in 0..n_seeds as u64 {
            let ((e, e_se), (m2, m2_se), (c13, c13_se)) = run_sw_xy_ring(beta, 0x350 + seed);
            let zs = [
                (e - exact.0) / e_se.max(1e-12),
                (m2 - exact.1) / m2_se.max(1e-12),
                (c13 - exact.2) / c13_se.max(1e-12),
            ];
            for (index, z) in zs.iter().enumerate() {
                z_max = z_max.max(z.abs());
                z_sum[index] += z;
            }
        }
        let n = n_seeds as f64;
        let mean_z: Vec<f64> = z_sum.iter().map(|z| z / n).collect();
        eprintln!(
            "[sw-continuous β={beta}] exact ({:.5}, {:.5}, {:.5}) | max|z| = {z_max:.2}, \
             z̄ = {:.2}/{:.2}/{:.2}",
            exact.0, exact.1, exact.2, mean_z[0], mean_z[1], mean_z[2]
        );
        assert!(
            z_max < 4.0,
            "β={beta}: SW continuous max |z| = {z_max:.2} vs exact quadrature"
        );
        assert!(
            mean_z.iter().all(|z| z.abs() < 2.0),
            "β={beta}: SW continuous mean z = {mean_z:?}"
        );
    }
}

/// SW run on the O(3) 4-ring: ((E, se), (m², se)).
fn run_sw_o3_ring(beta: f64, seed: u64, sweeps: usize, binsize: usize) -> (Estimate, Estimate) {
    let model = ONModel::<3>::new(COUPLING);
    let mut system = System::new(build_chain(4, true), 3, 0.0, beta);
    let mut rng = Rng::seed_from_u64(seed);
    for site in 0..system.n_sites() {
        let spin = model.random_spin(&mut rng);
        system.spin_at_mut(site, 3).copy_from_slice(&spin);
    }
    system.recompute_energy(&model);

    let mut kernel = SWCore::new();
    for _ in 0..THERM {
        kernel.sweep_with_phase(&mut system, &model, &mut rng, SimulationPhase::Measurement);
    }
    let mut energies = Vec::with_capacity(sweeps);
    let mut m2 = Vec::with_capacity(sweeps);
    for _ in 0..sweeps {
        kernel.sweep_with_phase(&mut system, &model, &mut rng, SimulationPhase::Measurement);
        energies.push(system.energy);
        let magnetization = model.magnetization(&system.spins);
        m2.push(magnetization * magnetization);
    }
    (binned_stats(&energies, binsize), binned_stats(&m2, binsize))
}

#[test]
fn sw_continuous_o3_analytic_limits() {
    // β → 0: exactly ⟨E⟩ = 0 and ⟨m²⟩ = 1/N for independent uniform spins.
    let ((e, e_se), (m2, m2_se)) = run_sw_o3_ring(0.001, 0x4E7, 20_000, 200);
    assert!(
        e.abs() < 4.0 * e_se.max(0.02) + 0.02,
        "β→0 SW O(3) ⟨E⟩ = {e:.4} ± {e_se:.4}, expected 0"
    );
    assert!(
        (m2 - 0.25).abs() < 4.0 * m2_se.max(0.005) + 0.01,
        "β→0 SW O(3) ⟨m²⟩ = {m2:.4} ± {m2_se:.4}, expected 1/N = 0.25"
    );

    // Strong coupling β = 8: harmonic (spin-wave) result. The O(3) ring has
    // 2 transverse modes per site minus 2 global-rotation zero modes ⇒
    // 2N−2 quadratic modes, each contributing kT/2:
    // ⟨E⟩ = −J·N + (N−1)/β; m² → 1. Anharmonic remainder is O(1/β²).
    let beta = 8.0;
    let n_sites = 4.0;
    let harmonic_energy = -COUPLING * n_sites + (n_sites - 1.0) / beta;
    let ((e, _), (m2, _)) = run_sw_o3_ring(beta, 0x4E8, 20_000, 200);
    assert!(
        (e - harmonic_energy).abs() < 0.05,
        "β=8 SW O(3) ⟨E⟩ = {e:.4}, harmonic {harmonic_energy:.4}"
    );
    assert!(m2 > 0.9, "β=8 SW O(3) ⟨m²⟩ = {m2:.4}, expected → 1");
}

// ── Cross-solver vs Wolff on 8×8 ───────────────────────────────────────────

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
        .run_one::<ClassicalMC<ONModel<D>, A>>(&params);
    let energy = results.get("Energy").expect("Energy observable");
    let m2 = results.get("M2").expect("M2 observable");
    ((energy.mean, energy.stderr), (m2.mean, m2.stderr))
}

fn assert_cross_solver<const D: usize>(beta: f64, label: &str) {
    let n_seeds = zscore_seed_count(8);
    let sw: Vec<_> = (0..n_seeds as u64)
        .map(|seed| run_8x8::<D, SWCore>(beta, 0x5B00 + seed))
        .collect();
    let wolff: Vec<_> = (0..n_seeds as u64)
        .map(|seed| run_8x8::<D, WolffCore>(beta, 0x5B10 + seed))
        .collect();
    let pool = |values: &[(f64, f64)]| {
        let n = values.len() as f64;
        let mean = values.iter().map(|(m, _)| m).sum::<f64>() / n;
        let sem = values.iter().map(|(_, s)| s * s).sum::<f64>().sqrt() / n;
        (mean, sem)
    };
    let sw_e: Vec<_> = sw.iter().map(|(e, _)| *e).collect();
    let wolff_e: Vec<_> = wolff.iter().map(|(e, _)| *e).collect();
    let sw_m2: Vec<_> = sw.iter().map(|(_, m)| *m).collect();
    let wolff_m2: Vec<_> = wolff.iter().map(|(_, m)| *m).collect();
    let (se, se_se) = pool(&sw_e);
    let (we, we_se) = pool(&wolff_e);
    let (sm, sm_se) = pool(&sw_m2);
    let (wm, wm_se) = pool(&wolff_m2);
    let z_e = (se - we) / (se_se * se_se + we_se * we_se).sqrt();
    let z_m2 = (sm - wm) / (sm_se * sm_se + wm_se * wm_se).sqrt();
    eprintln!(
        "[sw-continuous cross {label}] E: sw {se:.4}±{se_se:.4} wolff {we:.4}±{we_se:.4} \
         z={z_e:.2} | m²: sw {sm:.4}±{sm_se:.4} wolff {wm:.4}±{wm_se:.4} z={z_m2:.2}"
    );
    assert!(z_e.abs() < 4.0, "{label}: ⟨E⟩ pooled-z = {z_e:.2}");
    assert!(z_m2.abs() < 4.0, "{label}: ⟨m²⟩ pooled-z = {z_m2:.2}");
}

#[test]
fn sw_continuous_matches_wolff_8x8_xy() {
    assert_cross_solver::<2>(0.90, "XY β=0.90");
}

#[test]
fn sw_continuous_matches_wolff_8x8_heisenberg() {
    assert_cross_solver::<3>(0.90, "O(3) β=0.90");
}
