//! Z-score and ergodicity validation for Particle NVT (Lennard-Jones).
//!
//! LennardJonesNvt<3> is the canonical ensemble for continuous particles
//! with the 12-6 Lennard-Jones pair potential. At low density (ρ = 0.1)
//! and moderate temperature (β = 0.5), the system is nearly ideal gas —
//! the equation of state is dominated by the first virial correction.
//!
//! ## Z-score framework
//!
//! Since the Lennard-Jones fluid has no closed-form exact solution, the
//! scatter z-score approach is used. For 16 independent seeds:
//! ```text
//! pooled_mean = (1/N) Σ mean_i
//! pooled_var  = (1/(N-1)) Σ (mean_i − pooled_mean)²
//! z_i = (mean_i − pooled_mean) / √(pooled_var/N + stderr_i²)
//! ```
//! Assert: |z_i| < 4 per seed, mean |z| < 2.
//!
//! ## Ergodicity test
//!
//! 4-seed pairwise convergence: ⟨Energy⟩ and ⟨EnergyPerParticle⟩ must agree
//! within 4 combined standard errors across independent seeds.

use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use cmc_rs::LennardJonesNvt;

const N_SEEDS: usize = 16;
const N_PARTICLES: usize = 64;
const DENSITY: f64 = 0.1;
const BETA: f64 = 0.5;
const CUTOFF: f64 = 2.5;
const MAX_DISPLACEMENT: f64 = 0.2;
const SIGMA: f64 = 1.0;
const EPSILON: f64 = 1.0;
const THERM_SWEEPS: u64 = 2000;
const MEAS_SWEEPS: u64 = 8000;
const BINSIZE: usize = 200;

// ── Shared parameter construction ─────────────────────────────

fn make_params(seed: u64) -> (Params, RunConfig) {
    let mut params = Params::new();
    params.set("n_particles", N_PARTICLES);
    params.set("density", DENSITY);
    params.set("beta", BETA);
    params.set("cutoff", CUTOFF);
    params.set("max_displacement", MAX_DISPLACEMENT);
    params.set("sigma", SIGMA);
    params.set("epsilon", EPSILON);
    let config = RunConfig {
        thermalization_sweeps: THERM_SWEEPS,
        measurement_sweeps: MEAS_SWEEPS,
        binsize: BINSIZE,
        base_seed: seed,
        ..Default::default()
    };
    (params, config)
}

/// Run one NVT simulation and return (energy_per_particle_mean, stderr).
fn run_particle_nvt_epp(seed: u64) -> (f64, f64) {
    let (params, config) = make_params(seed);
    let r = Scheduler::new(RayonBackend::new(1), config).run_one::<LennardJonesNvt<3>>(&params);
    let epp = r
        .get("EnergyPerParticle")
        .expect("EnergyPerParticle missing");
    (epp.mean, epp.stderr)
}

// ── Scatter z-score analysis ──────────────────────────────────

/// Compute scatter z-scores with per-seed stderr in the denominator.
///
/// ```text
/// pooled_mean = (1/N) Σ mean_i
/// pooled_var  = (1/(N-1)) Σ (mean_i − pooled_mean)²
/// z_i = (mean_i − pooled_mean) / √(pooled_var/N + stderr_i²)
/// ```
fn particle_scatter_z(results: &[(f64, f64)]) -> (Vec<f64>, f64, f64) {
    let n = results.len() as f64;
    let means: Vec<f64> = results.iter().map(|(m, _)| *m).collect();
    let pooled_mean = means.iter().sum::<f64>() / n;
    let pooled_var = means.iter().map(|m| (m - pooled_mean).powi(2)).sum::<f64>() / (n - 1.0);
    let z_scores: Vec<f64> = results
        .iter()
        .map(|(m, s)| {
            let denom = (pooled_var / n + s.powi(2)).sqrt();
            (m - pooled_mean) / denom.max(1e-15)
        })
        .collect();
    let max_abs_z = z_scores.iter().map(|z| z.abs()).fold(0.0_f64, f64::max);
    let mean_abs_z = z_scores.iter().map(|z| z.abs()).sum::<f64>() / n;
    (z_scores, max_abs_z, mean_abs_z)
}

/// Assert scatter z-score criteria: |z_i| < 4 per seed, mean |z| < 2.
fn assert_particle_z(results: &[(f64, f64)], label: &str) {
    let (zs, max_z, mean_abs_z) = particle_scatter_z(results);
    for (i, &z) in zs.iter().enumerate() {
        assert!(
            z.abs() < 4.0,
            "{label}: seed {i} |z| = {:.2} should be < 4, z = {z:.3}",
            z.abs()
        );
    }
    assert!(max_z < 4.0, "{label}: max |z| = {max_z:.2} should be < 4");
    assert!(
        mean_abs_z < 2.0,
        "{label}: mean |z| = {mean_abs_z:.2} should be < 2"
    );
}

// ── Z-score test ──────────────────────────────────────────────

/// 16-seed scatter z-score on EnergyPerParticle for LennardJonesNvt<3>.
///
/// At ρ = 0.1 and β = 0.5, the LJ fluid is nearly ideal gas — pair
/// interactions are weak perturbations on the dominant kinetic
/// contribution (not measured). The scatter approach tests cross-seed
/// self-consistency without requiring an exact reference value.
#[test]
#[ignore = "~100s — 16 independent seeds × 10000 sweeps each; run via --ignored"]
fn particle_nvt_energy_per_particle_zscore_16_seeds() {
    let results: Vec<(f64, f64)> = (0..N_SEEDS as u64).map(run_particle_nvt_epp).collect();
    assert_particle_z(&results, "ParticleNVT EnergyPerParticle");
}

// ── Ergodicity test ───────────────────────────────────────────

/// Check that two estimates agree within `n` combined standard errors.
fn within_n_sigma(n: f64, a_mean: f64, a_err: f64, b_mean: f64, b_err: f64) -> bool {
    let combined = (a_err * a_err + b_err * b_err).sqrt();
    if combined == 0.0 {
        (a_mean - b_mean).abs() < 1e-10
    } else {
        (a_mean - b_mean).abs() < n * combined
    }
}

/// Run one NVT simulation and return (energy, energy_stderr, epp_mean, epp_stderr).
fn run_particle_nvt_full(seed: u64) -> (f64, f64, f64, f64) {
    let (params, config) = make_params(seed);
    let r = Scheduler::new(RayonBackend::new(1), config).run_one::<LennardJonesNvt<3>>(&params);
    let e = r.get("Energy").expect("Energy missing");
    let epp = r
        .get("EnergyPerParticle")
        .expect("EnergyPerParticle missing");
    (e.mean, e.stderr, epp.mean, epp.stderr)
}

/// 4-seed ergodicity check: ⟨Energy⟩ and ⟨EnergyPerParticle⟩ must agree
/// within 4σ across independent seeds.
#[test]
fn particle_nvt_ergodicity_4_seeds() {
    let seeds = [42u64, 999, 7, 314];
    let estimates: Vec<_> = seeds.iter().map(|&s| run_particle_nvt_full(s)).collect();

    for i in 0..estimates.len() {
        for j in (i + 1)..estimates.len() {
            let (ei, ei_err, eppi, eppi_err) = estimates[i];
            let (ej, ej_err, eppj, eppj_err) = estimates[j];
            assert!(
                within_n_sigma(4.0, ei, ei_err, ej, ej_err),
                "ParticleNVT ⟨E⟩ disagree: seeds {} vs {}: {:.4}±{:.4} vs {:.4}±{:.4}",
                seeds[i],
                seeds[j],
                ei,
                ei_err,
                ej,
                ej_err
            );
            assert!(
                within_n_sigma(4.0, eppi, eppi_err, eppj, eppj_err),
                "ParticleNVT ⟨E/N⟩ disagree: seeds {} vs {}: {:.4}±{:.4} vs {:.4}±{:.4}",
                seeds[i],
                seeds[j],
                eppi,
                eppi_err,
                eppj,
                eppj_err
            );
        }
    }
}
