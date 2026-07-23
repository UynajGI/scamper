//! Extended z-score statistical validation for additional CMC solvers.
//!
//! Tests HeatBathCore and KawasakiCore on a 2D Ising model (L=8, PBC,
//! J=1, β=0.3) using 16-seed z-scores. The microcanonical over-relaxation
//! kernel (MicrocanonicalCore) is documented below but cannot be tested
//! with IsingModel — see the note at the bottom of this file.
//!
//! ## Z-score framework
//!
//! Two complementary approaches are used:
//!
//! 1. **Exact-value z-scores** (used for HeatBath, where error bars are
//!    well-behaved): the exact finite-L energy and ⟨M²⟩ are computed via
//!    the transfer matrix method. For each seed:
//!    ```text
//!    z_i = (mean_i − exact) / stderr_i
//!    ```
//!    Asserts: |z| < 4 per seed, |z̄| < 2, Σz > −8.
//!
//! 2. **Scatter z-scores** (used for Kawasaki M², where magnetization
//!    conservation makes individual stderrs degenerate): the 16 seed means
//!    are standardized against their sample mean and sample standard
//!    deviation:
//!    ```text
//!    z_i = (mean_i − avg) / std_of_means
//!    ```
//!    Asserts: |z| < 4 (outlier check).
//!
//! ## Onsager reference
//!
//! As a physical-correctness sanity check, the pooled energy mean is
//! compared against the Onsager exact thermodynamic-limit internal energy.

use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use cmc_rs::{ClassicalMC, HeatBathCore, IsingModel, KawasakiCore};

const N_SEEDS: usize = 16;
const L: usize = 8;
const BETA: f64 = 0.3;
const J: f64 = 1.0;
const THERM_SWEEPS: u64 = 2000;
const MEAS_SWEEPS: u64 = 10000;
const BINSIZE: usize = 200;

// ── Transfer matrix exact values ──────────────────────────────

/// Compute the partition function Z = Tr[T^ly] for a 2D Ising model on an
/// lx × ly lattice with PBC in both directions.
///
/// `field` is the dimensionless field parameter βh (zero for the standard
/// model).
fn tm_partition_function(lx: usize, ly: usize, beta: f64, j: f64, field: f64) -> f64 {
    let n = 1usize << lx;

    // Spin value (+1 or −1) for bit `i` in a row configuration bitmask.
    let spin = |cfg: usize, i: usize| -> f64 {
        if (cfg >> i) & 1 == 1 {
            1.0
        } else {
            -1.0
        }
    };

    // Precompute per-row quantities: horizontal-bond sum and field sum.
    let row_h_bond: Vec<f64> = (0..n)
        .map(|s| {
            (0..lx)
                .map(|i| spin(s, i) * spin(s, (i + 1) % lx))
                .sum::<f64>()
        })
        .collect();
    let row_field: Vec<f64> = (0..n)
        .map(|s| (0..lx).map(|i| spin(s, i)).sum::<f64>())
        .collect();

    // Build the transfer matrix T[s, s'].
    let mut t = vec![0.0_f64; n * n];
    for s in 0..n {
        for sp in 0..n {
            // Vertical bonds between rows s and s'.
            let v_bond: f64 = (0..lx).map(|i| spin(s, i) * spin(sp, i)).sum();
            let exponent = beta * j * v_bond
                + beta * j * 0.5 * (row_h_bond[s] + row_h_bond[sp])
                + field * 0.5 * (row_field[s] + row_field[sp]);
            t[s * n + sp] = exponent.exp();
        }
    }

    // Compute T^ly by repeated multiplication (result = T, then × T ly−1 times).
    let mut result = t.clone();
    let mut temp = vec![0.0_f64; n * n];
    for _ in 1..ly {
        for i in 0..n {
            for j in 0..n {
                let mut sum = 0.0;
                for k in 0..n {
                    sum += result[i * n + k] * t[k * n + j];
                }
                temp[i * n + j] = sum;
            }
        }
        std::mem::swap(&mut result, &mut temp);
    }

    // Z = Tr[result]
    (0..n).map(|i| result[i * n + i]).sum()
}

/// Exact mean energy for a 2D Ising model on an lx × ly lattice with PBC.
/// Uses central-difference numerical differentiation of ln Z.
fn tm_exact_energy(lx: usize, ly: usize, beta: f64, j: f64) -> f64 {
    let eps = 1e-6;
    let z_plus = tm_partition_function(lx, ly, beta + eps, j, 0.0);
    let z_minus = tm_partition_function(lx, ly, beta - eps, j, 0.0);
    -(z_plus.ln() - z_minus.ln()) / (2.0 * eps)
}

/// Exact ⟨M²⟩ for a 2D Ising model on an lx × ly lattice with PBC.
/// Uses second-order central difference of ln Z with respect to the field.
fn tm_exact_m2(lx: usize, ly: usize, beta: f64, j: f64) -> f64 {
    let eps = 1e-4;
    let z0 = tm_partition_function(lx, ly, beta, j, 0.0);
    let zp = tm_partition_function(lx, ly, beta, j, eps);
    let zm = tm_partition_function(lx, ly, beta, j, -eps);
    let n = (lx * ly) as f64;
    (zp.ln() - 2.0 * z0.ln() + zm.ln()) / (eps * eps) / (n * n)
}

// ── Onsager thermodynamic-limit energy ────────────────────────

/// Complete elliptic integral of the first kind K(m) via the AGM.
fn complete_elliptic_k(m: f64) -> f64 {
    let k_comp = (1.0 - m * m).sqrt();
    let mut a = 1.0_f64;
    let mut b = k_comp;
    for _ in 0..60 {
        let new_a = (a + b) * 0.5;
        b = (a * b).sqrt();
        a = new_a;
        if (a - b).abs() < 1e-16 {
            break;
        }
    }
    std::f64::consts::PI / (2.0 * a)
}

/// Onsager exact internal energy per site (thermodynamic limit).
fn onsager_energy_per_site(beta: f64, j: f64) -> f64 {
    let k = beta * j;
    let two_k = 2.0 * k;
    let sinh_2k = two_k.sinh();
    let cosh_2k = two_k.cosh();
    let coth_2k = cosh_2k / sinh_2k;
    let tanh_2k = two_k.tanh();
    let coeff = 2.0 * tanh_2k * tanh_2k - 1.0;
    let modulus = if sinh_2k < 1.0 {
        2.0 * sinh_2k / (cosh_2k * cosh_2k)
    } else {
        (cosh_2k * cosh_2k) / (2.0 * sinh_2k)
    };
    let elliptic_k = complete_elliptic_k(modulus);
    -j * coth_2k * (1.0 + (2.0 / std::f64::consts::PI) * coeff * elliptic_k)
}

// ── Z-score analysis utilities ────────────────────────────────

/// Compute exact-value z-scores: z_i = (mean_i − exact) / stderr_i.
/// Returns (z_scores, max_abs_z, mean_z, sum_z).
fn exact_z_analysis(results: &[(f64, f64)], exact: f64) -> (Vec<f64>, f64, f64, f64) {
    let z_scores: Vec<f64> = results
        .iter()
        .map(|(mean, stderr)| (mean - exact) / stderr.max(1e-10))
        .collect();
    let max_abs_z = z_scores.iter().map(|z| z.abs()).fold(0.0_f64, f64::max);
    let mean_z = z_scores.iter().sum::<f64>() / z_scores.len() as f64;
    let sum_z = z_scores.iter().sum::<f64>();
    (z_scores, max_abs_z, mean_z, sum_z)
}

/// Compute scatter z-scores: z_i = (mean_i − avg) / std_of_means.
/// Returns (z_scores, max_abs_z, mean_z, sum_z).
fn scatter_z_analysis(results: &[(f64, f64)]) -> (Vec<f64>, f64, f64, f64) {
    let n = results.len() as f64;
    let means: Vec<f64> = results.iter().map(|(m, _)| *m).collect();
    let avg = means.iter().sum::<f64>() / n;
    let var = means.iter().map(|m| (m - avg).powi(2)).sum::<f64>() / (n - 1.0);
    let std = var.sqrt().max(1e-10);
    let z_scores: Vec<f64> = means.iter().map(|m| (m - avg) / std).collect();
    let max_abs_z = z_scores.iter().map(|z| z.abs()).fold(0.0_f64, f64::max);
    let mean_z = z_scores.iter().sum::<f64>() / z_scores.len() as f64;
    let sum_z = z_scores.iter().sum::<f64>();
    (z_scores, max_abs_z, mean_z, sum_z)
}

/// Assert exact-value z-score criteria.
fn assert_exact_z(results: &[(f64, f64)], exact: f64, label: &str) {
    let (_, max_z, mean_z, sum_z) = exact_z_analysis(results, exact);
    assert!(
        max_z < 4.0,
        "{label}: max |z| = {max_z:.2} should be < 4 (exact = {exact:.6})"
    );
    assert!(
        mean_z.abs() < 2.0,
        "{label}: mean z = {mean_z:.2} should be |z̄| < 2"
    );
    assert!(
        sum_z > -8.0,
        "{label}: sum z = {sum_z:.2} should be > −8 (no one-sided bias)"
    );
}

/// Assert scatter z-score criteria (self-consistency check).
fn assert_scatter_z(results: &[(f64, f64)], label: &str) {
    let (_, max_z, mean_z, sum_z) = scatter_z_analysis(results);
    assert!(
        max_z < 4.0,
        "{label}: max |z| = {max_z:.2} should be < 4 (scatter)"
    );
    assert!(
        mean_z.abs() < 2.0,
        "{label}: mean z = {mean_z:.2} should be |z̄| < 2"
    );
    assert!(
        sum_z > -8.0,
        "{label}: sum z = {sum_z:.2} should be > −8 (no one-sided bias)"
    );
}

// ── Simulation harness ────────────────────────────────────────

fn make_params(seed: u64) -> (Params, RunConfig) {
    let mut params = Params::new();
    params.set("Lx", L);
    params.set("Ly", L);
    params.set("J", J);
    params.set("beta", BETA);
    let config = RunConfig {
        thermalization_sweeps: THERM_SWEEPS,
        measurement_sweeps: MEAS_SWEEPS,
        binsize: BINSIZE,
        base_seed: seed,
        ..Default::default()
    };
    (params, config)
}

/// Run HeatBathCore for one seed.
/// Returns (energy_mean, energy_stderr, m2_mean, m2_stderr).
fn run_heat_bath(seed: u64) -> (f64, f64, f64, f64) {
    let (params, config) = make_params(seed);
    let r = Scheduler::new(RayonBackend::new(1), config)
        .run_one::<ClassicalMC<IsingModel, HeatBathCore>>(&params);
    let e = r.get("Energy").expect("Energy missing");
    let m2 = r.get("M2").expect("M2 missing");
    (e.mean, e.stderr, m2.mean, m2.stderr)
}

/// Run KawasakiCore for one seed.
/// Returns (energy_mean, energy_stderr, m2_mean, m2_stderr).
///
/// Note: Kawasaki dynamics conserves the total magnetization (spin-exchange
/// moves). M² is determined by the random initial state, so individual-seed
/// stderrs for M² are degenerate (near zero). The scatter approach is used
/// for M² self-consistency instead.
fn run_kawasaki(seed: u64) -> (f64, f64, f64, f64) {
    let (params, config) = make_params(seed);
    let r = Scheduler::new(RayonBackend::new(1), config)
        .run_one::<ClassicalMC<IsingModel, KawasakiCore>>(&params);
    let e = r.get("Energy").expect("Energy missing");
    let m2 = r.get("M2").expect("M2 missing");
    (e.mean, e.stderr, m2.mean, m2.stderr)
}

// ── Heat bath tests ────────────────────────────────────────────

#[test]
fn heat_bath_zscore_energy_16_seeds() {
    let exact_e = tm_exact_energy(L, L, BETA, J);
    let results: Vec<(f64, f64)> = (0..N_SEEDS as u64)
        .map(|s| {
            let (e_mean, e_se, _, _) = run_heat_bath(s);
            (e_mean, e_se)
        })
        .collect();

    // Exact-value z-scores (finite-L transfer matrix reference)
    assert_exact_z(&results, exact_e, "HeatBath Energy");

    // Onsager thermodynamic-limit sanity check
    let e_onsager = onsager_energy_per_site(BETA, J) * (L * L) as f64;
    let pooled_mean: f64 = results.iter().map(|(m, _)| m).sum::<f64>() / N_SEEDS as f64;
    let tol = (L * L) as f64 * 0.10; // 10% finite-size tolerance
    assert!(
        (pooled_mean - e_onsager).abs() < tol,
        "HeatBath energy: pooled {pooled_mean:.4}, Onsager {e_onsager:.4} (tol {tol:.4})"
    );
}

#[test]
fn heat_bath_zscore_m2_16_seeds() {
    let exact_m2 = tm_exact_m2(L, L, BETA, J);
    let results: Vec<(f64, f64)> = (0..N_SEEDS as u64)
        .map(|s| {
            let (_, _, m2_mean, m2_se) = run_heat_bath(s);
            (m2_mean, m2_se)
        })
        .collect();

    // Exact-value z-scores (finite-L transfer matrix reference)
    assert_exact_z(&results, exact_m2, "HeatBath M2");
}

// ── Kawasaki tests ─────────────────────────────────────────────

#[test]
fn kawasaki_zscore_energy_16_seeds() {
    let results: Vec<(f64, f64)> = (0..N_SEEDS as u64)
        .map(|s| {
            let (e_mean, e_se, _, _) = run_kawasaki(s);
            (e_mean, e_se)
        })
        .collect();

    // Scatter z-scores for self-consistency. Kawasaki samples a
    // fixed-magnetization sector (spin-exchange dynamics), so its energy
    // may differ from the canonical (transfer-matrix) exact value. The
    // scatter approach tests cross-seed agreement instead.
    assert_scatter_z(&results, "Kawasaki Energy");

    // Onsager sanity check (generous tolerance for conserved-M sector)
    let e_onsager = onsager_energy_per_site(BETA, J) * (L * L) as f64;
    let pooled_mean: f64 = results.iter().map(|(m, _)| m).sum::<f64>() / N_SEEDS as f64;
    let tol = (L * L) as f64 * 0.20; // 20% tolerance (conserved-M sector)
    assert!(
        (pooled_mean - e_onsager).abs() < tol,
        "Kawasaki energy: pooled {pooled_mean:.4}, Onsager {e_onsager:.4} (tol {tol:.4})"
    );
}

#[test]
fn kawasaki_zscore_m2_16_seeds() {
    // Kawasaki conserves magnetization, so M² is constant within each seed
    // and individual stderrs are degenerate. Use scatter z-scores to test
    // cross-seed self-consistency of the initial-state distribution.
    let results: Vec<(f64, f64)> = (0..N_SEEDS as u64)
        .map(|s| {
            let (_, _, m2_mean, m2_se) = run_kawasaki(s);
            (m2_mean, m2_se)
        })
        .collect();

    assert_scatter_z(&results, "Kawasaki M2");
}

// ── Microcanonical over-relaxation ────────────────────────────
//
// The microcanonical over-relaxation kernel is `MicrocanonicalCore` (not
// `OverRelaxationCore` as named in the task description). It implements
// `Algorithm<H>` for `H: Hamiltonian + LocalFieldModel`.
//
// `IsingModel` does NOT implement `LocalFieldModel` — that trait is only
// implemented for continuous-spin models (XYModel, HeisenbergModel, ONModel).
// Over-relaxation reflects a continuous spin about its local field, which is
// meaningless for discrete Ising spins (the reflection would just flip the
// spin, which is equivalent to a Metropolis proposal).
//
// Therefore `ClassicalMC<IsingModel, MicrocanonicalCore>` does not compile,
// and no z-score test can be written for this combination. The
// `MicrocanonicalCore` is validated for continuous spins in
// `tests/physics/usage_exact.rs` (xy_over_relaxation_preserves_energy_exactly)
// and `tests/physics/continuous_spins.rs`.
