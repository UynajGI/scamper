//! Exact-distribution validation for the two heat-bath kernels.
//!
//! Closes two VALIDATION.md gaps:
//!
//! - Discrete `HeatBathCore` — "NOT validated: Exact energy comparison":
//!   ⟨E⟩ and ⟨m²⟩ multi-seed z-scores against FULL exact enumeration
//!   (`common::exact_ising_moments`) on N=4 and N=8 Ising PBC chains at
//!   β ∈ {0.3, 0.8}. (The Potts part of that gap stays open.)
//!
//! - Continuous `ContinuousHeatBathCore` — "NOT validated: Finite-T
//!   distribution, XY (O(2))": finite-field conditional moments on S¹ and
//!   S², including the XY/O(2) case.
//!
//! ## Field mapping (documented, per the core source)
//!
//! `ContinuousHeatBathCore` has no external-field API: a site's local field
//! is built only from its neighbor bonds (`ONModel::local_field` sums
//! J·w·s_neighbor). A single spin in an external field h with weight
//! e^(β h·s) is therefore realized as a TWO-site lattice with ONE bond
//! (J=1, weight 1) — each site's local field is exactly its partner, so the
//! concentration is κ = β·J·w = x. The core updates both sites, but by O(N)
//! rotational invariance the relative-angle distribution after every sweep
//! is EXACTLY the single-spin-in-field conditional p(cosθ) ∝ e^(x cosθ):
//! the partner is the (instantaneously fixed) field direction. In fact the
//! last spin updated in each sweep is a fresh exact conditional draw, so
//! per-sweep cosθ samples are effectively independent. The measured
//! ⟨cosθ⟩ ≡ ⟨s_z⟩ and ⟨cos²θ⟩ ≡ ⟨s_z²⟩ (field frame z = partner direction)
//! are therefore the field-conditioned moments, and the O(2) run doubles as
//! the XY pair-equilibrium exact reference (two-site partition function
//! Z ∝ I₀(x)), which covers the "XY pair/lattice equilibrium" item.
//!
//! ## Exact references (in-code, no special-function dependency)
//!
//! - O(2)/S¹ (von Mises): ⟨cosθ⟩ = I₁(x)/I₀(x), ⟨cos²θ⟩ = (1+I₂(x)/I₀(x))/2
//!   with modified Bessel functions I_ν from the power series.
//! - O(3)/S² (Langevin): ⟨cosθ⟩ = coth(x) − 1/x, ⟨cos²θ⟩ = 1 − 2L(x)/x.
//! - Limits: ⟨cosθ⟩ → x/D for x→0 (linear response; x/3 for O(3), x/2 for
//!   O(2)) and → 1 for x→∞ (saturation). Asserted for the analytic
//!   references AND at the MC grid endpoints (x = 0.2 and x = 5) through a
//!   window = analytic truncation gap + 4 pooled seed SEs.

use super::common::{exact_ising_moments, zscore_seed_count};
use cmc_rs::{
    build_chain, Algorithm, Bond, BondType, ContinuousHeatBathCore, ContinuousHeatBathable,
    CsrLattice, Hamiltonian, HeatBathCore, IsingModel, ONModel, System,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

const J: f64 = 1.0;
const N_SEEDS: usize = 16;
const THERM_SWEEPS: usize = 1000;
const MEAS_SWEEPS: usize = 12_000;
const BINSIZE: usize = 200;
/// Field strengths x = βh = βJ spanning 0.2–5 (linear to near-saturated).
const X_GRID: [f64; 5] = [0.2, 0.5, 1.5, 3.0, 5.0];

// ── Statistics helpers (established multi-seed z pattern) ─────────────

/// Mean and stderr of per-sweep values from bin means (bin stderr / √bins).
fn mean_stderr_binned(values: &[f64], binsize: usize) -> (f64, f64) {
    let n_bins = values.len() / binsize;
    assert!(n_bins >= 2, "need at least 2 bins, got {n_bins}");
    let bin_means: Vec<f64> = (0..n_bins)
        .map(|b| {
            let start = b * binsize;
            values[start..start + binsize].iter().sum::<f64>() / binsize as f64
        })
        .collect();
    let mean = bin_means.iter().sum::<f64>() / n_bins as f64;
    let variance = bin_means
        .iter()
        .map(|bin| (bin - mean) * (bin - mean))
        .sum::<f64>()
        / (n_bins - 1) as f64;
    (mean, (variance / n_bins as f64).sqrt())
}

/// Per-seed z criteria (mirrors `zscore_validation`/`zscore_extended`):
/// |z| < 4 per seed, |z̄| < 2, Σz > −2√n (scale-invariant one-sided bound).
/// Returns max |z| for progress logging.
fn assert_exact_z(results: &[(f64, f64)], exact: f64, label: &str) -> f64 {
    let z_scores: Vec<f64> = results
        .iter()
        .map(|(mean, stderr)| (mean - exact) / stderr.max(1e-10))
        .collect();
    let max_abs_z = z_scores.iter().fold(0.0_f64, |acc, z| acc.max(z.abs()));
    let sum_z = z_scores.iter().sum::<f64>();
    let mean_z = sum_z / z_scores.len() as f64;
    let sum_floor = -2.0 * (z_scores.len() as f64).sqrt();
    assert!(
        max_abs_z < 4.0,
        "{label}: max |z| = {max_abs_z:.2} (exact = {exact:.6})"
    );
    assert!(
        mean_z.abs() < 2.0,
        "{label}: mean z = {mean_z:.2} (exact = {exact:.6})"
    );
    assert!(
        sum_z > sum_floor,
        "{label}: Σz = {sum_z:.2} should be > {sum_floor:.2} (one-sided bias)"
    );
    max_abs_z
}

/// Pooled mean over independent seeds and its stderr (√(Σ seᵢ²)/n).
fn pooled_mean_stderr(results: &[(f64, f64)]) -> (f64, f64) {
    let n = results.len() as f64;
    let mean = results.iter().map(|(m, _)| m).sum::<f64>() / n;
    let variance = results.iter().map(|(_, se)| se * se).sum::<f64>() / (n * n);
    (mean, variance.sqrt())
}

// ── Part A: discrete heat bath vs full exact enumeration ─────────────

/// One `HeatBathCore` seed on an Ising chain; returns binned ⟨E⟩ and ⟨m²⟩.
fn run_heat_bath_ising(lattice: &CsrLattice, beta: f64, seed: u64) -> ((f64, f64), (f64, f64)) {
    let model = IsingModel::new(J);
    let n_sites = lattice.n_sites;
    let mut system = System::new(lattice.clone(), 1, 1.0, beta);
    system.recompute_energy(&model);
    let mut core = HeatBathCore::new();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    for _ in 0..THERM_SWEEPS {
        core.sweep(&mut system, &model, &mut rng);
    }
    let mut energies = Vec::with_capacity(MEAS_SWEEPS);
    let mut m2 = Vec::with_capacity(MEAS_SWEEPS);
    for _ in 0..MEAS_SWEEPS {
        core.sweep(&mut system, &model, &mut rng);
        energies.push(system.energy);
        let magnetization = system.spins.iter().sum::<f64>() / n_sites as f64;
        m2.push(magnetization * magnetization);
    }
    assert!(
        system.energy_error(&model).abs() < 1e-9,
        "heat-bath energy cache drifted by {}",
        system.energy_error(&model)
    );
    (
        mean_stderr_binned(&energies, BINSIZE),
        mean_stderr_binned(&m2, BINSIZE),
    )
}

fn check_discrete_heat_bath(n_sites: usize, label: &str) {
    let lattice = build_chain(n_sites, true);
    let n_seeds = zscore_seed_count(N_SEEDS);
    for (beta_index, beta) in [0.3, 0.8].iter().enumerate() {
        let (_, exact_e, _, exact_m2) = exact_ising_moments(&lattice, J, *beta);
        let mut energy_runs = Vec::with_capacity(n_seeds);
        let mut m2_runs = Vec::with_capacity(n_seeds);
        for seed in 0..n_seeds as u64 {
            let ((e_mean, e_stderr), (m2_mean, m2_stderr)) =
                run_heat_bath_ising(&lattice, *beta, 0x8A7B + 0x100 * beta_index as u64 + seed);
            energy_runs.push((e_mean, e_stderr));
            m2_runs.push((m2_mean, m2_stderr));
        }
        let max_e_z = assert_exact_z(&energy_runs, exact_e, &format!("{label} β={beta} ⟨E⟩"));
        let max_m2_z = assert_exact_z(&m2_runs, exact_m2, &format!("{label} β={beta} ⟨m²⟩"));
        eprintln!(
            "[heat-bath-exact] {label} β={beta}: ⟨E⟩ exact {exact_e:.4}, max|z| {max_e_z:.2}; \
             ⟨m²⟩ exact {exact_m2:.4}, max|z| {max_m2_z:.2}"
        );
    }
}

#[test]
fn heat_bath_discrete_n4_energy_and_m2_match_exact_enumeration() {
    check_discrete_heat_bath(4, "N=4 chain PBC");
}

#[test]
fn heat_bath_discrete_n8_energy_and_m2_match_exact_enumeration() {
    check_discrete_heat_bath(8, "N=8 chain PBC");
}

// ── Analytic references ───────────────────────────────────────────────

/// Modified Bessel function I_ν(x), ν ∈ {0,1,2}, via the power series
/// I_ν(x) = Σ_k (x/2)^(2k+ν) / (k! (k+ν)!). Converges in ~30 terms for
/// x ≤ 5 (the grid used here).
fn bessel_i(order: usize, x: f64) -> f64 {
    let half = 0.5 * x;
    let mut term = 1.0;
    for k in 1..=order {
        term *= half / k as f64;
    }
    let mut sum = term;
    let squared = half * half;
    let mut k = 0.0_f64;
    loop {
        k += 1.0;
        term *= squared / (k * (k + order as f64));
        sum += term;
        if term <= 1e-17 * sum {
            break;
        }
    }
    sum
}

/// Langevin function L(x) = coth(x) − 1/x = ⟨cosθ⟩ on S².
fn langevin(x: f64) -> f64 {
    1.0 / x.tanh() - 1.0 / x
}

fn o2_exact_cos(x: f64) -> f64 {
    bessel_i(1, x) / bessel_i(0, x)
}

fn o2_exact_cos2(x: f64) -> f64 {
    0.5 * (1.0 + bessel_i(2, x) / bessel_i(0, x))
}

fn o3_exact_cos(x: f64) -> f64 {
    langevin(x)
}

fn o3_exact_cos2(x: f64) -> f64 {
    1.0 - 2.0 * langevin(x) / x
}

#[test]
fn bessel_and_langevin_references_match_known_values_and_limits() {
    // Independent anchors: I_0(1) and I_1(1) (standard reference values).
    assert!((bessel_i(0, 1.0) - 1.2660658777520084).abs() < 1e-12);
    assert!((bessel_i(1, 1.0) - 0.565_159_103_992_485).abs() < 1e-12);

    // Series self-consistency via the recurrence I_0 − I_2 = 2 I_1 / x.
    for &x in &X_GRID {
        let recurrence = bessel_i(0, x) - bessel_i(2, x) - 2.0 * bessel_i(1, x) / x;
        assert!(recurrence.abs() < 1e-12 * bessel_i(0, x), "at x={x}");
        // Jensen sanity for the second-moment formulas.
        assert!(
            o2_exact_cos2(x) >= o2_exact_cos(x) * o2_exact_cos(x),
            "at x={x}"
        );
        assert!(
            o3_exact_cos2(x) >= o3_exact_cos(x) * o3_exact_cos(x),
            "at x={x}"
        );
    }

    // Small-x linear response: L(x) = x/3 − x³/45 + …, I₁/I₀(x) = x/2 − x³/16 + …
    // Tolerances are the analytic truncations: (0.2)³/45 ≈ 1.8e-4 and
    // (0.2)³/16 = 5.0e-4.
    assert!((langevin(0.2) - 0.2 / 3.0).abs() < 2.0e-4);
    assert!((o2_exact_cos(0.2) - 0.1).abs() < 6.0e-4);

    // Large-x saturation: L(x) = 1 − 1/x + O(e^-x) (residue 9.1e-5 at x=5);
    // I₁/I₀(x) = 1 − 1/(2x) + … (asymptotic residue ≈ 6.6e-3 at x=5).
    assert!((langevin(5.0) - 0.8).abs() < 2.0e-4);
    assert!((o2_exact_cos(5.0) - 0.9).abs() < 1.0e-2);
}

// ── Part B: continuous O(N) heat bath, finite-field conditional ───────

/// Two-site lattice with a single unit bond: the field of each site is its
/// partner (κ = βJ = x).
fn pair_lattice() -> CsrLattice {
    CsrLattice::from_edges(2, vec![Bond::new(0, 1, BondType::Generic, 1.0)])
}

/// One `ContinuousHeatBathCore` seed on the pair at inverse temperature
/// `beta = x`. Returns binned (⟨cosθ⟩, ⟨cos²θ⟩) and the worst |‖s‖²−1| seen.
fn run_continuous_pair<H>(model: &H, beta: f64, seed: u64) -> ((f64, f64), (f64, f64), f64)
where
    H: Hamiltonian + ContinuousHeatBathable,
{
    let spin_dim = model.spin_dim();
    let mut system = System::new(pair_lattice(), spin_dim, 0.0, beta);
    // Generic unit-norm start (spin_dim >= 2): site 0 = e_0, site 1 = e_1.
    system.spins[0] = 1.0;
    system.spins[spin_dim + 1] = 1.0;
    system.recompute_energy(model);
    let mut core = ContinuousHeatBathCore::new();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    for _ in 0..THERM_SWEEPS {
        core.sweep(&mut system, model, &mut rng);
    }
    let mut cosines = Vec::with_capacity(MEAS_SWEEPS);
    let mut max_norm_error = 0.0_f64;
    for _ in 0..MEAS_SWEEPS {
        core.sweep(&mut system, model, &mut rng);
        let (left, right) = system.spins.split_at(spin_dim);
        let cosine = left.iter().zip(right).map(|(a, b)| a * b).sum::<f64>();
        cosines.push(cosine);
        for spin in [left, right] {
            let norm_squared = spin.iter().map(|v| v * v).sum::<f64>();
            max_norm_error = max_norm_error.max((norm_squared - 1.0).abs());
        }
    }
    // Both samplers emit normalized trigonometric constructions, so the
    // unit-norm invariant only ever deviates by rounding (~1e-16).
    assert!(
        max_norm_error < 1e-12,
        "heat-bath spin norm deviated by {max_norm_error:.3e}"
    );
    assert!(
        system.energy_error(model).abs() < 1e-9,
        "continuous heat-bath energy cache drifted by {}",
        system.energy_error(model)
    );
    let cosines_squared: Vec<f64> = cosines.iter().map(|c| c * c).collect();
    (
        mean_stderr_binned(&cosines, BINSIZE),
        mean_stderr_binned(&cosines_squared, BINSIZE),
        max_norm_error,
    )
}

/// Shared per-model check over the x grid.
///
/// `small_limit`/`large_limit` are the linear-response and saturation forms
/// checked at the grid endpoints through a window = analytic truncation gap
/// plus 4 pooled seed SEs (the z-tests pin the MC mean to the exact curve;
/// the window then pins it to the curve's limit).
fn check_continuous_pair<H>(
    model: &H,
    label: &str,
    seed_base: u64,
    exact_cos: fn(f64) -> f64,
    exact_cos2: fn(f64) -> f64,
    small_limit: fn(f64) -> f64,
    large_limit: fn(f64) -> f64,
) where
    H: Hamiltonian + ContinuousHeatBathable,
{
    let n_seeds = zscore_seed_count(N_SEEDS);
    for (x_index, &x) in X_GRID.iter().enumerate() {
        let exact_c = exact_cos(x);
        let exact_c2 = exact_cos2(x);
        let mut cos_runs = Vec::with_capacity(n_seeds);
        let mut cos2_runs = Vec::with_capacity(n_seeds);
        for seed in 0..n_seeds as u64 {
            let ((c_mean, c_stderr), (c2_mean, c2_stderr), _) =
                run_continuous_pair(model, x, seed_base + 0x100 * x_index as u64 + seed);
            cos_runs.push((c_mean, c_stderr));
            cos2_runs.push((c2_mean, c2_stderr));
        }
        let max_c_z = assert_exact_z(&cos_runs, exact_c, &format!("{label} x={x} ⟨cosθ⟩"));
        let max_c2_z = assert_exact_z(&cos2_runs, exact_c2, &format!("{label} x={x} ⟨cos²θ⟩"));
        eprintln!(
            "[heat-bath-exact] {label} x={x}: ⟨cosθ⟩ exact {exact_c:.5}, max|z| {max_c_z:.2}; \
             ⟨cos²θ⟩ exact {exact_c2:.5}, max|z| {max_c2_z:.2}"
        );

        if x == X_GRID[0] || x == X_GRID[X_GRID.len() - 1] {
            let limit = if x == X_GRID[0] {
                small_limit(x)
            } else {
                large_limit(x)
            };
            let (pooled, pooled_stderr) = pooled_mean_stderr(&cos_runs);
            let window = (exact_c - limit).abs() + 4.0 * pooled_stderr;
            assert!(
                (pooled - limit).abs() <= window,
                "{label} x={x}: pooled ⟨cosθ⟩ {pooled:.5} outside the limit window \
                 [{:.5}, {:.5}] around {limit:.5}",
                limit - window,
                limit + window
            );
        }
    }
}

#[test]
fn continuous_heat_bath_o2_xy_pair_matches_von_mises_moments() {
    check_continuous_pair(
        &ONModel::<2>::new(J),
        "O(2) XY",
        0x0B52,
        o2_exact_cos,
        o2_exact_cos2,
        |x| x / 2.0,
        |x| 1.0 - 1.0 / (2.0 * x),
    );
}

#[test]
fn continuous_heat_bath_o3_pair_matches_langevin_moments() {
    check_continuous_pair(
        &ONModel::<3>::new(J),
        "O(3)",
        0x03B5,
        o3_exact_cos,
        o3_exact_cos2,
        |x| x / 3.0,
        |x| 1.0 - 1.0 / x,
    );
}
