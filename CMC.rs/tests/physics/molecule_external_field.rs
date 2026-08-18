//! One-body external-field validation for `MolecularMetropolisCore`.
//!
//! Closes the VALIDATION.md molecule residue: "the API has no one-body
//! external term (blocks the true Langevin dipole case)". The new
//! `DipolarExternalField` couples a molecular dipole to a uniform field E;
//! this file samples the *free rotor* (one rigid dumbbell, no pair
//! interactions) through the real kernel (translation + rotation moves) and
//! checks the orientation distribution against the exact one-rotor Boltzmann
//! answer, which is the Langevin-dipole problem:
//!
//! - **2D** rotor (field in plane): von Mises distribution
//!   ⟨cosθ⟩ = I₁(x)/I₀(x), ⟨cos²θ⟩ = (1 + I₂(x)/I₀(x))/2, x = βpE.
//! - **3D** rotor: the Langevin function itself,
//!   ⟨cosθ⟩ = L(x) = coth x − 1/x, ⟨cos²θ⟩ = 1 − 2L(x)/x.
//!
//! The in-code Bessel/Langevin references are anchored against known values,
//! recurrence identities and small/large-x limits separately. A
//! machine-precision identity additionally pins the kernel's
//! `external_field_energy` to −E·μ for every sampled configuration, and the
//! constructor's neutrality/charge-table rejections are tested loudly.
//!
//! Statistical standard: the established multi-seed z pattern — per-seed
//! |z| < 4, |z̄| < 2, Σz > −2√n; seed count scales via `SCUTTLE_ZSCORE_SEEDS`.

use super::common::zscore_seed_count;
use cmc_rs::{
    DipolarExternalField, MolecularMetropolisCore, MoleculeTopology, OrthorhombicCell,
    PairPotential, ParticleAlgorithm, ParticleConfiguration, ParticleError, ParticleSystem,
    SimulationCell, SimulationPhase,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

const BETA: f64 = 1.0;
/// Dumbbell bond length; the two ±CHARGE atoms give dipole p = CHARGE·BOND.
const BOND: f64 = 1.0;
const CHARGE: f64 = 0.5;
const DIPOLE: f64 = CHARGE * BOND;
const CELL_SIDE: f64 = 8.0;
const N_SEEDS: usize = 8;
const THERM_SWEEPS: usize = 3_000;
const MEAS_SWEEPS: usize = 24_000;
/// Bin length chosen well above the orientation autocorrelation time at the
/// weakest field (≈60 sweeps in 3D), so binned stderrs are unbiased.
const BINSIZE: usize = 600;
/// Field strengths x = βpE spanning linear response to near-saturation.
const X_GRID: [f64; 4] = [0.5, 1.5, 3.0, 5.0];

// ── Free rotor: no pair interactions ───────────────────────────────────

/// Interaction-free potential (the rotor is a single molecule in vacuum).
#[derive(Debug, Clone, Copy)]
struct IdealRotor;

impl PairPotential for IdealRotor {
    fn cutoff_squared(&self) -> f64 {
        1.0
    }

    fn energy(&self, _species_i: u16, _species_j: u16, _distance_squared: f64) -> f64 {
        0.0
    }
}

// ── Statistics helpers (established multi-seed z pattern) ──────────────

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

/// Per-config gates: every seed |z| < 4 and |mean z| < 2. Returns the
/// per-seed z-scores; the one-sided Σz criterion is pooled over the whole
/// x-grid per dimension (see [`assert_no_one_sided_bias`]) so the ~2% tail
/// of a single configuration cannot fail the test by multiplicity.
fn assert_per_seed_z(results: &[(f64, f64)], exact: f64, label: &str) -> Vec<f64> {
    let z_scores: Vec<f64> = results
        .iter()
        .map(|(mean, stderr)| (mean - exact) / stderr.max(1e-10))
        .collect();
    let max_abs_z = z_scores.iter().fold(0.0_f64, |acc, z| acc.max(z.abs()));
    let mean_z = z_scores.iter().sum::<f64>() / z_scores.len() as f64;
    assert!(
        max_abs_z < 4.0,
        "{label}: max |z| = {max_abs_z:.2} (exact = {exact:.6})"
    );
    assert!(
        mean_z.abs() < 2.0,
        "{label}: mean z = {mean_z:.2} (exact = {exact:.6})"
    );
    z_scores
}

/// Pooled one-sided-bias gate over every z-score one test produced.
fn assert_no_one_sided_bias(z_scores: &[f64], label: &str) {
    let sum_z = z_scores.iter().sum::<f64>();
    let sum_floor = -2.0 * (z_scores.len() as f64).sqrt();
    assert!(
        sum_z > sum_floor,
        "{label}: pooled Σz = {sum_z:.2} over {} scores should be > {sum_floor:.2}",
        z_scores.len()
    );
}

// ── Exact references ───────────────────────────────────────────────────

/// Modified Bessel function I_ν(x), ν ∈ {0,1,2}, via the power series.
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

/// Langevin function L(x) = coth(x) − 1/x.
fn langevin(x: f64) -> f64 {
    1.0 / x.tanh() - 1.0 / x
}

fn rotor_2d_cos(x: f64) -> f64 {
    bessel_i(1, x) / bessel_i(0, x)
}

fn rotor_2d_cos2(x: f64) -> f64 {
    0.5 * (1.0 + bessel_i(2, x) / bessel_i(0, x))
}

fn rotor_3d_cos(x: f64) -> f64 {
    langevin(x)
}

fn rotor_3d_cos2(x: f64) -> f64 {
    1.0 - 2.0 * langevin(x) / x
}

#[test]
fn langevin_and_bessel_rotor_references_match_known_values_and_limits() {
    // Independent anchors: I_0(1), I_1(1) and L(1) (standard values).
    assert!((bessel_i(0, 1.0) - 1.2660658777520084).abs() < 1e-12);
    assert!((bessel_i(1, 1.0) - 0.565_159_103_992_485).abs() < 1e-12);
    assert!((langevin(1.0) - 0.31303528549933135).abs() < 1e-12);

    // Recurrence I_0 − I_2 = 2 I_1 / x ties the three orders together.
    for &x in &X_GRID {
        let recurrence = bessel_i(0, x) - bessel_i(2, x) - 2.0 * bessel_i(1, x) / x;
        assert!(recurrence.abs() < 1e-12 * bessel_i(0, x), "at x={x}");
        assert!(rotor_2d_cos2(x) >= rotor_2d_cos(x).powi(2), "at x={x}");
        assert!(rotor_3d_cos2(x) >= rotor_3d_cos(x).powi(2), "at x={x}");
    }

    // Small-x linear response: L(x) = x/3 − x³/45 + …, I₁/I₀(x) = x/2 − x³/16 + …
    // (tolerances are the analytic truncations at x = 0.5).
    assert!((langevin(0.5) - 0.5 / 3.0).abs() < 3.0e-3);
    assert!((rotor_2d_cos(0.5) - 0.25).abs() < 8.0e-3);

    // Large-x saturation: L(x) = 1 − 1/x + O(e^-x), I₁/I₀(x) = 1 − 1/(2x) + …
    assert!((langevin(5.0) - 0.8).abs() < 2.0e-4);
    assert!((rotor_2d_cos(5.0) - 0.9).abs() < 1.0e-2);
}

// ── Rotor chain through the real kernel ────────────────────────────────

/// One free-rotor seed at field strength x; returns binned (⟨cosθ⟩, ⟨cos²θ⟩).
fn run_free_rotor<const D: usize>(x: f64, seed: u64) -> ((f64, f64), (f64, f64)) {
    let cell = OrthorhombicCell::new([CELL_SIDE; D]).expect("valid cell lengths");
    // Start the bond along +y so the chain is not initialized aligned with
    // the +x field; the COM sits at the cell center.
    let mut axis = [0.0; D];
    axis[D - 1] = 1.0;
    let center = [CELL_SIDE / 2.0; D];
    let half = 0.5 * BOND;
    let mut positions = vec![[0.0; D]; 2];
    for component in 0..D {
        positions[0][component] = center[component] - half * axis[component];
        positions[1][component] = center[component] + half * axis[component];
    }
    let configuration =
        ParticleConfiguration::new(positions, vec![0, 0], cell).expect("valid configuration");
    let potential = IdealRotor;
    let mut system =
        ParticleSystem::new(configuration, &potential, BETA).expect("valid rotor system");
    let topology = MoleculeTopology::new(2, vec![vec![0, 1]]).expect("dumbbell topology");
    let mut field_vector = [0.0; D];
    field_vector[0] = x / (BETA * DIPOLE);
    let field = DipolarExternalField::new(field_vector, vec![CHARGE, -CHARGE])
        .expect("neutral dumbbell field is valid");
    let mut kernel = MolecularMetropolisCore::new(topology, 0.35, 0.4)
        .expect("valid kernel scales")
        .with_external_field(field)
        .expect("neutral molecule accepts the field");
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    for _ in 0..THERM_SWEEPS {
        kernel.sweep_with_phase(
            &mut system,
            &potential,
            &mut rng,
            SimulationPhase::Thermalization,
        );
    }
    let mut cosines = Vec::with_capacity(MEAS_SWEEPS);
    let mut bond_error = 0.0_f64;
    let mut field_energy_error = 0.0_f64;
    for _ in 0..MEAS_SWEEPS {
        kernel.sweep_with_phase(
            &mut system,
            &potential,
            &mut rng,
            SimulationPhase::Measurement,
        );
        let positions = system.configuration().positions();
        // μ = Σ_a q_a (r_a − r_0) = −CHARGE · bond ⇒ μ̂ = −bond/|bond|.
        let bond = cell.displacement(&positions[0], &positions[1]);
        let norm = bond.iter().map(|v| v * v).sum::<f64>().sqrt();
        bond_error = bond_error.max((norm - BOND).abs());
        cosines.push(-bond[0] / norm);
        // Machine-precision identity: E_field = −E·μ = −x/β · cosθ (β = 1).
        let field_energy = kernel
            .external_field_energy(system.configuration())
            .expect("field is attached");
        field_energy_error =
            field_energy_error.max((field_energy + x * cosines[cosines.len() - 1]).abs());
    }
    assert!(
        bond_error < 1e-12,
        "rigid bond length drifted by {bond_error:.3e}"
    );
    assert!(
        field_energy_error < 1e-12,
        "external-field energy deviated from −E·μ by {field_energy_error:.3e}"
    );
    let cosines_squared: Vec<f64> = cosines.iter().map(|c| c * c).collect();
    (
        mean_stderr_binned(&cosines, BINSIZE),
        mean_stderr_binned(&cosines_squared, BINSIZE),
    )
}

fn check_free_rotor<const D: usize>(
    label: &str,
    seed_base: u64,
    exact_cos: fn(f64) -> f64,
    exact_cos2: fn(f64) -> f64,
) {
    let n_seeds = zscore_seed_count(N_SEEDS);
    let mut pooled_z = Vec::new();
    for (x_index, &x) in X_GRID.iter().enumerate() {
        let mut cos_runs = Vec::with_capacity(n_seeds);
        let mut cos2_runs = Vec::with_capacity(n_seeds);
        for seed in 0..n_seeds as u64 {
            let ((c_mean, c_stderr), (c2_mean, c2_stderr)) =
                run_free_rotor::<D>(x, seed_base + 0x100 * x_index as u64 + seed);
            cos_runs.push((c_mean, c_stderr));
            cos2_runs.push((c2_mean, c2_stderr));
        }
        let cos_z = assert_per_seed_z(&cos_runs, exact_cos(x), &format!("{label} x={x} ⟨cosθ⟩"));
        let cos2_z =
            assert_per_seed_z(&cos2_runs, exact_cos2(x), &format!("{label} x={x} ⟨cos²θ⟩"));
        eprintln!(
            "[molecule-field] {label} x={x}: ⟨cosθ⟩ exact {:+.5}, max|z| {:.2}; \
             ⟨cos²θ⟩ exact {:+.5}, max|z| {:.2}",
            exact_cos(x),
            cos_z.iter().fold(0.0_f64, |acc, z| acc.max(z.abs())),
            exact_cos2(x),
            cos2_z.iter().fold(0.0_f64, |acc, z| acc.max(z.abs())),
        );
        pooled_z.extend(cos_z);
        pooled_z.extend(cos2_z);
    }
    assert_no_one_sided_bias(&pooled_z, label);
}

#[test]
fn external_field_free_rotor_2d_matches_von_mises_moments() {
    check_free_rotor::<2>("2D rotor", 0xD1F_0001, rotor_2d_cos, rotor_2d_cos2);
}

#[test]
fn external_field_free_rotor_3d_matches_langevin_moments() {
    check_free_rotor::<3>("3D rotor", 0xD1F_0003, rotor_3d_cos, rotor_3d_cos2);
}

// ── Loud input rejection (criterion G) ─────────────────────────────────

#[test]
fn external_field_rejects_invalid_charge_tables_loudly() {
    let topology = MoleculeTopology::new(2, vec![vec![0, 1]]).expect("dumbbell topology");

    // Non-neutral molecule: net charge couples to the wrapped position and
    // would silently break periodicity/detailed balance.
    let charged = DipolarExternalField::new([1.0, 0.0], vec![1.0, -0.5]).unwrap();
    let error = MolecularMetropolisCore::new(topology.clone(), 0.3, 0.3)
        .unwrap()
        .with_external_field(charged)
        .unwrap_err();
    assert!(
        matches!(error, ParticleError::InvalidPotential(_)),
        "expected InvalidPotential, got {error:?}"
    );

    // Charge table must cover exactly the topology's particles.
    let short_table = DipolarExternalField::new([1.0, 0.0], vec![0.5]).unwrap();
    let error = MolecularMetropolisCore::new(topology.clone(), 0.3, 0.3)
        .unwrap()
        .with_external_field(short_table)
        .unwrap_err();
    assert!(
        matches!(error, ParticleError::InvalidPotential(_)),
        "expected InvalidPotential, got {error:?}"
    );

    // Non-finite field and charge values are rejected at field construction.
    assert!(DipolarExternalField::new([f64::NAN, 0.0], vec![0.5, -0.5]).is_err());
    assert!(DipolarExternalField::new([1.0, 0.0], vec![f64::INFINITY, -1.0]).is_err());
    assert!(DipolarExternalField::new([1.0, 0.0], Vec::new()).is_err());

    // A neutral dumbbell with a full charge table is accepted, and the
    // kernel without a field reports no one-body energy.
    let valid = DipolarExternalField::new([2.0, 0.0], vec![0.5, -0.5]).unwrap();
    let kernel = MolecularMetropolisCore::new(topology, 0.3, 0.3)
        .unwrap()
        .with_external_field(valid)
        .unwrap();
    assert!(kernel.external_field().is_some());

    let bare = MolecularMetropolisCore::<2>::new(
        MoleculeTopology::new(2, vec![vec![0, 1]]).unwrap(),
        0.3,
        0.3,
    )
    .unwrap();
    assert!(bare.external_field().is_none());
}
