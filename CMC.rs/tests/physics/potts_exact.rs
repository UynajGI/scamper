//! Exact-distribution validation for q-state Potts (q > 2).
//!
//! Closes the "Potts q>2 unvalidated" residue from MATURITY_ASSESSMENT.md
//! item 19 for every Potts-capable kernel: `HeatBathCore`, `SWCore` and
//! `WolffCore` (all three implement the required capability traits for
//! `PottsModel`), with `MetropolisCore` joining for cross-solver checks.
//!
//! - **Exact enumeration** (the repository standard): full q^N enumeration on
//!   the 2×2 PBC square (N=4) and the 8-site PBC chain (N=8), q ∈ {3, 4},
//!   β ∈ {0.3, 0.8}; observables ⟨E⟩, ⟨m²⟩ and the specific heat
//!   C = β²(⟨E²⟩ − ⟨E⟩²), each with multi-seed z-scores (|z| < 4 per seed).
//! - **Enumeration anchors** (independent of the Boltzmann loop): q=1 is
//!   exactly degenerate at E = −J Σ w at any β, and q=2 maps onto the
//!   independently written Ising enumeration with E_potts(J) =
//!   E_ising(J/2) − (J/2) Σ w and shifted second moments.
//! - **Analytic limits**: β=0 gives the exactly uniform state distribution
//!   (⟨E⟩ = −J Σ w / q and the uniform ⟨m²⟩ from unweighted enumeration);
//!   strong coupling (β=8) freezes the ground state (E = −J Σ w, m = 1).
//! - **Cross-solver**: Metropolis vs heat bath vs SW vs Wolff on the 8×8
//!   q=3 lattice near βc = ln(1+√3); all six solver pairs agree on ⟨E⟩ and
//!   ⟨m²⟩ within 4 pooled σ.
//! - **Critical-coupling anchor** (cheap directional check): q=4 order
//!   parameter rises monotonically across βc = ln 3 on the 8×8 lattice.
//!
//! Setting `SCUTTLE_ZSCORE_SEEDS=<n>` raises the seed count (nightly).

use super::common::zscore_seed_count;
use cmc_rs::{
    build_chain, build_square, Algorithm, CsrLattice, HeatBathCore, Initializable, Measurable,
    MetropolisCore, PottsModel, SWCore, System, WolffCore,
};
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

const J: f64 = 1.0;
const N_SEEDS: usize = 16;
const THERM_SWEEPS: usize = 1_000;
const MEAS_SWEEPS: usize = 24_000;
/// Bin length above the longest observable autocorrelation time of the four
/// kernels near β = 0.8 (cluster sweeps decorrelate in a few sweeps, but the
/// m² time series of the small lattices relaxes over tens of sweeps).
const BINSIZE: usize = 400;
const BETAS: [f64; 2] = [0.3, 0.8];

// ── Exact enumeration ──────────────────────────────────────────────────

fn enumerate_potts(n_sites: usize, q: usize) -> Vec<Vec<f64>> {
    let mut states = vec![vec![0.0; n_sites]; q.pow(n_sites as u32)];
    for (index, state) in states.iter_mut().enumerate() {
        let mut remainder = index;
        for component in state.iter_mut() {
            *component = (remainder % q) as f64;
            remainder /= q;
        }
    }
    states
}

fn potts_energy(spins: &[f64], lattice: &CsrLattice, coupling: f64) -> f64 {
    lattice
        .edges
        .iter()
        .map(|edge| {
            if spins[edge.source] == spins[edge.target] {
                -coupling * edge.weight
            } else {
                0.0
            }
        })
        .sum()
}

fn total_weight(lattice: &CsrLattice) -> f64 {
    lattice.edges.iter().map(|edge| edge.weight).sum()
}

/// Potts order parameter (the `Measurable` convention), valid for q >= 1:
/// m = (q·max_a n_a − N) / (N(q−1)), and m ≡ 1 for q = 1.
fn potts_magnetization(spins: &[f64], q: usize) -> f64 {
    let n = spins.len();
    if q == 1 || n == 0 {
        return 1.0;
    }
    let mut counts = vec![0usize; q];
    for &spin in spins {
        counts[spin as usize] += 1;
    }
    let largest = counts.into_iter().max().unwrap_or(0);
    (q as f64 * largest as f64 - n as f64) / (n as f64 * (q - 1) as f64)
}

/// Exact (⟨E⟩, ⟨E²⟩, ⟨m²⟩) by full enumeration at inverse temperature β.
fn exact_potts_moments(
    lattice: &CsrLattice,
    q: usize,
    coupling: f64,
    beta: f64,
) -> (f64, f64, f64) {
    let mut z = 0.0;
    let mut e = 0.0;
    let mut e2 = 0.0;
    let mut m2 = 0.0;
    for spins in enumerate_potts(lattice.n_sites, q) {
        let energy = potts_energy(&spins, lattice, coupling);
        let weight = (-beta * energy).exp();
        let magnetization = potts_magnetization(&spins, q);
        z += weight;
        e += weight * energy;
        e2 += weight * energy * energy;
        m2 += weight * magnetization * magnetization;
    }
    (e / z, e2 / z, m2 / z)
}

#[test]
fn potts_enumeration_anchors_match_q1_and_ising_limits() {
    let lattice = build_chain(8, true);

    // The library `Measurable` convention agrees with the test-local one on a
    // sample of q=3 states (the exact references below use the local copy).
    let model = PottsModel::new(J, 3);
    for (index, spins) in enumerate_potts(4, 3).iter().enumerate() {
        if index % 17 != 0 {
            continue;
        }
        super::common::assert_close(
            model.magnetization(spins),
            potts_magnetization(spins, 3),
            1e-15,
        );
    }

    // q=1 is exactly degenerate: every bond satisfied at any temperature.
    for &beta in &BETAS {
        let (energy, energy2, m2) = exact_potts_moments(&lattice, 1, J, beta);
        super::common::assert_close(energy, -J * total_weight(&lattice), 1e-12);
        super::common::assert_close(energy2, (J * total_weight(&lattice)).powi(2), 1e-9);
        super::common::assert_close(m2, 1.0, 1e-12);
    }

    // q=2 maps onto the independent Ising enumeration:
    // σ = 2s − 1 ⇒ δ(s,s') = (1 + σσ')/2 ⇒
    // E_potts(J) = E_ising(J/2) − (J/2) Σ w,  E_ising = E_potts + (J/2) Σ w.
    for &beta in &BETAS {
        let (_, ising_e, ising_e2, _) = super::common::exact_ising_moments(&lattice, J / 2.0, beta);
        let (potts_e, potts_e2, _) = exact_potts_moments(&lattice, 2, J, beta);
        let shift = J * total_weight(&lattice) / 2.0;
        super::common::assert_close(potts_e, ising_e - shift, 1e-10);
        // ⟨(E_ising − shift)²⟩ = ⟨E_ising²⟩ − 2 shift ⟨E_ising⟩ + shift².
        super::common::assert_close(
            potts_e2,
            ising_e2 - 2.0 * shift * ising_e + shift * shift,
            1e-9,
        );
    }
}

// ── Multi-seed z statistics (established pattern) ──────────────────────

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
/// per-seed z-scores; the one-sided Σz criterion is applied once per solver
/// over the pooled scores (see [`assert_no_one_sided_bias`]) so that the
/// ~2% tail of a single configuration cannot fail the suite by multiplicity.
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

/// Pooled one-sided-bias gate over every z-score a solver produced.
fn assert_no_one_sided_bias(z_scores: &[f64], label: &str) {
    let sum_z = z_scores.iter().sum::<f64>();
    let sum_floor = -2.0 * (z_scores.len() as f64).sqrt();
    assert!(
        sum_z > sum_floor,
        "{label}: pooled Σz = {sum_z:.2} over {} scores should be > {sum_floor:.2}",
        z_scores.len()
    );
}

/// Specific heat C = β²(⟨E²⟩ − ⟨E⟩²) per system with a jackknife-over-bins
/// standard error (the fluctuation estimator is a nonlinear functional, so
/// its error bar needs jackknife rather than a plain bin mean).
fn specific_heat_binned(energies: &[f64], binsize: usize, beta: f64) -> (f64, f64) {
    let n_bins = energies.len() / binsize;
    assert!(n_bins >= 2, "need at least 2 bins, got {n_bins}");
    let bin_energy: Vec<f64> = (0..n_bins)
        .map(|b| {
            let start = b * binsize;
            energies[start..start + binsize].iter().sum::<f64>() / binsize as f64
        })
        .collect();
    let bin_energy_squared: Vec<f64> = (0..n_bins)
        .map(|b| {
            let start = b * binsize;
            energies[start..start + binsize]
                .iter()
                .map(|e| e * e)
                .sum::<f64>()
                / binsize as f64
        })
        .collect();
    let heat_of =
        |energy: f64, energy_squared: f64| beta * beta * (energy_squared - energy * energy);
    let mean = |values: &[f64]| values.iter().sum::<f64>() / values.len() as f64;
    let full = heat_of(mean(&bin_energy), mean(&bin_energy_squared));
    let mut jackknife_sum = 0.0;
    let mut jackknife_sum_squared = 0.0;
    for dropped in 0..n_bins {
        let kept_energy: Vec<f64> = bin_energy
            .iter()
            .enumerate()
            .filter(|(index, _)| *index != dropped)
            .map(|(_, &value)| value)
            .collect();
        let kept_squared: Vec<f64> = bin_energy_squared
            .iter()
            .enumerate()
            .filter(|(index, _)| *index != dropped)
            .map(|(_, &value)| value)
            .collect();
        let estimate = heat_of(mean(&kept_energy), mean(&kept_squared));
        jackknife_sum += estimate;
        jackknife_sum_squared += estimate * estimate;
    }
    let count = n_bins as f64;
    let variance =
        (jackknife_sum_squared - jackknife_sum * jackknife_sum / count) * (count - 1.0) / count;
    (full, variance.sqrt())
}

// ── Solver harness ─────────────────────────────────────────────────────

/// The four Potts-capable update kernels behind one interface.
enum PottsCore {
    Metropolis(MetropolisCore),
    HeatBath(HeatBathCore),
    SwendsenWang(SWCore),
    Wolff(WolffCore),
}

impl PottsCore {
    fn new(kind: &str) -> Self {
        match kind {
            "metropolis" => Self::Metropolis(MetropolisCore::new()),
            "heat-bath" => Self::HeatBath(HeatBathCore::new()),
            "sw" => Self::SwendsenWang(SWCore::new()),
            "wolff" => Self::Wolff(WolffCore::new()),
            _ => panic!("unknown Potts core `{kind}`"),
        }
    }

    fn sweep(&mut self, system: &mut System, model: &PottsModel, rng: &mut impl Rng) {
        match self {
            Self::Metropolis(core) => core.sweep(system, model, rng),
            Self::HeatBath(core) => core.sweep(system, model, rng),
            Self::SwendsenWang(core) => core.sweep(system, model, rng),
            Self::Wolff(core) => core.sweep(system, model, rng),
        }
    }
}

/// One chain of a named kernel; returns binned ⟨E⟩, ⟨m²⟩ and C.
fn run_potts_chain(
    kind: &str,
    lattice: &CsrLattice,
    q: usize,
    beta: f64,
    seed: u64,
) -> ((f64, f64), (f64, f64), (f64, f64)) {
    let model = PottsModel::new(J, q);
    let n_sites = lattice.n_sites;
    let mut system = System::new(lattice.clone(), 1, 0.0, beta);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    for site in 0..n_sites {
        let spin = model.random_spin(&mut rng);
        system.spin_at_mut(site, 1).copy_from_slice(&spin);
    }
    system.recompute_energy(&model);
    let mut core = PottsCore::new(kind);
    for _ in 0..THERM_SWEEPS {
        core.sweep(&mut system, &model, &mut rng);
    }
    let mut energies = Vec::with_capacity(MEAS_SWEEPS);
    let mut m2 = Vec::with_capacity(MEAS_SWEEPS);
    for _ in 0..MEAS_SWEEPS {
        core.sweep(&mut system, &model, &mut rng);
        energies.push(system.energy);
        let magnetization = model.magnetization(&system.spins);
        m2.push(magnetization * magnetization);
    }
    assert!(
        system.energy_error(&model).abs() < 1e-9,
        "{kind} q={q} energy cache drifted by {}",
        system.energy_error(&model)
    );
    (
        mean_stderr_binned(&energies, BINSIZE),
        mean_stderr_binned(&m2, BINSIZE),
        specific_heat_binned(&energies, BINSIZE, beta),
    )
}

/// Full exact-enumeration check of one kernel on both lattices, q ∈ {3,4}.
fn check_kernel_against_exact(kind: &str, seed_base: u64) {
    let n_seeds = zscore_seed_count(N_SEEDS);
    let lattices = [
        (build_square(2, 2, true), "2x2 square PBC"),
        (build_chain(8, true), "N=8 chain PBC"),
    ];
    let mut pooled_z = Vec::new();
    for (lattice, lattice_label) in &lattices {
        for q in [3usize, 4] {
            for (beta_index, &beta) in BETAS.iter().enumerate() {
                let (exact_e, exact_e2, exact_m2) = exact_potts_moments(lattice, q, J, beta);
                let exact_c = beta * beta * (exact_e2 - exact_e * exact_e);
                let mut energy_runs = Vec::with_capacity(n_seeds);
                let mut m2_runs = Vec::with_capacity(n_seeds);
                let mut heat_runs = Vec::with_capacity(n_seeds);
                for seed in 0..n_seeds as u64 {
                    let (energy, m2, heat) = run_potts_chain(
                        kind,
                        lattice,
                        q,
                        beta,
                        seed_base + 0x1000 * beta_index as u64 + 0x10 * q as u64 + seed,
                    );
                    energy_runs.push(energy);
                    m2_runs.push(m2);
                    heat_runs.push(heat);
                }
                let label = |observable: &str| {
                    format!("{kind} {lattice_label} q={q} β={beta} {observable}")
                };
                let energy_z = assert_per_seed_z(&energy_runs, exact_e, &label("⟨E⟩"));
                let m2_z = assert_per_seed_z(&m2_runs, exact_m2, &label("⟨m²⟩"));
                let heat_z = assert_per_seed_z(&heat_runs, exact_c, &label("C"));
                eprintln!(
                    "[potts-exact] {} q={q} β={beta}: ⟨E⟩ exact {exact_e:+.4}, max|z| {:.2}; \
                     ⟨m²⟩ exact {exact_m2:.4}, max|z| {:.2}; \
                     C exact {exact_c:.4}, max|z| {:.2}",
                    lattice_label,
                    energy_z.iter().fold(0.0_f64, |acc, z| acc.max(z.abs())),
                    m2_z.iter().fold(0.0_f64, |acc, z| acc.max(z.abs())),
                    heat_z.iter().fold(0.0_f64, |acc, z| acc.max(z.abs())),
                );
                pooled_z.extend(energy_z);
                pooled_z.extend(m2_z);
                pooled_z.extend(heat_z);
            }
        }
    }
    assert_no_one_sided_bias(&pooled_z, kind);
}

#[test]
fn potts_heat_bath_q3_q4_match_exact_enumeration() {
    check_kernel_against_exact("heat-bath", 0x0B0);
}

#[test]
fn potts_swendsen_wang_q3_q4_match_exact_enumeration() {
    check_kernel_against_exact("sw", 0x050);
}

#[test]
fn potts_wolff_q3_q4_match_exact_enumeration() {
    check_kernel_against_exact("wolff", 0x0F0);
}

// ── Analytic limits ────────────────────────────────────────────────────

#[test]
fn potts_beta_zero_limit_is_exactly_uniform() {
    // At β = 0 every kernel must sample the uniform state distribution:
    // ⟨E⟩ = −J Σ w / q exactly, and ⟨m²⟩ equals the unweighted enumeration
    // average (which the Boltzmann loop at β = 0 computes with all weights 1).
    let lattices = [
        (build_square(2, 2, true), "2x2 square PBC"),
        (build_chain(8, true), "N=8 chain PBC"),
    ];
    let mut pooled_z = Vec::new();
    for (lattice, lattice_label) in &lattices {
        for q in [3usize, 4] {
            let (exact_e, _, exact_m2) = exact_potts_moments(lattice, q, J, 0.0);
            super::common::assert_close(exact_e, -J * total_weight(lattice) / q as f64, 1e-12);
            for kind in ["heat-bath", "sw", "wolff"] {
                let n_seeds = zscore_seed_count(8);
                let mut energy_runs = Vec::with_capacity(n_seeds);
                let mut m2_runs = Vec::with_capacity(n_seeds);
                for seed in 0..n_seeds as u64 {
                    let (energy, m2, _) =
                        run_potts_chain(kind, lattice, q, 0.0, seed + 0xB00 * q as u64);
                    energy_runs.push(energy);
                    m2_runs.push(m2);
                }
                pooled_z.extend(assert_per_seed_z(
                    &energy_runs,
                    exact_e,
                    &format!("{kind} {lattice_label} q={q} β=0 ⟨E⟩"),
                ));
                pooled_z.extend(assert_per_seed_z(
                    &m2_runs,
                    exact_m2,
                    &format!("{kind} {lattice_label} q={q} β=0 ⟨m²⟩"),
                ));
            }
        }
    }
    assert_no_one_sided_bias(&pooled_z, "β=0 uniform limit");
}

#[test]
fn potts_strong_coupling_limit_freezes_the_ground_state() {
    // β = 8: single-site excitations cost ΔE ≥ 2J (a 2×2 site has degree 4
    // with PBC double edges ⇒ ΔE = 4J), so every kernel sits in the ordered
    // ground state to within e^{-β ΔE} ~ 1e-14 corrections.
    let lattice = build_square(2, 2, true);
    let ground_energy = -J * total_weight(&lattice);
    let n_seeds = zscore_seed_count(4);
    for kind in ["heat-bath", "sw", "wolff"] {
        for q in [3usize, 4] {
            for seed in 0..n_seeds as u64 {
                let ((energy_mean, _), (m2_mean, _), _) =
                    run_potts_chain(kind, &lattice, q, 8.0, seed + 0x5C0 * q as u64);
                assert!(
                    (energy_mean - ground_energy).abs() < 1e-6,
                    "{kind} q={q}: ⟨E⟩ = {energy_mean} at β=8, ground state is {ground_energy}"
                );
                assert!(
                    m2_mean > 1.0 - 1e-6,
                    "{kind} q={q}: ⟨m²⟩ = {m2_mean} at β=8, expected a frozen ordered state"
                );
            }
        }
    }
}

// ── Cross-solver agreement (no exact reference at 8×8 q=3) ─────────────

/// Pooled observable of one solver: (kind, ⟨E⟩ (mean, stderr), ⟨m²⟩ (mean, stderr)).
struct SolverMoments {
    kind: &'static str,
    energy: (f64, f64),
    m2: (f64, f64),
}

#[test]
fn potts_cross_solver_agreement_8x8_q3_near_criticality() {
    let lattice = build_square(8, 8, true);
    let q = 3usize;
    // βc(q=3) = ln(1+√3) ≈ 1.0051 — the hardest regime to agree in.
    let beta = 1.0;
    let n_seeds = zscore_seed_count(8);
    let kinds = ["metropolis", "heat-bath", "sw", "wolff"];
    let seed_bases = [0x0E7A, 0x0E7B, 0x0E7C, 0x0E7D];
    let mut means: Vec<SolverMoments> = Vec::new();
    for (index, kind) in kinds.iter().enumerate() {
        let mut energy_runs = Vec::with_capacity(n_seeds);
        let mut m2_runs = Vec::with_capacity(n_seeds);
        for seed in 0..n_seeds as u64 {
            // Longer chains and wider bins: local Metropolis near criticality
            // has the longest autocorrelation of the four kernels.
            let model = PottsModel::new(J, q);
            let mut system = System::new(lattice.clone(), 1, 0.0, beta);
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed_bases[index] + seed);
            for site in 0..system.n_sites() {
                let spin = model.random_spin(&mut rng);
                system.spin_at_mut(site, 1).copy_from_slice(&spin);
            }
            system.recompute_energy(&model);
            let mut core = PottsCore::new(kind);
            for _ in 0..2_000 {
                core.sweep(&mut system, &model, &mut rng);
            }
            let mut energies = Vec::with_capacity(20_000);
            let mut m2 = Vec::with_capacity(20_000);
            for _ in 0..20_000 {
                core.sweep(&mut system, &model, &mut rng);
                energies.push(system.energy);
                let magnetization = model.magnetization(&system.spins);
                m2.push(magnetization * magnetization);
            }
            energy_runs.push(mean_stderr_binned(&energies, 1_000));
            m2_runs.push(mean_stderr_binned(&m2, 1_000));
        }
        let pool = |runs: &[(f64, f64)]| {
            let count = runs.len() as f64;
            let mean = runs.iter().map(|(m, _)| m).sum::<f64>() / count;
            let variance = runs.iter().map(|(_, se)| se * se).sum::<f64>() / (count * count);
            (mean, variance.sqrt())
        };
        means.push(SolverMoments {
            kind,
            energy: pool(&energy_runs),
            m2: pool(&m2_runs),
        });
    }
    for (index, left) in means.iter().enumerate() {
        for right in means.iter().skip(index + 1) {
            let (left_kind, right_kind) = (left.kind, right.kind);
            let (left_e, left_e_se) = left.energy;
            let (right_e, right_e_se) = right.energy;
            let (left_m2, left_m2_se) = left.m2;
            let (right_m2, right_m2_se) = right.m2;
            let energy_z =
                (left_e - right_e) / (left_e_se * left_e_se + right_e_se * right_e_se).sqrt();
            let m2_z =
                (left_m2 - right_m2) / (left_m2_se * left_m2_se + right_m2_se * right_m2_se).sqrt();
            eprintln!(
                "[potts-cross] {left_kind} vs {right_kind}: ⟨E⟩ z = {energy_z:+.2}, ⟨m²⟩ z = {m2_z:+.2}"
            );
            assert!(
                energy_z.abs() < 4.0,
                "{left_kind} vs {right_kind}: ⟨E⟩ {left_e:.4} vs {right_e:.4}, z = {energy_z:.2}"
            );
            assert!(
                m2_z.abs() < 4.0,
                "{left_kind} vs {right_kind}: ⟨m²⟩ {left_m2:.4} vs {right_m2:.4}, z = {m2_z:.2}"
            );
        }
    }
}

#[test]
fn potts_q4_order_parameter_rises_across_the_square_critical_point() {
    // Cheap directional anchor on the exact square-lattice critical coupling
    // βc = ln(1+√q): for q = 4, βc = ln 3 ≈ 1.0986. The order parameter on
    // the 8×8 lattice must rise monotonically across it.
    let lattice = build_square(8, 8, true);
    let q = 4usize;
    let beta_c = 3.0_f64.ln();
    let n_seeds = zscore_seed_count(4);
    let mut m2_by_beta = Vec::new();
    for &beta in &[0.8 * beta_c, beta_c, 1.25 * beta_c] {
        let mut m2_means = Vec::with_capacity(n_seeds);
        for seed in 0..n_seeds as u64 {
            let (_, (m2_mean, _), _) = run_potts_chain("wolff", &lattice, q, beta, seed + 0xC40);
            m2_means.push(m2_mean);
        }
        m2_by_beta.push(m2_means.iter().sum::<f64>() / m2_means.len() as f64);
    }
    eprintln!(
        "[potts-beta-c] q=4 8x8 ⟨m²⟩: {:+.3} (0.8βc), {:+.3} (βc), {:+.3} (1.25βc)",
        m2_by_beta[0], m2_by_beta[1], m2_by_beta[2]
    );
    assert!(
        m2_by_beta[0] < m2_by_beta[1] && m2_by_beta[1] < m2_by_beta[2],
        "q=4 order parameter must rise across βc = ln 3: {:?}",
        m2_by_beta
    );
    assert!(
        m2_by_beta[2] - m2_by_beta[0] > 0.2,
        "q=4 ordering across βc must be unambiguous: {:?}",
        m2_by_beta
    );
}
