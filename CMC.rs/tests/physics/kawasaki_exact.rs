//! Kawasaki (`KawasakiCore`) equilibrium vs exact fixed-sector enumeration.
//!
//! Closes the VALIDATION.md gap "Kawasaki dynamics — NOT validated:
//! Quantitative equilibrium distribution".
//!
//! Kawasaki exchanges the spins of a uniformly chosen bond, so the total
//! magnetization M is conserved exactly and the equilibrium distribution is
//! the canonical ensemble restricted to the fixed-M exchange class reachable
//! from the initial state — NOT the full canonical ensemble.
//!
//! Methodology (the repo's explicit-Markov-chain approach):
//! 1. Enumerate the exchange graph: BFS from the initial state over single
//!    unequal-spin bond swaps. This is exactly the proposal graph of
//!    `KawasakiCore` (uniform edge choice + symmetric swap), so
//!    Metropolis-Hastings acceptance leaves the Boltzmann weights e^(−βE)
//!    stationary on the reachable set.
//! 2. Exact reference: ⟨E⟩, ⟨E²⟩ (and the fixed ⟨m²⟩ = (M/N)²) of the
//!    Boltzmann weights restricted to the BFS-reachable states.
//! 3. Run `KawasakiCore` from the same initial state for
//!    `zscore_seed_count(8)` seeds, asserting exact M conservation after
//!    every sweep and energy-cache integrity, then z-test the per-seed
//!    binned ⟨E⟩ and ⟨E²⟩ means against the sector-exact values with the
//!    established multi-seed criteria: max |z| < 4, |z̄| < 2, and the
//!    scale-invariant one-sided bound Σz > −2√n.
//!
//! Sector-connectivity finding (asserted below): on PBC rings the
//! BFS-reachable set equals the FULL fixed-M sector — N=4 M=0: 6/6 states,
//! N=8 M=+2: 56/56. The hypothesized sublattice-imbalance invariant does
//! NOT survive on a periodic chain: exchanges let an up-spin wrap around
//! the ring and pass other up-spins, so the sublattice up-count is not
//! conserved. The reference below always uses the BFS-reachable set, so the
//! test remains correct on any lattice/sector where extra invariants do
//! exist.
//!
//! The sector reference is also contrasted with the full-canonical ⟨E⟩
//! (from `common::exact_ising_moments`): the two differ by far more than
//! the MC error bars, i.e. testing Kawasaki against the unrestricted
//! canonical answer would fail — these tests pin the sector-restricted
//! distribution (M=0 on N=4 and the M≠0 sector M=+2 on N=8).

use super::common::{direct_ising_energy, enumerate_ising, exact_ising_moments, zscore_seed_count};
use cmc_rs::{build_chain, Algorithm, CsrLattice, IsingModel, KawasakiCore, System};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use std::collections::VecDeque;

const J: f64 = 1.0;
const N_SEEDS: usize = 8;
const THERM_SWEEPS: usize = 1000;
const MEAS_SWEEPS: usize = 12_000;
const BINSIZE: usize = 200;

/// Distinct undirected exchange pairs (self-loops excluded) of a lattice.
fn exchange_pairs(lattice: &CsrLattice) -> Vec<(usize, usize)> {
    let mut pairs: Vec<(usize, usize)> = lattice
        .edges
        .iter()
        .filter(|edge| edge.source != edge.target)
        .map(|edge| (edge.source, edge.target))
        .collect();
    pairs.sort_unstable();
    pairs.dedup();
    pairs
}

/// States reachable from `initial` under single unequal-spin swaps along
/// lattice edges (the Kawasaki proposal graph). A state is a spin-up bitmask
/// with bit `i` = site `i`, matching the `common::enumerate_ising` index.
fn reachable_sector(n_sites: usize, pairs: &[(usize, usize)], initial: usize) -> Vec<usize> {
    let mut seen = vec![false; 1usize << n_sites];
    seen[initial] = true;
    let mut queue = VecDeque::from(vec![initial]);
    while let Some(state) = queue.pop_front() {
        for &(a, b) in pairs {
            if (state >> a) & 1 != (state >> b) & 1 {
                // Swapping two different bits is flipping both bits.
                let next = state ^ (1 << a) ^ (1 << b);
                if !seen[next] {
                    seen[next] = true;
                    queue.push_back(next);
                }
            }
        }
    }
    (0..1usize << n_sites)
        .filter(|&state| seen[state])
        .collect()
}

fn spins_of_mask(mask: usize, n_sites: usize) -> Vec<f64> {
    (0..n_sites)
        .map(|site| if (mask >> site) & 1 == 1 { 1.0 } else { -1.0 })
        .collect()
}

/// (⟨E⟩, ⟨E²⟩, ⟨m²⟩) of the canonical distribution restricted to `states`.
fn sector_moments(lattice: &CsrLattice, states: &[usize], beta: f64) -> (f64, f64, f64) {
    let all_states = enumerate_ising(lattice.n_sites);
    let mut partition = 0.0;
    let mut weighted_e = 0.0;
    let mut weighted_e2 = 0.0;
    let mut weighted_m2 = 0.0;
    for &state in states {
        let spins = &all_states[state];
        let energy = direct_ising_energy(spins, lattice, J);
        let magnetization = spins.iter().sum::<f64>() / lattice.n_sites as f64;
        let weight = (-beta * energy).exp();
        partition += weight;
        weighted_e += weight * energy;
        weighted_e2 += weight * energy * energy;
        weighted_m2 += weight * magnetization * magnetization;
    }
    (
        weighted_e / partition,
        weighted_e2 / partition,
        weighted_m2 / partition,
    )
}

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

/// Established multi-seed z criteria (mirrors `zscore_extended`): per-seed
/// |z| < 4, |z̄| < 2, and the seed-count-scaled one-sided bound Σz > −2√n.
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

/// One Kawasaki seed: return binned (⟨E⟩ ± SE, ⟨E²⟩ ± SE).
///
/// Asserts the defining invariant — total magnetization stays EXACTLY
/// constant after every sweep (Ising spins are ±1, so the sum is exact in
/// f64) — and that the incremental energy cache never drifts.
fn run_kawasaki_seed(
    lattice: &CsrLattice,
    initial_mask: usize,
    beta: f64,
    seed: u64,
) -> ((f64, f64), (f64, f64)) {
    let model = IsingModel::new(J);
    let initial_spins = spins_of_mask(initial_mask, lattice.n_sites);
    let target_magnetization = initial_spins.iter().sum::<f64>();
    let mut system = System::new(lattice.clone(), 1, 1.0, beta);
    system.spins.copy_from_slice(&initial_spins);
    system.recompute_energy(&model);

    let mut kernel = KawasakiCore::new(2 * lattice.n_sites);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    for _ in 0..THERM_SWEEPS {
        kernel.sweep(&mut system, &model, &mut rng);
        assert_eq!(
            system.spins.iter().sum::<f64>(),
            target_magnetization,
            "Kawasaki lost magnetization conservation at beta={beta}, seed={seed}"
        );
    }
    let mut energies = Vec::with_capacity(MEAS_SWEEPS);
    for _ in 0..MEAS_SWEEPS {
        kernel.sweep(&mut system, &model, &mut rng);
        assert_eq!(
            system.spins.iter().sum::<f64>(),
            target_magnetization,
            "Kawasaki lost magnetization conservation at beta={beta}, seed={seed}"
        );
        energies.push(system.energy);
    }
    assert!(
        system.energy_error(&model).abs() < 1e-9,
        "Kawasaki energy cache drifted by {}",
        system.energy_error(&model)
    );
    let energies_squared: Vec<f64> = energies.iter().map(|e| e * e).collect();
    (
        mean_stderr_binned(&energies, BINSIZE),
        mean_stderr_binned(&energies_squared, BINSIZE),
    )
}

/// Shared check: BFS sector → exact moments → multi-seed Kawasaki z-scores.
fn check_sector_equilibrium(lattice: &CsrLattice, initial_mask: usize, label: &str) {
    let n_sites = lattice.n_sites;
    let pairs = exchange_pairs(lattice);
    let states = reachable_sector(n_sites, &pairs, initial_mask);

    // Connectivity finding: the reachable set must be the full fixed-M
    // sector (no hidden Kawasaki invariant on PBC rings).
    let up_count = initial_mask.count_ones() as usize;
    let sector_size = (0..1usize << n_sites)
        .filter(|state| state.count_ones() as usize == up_count)
        .count();
    assert_eq!(
        states.len(),
        sector_size,
        "{label}: BFS reachable {} of {sector_size} fixed-M states \
         (extra Kawasaki invariant discovered)",
        states.len()
    );

    // The sector restriction is the point: ⟨m²⟩ is frozen at (M/N)².
    let magnetization = (2 * up_count) as f64 - n_sites as f64;
    for beta in [0.3, 0.8] {
        let (exact_e, exact_e2, sector_m2) = sector_moments(lattice, &states, beta);
        assert!(
            (sector_m2 - (magnetization / n_sites as f64).powi(2)).abs() < 1e-12,
            "{label}: sector ⟨m²⟩ should be frozen at (M/N)²"
        );
        let (_, canonical_e, _, _) = exact_ising_moments(lattice, J, beta);
        assert!(
            (exact_e - canonical_e).abs() > 0.1,
            "{label} β={beta}: sector ⟨E⟩={exact_e:.4} should differ from canonical \
             ⟨E⟩={canonical_e:.4} (the wrong reference would make this test vacuous)"
        );

        let n_seeds = zscore_seed_count(N_SEEDS);
        let mut energy_runs = Vec::with_capacity(n_seeds);
        let mut energy2_runs = Vec::with_capacity(n_seeds);
        for seed in 0..n_seeds as u64 {
            let ((e_mean, e_stderr), (e2_mean, e2_stderr)) =
                run_kawasaki_seed(lattice, initial_mask, beta, 0x4A41 + seed);
            energy_runs.push((e_mean, e_stderr));
            energy2_runs.push((e2_mean, e2_stderr));
        }
        let max_e_z = assert_exact_z(&energy_runs, exact_e, &format!("{label} β={beta} ⟨E⟩"));
        let max_e2_z = assert_exact_z(&energy2_runs, exact_e2, &format!("{label} β={beta} ⟨E²⟩"));
        eprintln!(
            "[kawasaki-exact] {label} β={beta}: sector ⟨E⟩={exact_e:.4} \
             (canonical {canonical_e:.4}), max|z| ⟨E⟩={max_e_z:.2}, ⟨E²⟩={max_e2_z:.2}"
        );
    }
}

#[test]
fn kawasaki_n4_m0_sector_equilibrium_matches_exact_enumeration() {
    let lattice = build_chain(4, true);
    // Sites 0,1 up; sites 2,3 down → M = 0.
    check_sector_equilibrium(&lattice, 0b0011, "N=4 chain PBC M=0");
}

#[test]
fn kawasaki_n8_m2_sector_equilibrium_matches_exact_enumeration() {
    let lattice = build_chain(8, true);
    // Sites 0..=4 up; sites 5..=7 down → M = +2 (M≠0 sector).
    check_sector_equilibrium(&lattice, 0b0001_1111, "N=8 chain PBC M=+2");
}
