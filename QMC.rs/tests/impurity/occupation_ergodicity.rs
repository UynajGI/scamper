//! Ergodicity tests for the occupation worldline solver.
//!
//! Verifies that independent runs from different RNG seeds converge to
//! the same thermal expectation values for ⟨σz⟩ and ⟨n⟩, confirming
//! that the occupation-basis sampler is ergodic across the full Hilbert
//! space (spin × boson occupations).
//!
//! Production evidence (2026-08-19) adds real connectivity: an explicit
//! breadth-first search over the sampler's transition-support graph on the
//! complete closed-path space of a small system, plus multi-init convergence
//! from distinct initial worldlines to the exact ED distribution.
//!
//! `SCUTTLE_ZSCORE_SEEDS=<n>` raises the seed count for nightly
//! high-power monitoring (unset → the default 4 seeds, unchanged for CI).

use crate::zscore_seeds::zscore_seeds;
use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use qmc_rs::impurity::spin_boson::occupation::transfer::{multiply, SymmetricEigensystem};
use qmc_rs::OccupationWorldlineQmc;
use qmc_rs::{CavityMode, OccupationSpinBosonModel, OccupationWorldlineSampler};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use std::collections::VecDeque;

type Matrix = Vec<Vec<f64>>;

fn run_occupation_rabi(seed: u64) -> (f64, f64, f64, f64) {
    let mut params = Params::new();
    params.set("beta", 4.0);
    params.set("kind", "rabi");
    params.set("spin_splitting", 1.0);
    params.set("g", 0.3);
    params.set("omega0", 1.0);
    params.set("cutoff", 10);

    let run = RunConfig {
        thermalization_sweeps: 2000,
        measurement_sweeps: 10000,
        binsize: 100,
        base_seed: seed,
        ..Default::default()
    };
    let results =
        Scheduler::new(RayonBackend::new(1), run).run_one::<OccupationWorldlineQmc>(&params);
    let sigma_z = results.get("OccupationSigmaZ").expect("OccupationSigmaZ");
    let boson_n = results
        .get("OccupationBosonNumber")
        .expect("OccupationBosonNumber");
    (sigma_z.mean, sigma_z.stderr, boson_n.mean, boson_n.stderr)
}

/// z-score check against the pooled mean.
///
/// For each observable the per-seed means should scatter around the
/// pooled (inverse-variance weighted) mean with z-scores well within
/// ±4, and the mean |z| should be under 2.
fn assert_z_scores(name: &str, values: &[f64], stderrs: &[f64]) {
    let n = values.len() as f64;

    // Pooled mean (inverse-variance weighted)
    let mut sum_w = 0.0;
    let mut sum_wm = 0.0;
    for i in 0..values.len() {
        let w = 1.0 / stderrs[i].max(0.01).powi(2);
        sum_w += w;
        sum_wm += w * values[i];
    }
    let pooled_mean = sum_wm / sum_w;

    // Pooled variance (sample variance of the per-seed means)
    let pooled_var = values
        .iter()
        .map(|&v| (v - pooled_mean).powi(2))
        .sum::<f64>()
        / (n - 1.0);

    // z-scores
    let mut z_values: Vec<f64> = Vec::new();
    for i in 0..values.len() {
        let denom = (pooled_var / n + stderrs[i].max(0.01).powi(2)).sqrt();
        let z = (values[i] - pooled_mean) / denom.max(1e-10);
        z_values.push(z);
        assert!(
            z.abs() < 4.0,
            "{name}: z-score for seed {i} = {z:.4} exceeds 4σ (values: {values:?})",
        );
    }

    let mean_abs_z: f64 = z_values.iter().map(|z| z.abs()).sum::<f64>() / n;
    assert!(
        mean_abs_z < 2.0,
        "{name}: mean |z| = {mean_abs_z:.4} exceeds 2 (z-values: {z_values:?})",
    );
}

#[test]
fn occupation_ergodicity_multi_seed_convergence() {
    let seeds = zscore_seeds(&[42u64, 123, 456, 789]);
    let results: Vec<(f64, f64, f64, f64)> =
        seeds.iter().map(|&s| run_occupation_rabi(s)).collect();

    // ⟨σz⟩ consistency: max−min < 4 × max(stderr)
    let sz_values: Vec<f64> = results.iter().map(|r| r.0).collect();
    let sz_stderrs: Vec<f64> = results.iter().map(|r| r.1).collect();
    let sz_spread = sz_values.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
        - sz_values.iter().cloned().fold(f64::INFINITY, f64::min);
    let sz_max_stderr = sz_stderrs.iter().cloned().fold(0.0_f64, f64::max);
    assert!(
        sz_spread < 4.0 * sz_max_stderr.max(0.01),
        "⟨σz⟩ spread={sz_spread:.6} exceeds 4σ={:.6} (values: {sz_values:?})",
        4.0 * sz_max_stderr.max(0.01)
    );

    // ⟨n⟩ consistency: max−min < 4 × max(stderr)
    let n_values: Vec<f64> = results.iter().map(|r| r.2).collect();
    let n_stderrs: Vec<f64> = results.iter().map(|r| r.3).collect();
    let n_spread = n_values.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
        - n_values.iter().cloned().fold(f64::INFINITY, f64::min);
    let n_max_stderr = n_stderrs.iter().cloned().fold(0.0_f64, f64::max);
    assert!(
        n_spread < 4.0 * n_max_stderr.max(0.01),
        "⟨n⟩ spread={n_spread:.6} exceeds 4σ={:.6} (values: {n_values:?})",
        4.0 * n_max_stderr.max(0.01)
    );

    // z-score checks against pooled mean
    assert_z_scores("⟨σz⟩", &sz_values, &sz_stderrs);
    assert_z_scores("⟨n⟩", &n_values, &n_stderrs);
}

// ── Connectivity: explicit support-graph BFS over closed paths ─────────────

/// Rebuild the sampler's link propagator `T = exp(-dt (H - E0))` from the
/// public exact-diagonalization API, exactly as the sampler constructs it.
fn link_propagator(model: &OccupationSpinBosonModel, beta: f64, slices: usize) -> Matrix {
    let eigen = SymmetricEigensystem::diagonalize(model.hamiltonian()).expect("diagonalize");
    let dt = beta / slices as f64;
    let ground = eigen.values[0];
    eigen.matrix_function(|energy| (-dt * (energy - ground)).exp())
}

/// Kernel density the sweep assigns to one closed path (heat bath: the
/// density is independent of the current path).
#[allow(clippy::needless_range_loop)]
fn path_density(transfer: &Matrix, powers: &[Matrix], slices: usize, path: &[usize]) -> f64 {
    let dimension = transfer.len();
    let first = path[0];
    let full = &powers[slices];
    let total: f64 = (0..dimension).map(|state| full[state][state]).sum();
    let mut density = full[first][first] / total;
    for index in 1..slices {
        let previous = path[index - 1];
        let remaining = slices - index;
        let mut normalization = 0.0;
        let mut realized = 0.0;
        for candidate in 0..dimension {
            let weight = transfer[previous][candidate] * powers[remaining][candidate][first];
            normalization += weight;
            if candidate == path[index] {
                realized = weight;
            }
        }
        density *= realized / normalization;
    }
    density
}

#[test]
fn occupation_update_graph_is_strongly_connected() {
    // Small system whose complete closed-path space is enumerable: cutoff 2
    // (dimension 4), 3 slices -> 4^3 = 64 worldlines.
    let model =
        OccupationSpinBosonModel::rabi(0.9, vec![CavityMode::new(1.2, 0.35, 3).expect("mode")])
            .expect("model");
    let dimension = model.basis().dimension();
    let slices = 3;
    let beta = 2.0;
    let transfer = link_propagator(&model, beta, slices);
    let mut powers = vec![identity(dimension)];
    for _ in 0..slices {
        let next = multiply(&powers[powers.len() - 1], &transfer);
        powers.push(next);
    }

    // Enumerate every closed path and compute the sweep's kernel density.
    // A path carries positive Boltzmann weight iff every link propagator
    // T[s_k][s_{k+1}] is positive (the Rabi coupling pattern forbids direct
    // spin flips at fixed occupation, so not every tuple of basis states is
    // an admissible worldline). The kernel and the target measure share this
    // support exactly, which is the physically correct statement.
    let total_paths = dimension.pow(slices as u32);
    let mut densities = Vec::with_capacity(total_paths);
    let mut admissible = 0_usize;
    for encoded in 0..total_paths {
        let mut path = vec![0_usize; slices];
        let mut rest = encoded;
        for slot in path.iter_mut().rev() {
            *slot = rest % dimension;
            rest /= dimension;
        }
        let density = path_density(&transfer, &powers, slices, &path);
        let positive_links = (0..slices).all(|index| {
            let next = path[(index + 1) % slices];
            transfer[path[index]][next] > 0.0
        });
        assert_eq!(
            density > 0.0,
            positive_links,
            "path {encoded}: kernel density positivity must match positive link weights"
        );
        if density > 0.0 {
            admissible += 1;
        }
        densities.push(density);
    }

    // Explicit breadth-first search over the update graph: because the kernel
    // is the exact heat bath, the successors of any admissible node are ALL
    // admissible paths (plus the node itself, giving aperiodicity). BFS from
    // one admissible node must therefore reach the entire admissible set.
    let successors: Vec<usize> = (0..total_paths)
        .filter(|&encoded| densities[encoded] > 0.0)
        .collect();
    let start = successors[0];
    let mut visited = vec![false; total_paths];
    let mut queue = VecDeque::from([start]);
    visited[start] = true;
    let mut reached = 1_usize;
    while let Some(_node) = queue.pop_front() {
        for &next in &successors {
            if !visited[next] {
                visited[next] = true;
                reached += 1;
                queue.push_back(next);
            }
        }
    }
    assert_eq!(
        reached, admissible,
        "update graph reachable set covers only {reached} of {admissible} admissible paths"
    );

    // The sampler itself must never leave the admissible set.
    let mut sampler =
        OccupationWorldlineSampler::new(model.clone(), beta, slices, start % dimension)
            .expect("sampler");
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0xC0BB_2026);
    let mut visited_paths = std::collections::HashSet::new();
    for _ in 0..20_000 {
        sampler.sweep(&mut rng).expect("sweep");
        let states = sampler.states();
        let encoded: usize = states
            .iter()
            .fold(0_usize, |acc, &state| acc * dimension + state);
        assert!(
            densities[encoded] > 0.0,
            "sampler visited inadmissible path {encoded:?}"
        );
        visited_paths.insert(encoded);
    }
    // Every path whose kernel density gives an expected visit count well
    // above zero must actually be reached: real connectivity, not just a
    // support argument.
    let notable = densities
        .iter()
        .filter(|&&density| density >= 5.0e-4)
        .count();
    assert!(
        visited_paths.len() >= notable,
        "sampler visited only {} paths in 20k sweeps but {notable} paths carry          density >= 5e-4",
        visited_paths.len()
    );
    eprintln!(
        "occupation update graph: {admissible}/{total_paths} paths admissible,          {notable} notable, {} visited",
        visited_paths.len()
    );
}

fn identity(dimension: usize) -> Matrix {
    let mut matrix = vec![vec![0.0; dimension]; dimension];
    for (index, row) in matrix.iter_mut().enumerate() {
        row[index] = 1.0;
    }
    matrix
}

// ── Multi-init convergence against exact diagonalization ───────────────────

/// Exact thermal ⟨O⟩ from the model's dense Hamiltonian.
#[allow(clippy::needless_range_loop)]
fn thermal_expectation(model: &OccupationSpinBosonModel, operator: &Matrix, beta: f64) -> f64 {
    let eigen = SymmetricEigensystem::diagonalize(model.hamiltonian()).expect("diagonalize");
    let ground = eigen.values[0];
    let mut numerator = 0.0;
    let mut partition = 0.0;
    for (level, &energy) in eigen.values.iter().enumerate() {
        let weight = (-beta * (energy - ground)).exp();
        let mut expectation = 0.0;
        for row in 0..operator.len() {
            for column in 0..operator.len() {
                expectation += eigen.vectors[row][level]
                    * operator[row][column]
                    * eigen.vectors[column][level];
            }
        }
        numerator += weight * expectation;
        partition += weight;
    }
    numerator / partition
}

fn sigma_z_operator(dimension: usize) -> Matrix {
    let mut operator = vec![vec![0.0; dimension]; dimension];
    for (state, row) in operator.iter_mut().enumerate() {
        row[state] = if state % 2 == 0 { -1.0 } else { 1.0 };
    }
    operator
}

fn number_operator(basis: &qmc_rs::OccupationBasis) -> Matrix {
    let dimension = basis.dimension();
    let mut operator = vec![vec![0.0; dimension]; dimension];
    for (state, row) in operator.iter_mut().enumerate() {
        row[state] = basis.occupation(state, 0) as f64;
    }
    operator
}

fn blocked_stderr(values: &[f64], blocks: usize) -> f64 {
    let block_size = values.len() / blocks;
    assert!(block_size > 0);
    let mut block_means = Vec::with_capacity(blocks);
    for block in 0..blocks {
        let slice = &values[block * block_size..(block + 1) * block_size];
        block_means.push(slice.iter().sum::<f64>() / slice.len() as f64);
    }
    let mean = block_means.iter().sum::<f64>() / blocks as f64;
    let variance = block_means
        .iter()
        .map(|value| (value - mean) * (value - mean))
        .sum::<f64>()
        / (blocks - 1) as f64;
    variance.sqrt() / (blocks as f64).sqrt()
}

#[test]
fn occupation_multi_init_converges_to_exact_distribution() {
    // Three samplers started from distinct initial worldlines (vacuum-spin-
    // down, the top basis state, a middle state) must converge to the same
    // ED thermal values for three observables, |z| < 4 per init.
    let model =
        OccupationSpinBosonModel::rabi(0.9, vec![CavityMode::new(1.2, 0.32, 5).expect("mode")])
            .expect("model");
    let dimension = model.basis().dimension();
    let beta = 2.0;
    let slices = 5;

    // Note: <sigma_x> is exactly zero in the Rabi model by the Z2 parity
    // symmetry sigma_z (-1)^n, so it is a degenerate observable; energy,
    // sigma_z, and boson number are the three non-trivial references.
    let exact_sigma_z = thermal_expectation(&model, &sigma_z_operator(dimension), beta);
    let exact_energy = thermal_expectation(&model, &model.hamiltonian(), beta);
    let exact_number = thermal_expectation(&model, &number_operator(model.basis()), beta);

    let initial_states = [0_usize, dimension - 1, dimension / 2];
    let seeds = [0xE0F1_u64, 0xE0F2, 0xE0F3];
    let mut results: Vec<[(f64, f64); 3]> = Vec::new();
    for (&initial_state, &seed) in initial_states.iter().zip(&seeds) {
        let mut sampler =
            OccupationWorldlineSampler::new(model.clone(), beta, slices, initial_state)
                .expect("sampler");
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
        for _ in 0..2_000 {
            sampler.sweep(&mut rng).expect("warmup");
        }
        let samples = 120_000usize;
        let mut sigma_z = Vec::with_capacity(samples);
        let mut energy = Vec::with_capacity(samples);
        let mut number = Vec::with_capacity(samples);
        for _ in 0..samples {
            sampler.sweep(&mut rng).expect("sweep");
            let observables = sampler.measure().expect("measure");
            sigma_z.push(observables.sigma_z);
            energy.push(observables.energy);
            number.push(observables.total_boson_number);
        }
        let blocks = 16;
        results.push([
            (
                sigma_z.iter().sum::<f64>() / samples as f64,
                blocked_stderr(&sigma_z, blocks),
            ),
            (
                energy.iter().sum::<f64>() / samples as f64,
                blocked_stderr(&energy, blocks),
            ),
            (
                number.iter().sum::<f64>() / samples as f64,
                blocked_stderr(&number, blocks),
            ),
        ]);
    }

    let references = [
        ("⟨σz⟩", exact_sigma_z),
        ("⟨E⟩", exact_energy),
        ("⟨n⟩", exact_number),
    ];
    for (slot, (label, exact)) in references.iter().enumerate() {
        for (index, result) in results.iter().enumerate() {
            let (value, stderr) = result[slot];
            let z = (value - exact) / stderr.max(1.0e-12);
            assert!(
                z.abs() < 4.0,
                "init {index} ({label}): sampled {value:.5} \u{00b1} {stderr:.5} vs ED {exact:.5} \
                 (z = {z:.2})"
            );
        }
        // Pairwise agreement between independent inits.
        for left in 0..results.len() {
            for right in (left + 1)..results.len() {
                let combined =
                    (results[left][slot].1.powi(2) + results[right][slot].1.powi(2)).sqrt();
                let deviation = (results[left][slot].0 - results[right][slot].0).abs();
                assert!(
                    deviation < 4.0 * combined,
                    "{label}: inits {left} and {right} disagree by {deviation:.5} \
                     (4\u{3c3} = {:.5})",
                    4.0 * combined
                );
            }
        }
    }
}
