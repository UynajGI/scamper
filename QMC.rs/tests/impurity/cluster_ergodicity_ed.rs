//! Ergodicity evidence for the longitudinal spin-boson cluster solver
//! (criterion E): multi-init convergence against exact diagonalization plus
//! explicit sector access.
//!
//! The cluster update redraws the Poisson cut set from the full rate
//! `Delta/2 > 0` each sweep and then re-orients every cluster independently,
//! so the transition density between any two worldlines is strictly
//! positive: the update graph is complete on the space of admissible
//! (even-kink) worldlines. The tests below verify the observable
//! consequences: runs started from constant `sigma_z = +1`, constant
//! `sigma_z = -1`, and a prepared many-kink worldline converge to the same
//! ED distribution, and a single chain demonstrably visits both spin
//! sectors, the zero-kink sector, and many distinct multi-kink sectors.

use qmc_rs::{
    Bath, ContinuousTimeClusterEngine, LongitudinalSpinBosonModel, LongitudinalWorldline,
    SingleModeBath,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

// ── Exact diagonalization for the single-mode longitudinal model ──────────

fn state_index(boson: usize, spin: usize) -> usize {
    2 * boson + spin
}

/// H = omega b†b + (epsilon/2) sigma_z - (Delta/2) sigma_x
///     + g_sigma sigma_z (b + b†),
/// the convention the cluster model uses with `lambda = g_sigma^2 / omega`.
#[allow(clippy::needless_range_loop)]
fn longitudinal_matrix(
    boson_states: usize,
    omega: f64,
    tunnelling: f64,
    bias: f64,
    coupling_sigma: f64,
) -> Vec<Vec<f64>> {
    let dimension = 2 * boson_states;
    let mut hamiltonian = vec![vec![0.0; dimension]; dimension];
    for boson in 0..boson_states {
        for spin in 0..2 {
            let index = state_index(boson, spin);
            let sigma_z = if spin == 0 { -1.0 } else { 1.0 };
            hamiltonian[index][index] += omega * boson as f64 + 0.5 * bias * sigma_z;
            let flipped = state_index(boson, 1 - spin);
            hamiltonian[index][flipped] += -0.5 * tunnelling;
            if boson + 1 < boson_states {
                let raised = state_index(boson + 1, spin);
                let element = coupling_sigma * sigma_z * (boson as f64 + 1.0).sqrt();
                hamiltonian[index][raised] += element;
                hamiltonian[raised][index] += element;
            }
        }
    }
    hamiltonian
}

fn sigma_z_operator(boson_states: usize) -> Vec<Vec<f64>> {
    let dimension = 2 * boson_states;
    let mut operator = vec![vec![0.0; dimension]; dimension];
    for boson in 0..boson_states {
        operator[state_index(boson, 0)][state_index(boson, 0)] = -1.0;
        operator[state_index(boson, 1)][state_index(boson, 1)] = 1.0;
    }
    operator
}

fn sigma_x_operator(boson_states: usize) -> Vec<Vec<f64>> {
    let dimension = 2 * boson_states;
    let mut operator = vec![vec![0.0; dimension]; dimension];
    for boson in 0..boson_states {
        operator[state_index(boson, 0)][state_index(boson, 1)] = 1.0;
        operator[state_index(boson, 1)][state_index(boson, 0)] = 1.0;
    }
    operator
}

#[allow(clippy::needless_range_loop)]
fn jacobi_eigensystem(mut matrix: Vec<Vec<f64>>) -> (Vec<f64>, Vec<Vec<f64>>) {
    let dimension = matrix.len();
    let mut vectors = vec![vec![0.0; dimension]; dimension];
    for (index, row) in vectors.iter_mut().enumerate() {
        row[index] = 1.0;
    }
    for _ in 0..(120 * dimension * dimension) {
        let mut left = 0usize;
        let mut right = 1usize;
        let mut largest = 0.0;
        for row in 0..dimension {
            for column in (row + 1)..dimension {
                if matrix[row][column].abs() > largest {
                    largest = matrix[row][column].abs();
                    left = row;
                    right = column;
                }
            }
        }
        if largest < 1e-13 {
            break;
        }
        let angle =
            0.5 * (2.0 * matrix[left][right]).atan2(matrix[right][right] - matrix[left][left]);
        let (sine, cosine) = angle.sin_cos();
        for index in 0..dimension {
            if index != left && index != right {
                let index_left = matrix[index][left];
                let index_right = matrix[index][right];
                matrix[index][left] = cosine * index_left - sine * index_right;
                matrix[left][index] = matrix[index][left];
                matrix[index][right] = sine * index_left + cosine * index_right;
                matrix[right][index] = matrix[index][right];
            }
        }
        let left_left = matrix[left][left];
        let right_right = matrix[right][right];
        let left_right = matrix[left][right];
        matrix[left][left] = cosine * cosine * left_left - 2.0 * sine * cosine * left_right
            + sine * sine * right_right;
        matrix[right][right] = sine * sine * left_left
            + 2.0 * sine * cosine * left_right
            + cosine * cosine * right_right;
        matrix[left][right] = 0.0;
        matrix[right][left] = 0.0;
        for row in &mut vectors {
            let vector_left = row[left];
            let vector_right = row[right];
            row[left] = cosine * vector_left - sine * vector_right;
            row[right] = sine * vector_left + cosine * vector_right;
        }
    }
    let mut order = (0..dimension).collect::<Vec<_>>();
    order.sort_by(|&left, &right| matrix[left][left].total_cmp(&matrix[right][right]));
    let energies = order.iter().map(|&index| matrix[index][index]).collect();
    let sorted_vectors = (0..dimension)
        .map(|row| order.iter().map(|&column| vectors[row][column]).collect())
        .collect();
    (energies, sorted_vectors)
}

#[allow(clippy::needless_range_loop)]
fn thermal_expectation(
    energies: &[f64],
    vectors: &[Vec<f64>],
    operator: &[Vec<f64>],
    beta: f64,
) -> f64 {
    let ground = energies[0];
    let mut numerator = 0.0;
    let mut partition = 0.0;
    for (state, &energy) in energies.iter().enumerate() {
        let weight = (-beta * (energy - ground)).exp();
        let mut expectation = 0.0;
        for row in 0..operator.len() {
            for column in 0..operator.len() {
                expectation += vectors[row][state] * operator[row][column] * vectors[column][state];
            }
        }
        numerator += weight * expectation;
        partition += weight;
    }
    numerator / partition
}

/// `C(beta/2) = Z^-1 sum_{k,l} |<k|O|l>|^2 e^{-beta(E_k+E_l)/2}`.
#[allow(clippy::needless_range_loop)]
fn correlation_half(
    energies: &[f64],
    vectors: &[Vec<f64>],
    operator: &[Vec<f64>],
    beta: f64,
) -> f64 {
    let dimension = energies.len();
    let ground = energies[0];
    let weights: Vec<f64> = energies
        .iter()
        .map(|&energy| (-0.5 * beta * (energy - ground)).exp())
        .collect();
    let partition: f64 = weights.iter().map(|weight| weight * weight).sum();
    let mut result = 0.0;
    for k in 0..dimension {
        for l in 0..dimension {
            let mut element = 0.0;
            for row in 0..dimension {
                for column in 0..dimension {
                    element += vectors[row][k] * operator[row][column] * vectors[column][l];
                }
            }
            result += element * element * weights[k] * weights[l];
        }
    }
    result / partition
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

// ── Multi-init convergence against ED ──────────────────────────────────────

#[test]
fn cluster_multi_init_converges_to_ed() {
    let beta = 2.4;
    let omega = 1.25;
    let tunnelling = 0.9;
    let bias = 0.35;
    let coupling_sigma = 0.31;
    let boson_states = 12;
    let lambda = coupling_sigma * coupling_sigma / omega;

    let hamiltonian = longitudinal_matrix(boson_states, omega, tunnelling, bias, coupling_sigma);
    let (energies, vectors) = jacobi_eigensystem(hamiltonian);
    let exact_sigma_z =
        thermal_expectation(&energies, &vectors, &sigma_z_operator(boson_states), beta);
    let exact_sigma_x =
        thermal_expectation(&energies, &vectors, &sigma_x_operator(boson_states), beta);
    let exact_correlation =
        correlation_half(&energies, &vectors, &sigma_z_operator(boson_states), beta);

    let model = LongitudinalSpinBosonModel::with_default_quadrature(
        Bath::SingleMode(SingleModeBath::new(omega).expect("mode")),
        lambda,
        tunnelling,
        bias,
    )
    .expect("model");

    // Three genuinely distinct initial worldlines.
    let inits: Vec<(&str, LongitudinalWorldline, u64)> = vec![
        (
            "constant +1",
            LongitudinalWorldline::new(beta, 1).expect("worldline"),
            0xCEA1_u64,
        ),
        (
            "constant -1",
            LongitudinalWorldline::new(beta, -1).expect("worldline"),
            0xCEA2,
        ),
        (
            "six prepared kinks",
            LongitudinalWorldline::from_kinks(beta, 1, vec![0.3, 0.8, 1.1, 1.5, 1.9, 2.2])
                .expect("worldline"),
            0xCEA3,
        ),
    ];

    let warmup = 5_000;
    let samples = 100_000usize;
    let mut results: Vec<[(f64, f64); 3]> = Vec::new();
    for (label, worldline, seed) in inits {
        let mut engine = ContinuousTimeClusterEngine::new(model.clone());
        let mut worldline = worldline;
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
        for _ in 0..warmup {
            engine.update(&mut worldline, &mut rng).expect("warmup");
        }
        let mut sigma_z = Vec::with_capacity(samples);
        let mut sigma_x = Vec::with_capacity(samples);
        let mut correlation = Vec::with_capacity(samples);
        for _ in 0..samples {
            let report = engine.update(&mut worldline, &mut rng).expect("update");
            sigma_z.push(report.improved_magnetization_sigma_z);
            sigma_x.push(2.0 * worldline.kink_count() as f64 / (beta * tunnelling));
            correlation.push(worldline.correlation_sigma_z(0.5 * beta));
        }
        let blocks = 16;
        results.push([
            (
                sigma_z.iter().sum::<f64>() / samples as f64,
                blocked_stderr(&sigma_z, blocks),
            ),
            (
                sigma_x.iter().sum::<f64>() / samples as f64,
                blocked_stderr(&sigma_x, blocks),
            ),
            (
                correlation.iter().sum::<f64>() / samples as f64,
                blocked_stderr(&correlation, blocks),
            ),
        ]);
        eprintln!(
            "init {label:18}: <sz> = {:.4} ± {:.4}, <sx> = {:.4} ± {:.4}, C(b/2) = {:.4} ± {:.4}",
            results.last().expect("result")[0].0,
            results.last().expect("result")[0].1,
            results.last().expect("result")[1].0,
            results.last().expect("result")[1].1,
            results.last().expect("result")[2].0,
            results.last().expect("result")[2].1,
        );
    }

    let references = [
        ("⟨σz⟩", exact_sigma_z),
        ("⟨σx⟩", exact_sigma_x),
        ("C(β/2)", exact_correlation),
    ];
    for (slot, (label, exact)) in references.iter().enumerate() {
        for (index, result) in results.iter().enumerate() {
            let (value, stderr) = result[slot];
            let z = (value - exact) / stderr.max(1.0e-12);
            assert!(
                z.abs() < 4.0,
                "init {index} ({label}): sampled {value:.5} ± {stderr:.5} vs ED {exact:.5} \
                 (z = {z:.2})"
            );
        }
        for left in 0..results.len() {
            for right in (left + 1)..results.len() {
                let combined =
                    (results[left][slot].1.powi(2) + results[right][slot].1.powi(2)).sqrt();
                let deviation = (results[left][slot].0 - results[right][slot].0).abs();
                assert!(
                    deviation < 4.0 * combined,
                    "{label}: inits {left} and {right} disagree by {deviation:.5} \
                     (4σ = {:.5})",
                    4.0 * combined
                );
            }
        }
    }
}

// ── Sector access from a single chain ──────────────────────────────────────

#[test]
fn cluster_chain_visits_both_spin_and_many_kink_sectors() {
    let beta = 3.0;
    let omega = 1.1;
    let coupling_sigma = 0.35;
    let lambda = coupling_sigma * coupling_sigma / omega;
    let model = LongitudinalSpinBosonModel::with_default_quadrature(
        Bath::SingleMode(SingleModeBath::new(omega).expect("mode")),
        lambda,
        0.8,
        0.15,
    )
    .expect("model");
    let mut engine = ContinuousTimeClusterEngine::new(model);
    let mut worldline = LongitudinalWorldline::new(beta, 1).expect("worldline");
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x5EC7_2026);

    let mut saw_plus = false;
    let mut saw_minus = false;
    let mut kink_sectors = std::collections::HashSet::new();
    let mut distinct_worldlines = std::collections::HashSet::new();
    let sweeps = 40_000usize;
    for _ in 0..sweeps {
        engine.update(&mut worldline, &mut rng).expect("update");
        assert_eq!(worldline.kink_count() % 2, 0, "periodic invariant broken");
        if worldline.spin_at_zero() == 1 {
            saw_plus = true;
        } else {
            saw_minus = true;
        }
        kink_sectors.insert(worldline.kink_count());
        distinct_worldlines.insert((
            worldline.spin_at_zero(),
            worldline.kink_count(),
            worldline
                .kinks()
                .iter()
                .map(|time| (time * 1000.0) as u64)
                .collect::<Vec<_>>(),
        ));
    }
    assert!(
        saw_plus && saw_minus,
        "chain never visited both spin-at-zero sectors"
    );
    assert!(
        kink_sectors.contains(&0),
        "chain never visited the zero-kink sector"
    );
    assert!(
        kink_sectors.iter().any(|&kinks| kinks >= 4),
        "chain never visited a >=4-kink sector (visited: {kink_sectors:?})"
    );
    assert!(
        kink_sectors.len() >= 5,
        "chain visited only {} kink sectors",
        kink_sectors.len()
    );
    assert!(
        distinct_worldlines.len() >= 1_000,
        "chain visited only {} distinct worldlines in {sweeps} sweeps",
        distinct_worldlines.len()
    );
}
