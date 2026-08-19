//! Multi-mode bath support for the longitudinal spin-boson cluster solver
//! (criterion H): the former "single-mode only" limitation closed by
//! validation, not documentation.
//!
//! The cluster engine consumes the bath only through its retarded kernel,
//! and `RetardedKernel` carries the discrete spectral measure of a
//! `Bath::Tabulated` multi-mode bath. Three pieces of evidence:
//!
//! 1. A machine-precision kernel identity: the tabulated multi-mode kernel
//!    equals the mass-weighted sum of single-mode kernels (criterion A).
//! 2. Interacting multi-mode MC-vs-ED on three observables
//!    (⟨σz⟩, ⟨σx⟩ via kinks, C(β/2)) with per-seed |z| < 4 (criterion B).
//! 3. A cross-solver comparison against the wormhole solver on the same
//!    multi-mode physical Hamiltonian (criterion F).
//!
//! The shared physical model is
//! `H = sum_k omega_k b_k†b_k + (epsilon/2) sigma_z - (Delta/2) sigma_x
//!      + sigma_z sum_k g_k (b_k + b_k†)`,
//! with the cluster convention `lambda = sum_k g_k^2 / omega_k`
//! (`g_k^2 = lambda w_k omega_k` for normalized tabulated masses `w_k`).

use crate::zscore_seeds::zscore_seeds;
use qmc_rs::impurity::ImpurityQmc;
use qmc_rs::{
    Bath, ContinuousTimeClusterEngine, LongitudinalSpinBosonModel, LongitudinalWorldline,
    RetardedKernel, SingleModeBath, TabulatedBath,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

const MODE_OMEGAS: [f64; 2] = [0.9, 1.7];
const MODE_WEIGHTS: [f64; 2] = [1.0, 0.6];
const LAMBDA: f64 = 0.24;
const TUNNELLING: f64 = 0.8;
const BIAS: f64 = 0.3;
const BETA: f64 = 2.5;
const BOSON_CUTOFFS: [usize; 2] = [7, 5];

fn multimode_bath() -> Bath {
    Bath::Tabulated(
        TabulatedBath::new(MODE_OMEGAS.to_vec(), MODE_WEIGHTS.to_vec()).expect("tabulated bath"),
    )
}

fn normalized_weights() -> [f64; 2] {
    let total: f64 = MODE_WEIGHTS.iter().sum();
    [MODE_WEIGHTS[0] / total, MODE_WEIGHTS[1] / total]
}

/// Per-mode coupling `g_k` of the shared physical Hamiltonian.
fn mode_couplings(lambda: f64) -> [f64; 2] {
    let weights = normalized_weights();
    [
        (lambda * weights[0] * MODE_OMEGAS[0]).sqrt(),
        (lambda * weights[1] * MODE_OMEGAS[1]).sqrt(),
    ]
}

// ── 1. Kernel identity: tabulated = mass-weighted sum of single modes ──────

#[test]
fn multimode_kernel_equals_mass_weighted_single_mode_sum() {
    let tabulated = RetardedKernel::with_default_quadrature(&multimode_bath(), LAMBDA)
        .expect("tabulated kernel");
    let weights = normalized_weights();
    let single: Vec<RetardedKernel> = MODE_OMEGAS
        .iter()
        .zip(weights)
        .map(|(&omega, weight)| {
            RetardedKernel::with_default_quadrature(
                &Bath::SingleMode(SingleModeBath::new(omega).expect("mode")),
                LAMBDA * weight,
            )
            .expect("single-mode kernel")
        })
        .collect();

    for beta in [1.5, 2.5, 6.0] {
        for delta in [0.0, 0.17, beta / 2.0, beta - 0.31] {
            let combined: f64 = single
                .iter()
                .map(|kernel| kernel.value(beta, delta).expect("value"))
                .sum();
            let direct = tabulated.value(beta, delta).expect("value");
            assert!(
                (combined - direct).abs() < 1.0e-12,
                "K_tabulated(beta={beta}, tau={delta:.4}) = {direct} vs single-mode sum \
                 {combined}"
            );
        }
    }

    // Same identity for the interval-integrated kernel used by the cluster
    // retarded bonds.
    use qmc_rs::impurity::spin_boson::cluster::segments::TimeInterval;
    let beta = BETA;
    let left = TimeInterval::new(0.2, 0.9, beta).expect("interval");
    let right = TimeInterval::new(1.4, 2.3, beta).expect("interval");
    let combined: f64 = single
        .iter()
        .map(|kernel| {
            kernel
                .integrated_intervals(beta, left, right)
                .expect("integral")
        })
        .sum();
    let direct = tabulated
        .integrated_intervals(beta, left, right)
        .expect("integral");
    assert!(
        (combined - direct).abs() < 1.0e-12,
        "integrated tabulated kernel {direct} vs single-mode sum {combined}"
    );
}

// ── 2. Interacting multi-mode MC vs ED ─────────────────────────────────────

/// Dense two-mode Hamiltonian in the |n1, n2, spin⟩ occupation basis.
#[allow(clippy::needless_range_loop)]
fn multimode_matrix_with(lambda: f64, tunnelling: f64, bias: f64) -> Vec<Vec<f64>> {
    let [c1, c2] = BOSON_CUTOFFS;
    let [g1, g2] = mode_couplings(lambda);
    let boson_dimension = c1 * c2;
    let dimension = 2 * boson_dimension;
    let mut hamiltonian = vec![vec![0.0; dimension]; dimension];
    let index = |n1: usize, n2: usize, spin: usize| 2 * (n1 * c2 + n2) + spin;
    for n1 in 0..c1 {
        for n2 in 0..c2 {
            for spin in 0..2 {
                let state = index(n1, n2, spin);
                let sigma_z = if spin == 0 { -1.0 } else { 1.0 };
                hamiltonian[state][state] +=
                    MODE_OMEGAS[0] * n1 as f64 + MODE_OMEGAS[1] * n2 as f64 + 0.5 * bias * sigma_z;
                // Tunnelling flips the spin at fixed occupations.
                let flipped = index(n1, n2, 1 - spin);
                hamiltonian[state][flipped] += -0.5 * tunnelling;
                // sigma_z (b_k + b_k†) per mode.
                for (mode, n) in [(0_usize, n1), (1, n2)] {
                    let cutoff = BOSON_CUTOFFS[mode];
                    let coupling = [g1, g2][mode];
                    if n + 1 < cutoff {
                        let raised = if mode == 0 {
                            index(n + 1, n2, spin)
                        } else {
                            index(n1, n + 1, spin)
                        };
                        let element = coupling * sigma_z * (n as f64 + 1.0).sqrt();
                        hamiltonian[state][raised] += element;
                        hamiltonian[raised][state] += element;
                    }
                }
            }
        }
    }
    hamiltonian
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

fn multimode_matrix(bias: f64) -> Vec<Vec<f64>> {
    multimode_matrix_with(LAMBDA, TUNNELLING, bias)
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
fn multimode_cluster_matches_ed_three_observables() {
    let dimension = 2 * BOSON_CUTOFFS[0] * BOSON_CUTOFFS[1];
    let mut sigma_z = vec![vec![0.0; dimension]; dimension];
    let mut sigma_x = vec![vec![0.0; dimension]; dimension];
    let encode = |n1: usize, n2: usize, spin: usize| 2 * (n1 * BOSON_CUTOFFS[1] + n2) + spin;
    for n1 in 0..BOSON_CUTOFFS[0] {
        for n2 in 0..BOSON_CUTOFFS[1] {
            sigma_z[encode(n1, n2, 0)][encode(n1, n2, 0)] = -1.0;
            sigma_z[encode(n1, n2, 1)][encode(n1, n2, 1)] = 1.0;
            sigma_x[encode(n1, n2, 0)][encode(n1, n2, 1)] = 1.0;
            sigma_x[encode(n1, n2, 1)][encode(n1, n2, 0)] = 1.0;
        }
    }
    let (energies, vectors) = jacobi_eigensystem(multimode_matrix(BIAS));
    let exact_sigma_z = thermal_expectation(&energies, &vectors, &sigma_z, BETA);
    let exact_sigma_x = thermal_expectation(&energies, &vectors, &sigma_x, BETA);
    let exact_correlation = correlation_half(&energies, &vectors, &sigma_z, BETA);
    eprintln!(
        "multi-mode ED (dim {dimension}): <sz>={exact_sigma_z:.5}, \
         <sx>={exact_sigma_x:.5}, C(b/2)={exact_correlation:.5}"
    );

    let model = LongitudinalSpinBosonModel::with_default_quadrature(
        multimode_bath(),
        LAMBDA,
        TUNNELLING,
        BIAS,
    )
    .expect("model");

    let seeds = zscore_seeds(&[0x40A1_u64, 0x40A2, 0x40A3, 0x40A4]);
    let warmup = 5_000;
    let samples = 100_000usize;
    for (run, &seed) in seeds.iter().enumerate() {
        let mut engine = ContinuousTimeClusterEngine::new(model.clone());
        let mut worldline = LongitudinalWorldline::new(BETA, 1).expect("worldline");
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
        for _ in 0..warmup {
            engine.update(&mut worldline, &mut rng).expect("warmup");
        }
        let mut sigma_z_samples = Vec::with_capacity(samples);
        let mut sigma_x_samples = Vec::with_capacity(samples);
        let mut correlation_samples = Vec::with_capacity(samples);
        for _ in 0..samples {
            let report = engine.update(&mut worldline, &mut rng).expect("update");
            sigma_z_samples.push(report.improved_magnetization_sigma_z);
            sigma_x_samples.push(2.0 * worldline.kink_count() as f64 / (BETA * TUNNELLING));
            correlation_samples.push(worldline.correlation_sigma_z(0.5 * BETA));
        }
        let blocks = 16;
        for (label, values, exact) in [
            ("⟨σz⟩", &sigma_z_samples, exact_sigma_z),
            ("⟨σx⟩", &sigma_x_samples, exact_sigma_x),
            ("C(β/2)", &correlation_samples, exact_correlation),
        ] {
            let mean = values.iter().sum::<f64>() / samples as f64;
            let stderr = blocked_stderr(values, blocks);
            let z = (mean - exact) / stderr.max(1.0e-12);
            eprintln!("seed {run}: {label} = {mean:.5} ± {stderr:.5} vs ED {exact:.5} (z={z:.2})");
            assert!(
                z.abs() < 4.0,
                "multi-mode cluster seed {run} ({label}): {mean:.5} ± {stderr:.5} vs ED \
                 {exact:.5} (z = {z:.2})"
            );
        }
    }
}

// ── 3. Cross-solver: wormhole on the same multi-mode model ─────────────────

#[test]
fn multimode_cluster_agrees_with_wormhole() {
    use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};

    // Two shared parameter points of the unbiased multi-mode model (the
    // wormhole Rabi catalog carries no bias). Both solvers are compared to
    // ED and to each other on <sigma_x>, which both measure through
    // configuration-level (already-validated) estimators.
    //
    // Note on the wormhole transverse improved estimator: its
    // `PhysicalCorrelationS*` outputs are NOT used here. A controlled
    // single-mode ED comparison (2026-08-19) shows that family does not
    // reproduce ED transverse correlators on interacting models; it is
    // recorded as an open item in VALIDATION.md.
    for (lambda, tunnelling, beta, seed) in [
        (LAMBDA, TUNNELLING, BETA, 0x40F2_u64),
        (0.35, 0.5, 3.0, 0x40F3),
    ] {
        let (energies, vectors) =
            jacobi_eigensystem(multimode_matrix_with(lambda, tunnelling, 0.0));
        let dimension = energies.len();
        let mut sigma_x = vec![vec![0.0; dimension]; dimension];
        let encode = |n1: usize, n2: usize, spin: usize| 2 * (n1 * BOSON_CUTOFFS[1] + n2) + spin;
        for n1 in 0..BOSON_CUTOFFS[0] {
            for n2 in 0..BOSON_CUTOFFS[1] {
                sigma_x[encode(n1, n2, 0)][encode(n1, n2, 1)] = 1.0;
                sigma_x[encode(n1, n2, 1)][encode(n1, n2, 0)] = 1.0;
            }
        }
        let exact_sigma_x = thermal_expectation(&energies, &vectors, &sigma_x, beta);

        // Cluster solver on the multi-mode tabulated bath.
        let model = LongitudinalSpinBosonModel::with_default_quadrature(
            multimode_bath(),
            lambda,
            tunnelling,
            0.0,
        )
        .expect("model");
        let mut engine = ContinuousTimeClusterEngine::new(model);
        let mut worldline = LongitudinalWorldline::new(beta, 1).expect("worldline");
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
        for _ in 0..5_000 {
            engine.update(&mut worldline, &mut rng).expect("warmup");
        }
        let samples = 100_000usize;
        let mut cluster_sigma_x = Vec::with_capacity(samples);
        for _ in 0..samples {
            engine.update(&mut worldline, &mut rng).expect("update");
            cluster_sigma_x.push(2.0 * worldline.kink_count() as f64 / (beta * tunnelling));
        }
        let cluster_mean = cluster_sigma_x.iter().sum::<f64>() / samples as f64;
        let cluster_err = blocked_stderr(&cluster_sigma_x, 16);

        // Wormhole solver on the same physical Hamiltonian. Its Rabi
        // convention parametrizes the coupling as `(g/2) sigma_z (b + b†)`
        // with `lambda_wh = g^2/omega`, i.e. four times the cluster mass
        // for the same physical couplings
        // (`cluster lambda = g_phys^2/omega`).
        let mut params = Params::new();
        params.set("beta", beta);
        params.set("model", "rabi");
        params.set("bath", "tabulated");
        params.set("bath_omegas", "0.9,1.7");
        params.set("bath_weights", "1.0,0.6");
        params.set("lambda", 4.0 * lambda);
        params.set("tunnelling", tunnelling);
        params.set("validate_each_sweep", true);
        let run = RunConfig {
            thermalization_sweeps: 5_000,
            measurement_sweeps: 40_000,
            binsize: 100,
            base_seed: seed,
            ..Default::default()
        };
        let results = Scheduler::new(RayonBackend::new(1), run).run_one::<ImpurityQmc>(&params);
        let wormhole = results
            .get("MagnetizationSigmaZ")
            .expect("wormhole MagnetizationSigmaZ (= physical <sigma_x>)");

        for (label, value, stderr) in [
            ("cluster ⟨σx⟩", cluster_mean, cluster_err),
            ("wormhole ⟨σx⟩", wormhole.mean, wormhole.stderr),
        ] {
            let z = (value - exact_sigma_x) / stderr.max(1.0e-12);
            eprintln!(
                "lambda={lambda} Delta={tunnelling} beta={beta}: {label} = {value:.5} ±                  {stderr:.5} vs ED {exact_sigma_x:.5} (z = {z:.2})"
            );
            assert!(
                z.abs() < 4.0,
                "{label}: {value:.5} ± {stderr:.5} vs ED {exact_sigma_x:.5} (z = {z:.2})"
            );
        }
        let combined = (cluster_err.powi(2) + wormhole.stderr.powi(2)).sqrt();
        assert!(
            (cluster_mean - wormhole.mean).abs() < 4.0 * combined,
            "cross-solver ⟨σx⟩: cluster {cluster_mean:.5} vs wormhole {:.5} (4σ = {:.5})",
            wormhole.mean,
            4.0 * combined
        );
    }
}
