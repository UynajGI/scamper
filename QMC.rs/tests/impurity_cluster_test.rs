use qmc_rs::{
    Bath, ContinuousTimeClusterEngine, LongitudinalSpinBosonModel, LongitudinalWorldline,
    PowerLawBath, SingleModeBath,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

fn free_model(tunnelling: f64, bias: f64) -> LongitudinalSpinBosonModel {
    LongitudinalSpinBosonModel::with_default_quadrature(
        Bath::SingleMode(SingleModeBath::new(1.0).expect("mode")),
        0.0,
        tunnelling,
        bias,
    )
    .expect("model")
}

#[test]
fn free_two_level_cluster_solver_matches_exact_thermodynamics() {
    let beta: f64 = 3.0;
    let tunnelling: f64 = 1.1;
    let bias: f64 = 0.7;
    let gamma: f64 = 0.5 * tunnelling;
    let field: f64 = 0.5 * bias;
    let gap: f64 = gamma.hypot(field);
    let exact_sigma_z: f64 = -field / gap * (beta * gap).tanh();
    let exact_kinks: f64 = beta * gamma * gamma / gap * (beta * gap).tanh();

    let mut engine = ContinuousTimeClusterEngine::new(free_model(tunnelling, bias));
    let mut worldline = LongitudinalWorldline::new(beta, 1).expect("worldline");
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x51A7_2026);
    for _ in 0..4_000 {
        engine.update(&mut worldline, &mut rng).expect("warmup");
    }

    let samples = 80_000usize;
    let mut magnetization = 0.0;
    let mut improved_magnetization = 0.0;
    let mut kinks = 0.0;
    for _ in 0..samples {
        let report = engine.update(&mut worldline, &mut rng).expect("update");
        magnetization += worldline.integrated_sigma_z();
        improved_magnetization += report.improved_magnetization_sigma_z;
        kinks += worldline.kink_count() as f64;
    }
    magnetization /= samples as f64;
    improved_magnetization /= samples as f64;
    kinks /= samples as f64;

    assert!(
        (magnetization - exact_sigma_z).abs() < 0.025,
        "sampled={magnetization}, exact={exact_sigma_z}"
    );
    assert!(
        (improved_magnetization - exact_sigma_z).abs() < 0.02,
        "improved={improved_magnetization}, exact={exact_sigma_z}"
    );
    assert!(
        (kinks - exact_kinks).abs() < 0.04,
        "sampled={kinks}, exact={exact_kinks}"
    );
}

#[test]
fn unbiased_free_correlation_matches_exact_imaginary_time_result() {
    let beta: f64 = 4.0;
    let tunnelling: f64 = 0.9;
    let gamma: f64 = 0.5 * tunnelling;
    let delta_tau: f64 = 0.5 * beta;
    let exact: f64 = ((beta - 2.0 * delta_tau) * gamma).cosh() / (beta * gamma).cosh();

    let mut engine = ContinuousTimeClusterEngine::new(free_model(tunnelling, 0.0));
    let mut worldline = LongitudinalWorldline::new(beta, -1).expect("worldline");
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x0C01_1EC7);
    for _ in 0..4_000 {
        engine.update(&mut worldline, &mut rng).expect("warmup");
    }

    let samples = 80_000usize;
    let correlation = (0..samples)
        .map(|_| {
            engine.update(&mut worldline, &mut rng).expect("update");
            worldline.correlation_sigma_z(delta_tau)
        })
        .sum::<f64>()
        / samples as f64;
    assert!(
        (correlation - exact).abs() < 0.025,
        "sampled={correlation}, exact={exact}"
    );
}

#[test]
fn sub_ohmic_retarded_kernel_is_positive_and_periodic() {
    let bath = Bath::PowerLaw(PowerLawBath::new(0.4, 2.0).expect("bath"));
    let model = LongitudinalSpinBosonModel::new(bath, 0.3, 0.2, 0.0, 128).expect("model");
    let beta = 8.0;
    for delta in [0.0, 0.1, 1.0, 3.5, 7.9] {
        let left = model.kernel().value(beta, delta).expect("kernel");
        let right = model.kernel().value(beta, beta - delta).expect("kernel");
        assert!(left.is_finite() && left > 0.0);
        assert!((left - right).abs() < 1e-12);
    }
}

#[test]
fn every_cluster_update_preserves_periodic_worldline_invariants() {
    let bath = Bath::SingleMode(SingleModeBath::new(1.3).expect("mode"));
    let model =
        LongitudinalSpinBosonModel::with_default_quadrature(bath, 0.8, 1.4, 0.2).expect("model");
    let mut engine = ContinuousTimeClusterEngine::new(model);
    engine.set_validate_each_sweep(true);
    let mut worldline = LongitudinalWorldline::new(6.0, 1).expect("worldline");
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(9921);
    for _ in 0..10_000 {
        engine.update(&mut worldline, &mut rng).expect("update");
        assert_eq!(worldline.kink_count() % 2, 0);
        worldline.validate().expect("valid worldline");
    }
}

#[test]
fn interacting_single_mode_cluster_matches_truncated_ed() {
    let beta = 2.4;
    let omega = 1.3;
    let tunnelling = 0.9;
    let bias = 0.35;
    let coupling_sigma = 0.28;
    let boson_states = 12;

    let hamiltonian =
        longitudinal_single_mode_matrix(boson_states, omega, tunnelling, bias, coupling_sigma);
    let sigma_z = pauli_operator(boson_states, false);
    let sigma_x = pauli_operator(boson_states, true);
    let (energies, vectors) = jacobi_eigensystem(hamiltonian);
    let exact_sigma_z = thermal_expectation(&energies, &vectors, &sigma_z, beta);
    let exact_sigma_x = thermal_expectation(&energies, &vectors, &sigma_x, beta);

    let bath = Bath::SingleMode(SingleModeBath::new(omega).expect("mode"));
    let lambda = coupling_sigma * coupling_sigma / omega;
    let model = LongitudinalSpinBosonModel::with_default_quadrature(bath, lambda, tunnelling, bias)
        .expect("model");
    let mut engine = ContinuousTimeClusterEngine::new(model);
    let mut worldline = LongitudinalWorldline::new(beta, 1).expect("worldline");
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0xEDC1_057E);
    for _ in 0..5_000 {
        engine.update(&mut worldline, &mut rng).expect("warmup");
    }

    let samples = 100_000usize;
    let mut sampled_sigma_z = 0.0;
    let mut sampled_kinks = 0.0;
    for _ in 0..samples {
        let report = engine.update(&mut worldline, &mut rng).expect("update");
        sampled_sigma_z += report.improved_magnetization_sigma_z;
        sampled_kinks += worldline.kink_count() as f64;
    }
    sampled_sigma_z /= samples as f64;
    let sampled_sigma_x = 2.0 * sampled_kinks / (samples as f64 * beta * tunnelling);

    assert!(
        (sampled_sigma_z - exact_sigma_z).abs() < 0.035,
        "cluster sigma_z={sampled_sigma_z}, ED={exact_sigma_z}"
    );
    assert!(
        (sampled_sigma_x - exact_sigma_x).abs() < 0.04,
        "cluster sigma_x={sampled_sigma_x}, ED={exact_sigma_x}"
    );
}

fn state_index(boson: usize, spin: usize) -> usize {
    2 * boson + spin
}

fn longitudinal_single_mode_matrix(
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
                let matrix_element = coupling_sigma * sigma_z * (boson as f64 + 1.0).sqrt();
                hamiltonian[index][raised] += matrix_element;
                hamiltonian[raised][index] += matrix_element;
            }
        }
    }
    hamiltonian
}

fn pauli_operator(boson_states: usize, transverse: bool) -> Vec<Vec<f64>> {
    let dimension = 2 * boson_states;
    let mut operator = vec![vec![0.0; dimension]; dimension];
    for boson in 0..boson_states {
        if transverse {
            operator[state_index(boson, 0)][state_index(boson, 1)] = 1.0;
            operator[state_index(boson, 1)][state_index(boson, 0)] = 1.0;
        } else {
            operator[state_index(boson, 0)][state_index(boson, 0)] = -1.0;
            operator[state_index(boson, 1)][state_index(boson, 1)] = 1.0;
        }
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
    for (state, energy) in energies.iter().enumerate() {
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
