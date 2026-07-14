//! Worldline and basis-aware estimators for spin-boson impurity models.

use carlo_rs::{CarloError, Evaluator};
use rand::Rng;
use rand::RngExt;

use crate::impurity::core::estimator::register_connected_susceptibility;
use crate::impurity::core::operators::{PhysicalAxis, SignedAxis};
use crate::impurity::spin_boson::model::ImpurityModel;
use crate::impurity::spin_boson::wormhole::configuration::WormholeConfiguration;
use crate::impurity::ImpurityError;

/// One set of scalar measurements from a wormhole configuration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImpurityObservables {
    /// Sample-basis imaginary-time averaged `sigma_z`.
    pub magnetization_sigma_z: f64,
    /// Sample-basis imaginary-time averaged `S_z = sigma_z/2`.
    pub magnetization_s_z: f64,
    pub magnetization_sigma_z_squared: f64,
    pub magnetization_sigma_z_fourth: f64,
    /// Raw sample `beta * m_Sz^2`; the connected response is registered as a
    /// Carlo.rs derived observable from `<m>` and `<m^2>`.
    pub susceptibility_z: f64,
    pub susceptibility_z_raw: f64,
    pub correlation_sigma_z_half: f64,
    pub correlation_s_z_half: f64,
    /// Physical signed axis represented by sampled `Z`.
    pub physical_axis_for_sampled_z: SignedAxis,
    /// Basis-corrected one-point value on that physical axis.
    pub physical_magnetization_s: f64,
    /// Same-axis correlation; the basis sign cancels.
    pub physical_correlation_s_half: f64,
    pub expansion_order: f64,
    pub diagonal_order: f64,
    pub offdiagonal_order: f64,
    /// Backward-compatible alias for `-n/beta` in the shifted expansion.
    pub shifted_interaction_energy: f64,
    /// Explicit name for the shifted expansion-order estimator.
    pub expansion_order_energy: f64,
    /// Shift-corrected spin-plus-coupling energy.
    pub spin_coupling_energy: f64,
}

/// Measure all built-in scalar observables.
pub fn measure_observables<R: Rng + ?Sized>(
    configuration: &WormholeConfiguration,
    model: &ImpurityModel,
    correlation_samples: usize,
    rng: &mut R,
) -> Result<ImpurityObservables, ImpurityError> {
    let magnetization_sigma_z = integrated_sigma_z(configuration, model)?;
    let magnetization_s_z = 0.5 * magnetization_sigma_z;
    let correlation_sigma_z_half = correlation_sigma_z(
        configuration,
        model,
        0.5 * configuration.beta(),
        correlation_samples,
        rng,
    )?;
    let expansion_order = configuration.expansion_order() as f64;
    let expansion_order_energy = -expansion_order / configuration.beta();
    let physical_axis_for_sampled_z = model
        .basis_transform()
        .physical_for_sampled(PhysicalAxis::Z);
    Ok(ImpurityObservables {
        magnetization_sigma_z,
        magnetization_s_z,
        magnetization_sigma_z_squared: magnetization_sigma_z * magnetization_sigma_z,
        magnetization_sigma_z_fourth: magnetization_sigma_z.powi(4),
        susceptibility_z: configuration.beta() * magnetization_s_z * magnetization_s_z,
        susceptibility_z_raw: configuration.beta() * magnetization_s_z * magnetization_s_z,
        correlation_sigma_z_half,
        correlation_s_z_half: 0.25 * correlation_sigma_z_half,
        physical_axis_for_sampled_z,
        physical_magnetization_s: f64::from(physical_axis_for_sampled_z.sign) * magnetization_s_z,
        physical_correlation_s_half: 0.25 * correlation_sigma_z_half,
        expansion_order,
        diagonal_order: configuration.diagonal_order() as f64,
        offdiagonal_order: configuration.offdiagonal_order() as f64,
        shifted_interaction_energy: expansion_order_energy,
        expansion_order_energy,
        spin_coupling_energy: model
            .corrected_spin_coupling_energy(expansion_order, configuration.beta()),
    })
}

/// Register connected susceptibilities after binned measurements have been
/// accumulated. This subtraction must be done at ensemble level, not on each
/// configuration separately.
pub fn register_impurity_evaluables(
    evaluator: &mut Evaluator,
    beta: f64,
    model: &ImpurityModel,
) -> Result<(), CarloError> {
    register_connected_susceptibility(
        evaluator,
        "ChiSampledZConnected",
        "SampledMagnetizationSz",
        "SampledM2Sz",
        beta,
    )?;
    // Preserve the historical identity-basis label when it is physically Z.
    let physical = model
        .basis_transform()
        .physical_for_sampled(PhysicalAxis::Z);
    let magnetization = format!("PhysicalMagnetizationS{}", physical.axis.label());
    let squared = format!("PhysicalM2S{}", physical.axis.label());
    let name = format!("ChiPhysical{}Connected", physical.axis.label());
    register_connected_susceptibility(evaluator, &name, &magnetization, &squared, beta)?;
    if physical.axis == PhysicalAxis::Z {
        register_connected_susceptibility(
            evaluator,
            "ChiZConnected",
            "MagnetizationSz",
            "M2Sz",
            beta,
        )?;
    }
    Ok(())
}

/// Imaginary-time average of sampled `sigma_z`.
pub fn integrated_sigma_z(
    configuration: &WormholeConfiguration,
    model: &ImpurityModel,
) -> Result<f64, ImpurityError> {
    if configuration.expansion_order() == 0 {
        return Ok(f64::from(configuration.empty_spin()));
    }
    let mut spin = configuration.spin_at(model, 0.0)?;
    let mut previous = 0.0;
    let mut total = 0.0;
    for (time, endpoint) in configuration.time_ordered_endpoints() {
        total += f64::from(spin) * (time - previous);
        spin = configuration.endpoint_outgoing_spin(endpoint, model)?;
        previous = time;
    }
    total += f64::from(spin) * (configuration.beta() - previous);
    Ok(total / configuration.beta())
}

/// Random-origin estimator of the sampled longitudinal correlation.
pub fn correlation_sigma_z<R: Rng + ?Sized>(
    configuration: &WormholeConfiguration,
    model: &ImpurityModel,
    delta_tau: f64,
    samples: usize,
    rng: &mut R,
) -> Result<f64, ImpurityError> {
    let sample_count = samples.max(1);
    let mut total = 0.0;
    for _ in 0..sample_count {
        let tau = rng.random::<f64>() * configuration.beta();
        let left = configuration.spin_at(model, tau)?;
        let right = configuration.spin_at(model, tau + delta_tau)?;
        total += f64::from(left * right);
    }
    Ok(total / sample_count as f64)
}

#[cfg(test)]
mod tests {
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256PlusPlus;

    use crate::impurity::spin_boson::bath::{Bath, SingleModeBath};

    use super::*;

    #[test]
    fn empty_worldline_estimators_are_exact() {
        let model = ImpurityModel::jaynes_cummings(
            Bath::SingleMode(SingleModeBath::new(1.0).expect("mode")),
            0.2,
            0.0,
            None,
        )
        .expect("model");
        let configuration = WormholeConfiguration::new(4.0, -1).expect("configuration");
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(5);
        let observables =
            measure_observables(&configuration, &model, 8, &mut rng).expect("observables");
        assert!((observables.magnetization_sigma_z + 1.0).abs() < f64::EPSILON);
        assert!((observables.correlation_sigma_z_half - 1.0).abs() < f64::EPSILON);
        assert!(observables.expansion_order.abs() < f64::EPSILON);
    }

    #[test]
    fn rotated_rabi_maps_sampled_z_to_physical_x() {
        let model = ImpurityModel::rotated_impurity(
            Bath::SingleMode(SingleModeBath::new(1.0).expect("mode")),
            0.2,
            0.4,
            Some(0.5),
        )
        .expect("model");
        let configuration = WormholeConfiguration::new(3.0, 1).expect("configuration");
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(2);
        let observables = measure_observables(&configuration, &model, 4, &mut rng).unwrap();
        assert_eq!(
            observables.physical_axis_for_sampled_z.axis,
            PhysicalAxis::X
        );
        assert!((observables.physical_magnetization_s - 0.5).abs() < 1e-14);
    }

    // Interacting exact-diagonalization comparison.  The original single-mode
    // Rabi Hamiltonian and the sampled rotated Hamiltonian are related by the
    // declared spin-basis rotation, so their spectra and mapped thermal
    // one-point functions must agree at every boson cutoff.
    #[test]
    fn interacting_single_mode_rabi_matches_rotated_ed() {
        let boson_states = 5;
        let omega = 1.3;
        let tunnelling = 0.7;
        let coupling = 0.45;
        let beta = 2.1;
        let original = rabi_matrix(boson_states, omega, tunnelling, coupling, false);
        let rotated = rabi_matrix(boson_states, omega, tunnelling, coupling, true);
        let original_sx = spin_operator(boson_states, PhysicalAxis::X);
        let rotated_sz = spin_operator(boson_states, PhysicalAxis::Z);
        let (e_original, v_original) = jacobi_eigensystem(original);
        let (e_rotated, v_rotated) = jacobi_eigensystem(rotated);
        for (left, right) in e_original.iter().zip(&e_rotated) {
            assert!((left - right).abs() < 1e-10);
        }
        let physical_sx = thermal_expectation(&e_original, &v_original, &original_sx, beta);
        let sampled_sz = thermal_expectation(&e_rotated, &v_rotated, &rotated_sz, beta);
        assert!((physical_sx - sampled_sz).abs() < 1e-10);
    }

    fn index(n: usize, spin: usize) -> usize {
        2 * n + spin
    }

    fn rabi_matrix(
        boson_states: usize,
        omega: f64,
        tunnelling: f64,
        coupling: f64,
        rotated: bool,
    ) -> Vec<Vec<f64>> {
        let dim = 2 * boson_states;
        let mut h = vec![vec![0.0; dim]; dim];
        for n in 0..boson_states {
            for spin in 0..2 {
                let i = index(n, spin);
                let sz = if spin == 0 { -0.5 } else { 0.5 };
                h[i][i] += omega * n as f64;
                if rotated {
                    h[i][i] += -tunnelling * sz;
                } else {
                    let j = index(n, 1 - spin);
                    h[i][j] += -0.5 * tunnelling;
                }
                if n + 1 < boson_states {
                    let amplitude = (n as f64 + 1.0).sqrt();
                    if rotated {
                        let j = index(n + 1, 1 - spin);
                        h[i][j] += -0.5 * coupling * amplitude;
                        h[j][i] += -0.5 * coupling * amplitude;
                    } else {
                        let j = index(n + 1, spin);
                        h[i][j] += coupling * sz * amplitude;
                        h[j][i] += coupling * sz * amplitude;
                    }
                }
            }
        }
        h
    }

    fn spin_operator(boson_states: usize, axis: PhysicalAxis) -> Vec<Vec<f64>> {
        let dim = 2 * boson_states;
        let mut operator = vec![vec![0.0; dim]; dim];
        for n in 0..boson_states {
            match axis {
                PhysicalAxis::Z => {
                    operator[index(n, 0)][index(n, 0)] = -0.5;
                    operator[index(n, 1)][index(n, 1)] = 0.5;
                }
                PhysicalAxis::X => {
                    operator[index(n, 0)][index(n, 1)] = 0.5;
                    operator[index(n, 1)][index(n, 0)] = 0.5;
                }
                PhysicalAxis::Y => unreachable!("real ED helper only needs X and Z"),
            }
        }
        operator
    }

    #[allow(clippy::needless_range_loop)]
    fn jacobi_eigensystem(mut matrix: Vec<Vec<f64>>) -> (Vec<f64>, Vec<Vec<f64>>) {
        let n = matrix.len();
        let mut vectors = vec![vec![0.0; n]; n];
        for (i, row) in vectors.iter_mut().enumerate() {
            row[i] = 1.0;
        }
        for _ in 0..(100 * n * n) {
            let mut p = 0;
            let mut q = 1;
            let mut largest = 0.0;
            for i in 0..n {
                for j in (i + 1)..n {
                    if matrix[i][j].abs() > largest {
                        largest = matrix[i][j].abs();
                        p = i;
                        q = j;
                    }
                }
            }
            if largest < 1e-13 {
                break;
            }
            let angle = 0.5 * (2.0 * matrix[p][q]).atan2(matrix[q][q] - matrix[p][p]);
            let (s, c) = angle.sin_cos();
            for k in 0..n {
                if k != p && k != q {
                    let mkp = matrix[k][p];
                    let mkq = matrix[k][q];
                    matrix[k][p] = c * mkp - s * mkq;
                    matrix[p][k] = matrix[k][p];
                    matrix[k][q] = s * mkp + c * mkq;
                    matrix[q][k] = matrix[k][q];
                }
            }
            let app = matrix[p][p];
            let aqq = matrix[q][q];
            let apq = matrix[p][q];
            matrix[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
            matrix[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
            matrix[p][q] = 0.0;
            matrix[q][p] = 0.0;
            for row in &mut vectors {
                let vkp = row[p];
                let vkq = row[q];
                row[p] = c * vkp - s * vkq;
                row[q] = s * vkp + c * vkq;
            }
        }
        let mut order: Vec<_> = (0..n).collect();
        order.sort_by(|&i, &j| matrix[i][i].total_cmp(&matrix[j][j]));
        let values = order.iter().map(|&i| matrix[i][i]).collect();
        let sorted_vectors = (0..n)
            .map(|row| order.iter().map(|&column| vectors[row][column]).collect())
            .collect();
        (values, sorted_vectors)
    }

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
            for i in 0..operator.len() {
                for j in 0..operator.len() {
                    expectation += vectors[i][state] * operator[i][j] * vectors[j][state];
                }
            }
            numerator += weight * expectation;
            partition += weight;
        }
        numerator / partition
    }
}
