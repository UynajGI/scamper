//! Explicit finite-occupation Jaynes-Cummings and Rabi Hamiltonians.

use crate::impurity::spin_boson::occupation::basis::OccupationBasis;
use crate::impurity::ImpurityError;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OccupationModelKind {
    JaynesCummings,
    Rabi,
}

/// One explicitly retained cavity/boson mode.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CavityMode {
    pub omega: f64,
    pub coupling: f64,
    pub cutoff: usize,
}

impl CavityMode {
    pub fn new(omega: f64, coupling: f64, cutoff: usize) -> Result<Self, ImpurityError> {
        if !omega.is_finite() || omega <= 0.0 {
            return Err(ImpurityError::parameter(
                "omega",
                format!("must be finite and positive, got {omega}"),
            ));
        }
        if !coupling.is_finite() {
            return Err(ImpurityError::parameter(
                "coupling",
                format!("must be finite, got {coupling}"),
            ));
        }
        if cutoff == 0 {
            return Err(ImpurityError::parameter(
                "cutoff",
                "must retain at least the vacuum state",
            ));
        }
        Ok(Self {
            omega,
            coupling,
            cutoff,
        })
    }
}

/// Finite-dimensional explicit-boson impurity Hamiltonian.
///
/// Rabi convention:
/// `H = omega_q/2 sigma_z + sum_m omega_m n_m
///      + sum_m g_m sigma_x (a_m + a_m^dagger)`.
///
/// Jaynes-Cummings convention:
/// `H = omega_q/2 sigma_z + sum_m omega_m n_m
///      + sum_m g_m (sigma_+ a_m + sigma_- a_m^dagger)`.
#[derive(Debug, Clone, PartialEq)]
pub struct OccupationSpinBosonModel {
    kind: OccupationModelKind,
    spin_splitting: f64,
    modes: Vec<CavityMode>,
    basis: OccupationBasis,
}

impl OccupationSpinBosonModel {
    pub fn new(
        kind: OccupationModelKind,
        spin_splitting: f64,
        modes: Vec<CavityMode>,
    ) -> Result<Self, ImpurityError> {
        if !spin_splitting.is_finite() {
            return Err(ImpurityError::parameter("spin_splitting", "must be finite"));
        }
        if modes.is_empty() {
            return Err(ImpurityError::parameter(
                "modes",
                "at least one cavity mode is required",
            ));
        }
        let basis = OccupationBasis::new(modes.iter().map(|mode| mode.cutoff).collect())?;
        Ok(Self {
            kind,
            spin_splitting,
            modes,
            basis,
        })
    }

    pub fn rabi(spin_splitting: f64, modes: Vec<CavityMode>) -> Result<Self, ImpurityError> {
        Self::new(OccupationModelKind::Rabi, spin_splitting, modes)
    }

    pub fn jaynes_cummings(
        spin_splitting: f64,
        modes: Vec<CavityMode>,
    ) -> Result<Self, ImpurityError> {
        Self::new(OccupationModelKind::JaynesCummings, spin_splitting, modes)
    }

    #[inline]
    pub const fn kind(&self) -> OccupationModelKind {
        self.kind
    }
    #[inline]
    pub const fn spin_splitting(&self) -> f64 {
        self.spin_splitting
    }
    #[inline]
    pub fn modes(&self) -> &[CavityMode] {
        &self.modes
    }
    #[inline]
    pub const fn basis(&self) -> &OccupationBasis {
        &self.basis
    }

    /// Dense real symmetric Hamiltonian in the sigma-z/occupation basis.
    pub fn hamiltonian(&self) -> Vec<Vec<f64>> {
        let dimension = self.basis.dimension();
        let mut matrix = vec![vec![0.0; dimension]; dimension];
        for state in 0..dimension {
            let sigma_z = self.basis.spin(state).sigma_z();
            matrix[state][state] += 0.5 * self.spin_splitting * sigma_z;
            for (mode_index, mode) in self.modes.iter().enumerate() {
                let occupation = self.basis.occupation(state, mode_index);
                matrix[state][state] += mode.omega * occupation as f64;
                match self.kind {
                    OccupationModelKind::Rabi => {
                        let flipped = self.basis.flipped_spin(state);
                        if let Some(raised) = self.basis.shifted_state(flipped, mode_index, 1) {
                            add_symmetric(
                                &mut matrix,
                                state,
                                raised,
                                -mode.coupling * (occupation as f64 + 1.0).sqrt(),
                            );
                        }
                    }
                    OccupationModelKind::JaynesCummings => {
                        let spin_up = self.basis.spin(state).sigma_z() > 0.0;
                        if spin_up {
                            if let Some(raised) = self.basis.shifted_state(
                                self.basis.flipped_spin(state),
                                mode_index,
                                1,
                            ) {
                                add_symmetric(
                                    &mut matrix,
                                    state,
                                    raised,
                                    -mode.coupling * (occupation as f64 + 1.0).sqrt(),
                                );
                            }
                        }
                    }
                }
            }
        }
        matrix
    }
}

fn add_symmetric(matrix: &mut [Vec<f64>], left: usize, right: usize, value: f64) {
    matrix[left][right] += value;
    matrix[right][left] += value;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::impurity::spin_boson::occupation::basis::SpinState;

    fn mode(omega: f64, coupling: f64, cutoff: usize) -> CavityMode {
        CavityMode::new(omega, coupling, cutoff).expect("mode")
    }

    #[test]
    fn cavity_mode_rejects_nonpositive_omega() {
        assert!(CavityMode::new(0.0, 1.0, 3).is_err());
        assert!(CavityMode::new(-1.0, 1.0, 3).is_err());
        assert!(CavityMode::new(f64::NAN, 1.0, 3).is_err());
    }

    #[test]
    fn cavity_mode_rejects_zero_cutoff() {
        assert!(CavityMode::new(1.0, 1.0, 0).is_err());
    }

    #[test]
    fn model_construction_propagates_basis_dimension() {
        let model =
            OccupationSpinBosonModel::rabi(1.0, vec![mode(1.0, 0.5, 4)]).expect("rabi model");
        assert_eq!(model.kind(), OccupationModelKind::Rabi);
        assert_eq!(model.basis().dimension(), 8);
        assert_eq!(model.modes().len(), 1);
        assert!((model.spin_splitting() - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn jaynes_cummings_constructor_sets_kind() {
        let model = OccupationSpinBosonModel::jaynes_cummings(0.5, vec![mode(2.0, 0.3, 3)])
            .expect("jc model");
        assert_eq!(model.kind(), OccupationModelKind::JaynesCummings);
    }

    #[test]
    fn model_rejects_empty_mode_list() {
        assert!(OccupationSpinBosonModel::rabi(1.0, vec![]).is_err());
    }

    #[test]
    #[allow(clippy::needless_range_loop)]
    fn uncoupled_hamiltonian_is_diagonal_with_analytical_elements() {
        // H = (omega_q/2) sigma_z + omega n with no coupling: purely diagonal.
        let model = OccupationSpinBosonModel::rabi(2.0, vec![mode(1.5, 0.0, 3)]).expect("model");
        let matrix = model.hamiltonian();
        let basis = model.basis();
        for state in 0..basis.dimension() {
            for other in 0..basis.dimension() {
                if state == other {
                    continue;
                }
                assert_eq!(matrix[state][other], 0.0, "unexpected off-diagonal entry");
            }
        }
        // |down, n> at omega_q=2, omega=1.5: E = -1 + 1.5 n.
        for n in 0..basis.cutoff(0) {
            let state = basis.encode(SpinState::Down, &[n]).expect("encode");
            assert!((matrix[state][state] - (-1.0 + 1.5 * n as f64)).abs() < 1e-14);
        }
        // |up, n>: E = +1 + 1.5 n.
        for n in 0..basis.cutoff(0) {
            let state = basis.encode(SpinState::Up, &[n]).expect("encode");
            assert!((matrix[state][state] - (1.0 + 1.5 * n as f64)).abs() < 1e-14);
        }
    }

    #[test]
    #[allow(clippy::needless_range_loop)]
    fn hamiltonian_is_symmetric() {
        let model = OccupationSpinBosonModel::rabi(0.7, vec![mode(1.1, 0.4, 4)]).expect("model");
        let matrix = model.hamiltonian();
        let n = matrix.len();
        for i in 0..n {
            for j in (i + 1)..n {
                assert!((matrix[i][j] - matrix[j][i]).abs() < 1e-14);
            }
        }
    }

    #[test]
    fn rabi_coupling_contributes_at_largest_off_diagonal_scale() {
        // For Rabi the coupling connects |s, n> <-> |-s, n+1> with matrix
        // element -g sqrt(n+1). Verify that the ground-state off-diagonal
        // matches this value analytically.
        let g = 0.6;
        let model = OccupationSpinBosonModel::rabi(0.0, vec![mode(1.0, g, 3)]).expect("model");
        let matrix = model.hamiltonian();
        let basis = model.basis();
        let down0 = basis.encode(SpinState::Down, &[0]).expect("encode");
        let up1 = basis.encode(SpinState::Up, &[1]).expect("encode");
        let expected = -g * (1.0_f64).sqrt();
        assert!((matrix[down0][up1] - expected).abs() < 1e-14);
        assert!((matrix[up1][down0] - expected).abs() < 1e-14);
    }

    #[test]
    #[allow(clippy::needless_range_loop)]
    fn jaynes_cummings_only_couples_within_excitation_sectors() {
        // JC conserves total excitation number N = n + (spin up ? 1 : 0).
        // Here we check that |down, 0> is completely decoupled (no off-diagonal).
        let model =
            OccupationSpinBosonModel::jaynes_cummings(0.0, vec![mode(1.0, 0.5, 3)]).expect("jc");
        let matrix = model.hamiltonian();
        let basis = model.basis();
        let down0 = basis.encode(SpinState::Down, &[0]).expect("encode");
        for other in 0..basis.dimension() {
            if other == down0 {
                continue;
            }
            assert_eq!(matrix[down0][other], 0.0, "JC ground state should be dark");
        }
    }

    #[test]
    fn jaynes_cummings_ground_state_energy_is_zero() {
        // With the parameters above the |down, 0> state has zero diagonal energy.
        let model =
            OccupationSpinBosonModel::jaynes_cummings(0.0, vec![mode(1.0, 0.5, 3)]).expect("jc");
        let matrix = model.hamiltonian();
        let basis = model.basis();
        let down0 = basis.encode(SpinState::Down, &[0]).expect("encode");
        assert!(matrix[down0][down0].abs() < 1e-14);
    }
}
