//! Exact transfer-matrix occupation worldline sampler.
//!
//! The link weight is the full finite-basis propagator
//! `T = exp(-beta H / slices)`. Therefore `Tr(T^slices)=Tr(exp(-beta H))`
//! exactly: the slice count controls autocorrelation and storage, not a
//! Trotter approximation.

use rand::Rng;
use rand::RngExt;

use crate::impurity::spin_boson::occupation::basis::OccupationBasis;
use crate::impurity::spin_boson::occupation::model::OccupationSpinBosonModel;
use crate::impurity::spin_boson::occupation::transfer::{multiply, SymmetricEigensystem};
use crate::impurity::ImpurityError;

#[derive(Debug, Clone, PartialEq)]
pub struct OccupationObservables {
    pub energy: f64,
    pub sigma_z: f64,
    pub sigma_x: f64,
    pub total_boson_number: f64,
    pub mode_occupations: Vec<f64>,
    pub mode_number_squared: Vec<f64>,
    /// Per-mode equal-time factorial moment `<n(n-1)>`.
    pub mode_factorial_moments: Vec<f64>,
    /// Per-mode equal-time photon coherence `g2(0)`; zero when `<n>=0`.
    pub mode_g2_zero: Vec<f64>,
    pub mode_cross_correlations: Vec<f64>,
    pub spin_boson_covariance_z_n: Vec<f64>,
    pub parity: f64,
    /// Purity of the exact finite-cutoff thermal reduced spin density matrix.
    pub reduced_spin_purity: f64,
    /// Variance of the bosonic quadrature x = a + a†, i.e. ⟨x²⟩_thermal.
    pub quadrature_variance: f64,
    /// Fourth moment of the bosonic quadrature, i.e. ⟨x⁴⟩_thermal.
    pub quadrature_fourth_moment: f64,
}

#[derive(Debug, Clone)]
pub struct OccupationWorldlineSampler {
    model: OccupationSpinBosonModel,
    beta: f64,
    slices: usize,
    states: Vec<usize>,
    transfer: Vec<Vec<f64>>,
    h_transfer: Vec<Vec<f64>>,
    sx_transfer: Vec<Vec<f64>>,
    x2_transfer: Vec<Vec<f64>>,
    x4_transfer: Vec<Vec<f64>>,
    eigensystem: SymmetricEigensystem,
    transfer_powers: Vec<Vec<Vec<f64>>>,
    reduced_spin_purity: f64,
    sweep_attempts: u64,
    sweep_changes: u64,
}

impl OccupationWorldlineSampler {
    const MAX_DENSE_DIMENSION: usize = 512;

    pub fn new(
        model: OccupationSpinBosonModel,
        beta: f64,
        slices: usize,
        initial_state: usize,
    ) -> Result<Self, ImpurityError> {
        if !beta.is_finite() || beta <= 0.0 {
            return Err(ImpurityError::parameter(
                "beta",
                "must be finite and positive",
            ));
        }
        if slices < 2 {
            return Err(ImpurityError::parameter("slices", "must be at least 2"));
        }
        if initial_state >= model.basis().dimension() {
            return Err(ImpurityError::InvalidConfiguration(
                "initial occupation state is outside the basis".into(),
            ));
        }
        if model.basis().dimension() > Self::MAX_DENSE_DIMENSION {
            return Err(ImpurityError::parameter(
                "cutoffs",
                format!("dense occupation solver dimension {} exceeds safety limit {}; use smaller cutoffs/fewer modes or a future sparse SSE backend", model.basis().dimension(), Self::MAX_DENSE_DIMENSION),
            ));
        }
        let hamiltonian = model.hamiltonian();
        let eigensystem = SymmetricEigensystem::diagonalize(hamiltonian.clone())?;
        let dt = beta / slices as f64;
        let ground = eigensystem.values[0];
        let transfer = eigensystem.matrix_function(|energy| (-dt * (energy - ground)).exp());
        for row in &transfer {
            for &weight in row {
                positive_weight(weight)?;
            }
        }
        let mut sx = vec![vec![0.0; model.basis().dimension()]; model.basis().dimension()];
        for state in 0..model.basis().dimension() {
            sx[state][model.basis().flipped_spin(state)] = 1.0;
        }
        // x² = (a+a†)² in the occupation basis (spin-diagonal).
        // x²|n,s⟩ = (2n+1)|n,s⟩ + √((n+1)(n+2))|n+2,s⟩ + √(n(n-1))|n-2,s⟩.
        let dim = model.basis().dimension();
        let cutoff = model.basis().boson_dimension();
        let mut x2 = vec![vec![0.0; dim]; dim];
        #[allow(clippy::needless_range_loop)]
        for state in 0..dim {
            let n = model.basis().occupation(state, 0) as f64;
            let s = state & 1;
            x2[state][state] = 2.0 * n + 1.0;
            let np2 = n as usize + 2;
            if np2 < cutoff {
                let target = 2 * np2 + s;
                let amp = ((n + 1.0) * (n + 2.0)).sqrt();
                x2[state][target] = amp;
                x2[target][state] = amp;
            }
        }
        let h_transfer = multiply(&hamiltonian, &transfer);
        let sx_transfer = multiply(&sx, &transfer);
        let x2_transfer = multiply(&x2, &transfer);
        let x4 = multiply(&x2, &x2);
        let x4_transfer = multiply(&x4, &transfer);
        let rho = eigensystem.thermal_density_matrix(beta);
        let reduced_spin_purity = reduced_spin_purity(model.basis(), &rho);
        let mut transfer_powers = Vec::with_capacity(slices + 1);
        transfer_powers.push(identity(model.basis().dimension()));
        for power in 1..=slices {
            transfer_powers.push(multiply(&transfer_powers[power - 1], &transfer));
        }
        Ok(Self {
            model,
            beta,
            slices,
            states: vec![initial_state; slices],
            transfer,
            h_transfer,
            sx_transfer,
            x2_transfer,
            x4_transfer,
            eigensystem,
            transfer_powers,
            reduced_spin_purity,
            sweep_attempts: 0,
            sweep_changes: 0,
        })
    }

    #[inline]
    pub const fn model(&self) -> &OccupationSpinBosonModel {
        &self.model
    }
    #[inline]
    pub const fn beta(&self) -> f64 {
        self.beta
    }
    #[inline]
    pub const fn slices(&self) -> usize {
        self.slices
    }
    #[inline]
    pub fn states(&self) -> &[usize] {
        &self.states
    }
    pub fn acceptance_fraction(&self) -> f64 {
        if self.sweep_attempts == 0 {
            0.0
        } else {
            self.sweep_changes as f64 / self.sweep_attempts as f64
        }
    }

    /// Draw a complete closed occupation worldline from its exact conditional
    /// distribution. This bridge heat bath can move between conserved JC
    /// excitation sectors and Rabi parity sectors, unlike a purely local
    /// single-slice update.
    pub fn sweep<R: Rng + ?Sized>(&mut self, rng: &mut R) -> Result<(), ImpurityError> {
        let dimension = self.model.basis().dimension();
        let full = &self.transfer_powers[self.slices];
        let total = (0..dimension).try_fold(0.0, |sum, state| {
            Ok::<_, ImpurityError>(sum + positive_weight(full[state][state])?)
        })?;
        if total <= 0.0 {
            return Err(ImpurityError::InvalidConfiguration(
                "zero finite-cutoff partition weight".into(),
            ));
        }
        let mut draw = rng.random::<f64>() * total;
        let mut first = 0usize;
        for (state, full_row) in full.iter().enumerate().take(dimension) {
            draw -= positive_weight(full_row[state])?;
            if draw <= 0.0 {
                first = state;
                break;
            }
        }
        let old = self.states.clone();
        self.states[0] = first;
        for index in 1..self.slices {
            let previous = self.states[index - 1];
            let remaining = self.slices - index;
            let mut normalization = 0.0;
            for candidate in 0..dimension {
                normalization += positive_weight(
                    self.transfer[previous][candidate]
                        * self.transfer_powers[remaining][candidate][first],
                )?;
            }
            if normalization <= 0.0 {
                return Err(ImpurityError::InvalidConfiguration(
                    "zero bridge conditional weight".into(),
                ));
            }
            let mut bridge_draw = rng.random::<f64>() * normalization;
            let mut chosen = 0usize;
            for candidate in 0..dimension {
                bridge_draw -= positive_weight(
                    self.transfer[previous][candidate]
                        * self.transfer_powers[remaining][candidate][first],
                )?;
                if bridge_draw <= 0.0 {
                    chosen = candidate;
                    break;
                }
            }
            self.states[index] = chosen;
        }
        self.sweep_attempts = self.sweep_attempts.saturating_add(self.slices as u64);
        self.sweep_changes = self
            .sweep_changes
            .saturating_add(old.iter().zip(&self.states).filter(|(a, b)| a != b).count() as u64);
        Ok(())
    }

    pub fn measure(&self) -> Result<OccupationObservables, ImpurityError> {
        let basis = self.model.basis();
        let mut sigma_z = 0.0;
        let mut total_boson = 0.0;
        let mut mode_occupations = vec![0.0; basis.modes()];
        let mut mode_n2 = vec![0.0; basis.modes()];
        let mut mode_factorial = vec![0.0; basis.modes()];
        let mut cross = vec![0.0; basis.modes() * basis.modes()];
        let mut zn = vec![0.0; basis.modes()];
        let mut parity = 0.0;
        let mut energy = 0.0;
        let mut sigma_x = 0.0;
        let mut quadrature_variance = 0.0;
        let mut quadrature_fourth = 0.0;
        for index in 0..self.slices {
            let state = self.states[index];
            let next = self.states[(index + 1) % self.slices];
            let z = basis.spin(state).sigma_z();
            sigma_z += z;
            let denominator = self.transfer[state][next];
            if denominator <= 0.0 {
                return Err(ImpurityError::InvalidConfiguration(
                    "sampled a non-positive transfer link".into(),
                ));
            }
            energy += self.h_transfer[state][next] / denominator;
            sigma_x += self.sx_transfer[state][next] / denominator;
            quadrature_variance += self.x2_transfer[state][next] / denominator;
            quadrature_fourth += self.x4_transfer[state][next] / denominator;
            let mut total_n = 0usize;
            for mode in 0..basis.modes() {
                let n = basis.occupation(state, mode) as f64;
                mode_occupations[mode] += n;
                mode_n2[mode] += n * n;
                mode_factorial[mode] += n * (n - 1.0);
                zn[mode] += z * n;
                total_n += n as usize;
            }
            for left in 0..basis.modes() {
                for right in 0..basis.modes() {
                    cross[left * basis.modes() + right] += basis.occupation(state, left) as f64
                        * basis.occupation(state, right) as f64;
                }
            }
            total_boson += total_n as f64;
            parity += z * if total_n.is_multiple_of(2) { 1.0 } else { -1.0 };
        }
        let inv = 1.0 / self.slices as f64;
        sigma_z *= inv;
        total_boson *= inv;
        energy *= inv;
        sigma_x *= inv;
        quadrature_variance *= inv;
        quadrature_fourth *= inv;
        parity *= inv;
        for value in &mut mode_occupations {
            *value *= inv;
        }
        for value in &mut mode_n2 {
            *value *= inv;
        }
        for value in &mut mode_factorial {
            *value *= inv;
        }
        for value in &mut cross {
            *value *= inv;
        }
        // Connected spin-boson covariance.
        for mode in 0..basis.modes() {
            let raw = (0..self.slices)
                .map(|i| {
                    basis.spin(self.states[i]).sigma_z()
                        * basis.occupation(self.states[i], mode) as f64
                })
                .sum::<f64>()
                * inv;
            zn[mode] = raw - sigma_z * mode_occupations[mode];
        }
        let mode_g2_zero = mode_factorial
            .iter()
            .zip(&mode_occupations)
            .map(|(&factorial, &mean)| {
                if mean > 0.0 {
                    factorial / (mean * mean)
                } else {
                    0.0
                }
            })
            .collect();
        Ok(OccupationObservables {
            energy,
            sigma_z,
            sigma_x,
            total_boson_number: total_boson,
            mode_occupations,
            mode_number_squared: mode_n2,
            mode_factorial_moments: mode_factorial,
            mode_g2_zero,
            mode_cross_correlations: cross,
            spin_boson_covariance_z_n: zn,
            parity,
            quadrature_variance,
            quadrature_fourth_moment: quadrature_fourth,
            reduced_spin_purity: self.reduced_spin_purity,
        })
    }

    pub fn exact_partition_function(&self) -> f64 {
        let ground = self.eigensystem.values[0];
        self.eigensystem
            .values
            .iter()
            .map(|&e| (-self.beta * (e - ground)).exp())
            .sum()
    }
}

fn identity(dimension: usize) -> Vec<Vec<f64>> {
    let mut matrix = vec![vec![0.0; dimension]; dimension];
    for (index, row) in matrix.iter_mut().enumerate() {
        row[index] = 1.0;
    }
    matrix
}

fn positive_weight(value: f64) -> Result<f64, ImpurityError> {
    let tolerance = 1e-12;
    if value < -tolerance {
        Err(ImpurityError::InvalidConfiguration(format!("non-stoquastic transfer weight {value}; choose a sign gauge with non-positive off-diagonal Hamiltonian elements")))
    } else {
        Ok(value.max(0.0))
    }
}

fn reduced_spin_purity(basis: &OccupationBasis, rho: &[Vec<f64>]) -> f64 {
    let mut reduced = [[0.0; 2]; 2];
    for boson in 0..basis.boson_dimension() {
        for left in 0..2 {
            for right in 0..2 {
                reduced[left][right] += rho[2 * boson + left][2 * boson + right];
            }
        }
    }
    reduced[0][0] * reduced[0][0]
        + reduced[1][1] * reduced[1][1]
        + 2.0 * reduced[0][1] * reduced[1][0]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::impurity::spin_boson::occupation::model::{
        CavityMode, OccupationModelKind, OccupationSpinBosonModel,
    };
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256PlusPlus;

    fn mode(omega: f64, coupling: f64, cutoff: usize) -> CavityMode {
        CavityMode::new(omega, coupling, cutoff).expect("mode")
    }

    /// Build a sampler, warm it up, and return it ready for production sampling.
    fn warmed(
        model: OccupationSpinBosonModel,
        beta: f64,
        slices: usize,
        seed: u64,
    ) -> OccupationWorldlineSampler {
        let dimension = model.basis().dimension();
        let mut sampler = OccupationWorldlineSampler::new(model, beta, slices, 0).expect("sampler");
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
        for _ in 0..2000 {
            sampler.sweep(&mut rng).expect("warmup");
        }
        // Sampler must only visit states inside the basis.
        for &state in sampler.states() {
            assert!(state < dimension, "sampler visited out-of-basis state");
        }
        sampler
    }

    /// Exact thermal <O> = Tr(O exp(-beta H)) / Tr(exp(-beta H)) via the
    /// model's own dense Hamiltonian and a Jacobi eigensystem.
    #[allow(clippy::needless_range_loop)]
    fn thermal_expectation(
        model: &OccupationSpinBosonModel,
        operator: &[Vec<f64>],
        beta: f64,
    ) -> f64 {
        let hamiltonian = model.hamiltonian();
        let eigen = SymmetricEigensystem::diagonalize(hamiltonian).expect("diagonalize");
        let ground = eigen.values[0];
        let mut numerator = 0.0;
        let mut partition = 0.0;
        for (state, &energy) in eigen.values.iter().enumerate() {
            let weight = (-beta * (energy - ground)).exp();
            let mut expectation = 0.0;
            for row in 0..operator.len() {
                for col in 0..operator.len() {
                    expectation +=
                        eigen.vectors[row][state] * operator[row][col] * eigen.vectors[col][state];
                }
            }
            numerator += weight * expectation;
            partition += weight;
        }
        numerator / partition
    }

    fn sigma_z_operator(dimension: usize) -> Vec<Vec<f64>> {
        let mut operator = vec![vec![0.0; dimension]; dimension];
        for (state, row) in operator.iter_mut().enumerate() {
            row[state] = if state & 1 == 0 { -1.0 } else { 1.0 };
        }
        operator
    }

    fn sigma_x_operator(basis: &OccupationBasis) -> Vec<Vec<f64>> {
        let dimension = basis.dimension();
        let mut operator = vec![vec![0.0; dimension]; dimension];
        for state in 0..dimension {
            operator[state][basis.flipped_spin(state)] = 1.0;
        }
        operator
    }

    #[test]
    fn rejects_non_positive_beta() {
        let model = OccupationSpinBosonModel::rabi(1.0, vec![mode(1.0, 0.0, 2)]).expect("model");
        assert!(OccupationWorldlineSampler::new(model, 0.0, 4, 0).is_err());
        assert!(OccupationWorldlineSampler::new(
            OccupationSpinBosonModel::rabi(1.0, vec![mode(1.0, 0.0, 2)]).expect("model"),
            -1.0,
            4,
            0,
        )
        .is_err());
    }

    #[test]
    fn rejects_too_few_slices() {
        let model = OccupationSpinBosonModel::rabi(1.0, vec![mode(1.0, 0.0, 2)]).expect("model");
        assert!(OccupationWorldlineSampler::new(model, 1.0, 1, 0).is_err());
    }

    #[test]
    fn rejects_out_of_basis_initial_state() {
        let model = OccupationSpinBosonModel::rabi(1.0, vec![mode(1.0, 0.0, 2)]).expect("model");
        // dimension is 4; state 4 is out of range.
        assert!(OccupationWorldlineSampler::new(model, 1.0, 4, 4).is_err());
    }

    #[test]
    fn exact_partition_function_matches_analytic_two_level_system() {
        // No coupling: H = (omega_q/2) sigma_z + omega n. With omega_q = 2 and
        // omega = 1, cutoff 1 (only vacuum): the two spin states have energies
        // {-1, +1}, so the ground-referenced partition function is
        // 1 + exp(-beta * gap) = 1 + exp(-2 beta).
        let model = OccupationSpinBosonModel::rabi(2.0, vec![mode(1.0, 0.0, 1)]).expect("model");
        let beta = 1.5;
        let sampler = OccupationWorldlineSampler::new(model, beta, 4, 0).expect("sampler");
        let exact = 1.0 + (-2.0 * beta).exp();
        assert!(
            (sampler.exact_partition_function() - exact).abs() < 1e-10,
            "Z = {}, expected {}",
            sampler.exact_partition_function(),
            exact
        );
    }

    #[test]
    fn free_spin_magnetization_matches_exact_tanh() {
        // Two-level system with omega_q = 2 epsilon: <sigma_z> = -tanh(beta epsilon).
        let epsilon = 0.5;
        let beta = 2.0;
        let model =
            OccupationSpinBosonModel::rabi(2.0 * epsilon, vec![mode(1.0, 0.0, 1)]).expect("model");
        let mut sampler = warmed(model, beta, 6, 123457);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(123457);
        let samples = 40_000usize;
        let mut mean_sigma_z = 0.0;
        for _ in 0..samples {
            sampler.sweep(&mut rng).expect("sweep");
            mean_sigma_z += sampler.measure().expect("measure").sigma_z;
        }
        mean_sigma_z /= samples as f64;
        let exact = -(beta * epsilon).tanh();
        assert!(
            (mean_sigma_z - exact).abs() < 0.02,
            "sampled sigma_z = {mean_sigma_z}, exact = {exact}"
        );
    }

    #[test]
    fn rabi_thermal_observables_match_exact_diagonalization() {
        let omega_q = 0.8;
        let omega = 1.3;
        let g = 0.45;
        let beta = 2.5;
        let cutoff = 6;
        let model =
            OccupationSpinBosonModel::rabi(omega_q, vec![mode(omega, g, cutoff)]).expect("model");

        let exact_sigma_z =
            thermal_expectation(&model, &sigma_z_operator(model.basis().dimension()), beta);
        let exact_sigma_x = thermal_expectation(&model, &sigma_x_operator(model.basis()), beta);

        let mut sampler = warmed(model, beta, 8, 987654);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(987654);
        let samples = 60_000usize;
        let mut sampled_sigma_z = 0.0;
        let mut sampled_sigma_x = 0.0;
        let mut sampled_energy = 0.0;
        for _ in 0..samples {
            sampler.sweep(&mut rng).expect("sweep");
            let obs = sampler.measure().expect("measure");
            sampled_sigma_z += obs.sigma_z;
            sampled_sigma_x += obs.sigma_x;
            sampled_energy += obs.energy;
        }
        sampled_sigma_z /= samples as f64;
        sampled_sigma_x /= samples as f64;
        sampled_energy /= samples as f64;

        assert!(
            (sampled_sigma_z - exact_sigma_z).abs() < 0.02,
            "sigma_z: sampled = {sampled_sigma_z}, exact = {exact_sigma_z}"
        );
        assert!(
            (sampled_sigma_x - exact_sigma_x).abs() < 0.03,
            "sigma_x: sampled = {sampled_sigma_x}, exact = {exact_sigma_x}"
        );
        // Energy from <H> via exact diagonalization for cross-check.
        let hamiltonian = sampler.model().hamiltonian();
        let exact_energy = thermal_expectation(sampler.model(), &hamiltonian, beta);
        assert!(
            (sampled_energy - exact_energy).abs() < 0.03,
            "energy: sampled = {sampled_energy}, exact = {exact_energy}"
        );
    }

    #[test]
    fn jaynes_cummings_thermal_observables_match_exact_diagonalization() {
        let omega_q = 1.0;
        let omega = 1.0;
        let g = 0.3;
        let beta = 3.0;
        let cutoff = 6;
        let model =
            OccupationSpinBosonModel::jaynes_cummings(omega_q, vec![mode(omega, g, cutoff)])
                .expect("model");

        let exact_sigma_z =
            thermal_expectation(&model, &sigma_z_operator(model.basis().dimension()), beta);
        let exact_n = {
            let dimension = model.basis().dimension();
            let mut n_op = vec![vec![0.0; dimension]; dimension];
            for (state, row) in n_op.iter_mut().enumerate() {
                row[state] = model.basis().occupation(state, 0) as f64;
            }
            thermal_expectation(&model, &n_op, beta)
        };

        let mut sampler = warmed(model, beta, 8, 246813);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(246813);
        let samples = 60_000usize;
        let mut sampled_sigma_z = 0.0;
        let mut sampled_n = 0.0;
        for _ in 0..samples {
            sampler.sweep(&mut rng).expect("sweep");
            let obs = sampler.measure().expect("measure");
            sampled_sigma_z += obs.sigma_z;
            sampled_n += obs.total_boson_number;
        }
        sampled_sigma_z /= samples as f64;
        sampled_n /= samples as f64;
        assert!(
            (sampled_sigma_z - exact_sigma_z).abs() < 0.02,
            "JC sigma_z: sampled = {sampled_sigma_z}, exact = {exact_sigma_z}"
        );
        assert!(
            (sampled_n - exact_n).abs() < 0.02,
            "JC <n>: sampled = {sampled_n}, exact = {exact_n}"
        );
    }

    #[test]
    fn sweep_kernel_is_exact_heat_bath_on_closed_paths() {
        // Per-update detailed balance (criterion D): one `sweep` must apply
        // the exact heat-bath kernel over closed worldlines,
        //   P(x -> x') = w(x') / Z   with   w(x) = prod_links T[s_k][s_{k+1}],
        // which satisfies w(x) P(x'|x) = w(x') P(x|x') identically. Verify
        // that the sweep's own bridge recipe (first-slice marginal followed
        // by exact bridge conditionals, exactly as `sweep` evaluates them)
        // reproduces that density for the realized path at machine precision.
        let model = OccupationSpinBosonModel::rabi(0.9, vec![mode(1.15, 0.35, 4)]).expect("model");
        let beta = 1.7;
        let slices = 4;
        let mut sampler = OccupationWorldlineSampler::new(model, beta, slices, 0).expect("sampler");
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0xDBA1_2026);
        for _ in 0..50 {
            sampler.sweep(&mut rng).expect("sweep");
        }
        let dimension = sampler.model.basis().dimension();
        let states = sampler.states.clone();
        let first = states[0];
        let full = &sampler.transfer_powers[slices];
        let total: f64 = (0..dimension).map(|state| full[state][state]).sum();

        // Log density of the realized path under the sweep's sampling recipe.
        let mut log_probability = (full[first][first] / total).ln();
        for index in 1..slices {
            let previous = states[index - 1];
            let remaining = slices - index;
            let mut normalization = 0.0;
            let mut realized = 0.0;
            for candidate in 0..dimension {
                let weight = sampler.transfer[previous][candidate]
                    * sampler.transfer_powers[remaining][candidate][first];
                normalization += weight;
                if candidate == states[index] {
                    realized = weight;
                }
            }
            log_probability += (realized / normalization).ln();
        }

        // Exact closed-path weight: product of link propagators over Z.
        let mut log_weight = 0.0;
        for index in 0..slices {
            log_weight += sampler.transfer[states[index]][states[(index + 1) % slices]].ln();
        }
        let log_exact = log_weight - total.ln();
        assert!(
            (log_probability - log_exact).abs() < 1.0e-12,
            "bridge recipe log-density {log_probability} differs from the exact \
             heat-bath density {log_exact}"
        );
    }

    #[test]
    fn sampler_only_visits_states_within_basis() {
        let model = OccupationSpinBosonModel::rabi(0.5, vec![mode(1.0, 0.3, 4)]).expect("model");
        let dimension = model.basis().dimension();
        let mut sampler = OccupationWorldlineSampler::new(model, 2.0, 6, 0).expect("sampler");
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(123);
        for _ in 0..500 {
            sampler.sweep(&mut rng).expect("sweep");
            for &state in sampler.states() {
                assert!(state < dimension, "state {state} >= dimension {dimension}");
            }
        }
    }

    #[test]
    fn reduced_spin_purity_is_bounded_between_zero_and_one() {
        // For a pure reduced density matrix the purity is 1; for the maximally
        // mixed spin it is 1/2. Any physical finite-temperature state lies in
        // [1/2, 1] for a spin-1/2.
        let model = OccupationSpinBosonModel::rabi(0.5, vec![mode(1.0, 0.4, 5)]).expect("model");
        let sampler = OccupationWorldlineSampler::new(model, 1.5, 6, 0).expect("sampler");
        let purity = sampler.measure().expect("measure").reduced_spin_purity;
        assert!(
            (0.5..=1.0).contains(&purity),
            "spin purity {purity} outside [0.5, 1.0]"
        );
        assert!(purity.is_finite());
    }

    #[test]
    fn acceptance_fraction_starts_at_zero_and_updates() {
        let model = OccupationSpinBosonModel::rabi(1.0, vec![mode(1.0, 0.0, 2)]).expect("model");
        let mut sampler = OccupationWorldlineSampler::new(model, 1.0, 4, 0).expect("sampler");
        assert_eq!(sampler.acceptance_fraction(), 0.0);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(7);
        sampler.sweep(&mut rng).expect("sweep");
        // After one sweep the fraction is well-defined and in [0, 1].
        let fraction = sampler.acceptance_fraction();
        assert!((0.0..=1.0).contains(&fraction));
    }

    #[test]
    fn g2_zero_of_vacuum_dominated_state_is_finite() {
        // At strong splitting and weak coupling the bosons mostly sit in the
        // vacuum; g2(0) must still be a finite number (zero mean -> zero g2).
        let model = OccupationSpinBosonModel::rabi(5.0, vec![mode(1.0, 0.01, 4)]).expect("model");
        let sampler = OccupationWorldlineSampler::new(model, 1.0, 6, 0).expect("sampler");
        let obs = sampler.measure().expect("measure");
        assert!(obs.mode_g2_zero[0].is_finite());
    }

    #[test]
    fn non_stoquastic_hamiltonian_is_rejected() {
        // sigma_x coupling in the Rabi convention enters with a minus sign in
        // the off-diagonal, which is stoquastic. To force a non-stoquastic
        // matrix we use a negative coupling and check the constructor rejects
        // the resulting positive off-diagonal transfer matrix entries.
        // With g < 0 the off-diagonal H element is -g sqrt(n+1) > 0, so
        // exp(-dt H) can develop negative off-diagonal entries.
        let model = OccupationSpinBosonModel::new(
            OccupationModelKind::Rabi,
            0.0,
            vec![CavityMode::new(1.0, -0.5, 3).expect("mode")],
        )
        .expect("model");
        let result = OccupationWorldlineSampler::new(model, 1.0, 4, 0);
        assert!(
            result.is_err(),
            "non-stoquastic Hamiltonian should be rejected by the sampler constructor"
        );
    }
}
