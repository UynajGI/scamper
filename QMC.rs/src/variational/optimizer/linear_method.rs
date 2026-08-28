//! The linear method (Umrigar–Nightingale 2001).
//!
//! Work in the linearized parameter-displacement basis
//! `{ψ, ψ δ_1, …, ψ δ_P}` with `ψ δ_k` denoting `ψ·(1 + s_k Ȯ_k)` to
//! first order in the per-parameter shift `s_k`. The generalized
//! eigenproblem `H c = μ S c` over that basis, with the matrix elements
//! estimated from one block's samples,
//!
//! ```text
//! S_00 = 1,      S_0k = 0,      S_kl = s_k s_l G_kl
//! H_00 = E,      H_0k = s_k S_k,
//! H_kl = E + s_k S_k + s_l S_l + s_k s_l T_kl
//! ```
//!
//! (derived by expanding `⟨(1+δ_i·Ȯ) E_L (1+δ_j·Ȯ)⟩` with centered
//! `Ȯ`; `G`, `S`, `T` are the [`BlockStats`] moments), has eigenvectors
//! whose lowest eigenvalue `μ` estimates the best linear-combination
//! energy. Enforcing the `c_0 = 1` gauge (the combination must remain
//! normalized to leading order) turns that eigenvector into the update
//! `Δp_k = c_k s_k`.
//!
//! The metric side is regularized by the same diagonal-shift trust region
//! as SR: escalates on rejected steps, relaxes on accepted ones. This
//! makes the linear method's "step quality" comparable to SR on equal
//! footing — the L2 gate asks that it never be slower to converge.

use nalgebra::{DMatrix, SymmetricEigen};

use super::VariationalError;
use super::{require_finite_positive, BlockStats, Optimizer, ParamUpdate};

/// Trust-region linear method on the linearized displacement basis.
#[derive(Debug, Clone)]
pub struct LinearMethod {
    /// Per-parameter displacement scale `s` of the linearized basis
    /// (recommended O(0.1..=1); it is a real scale, not a gauge: tiny
    /// shifts collapse the displacement block of the overlap and degrade
    /// step quality, large shifts push the linearization outside its
    /// validity — on the deterministic Gaussian toy, s = 0.05 needs
    /// 500+ iterations where s = 0.5 converges in ~17).
    shift: f64,
    metric_shift: f64,
    force_tol: f64,
    patience: usize,
    consecutive_below: usize,
    converged: bool,
}

impl LinearMethod {
    /// Trust-region linear method with explicit hyperparameters.
    pub fn new(
        displacement_shift: f64,
        metric_shift: f64,
        force_tol: f64,
        patience: usize,
    ) -> Result<Self, VariationalError> {
        require_finite_positive("displacement_shift", displacement_shift)?;
        require_finite_positive("metric_shift", metric_shift)?;
        require_finite_positive("force_tol", force_tol)?;
        if patience == 0 {
            return Err(VariationalError::invalid(
                "patience",
                "must be at least 1 block",
            ));
        }
        Ok(Self {
            shift: displacement_shift,
            metric_shift,
            force_tol,
            patience,
            consecutive_below: 0,
            converged: false,
        })
    }

    /// The current metric trust-region shift (checkpointed driver state).
    pub const fn metric_shift(&self) -> f64 {
        self.metric_shift
    }
}

impl Optimizer for LinearMethod {
    fn propose(&mut self, stats: &BlockStats) -> Result<ParamUpdate, VariationalError> {
        if stats.n_samples() == 0 {
            return Err(VariationalError::invalid(
                "block",
                "statistics contain no samples",
            ));
        }
        let norm = super::natural_force_norm(stats, self.metric_shift);
        if norm < self.force_tol {
            self.consecutive_below += 1;
            if self.consecutive_below >= self.patience {
                self.converged = true;
            }
            return Ok(vec![0.0; stats.n_params()]);
        }
        self.consecutive_below = 0;

        let p = stats.n_params();
        let energy = stats.energy();
        let force = stats.force();
        let metric = stats.metric();
        let three_point = stats.metric_energy();

        // Overlap matrix on the {0, e_1..e_P} basis with the diagonal
        // trust-region shift (keeps the Cholesky transform well-defined).
        let mut overlap = DMatrix::identity(p + 1, p + 1);
        for k in 0..p {
            for l in 0..p {
                overlap[(k + 1, l + 1)] = self.shift * self.shift * metric[(k, l)]
                    + self.metric_shift * metric[(k, k)].max(1e-300) * (k == l) as i32 as f64;
            }
        }
        // Hamiltonian in the same basis (module docs equations).
        let mut hamiltonian = DMatrix::zeros(p + 1, p + 1);
        hamiltonian[(0, 0)] = energy;
        for k in 0..p {
            hamiltonian[(0, k + 1)] = self.shift * force[k];
            hamiltonian[(k + 1, 0)] = self.shift * force[k];
            for l in 0..p {
                hamiltonian[(k + 1, l + 1)] = energy
                    + self.shift * force[k]
                    + self.shift * force[l]
                    + self.shift * self.shift * three_point[(k, l)];
            }
        }

        // Symmetric reduction H c = μ S c with S = L Lᵀ: substitute
        // c = L^{-T} y to get (L^{-1} H L^{-T}) y = μ y, and recover
        // c = L^{-T} y from the lowest eigenvector.
        let cholesky = overlap.cholesky().ok_or_else(|| {
            VariationalError::invalid(
                "overlap",
                "linear-method overlap is not positive definite under the shift",
            )
        })?;
        let l_inverse = cholesky.inverse();
        let reduced = &l_inverse * &hamiltonian * l_inverse.transpose();
        let eigen = SymmetricEigen::new(reduced);
        // Lowest eigenvalue's eigenvector, back in the c = L^{-T} y form.
        let (mut best_index, mut best_value) = (0usize, f64::INFINITY);
        for (index, &value) in eigen.eigenvalues.iter().enumerate() {
            if value < best_value {
                best_value = value;
                best_index = index;
            }
        }
        let mut coefficients = l_inverse.transpose() * eigen.eigenvectors.column(best_index);
        // c_0 = 1 gauge (normalize the sign; a near-zero c_0 means the
        // eigenvector left the physical sector — treat as a rejected step).
        let c0 = coefficients[0];
        if c0.abs() < 1e-10 {
            return Ok(vec![0.0; p]);
        }
        coefficients *= 1.0 / c0;
        Ok((1..=p).map(|k| coefficients[k] * self.shift).collect())
    }

    fn feedback(&mut self, accepted: bool) {
        if accepted {
            self.metric_shift = (self.metric_shift * 0.5).max(1e-8);
        } else {
            self.metric_shift = (self.metric_shift * 2.0).min(1e8);
        }
    }

    fn converged(&self) -> bool {
        self.converged
    }
}

/// Expose the raw (unregularized) predicted energy of the last
/// eigenvector for diagnostics: the lowest generalized eigenvalue is the
/// linear method's estimate of the post-update energy. Kept out of the
/// trait because only drivers optimizing with LM can use it.
impl LinearMethod {
    /// Predicted post-update energy from the block statistics alone
    /// (solves the unshifted problem; falls back to the current block
    /// energy when the unshifted overlap is singular).
    pub fn predicted_energy(stats: &BlockStats) -> f64 {
        let p = stats.n_params();
        let energy = stats.energy();
        let force = stats.force();
        let metric = stats.metric();
        let three_point = stats.metric_energy();
        let mut overlap = DMatrix::identity(p + 1, p + 1);
        for k in 0..p {
            for l in 0..p {
                overlap[(k + 1, l + 1)] = metric[(k, l)];
            }
        }
        let mut hamiltonian = DMatrix::zeros(p + 1, p + 1);
        hamiltonian[(0, 0)] = energy;
        for k in 0..p {
            hamiltonian[(0, k + 1)] = force[k];
            hamiltonian[(k + 1, 0)] = force[k];
            for l in 0..p {
                hamiltonian[(k + 1, l + 1)] = energy + force[k] + force[l] + three_point[(k, l)];
            }
        }
        let Some(cholesky) = overlap.cholesky() else {
            return energy;
        };
        let l_inverse = cholesky.inverse();
        let reduced = &l_inverse * &hamiltonian * l_inverse.transpose();
        SymmetricEigen::new(reduced)
            .eigenvalues
            .iter()
            .copied()
            .fold(f64::INFINITY, f64::min)
    }
}
