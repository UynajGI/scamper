//! Stochastic reconfiguration / natural gradient (Sorella 1998).
//!
//! The covariant energy gradient is `∇_k E = 2 S_k` with
//! `S = ⟨(E_L − E)(O − ⟨O⟩)⟩` (integration by parts for a real trial
//! state sampled at `|ψ_T|²`); the SR metric `G = ⟨Ȯ Ȯᵀ⟩` is the
//! Fisher information of the parameterized family. The update solves
//!
//! ```text
//! (G + λ diag(G)) Δp = − ε S
//! ```
//!
//! with the Toulouse–Umrigar-style diagonal shift `λ` as the trust
//! region: a step that worsens the block energy escalates `λ` (shrinking
//! the step toward steepest descent in the G-norm) and is retried from
//! the old parameters; an accepted step relaxes `λ` toward the minimum,
//! restoring the full natural gradient. `ε` is the dimensionless learning
//! rate (≈ the fraction of the predicted energy gain to take).

use nalgebra::{DMatrix, DVector};

use super::VariationalError;
use super::{require_finite_positive, BlockStats, Optimizer, ParamUpdate};

/// Trust-region stochastic reconfiguration.
#[derive(Debug, Clone)]
pub struct StochasticReconfiguration {
    learning_rate: f64,
    shift: f64,
    shift_min: f64,
    shift_max: f64,
    /// Convergence: natural force norm below `force_tol` for
    /// `patience` consecutive blocks.
    force_tol: f64,
    patience: usize,
    consecutive_below: usize,
    converged: bool,
}

impl StochasticReconfiguration {
    /// Trust-region SR with explicit hyperparameters.
    pub fn new(
        learning_rate: f64,
        initial_shift: f64,
        force_tol: f64,
        patience: usize,
    ) -> Result<Self, VariationalError> {
        require_finite_positive("learning_rate", learning_rate)?;
        require_finite_positive("initial_shift", initial_shift)?;
        require_finite_positive("force_tol", force_tol)?;
        if patience == 0 {
            return Err(VariationalError::invalid(
                "patience",
                "must be at least 1 block",
            ));
        }
        Ok(Self {
            learning_rate,
            shift: initial_shift,
            shift_min: 1e-8,
            shift_max: 1e8,
            force_tol,
            patience,
            consecutive_below: 0,
            converged: false,
        })
    }

    /// The current trust-region shift (checkpointed driver state).
    pub const fn shift(&self) -> f64 {
        self.shift
    }
}

impl Optimizer for StochasticReconfiguration {
    fn propose(&mut self, stats: &BlockStats) -> Result<ParamUpdate, VariationalError> {
        if stats.n_samples() == 0 {
            return Err(VariationalError::invalid(
                "block",
                "statistics contain no samples",
            ));
        }
        let norm = super::natural_force_norm(stats, self.shift);
        if norm < self.force_tol {
            self.consecutive_below += 1;
            if self.consecutive_below >= self.patience {
                self.converged = true;
            }
            return Ok(vec![0.0; stats.n_params()]);
        }
        self.consecutive_below = 0;

        let metric = stats.metric();
        let force = DVector::from_vec(stats.force());
        let regularized = &metric
            + DMatrix::from_diagonal(&DVector::from_fn(metric.nrows(), |k, _| {
                self.shift * metric[(k, k)].max(1e-300)
            }));
        let cholesky = regularized.cholesky().ok_or_else(|| {
            VariationalError::invalid(
                "metric",
                "SR metric is not positive definite even under the trust-region shift",
            )
        })?;
        let direction: DVector<f64> = cholesky.solve(&force) * -self.learning_rate;
        Ok(direction.as_slice().to_vec())
    }

    fn feedback(&mut self, accepted: bool) {
        if accepted {
            self.shift = (self.shift * 0.5).max(self.shift_min);
        } else {
            self.shift = (self.shift * 2.0).min(self.shift_max);
        }
    }

    fn converged(&self) -> bool {
        self.converged
    }
}
