//! Parameter optimizers for the variational family (L2).
//!
//! The layer contract (DESIGN.md §3): optimizers are an OUTER loop that
//! never runs inside `sweep`. A driver alternating measurement blocks and
//! parameter updates feeds each block's statistics to
//! [`Optimizer::propose`], applies the returned delta through
//! `WaveFunction::update_params`, and reports the energetic outcome via
//! [`Optimizer::feedback`]. Optimizer state is checkpointed alongside the
//! walkers by the driver.
//!
//! Physics split (mature-crate policy): the linear algebra — Cholesky
//! solves for stochastic reconfiguration, the symmetric eigenproblem for
//! the linear method — is `nalgebra`; the hand-derived content is the
//! statistics below and the update rules of Sorella (1998) and
//! Umrigar–Nightingale (2001).
//!
//! # Block statistics
//!
//! Per measurement sample (one walker configuration) we record the local
//! energy `E_L` and the parameter gradient of the log-amplitude
//! `O_k = ∂_k ln|ψ_T|`. With centered variables `Ȯ_k = O_k − ⟨O_k⟩` the
//! optimization-relevant moments are (all derived from raw moments by
//! exact algebra in [`BlockStats`]):
//!
//! ```text
//! E       = ⟨E_L⟩,                σ² = ⟨E_L²⟩ − E²
//! S_k     = ⟨Ȯ_k E_L⟩             (force; covariant energy gradient
//!                                  ∇_k E = 2 S_k — zero for an exact
//!                                  ψ_T, since E_L is then constant)
//! G_kl    = ⟨Ȯ_k Ȯ_l⟩             (stochastic-reconfiguration metric)
//! T_kl    = ⟨Ȯ_k Ȯ_l E_L⟩         (three-point moment, linear method)
//! ```
//!
//! L2 entry points: [`StochasticReconfiguration`] (natural gradient with a
//! trust-region shift), [`LinearMethod`] (Umrigar–Nightingale generalized
//! eigenproblem on the linearized basis), and [`VarianceMinimization`]
//! (correlated-sampling variance minimization on `argmin`'s Nelder–Mead).

pub mod linear_method;
pub mod sr;
pub mod variance;

pub use linear_method::LinearMethod;
pub use sr::StochasticReconfiguration;
pub use variance::{
    ReferenceSample, VarianceMinimization, VarianceMinimizationResult, VarianceObjective,
};

use nalgebra::{DMatrix, DVector};

pub(crate) use super::error::VariationalError;

/// Raw-moment accumulator over one measurement block; see the module docs
/// for the derived (centered) quantities and their physical meaning.
#[derive(Debug, Clone)]
pub struct BlockStats {
    n_samples: usize,
    weight_total: f64,
    n_params: usize,
    sum_e: f64,
    sum_e2: f64,
    /// `sum_o[k] = Σ O_k`.
    sum_o: Vec<f64>,
    /// Symmetric `P x P` raw moment `Σ O_k O_l` (row-major upper-triangular
    /// storage is unnecessary at these sizes; a dense `DMatrix` keeps the
    /// nalgebra path direct).
    sum_oo: DMatrix<f64>,
    sum_oe: Vec<f64>,
    /// `Σ O_k O_l E_L` (symmetric).
    sum_ooe: DMatrix<f64>,
}

impl BlockStats {
    /// An empty accumulator for `n_params` variational parameters.
    pub fn new(n_params: usize) -> Self {
        Self {
            n_samples: 0,
            weight_total: 0.0,
            n_params,
            sum_e: 0.0,
            sum_e2: 0.0,
            sum_o: vec![0.0; n_params],
            sum_oo: DMatrix::zeros(n_params, n_params),
            sum_oe: vec![0.0; n_params],
            sum_ooe: DMatrix::zeros(n_params, n_params),
        }
    }

    /// Number of variational parameters the statistics describe.
    pub const fn n_params(&self) -> usize {
        self.n_params
    }

    /// Samples accumulated so far.
    pub const fn n_samples(&self) -> usize {
        self.n_samples
    }

    /// Fold one walker measurement (local energy and its parameter
    /// gradient) into the block.
    pub fn push(&mut self, local_energy: f64, o: &[f64]) {
        self.push_weighted(local_energy, o, 1.0);
    }

    /// Fold one importance-weighted sample. `n_samples` counts the raw
    /// pushes (the statistical-error estimator stays honest for equal
    /// walker weights); every moment scales by `weight`, so deterministic
    /// quadrature over `|ψ_T|²`-weighted grids produces the exact same
    /// centered statistics as MC sampling would in expectation.
    pub fn push_weighted(&mut self, local_energy: f64, o: &[f64], weight: f64) {
        assert_eq!(
            o.len(),
            self.n_params,
            "parameter-gradient length must match the block's parameter count"
        );
        assert!(
            weight.is_finite() && weight >= 0.0,
            "sample weights must be finite and non-negative"
        );
        self.n_samples += 1;
        self.weight_total += weight;
        self.sum_e += weight * local_energy;
        self.sum_e2 += weight * local_energy * local_energy;
        for (k, &o_k) in o.iter().enumerate() {
            self.sum_o[k] += weight * o_k;
            self.sum_oe[k] += weight * o_k * local_energy;
        }
        for k in 0..self.n_params {
            for l in k..self.n_params {
                let moment = o[k] * o[l];
                let weighted = moment * local_energy;
                self.sum_oo[(k, l)] += weight * moment;
                self.sum_oo[(l, k)] = self.sum_oo[(k, l)];
                self.sum_ooe[(k, l)] += weight * weighted;
                self.sum_ooe[(l, k)] = self.sum_ooe[(k, l)];
            }
        }
    }

    /// Reset for the next block (keeps the allocation).
    pub fn reset(&mut self) {
        let n = self.n_params;
        *self = Self::new(n);
    }

    /// Block mean local energy (weighted; the weight sum normalizes).
    pub fn energy(&self) -> f64 {
        self.sum_e / self.weight_total
    }

    /// Block variance of the local energy (the zero-variance principle's
    /// diagnostic).
    pub fn energy_variance(&self) -> f64 {
        let mean = self.energy();
        self.sum_e2 / self.weight_total - mean * mean
    }

    /// Statistical error of the mean via the naive i.i.d. estimator (the
    /// driver's blocks are short; correlated-time corrections are the
    /// measurement layer's business, not the optimizer's).
    pub fn energy_stderr_naive(&self) -> f64 {
        (self.energy_variance() / self.n_samples.max(1) as f64).sqrt()
    }

    /// Centered force `S_k = ⟨Ȯ_k E_L⟩` (covariant gradient / 2). Vanishes
    /// with the variance for an exact trial state.
    pub fn force(&self) -> Vec<f64> {
        let n = self.weight_total;
        let energy = self.energy();
        (0..self.n_params)
            .map(|k| self.sum_oe[k] / n - (self.sum_o[k] / n) * energy)
            .collect()
    }

    /// Centered stochastic-reconfiguration metric `G_kl = ⟨Ȯ_k Ȯ_l⟩`.
    pub fn metric(&self) -> DMatrix<f64> {
        let n = self.weight_total;
        let mut g = self.sum_oo.clone();
        g *= 1.0 / n;
        for k in 0..self.n_params {
            let m_k = self.sum_o[k] / n;
            for l in 0..self.n_params {
                g[(k, l)] -= m_k * (self.sum_o[l] / n);
            }
        }
        g
    }

    /// Centered three-point moment `T_kl = ⟨Ȯ_k Ȯ_l E_L⟩`, obtained from
    /// the raw moments by the exact multivariate centering identity
    /// `⟨(X−x̄)(Y−ȳ)(Z−z̄)⟩ = ⟨XYZ⟩ − x̄⟨YZ⟩ − ȳ⟨XZ⟩ − z̄⟨XY⟩
    /// + 2 x̄ ȳ z̄` with `Z = E_L`.
    pub fn metric_energy(&self) -> DMatrix<f64> {
        let n = self.weight_total;
        let energy = self.energy();
        let mean_o: Vec<f64> = (0..self.n_params).map(|k| self.sum_o[k] / n).collect();
        let mut t = DMatrix::zeros(self.n_params, self.n_params);
        for k in 0..self.n_params {
            for l in 0..self.n_params {
                let raw_oo = self.sum_oo[(k, l)] / n;
                t[(k, l)] = self.sum_ooe[(k, l)] / n
                    - mean_o[k] * (self.sum_oe[l] / n)
                    - mean_o[l] * (self.sum_oe[k] / n)
                    - raw_oo * energy
                    + 2.0 * mean_o[k] * mean_o[l] * energy;
            }
        }
        t
    }
}

/// A proposed parameter delta to feed `WaveFunction::update_params`.
pub type ParamUpdate = Vec<f64>;

/// Outer-loop optimizer contract; see the module docs.
pub trait Optimizer {
    /// Propose the next parameter update from one block's statistics.
    fn propose(&mut self, stats: &BlockStats) -> Result<ParamUpdate, VariationalError>;

    /// Report whether the last applied update improved the block energy
    /// (the driver's call, right after measuring the post-update block).
    fn feedback(&mut self, accepted: bool);

    /// Converged when the normalized force has stayed below tolerance for
    /// the required number of consecutive blocks.
    fn converged(&self) -> bool;
}

/// Shared convergence bookkeeping for the concrete optimizers: the
/// force norm `sqrt(SᵀG⁻¹S)` is the natural dimensionless scale (it is the
/// predicted energy gain of the natural-gradient step); `sqrt` guards the
/// degenerate zero-variance/exact-state case where both S and G vanish.
pub(crate) fn natural_force_norm(stats: &BlockStats, shift: f64) -> f64 {
    let force = stats.force();
    let metric = stats.metric();
    let regularized = &metric
        + DMatrix::from_diagonal(&nalgebra::DVector::from_fn(metric.nrows(), |k, _| {
            shift.max(1e-12) * metric[(k, k)].max(1e-12)
        }));
    match regularized.cholesky() {
        Some(cholesky) => {
            let force_vector = DVector::from_vec(force);
            let solved = cholesky.solve(&force_vector);
            force_vector.dot(&solved).max(0.0).sqrt()
        }
        None => force.iter().map(|s| s.abs()).sum(),
    }
}

/// Validate common optimizer hyperparameters (criterion G).
pub(crate) fn require_finite_positive(label: &str, value: f64) -> Result<(), VariationalError> {
    if !value.is_finite() || value <= 0.0 {
        return Err(VariationalError::invalid(
            label,
            format!("must be finite and positive, got {value}"),
        ));
    }
    Ok(())
}
