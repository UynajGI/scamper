//! Correlated-sampling variance minimization (Umrigar-style) on `argmin`.
//!
//! The third L2 entry point. Unlike [`StochasticReconfiguration`] and
//! [`LinearMethod`](super::LinearMethod) it is not block-adaptive: it is a
//! one-shot batch optimization over a FIXED reference sample, so it does
//! not implement the [`Optimizer`](super::Optimizer) trait (whose
//! propose/feedback cycle assumes fresh Markov blocks).
//!
//! # Method
//!
//! Reference configurations `{R_i}` are drawn once from
//! `|psi_{p0}|^2` (a kernel's walker snapshots). For a candidate `p` the
//! correlated-sampling reweighting gives
//!
//! ```text
//! w_i(p)  = |psi_p(R_i) / psi_{p0}(R_i)|^2
//! E(p)    = sum_i w_i E_L,p(R_i) / sum_i w_i
//! Var(p)  = sum_i w_i (E_L,p(R_i) - E(p))^2 / sum_i w_i
//! ```
//!
//! and `Var(p)` — Umrigar's variance-minimization objective — is handed
//! to `argmin`'s Nelder-Mead simplex (mature-crate policy: the search is
//! the library's business; the physics — weights and the estimator — is
//! ours). The objective is minimized to zero variance at an exact trial
//! state, where it is stationary in a particularly strong sense: every
//! `E_L` sample equals the exact energy.
//!
//! The gradient of `Var` w.r.t. `p` would require third-order parameter
//! derivative chains (`d/dp of grad ln psi` and of `lap ln psi`) that the
//! `WaveFunction` trait deliberately does not carry, which is why the
//! solver is derivative-free.
//!
//! # Determinism and validation
//!
//! With a uniform grid as the "sample" the reweighting reproduces exact
//! quadrature at any `p` (weight `exp(-2 p r^2) * exp(2(p - p0) r^2)`),
//! so the whole pipeline — weights, two-pass variance, simplex search —
//! is testable deterministically against the closed forms
//! `Var(alpha) = c(alpha)^2 * 3/(8 alpha^2)`, zero at `alpha = omega/2`.
//!
//! Out-of-domain candidates (ansatz rejects the parameters) cost a large
//! finite penalty instead of aborting, so the simplex walks around them.

use argmin::core::{CostFunction, Error as ArgminError, Executor, State};
use argmin::solver::neldermead::NelderMead;

use super::super::error::VariationalError;
use super::super::hamiltonian::ContinuumHamiltonian;
use super::super::wavefunction::{Positions, WaveFunctionParams};
use super::require_finite_positive;

/// One correlated-sampling reference point: a configuration, the
/// reference-state log-amplitude, and the log-density of the sampling
/// measure the configuration was drawn from (importance base weight).
#[derive(Debug, Clone)]
pub struct ReferenceSample {
    pub configuration: Positions,
    pub reference_log_psi: f64,
    /// `ln(base density)` of the sampling measure at this configuration,
    /// up to a sample-independent constant. The effective weight of the
    /// sample at candidate parameters `p` is
    /// `exp(base_log_weight + 2 (ln|psi_p| - ln|psi_ref|))`.
    pub base_log_weight: f64,
}

impl ReferenceSample {
    /// Snapshot a configuration sampled from `|psi_ref|^2` (the kernel
    /// case). For a deterministic uniform grid the same constructor is
    /// correct: the flat base measure is a constant that cancels in the
    /// normalized estimator, and the effective weight
    /// `exp(2 ln|psi_ref| + 2 (ln|psi_p| - ln|psi_ref|)) = |psi_p|^2`
    /// is exactly the target measure at any candidate `p`.
    pub fn new<W: WaveFunctionParams<Config = Positions>>(
        wave_function: &W,
        configuration: Positions,
    ) -> Self {
        let reference_log_psi = wave_function.log_psi(&configuration);
        Self {
            configuration,
            base_log_weight: 2.0 * reference_log_psi,
            reference_log_psi,
        }
    }
}

/// The correlated-sampling variance objective, `CostFunction`-shaped for
/// `argmin`.
pub struct VarianceObjective<W: WaveFunctionParams<Config = Positions> + Clone> {
    wave_function: W,
    hamiltonian: ContinuumHamiltonian,
    samples: Vec<ReferenceSample>,
    /// Cost returned for out-of-domain parameter vectors (finite, so the
    /// simplex can walk around infeasible vertices).
    penalty: f64,
}

impl<W: WaveFunctionParams<Config = Positions> + Clone> VarianceObjective<W> {
    /// Build the objective over reference samples (at the reference
    /// parameters of `wave_function`).
    pub fn new(
        wave_function: W,
        hamiltonian: ContinuumHamiltonian,
        samples: Vec<ReferenceSample>,
    ) -> Result<Self, VariationalError> {
        if samples.is_empty() {
            return Err(VariationalError::invalid(
                "samples",
                "at least one reference configuration is required",
            ));
        }
        require_finite_positive("penalty", 1e6)?;
        Ok(Self {
            wave_function,
            hamiltonian,
            samples,
            penalty: 1e6,
        })
    }

    /// The variance objective at an explicit parameter vector, without
    /// the argmin machinery (deterministic-test entry point).
    pub fn variance_at(&self, params: &[f64]) -> f64 {
        let mut candidate = self.wave_function.clone();
        if candidate.try_set_params(params).is_err() {
            return self.penalty;
        }
        self.weighted_variance(&candidate)
    }

    /// Two-pass weighted variance of the local energy under the
    /// correlated-sampling reweighting (two-pass: no cancellation floor —
    /// the samples are retained, unlike `BlockStats`' accumulated sums).
    fn weighted_variance(&self, candidate: &W) -> f64 {
        let Self {
            hamiltonian,
            samples,
            ..
        } = self;
        let mut grad_scratch =
            super::super::wavefunction::GradBuffer::new(samples[0].configuration.n_particles());
        let mut weight_sum = 0.0_f64;
        let mut mean = 0.0_f64;
        let mut weighted_energies = Vec::with_capacity(samples.len());
        for sample in samples {
            // Importance-sampling weight relative to the sampling measure:
            // w = base_density * |psi_new / psi_ref|^2
            //   = exp(base_log_weight + 2 (ln|psi_new| - ln|psi_ref|)),
            // which equals |psi_new|^2 up to the constant of the flat or
            // reference measure (absorbed by the normalization below).
            let delta = candidate.log_psi(&sample.configuration) - sample.reference_log_psi;
            let weight = (sample.base_log_weight + 2.0 * delta).exp();
            let local = super::super::estimators::local_energy(
                candidate,
                hamiltonian,
                &sample.configuration,
                &mut grad_scratch,
            )
            .value;
            weighted_energies.push((weight, local));
            weight_sum += weight;
            mean += weight * local;
        }
        if !weight_sum.is_finite() || weight_sum <= 0.0 {
            return f64::INFINITY;
        }
        mean /= weight_sum;
        let mut variance = 0.0_f64;
        for (weight, local) in weighted_energies {
            let deviation = local - mean;
            variance += weight * deviation * deviation;
        }
        variance / weight_sum
    }
}

impl<W: WaveFunctionParams<Config = Positions> + Clone> CostFunction for VarianceObjective<W> {
    type Param = Vec<f64>;
    type Output = f64;

    fn cost(&self, param: &Self::Param) -> Result<Self::Output, ArgminError> {
        Ok(self.variance_at(param))
    }
}

/// Outcome of a variance-minimization run.
#[derive(Debug, Clone)]
pub struct VarianceMinimizationResult {
    /// Optimized variational parameters.
    pub params: Vec<f64>,
    /// Objective value (correlated-sampling variance) at `params`.
    pub variance: f64,
    /// Simplex iterations consumed.
    pub iterations: u64,
}

/// One-shot correlated-sampling variance minimization driven by argmin's
/// Nelder-Mead simplex.
#[derive(Debug, Clone)]
pub struct VarianceMinimization {
    max_iterations: u64,
    sd_tolerance: f64,
    /// Initial-simplex edge scale per parameter, as a fraction of the
    /// reference parameter value (floored for near-zero parameters).
    simplex_scale: f64,
}

impl VarianceMinimization {
    /// Configure the search.
    pub fn new(
        max_iterations: u64,
        sd_tolerance: f64,
        simplex_scale: f64,
    ) -> Result<Self, VariationalError> {
        require_finite_positive("sd_tolerance", sd_tolerance)?;
        require_finite_positive("simplex_scale", simplex_scale)?;
        Ok(Self {
            max_iterations,
            sd_tolerance,
            simplex_scale,
        })
    }

    /// Minimize the correlated-sampling variance over `samples` starting
    /// from the wave function's current (reference) parameters.
    pub fn minimize<W: WaveFunctionParams<Config = Positions> + Clone>(
        &self,
        wave_function: W,
        hamiltonian: ContinuumHamiltonian,
        samples: Vec<ReferenceSample>,
    ) -> Result<VarianceMinimizationResult, VariationalError> {
        let initial = wave_function.param_values();
        let objective = VarianceObjective::new(wave_function, hamiltonian, samples)?;

        // Initial simplex: the reference vertex plus one displaced vertex
        // per parameter (scale = fraction of |p0| with a floor).
        let mut vertices = vec![initial.clone()];
        for (k, &p) in initial.iter().enumerate() {
            let mut displaced = initial.clone();
            displaced[k] += self.simplex_scale * p.abs().max(0.1);
            vertices.push(displaced);
        }
        let solver = NelderMead::new(vertices)
            .with_sd_tolerance(self.sd_tolerance)
            .map_err(|error| {
                VariationalError::invalid(
                    "nelder_mead",
                    format!("solver construction failed: {error}"),
                )
            })?;

        let result = Executor::new(objective, solver)
            .configure(|state| state.max_iters(self.max_iterations))
            .run()
            .map_err(|error| {
                VariationalError::invalid(
                    "nelder_mead",
                    format!("variance minimization failed: {error}"),
                )
            })?;
        let state = result.state();
        let params = state
            .get_best_param()
            .cloned()
            .ok_or_else(|| VariationalError::invalid("nelder_mead", "no best vertex"))?;
        let variance = state.get_best_cost();
        Ok(VarianceMinimizationResult {
            params,
            variance,
            iterations: state.get_iter(),
        })
    }
}
