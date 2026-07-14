use rand::{Rng, RngExt};
use serde::{Deserialize, Serialize};

use crate::adaptation::{
    find_reasonable_step_size, HmcWarmup, MetricAdaptation, MetricUpdate, StepSizeSearch,
    WarmupWindowConfig,
};
use crate::integrator::{LeapfrogIntegrator, PhasePoint};
use crate::metric::{DenseMetric, DiagonalMetric, Metric, MetricKind, UnitMetric};
use crate::{
    DifferentiableLogDensity, EuclideanState, McmcError, SamplingPhase, TransitionKernel,
    TransitionReport,
};

/// Fixed-trajectory Hamiltonian Monte Carlo with optional warmup adaptation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StaticHmc<M> {
    metric: M,
    step_size: f64,
    leapfrog_steps: usize,
    max_energy_error: f64,
    warmup: Option<HmcWarmup>,
    #[serde(default)]
    step_size_search: Option<StepSizeSearch>,
    #[serde(default = "step_size_search_complete_default")]
    step_size_search_complete: bool,
    // Per-transition workspaces are rebuilt from the accepted state. Skipping
    // them keeps checkpoints finite even after an invalid divergent path.
    #[serde(skip, default)]
    phase_point: PhasePoint,
    #[serde(skip, default)]
    current_gradient: Vec<f64>,
    #[serde(skip, default)]
    integrator: LeapfrogIntegrator,
}

impl<M> StaticHmc<M>
where
    M: Metric,
{
    pub fn new(metric: M, step_size: f64, leapfrog_steps: usize) -> Result<Self, McmcError> {
        validate_hmc_config(metric.dimension(), step_size, leapfrog_steps, 1_000.0)?;
        let dimension = metric.dimension();
        Ok(Self {
            metric,
            step_size,
            leapfrog_steps,
            max_energy_error: 1_000.0,
            warmup: None,
            step_size_search: None,
            step_size_search_complete: true,
            phase_point: PhasePoint::with_dimension(dimension),
            current_gradient: vec![0.0; dimension],
            integrator: LeapfrogIntegrator::with_dimension(dimension),
        })
    }

    pub fn with_max_energy_error(mut self, max_energy_error: f64) -> Result<Self, McmcError> {
        if !max_energy_error.is_finite() || max_energy_error <= 0.0 {
            return Err(McmcError::InvalidConfig(
                "maximum HMC energy error must be finite and positive".to_string(),
            ));
        }
        self.max_energy_error = max_energy_error;
        Ok(self)
    }

    pub fn with_dual_averaging(
        self,
        total_warmup: u64,
        target_acceptance: f64,
    ) -> Result<Self, McmcError> {
        let windows = WarmupWindowConfig::default_for(total_warmup)?;
        self.with_warmup_adaptation(
            total_warmup,
            target_acceptance,
            MetricAdaptation::None,
            windows,
        )
    }

    pub fn with_warmup_adaptation(
        mut self,
        total_warmup: u64,
        target_acceptance: f64,
        metric_adaptation: MetricAdaptation,
        windows: WarmupWindowConfig,
    ) -> Result<Self, McmcError> {
        validate_metric_adaptation(self.metric.kind(), metric_adaptation)?;
        self.warmup = Some(HmcWarmup::new(
            self.metric.dimension(),
            total_warmup,
            self.step_size,
            target_acceptance,
            metric_adaptation,
            windows,
        )?);
        Ok(self)
    }

    /// Enable a one-step warmup search for a reasonable initial step size.
    pub fn with_step_size_search(mut self, search: StepSizeSearch) -> Result<Self, McmcError> {
        self.step_size_search = Some(search.validate()?);
        self.step_size_search_complete = false;
        Ok(self)
    }

    pub const fn metric(&self) -> &M {
        &self.metric
    }

    pub fn metric_mut(&mut self) -> &mut M {
        &mut self.metric
    }

    pub const fn step_size(&self) -> f64 {
        self.step_size
    }

    pub const fn leapfrog_steps(&self) -> usize {
        self.leapfrog_steps
    }

    pub const fn max_energy_error(&self) -> f64 {
        self.max_energy_error
    }

    pub fn adaptation_is_frozen(&self) -> bool {
        self.warmup.as_ref().is_none_or(HmcWarmup::is_frozen)
    }

    fn ensure_workspace(&mut self, dimension: usize) -> Result<(), McmcError> {
        if self.metric.dimension() != dimension {
            return Err(McmcError::DimensionMismatch {
                expected: self.metric.dimension(),
                actual: dimension,
            });
        }
        if self.phase_point.position.is_empty()
            && self.phase_point.momentum.is_empty()
            && self.phase_point.gradient.is_empty()
        {
            self.phase_point = PhasePoint::with_dimension(dimension);
        }
        self.phase_point.validate_dimension(dimension)?;
        if self.current_gradient.is_empty() {
            self.current_gradient.resize(dimension, 0.0);
        } else if self.current_gradient.len() != dimension {
            return Err(McmcError::DimensionMismatch {
                expected: dimension,
                actual: self.current_gradient.len(),
            });
        }
        Ok(())
    }

    fn prepare_sampling(&mut self) -> Result<(), McmcError> {
        if let Some(warmup) = &mut self.warmup {
            self.step_size = warmup.finish()?;
        }
        Ok(())
    }

    fn apply_metric_update(&mut self, update: MetricUpdate) -> Result<(), McmcError> {
        match update {
            MetricUpdate::Diagonal(diagonal) => self.metric.set_diagonal_inverse_mass(&diagonal),
            MetricUpdate::Dense {
                dimension,
                covariance,
                jitter,
            } => self
                .metric
                .set_dense_inverse_mass(dimension, &covariance, jitter),
        }
    }
}

impl StaticHmc<UnitMetric> {
    pub fn unit(
        dimension: usize,
        step_size: f64,
        leapfrog_steps: usize,
    ) -> Result<Self, McmcError> {
        Self::new(UnitMetric::new(dimension)?, step_size, leapfrog_steps)
    }
}

impl StaticHmc<DiagonalMetric> {
    pub fn diagonal(
        inverse_mass: Vec<f64>,
        step_size: f64,
        leapfrog_steps: usize,
    ) -> Result<Self, McmcError> {
        Self::new(
            DiagonalMetric::new(inverse_mass)?,
            step_size,
            leapfrog_steps,
        )
    }

    pub fn with_diagonal_adaptation(
        self,
        total_warmup: u64,
        target_acceptance: f64,
        regularization: f64,
    ) -> Result<Self, McmcError> {
        let windows = WarmupWindowConfig::default_for(total_warmup)?;
        self.with_warmup_adaptation(
            total_warmup,
            target_acceptance,
            MetricAdaptation::Diagonal { regularization },
            windows,
        )
    }
}

impl StaticHmc<DenseMetric> {
    pub fn dense(
        dimension: usize,
        inverse_mass: &[f64],
        jitter: f64,
        step_size: f64,
        leapfrog_steps: usize,
    ) -> Result<Self, McmcError> {
        Self::new(
            DenseMetric::from_inverse_mass(dimension, inverse_mass, jitter)?,
            step_size,
            leapfrog_steps,
        )
    }

    pub fn with_dense_adaptation(
        self,
        total_warmup: u64,
        target_acceptance: f64,
        regularization: f64,
    ) -> Result<Self, McmcError> {
        let windows = WarmupWindowConfig::default_for(total_warmup)?;
        self.with_warmup_adaptation(
            total_warmup,
            target_acceptance,
            MetricAdaptation::Dense { regularization },
            windows,
        )
    }
}

impl<T, M> TransitionKernel<T> for StaticHmc<M>
where
    T: DifferentiableLogDensity + ?Sized,
    M: Metric,
{
    fn transition<R>(
        &mut self,
        target: &mut T,
        state: &mut EuclideanState,
        rng: &mut R,
        phase: SamplingPhase,
    ) -> Result<TransitionReport, McmcError>
    where
        R: Rng + ?Sized,
    {
        if phase == SamplingPhase::Sampling && !self.adaptation_is_frozen() {
            self.prepare_sampling()?;
        }
        state.validate()?;
        let dimension = state.dimension();
        self.ensure_workspace(dimension)?;

        let mut target_evaluations = 0_u32;
        let mut gradient_evaluations = 0_u32;
        if let Some(gradient) = state.cache().gradient() {
            self.current_gradient.copy_from_slice(gradient);
        } else {
            let log_density =
                target.log_density_and_gradient(state.position(), &mut self.current_gradient);
            target_evaluations = target_evaluations.saturating_add(1);
            gradient_evaluations = gradient_evaluations.saturating_add(1);
            if !log_density.is_finite()
                || self.current_gradient.iter().any(|value| !value.is_finite())
            {
                return Err(McmcError::InvalidConfig(
                    "differentiable target returned a non-finite accepted-state gradient"
                        .to_string(),
                ));
            }
            state.synchronize_gradient(log_density, &self.current_gradient)?;
        }

        if let Some(search) = self.step_size_search {
            if self.step_size_search_complete {
                // already searched; skip
            } else if phase != SamplingPhase::Warmup || self.warmup.is_none() {
                return Err(McmcError::InvalidConfig(
                    "automatic HMC step-size search requires configured warmup".to_string(),
                ));
            } else {
                let result = find_reasonable_step_size(
                    target,
                    &self.metric,
                    state.position(),
                    state.log_density(),
                    &self.current_gradient,
                    self.step_size,
                    search,
                    rng,
                    &mut self.integrator,
                )?;
                target_evaluations = target_evaluations.saturating_add(result.target_evaluations);
                gradient_evaluations =
                    gradient_evaluations.saturating_add(result.gradient_evaluations);
                self.step_size = result.step_size;
                self.warmup
                    .as_mut()
                    .expect("checked HMC warmup")
                    .restart_step_size(result.step_size)?;
                self.step_size_search_complete = true;
            }
        }
        let used_step_size = self.step_size;

        self.phase_point.position.copy_from_slice(state.position());
        self.phase_point
            .gradient
            .copy_from_slice(&self.current_gradient);
        self.phase_point.log_density = state.log_density();
        self.metric
            .sample_momentum(&mut self.phase_point.momentum, rng)?;
        let initial_kinetic = self.metric.kinetic_energy(&self.phase_point.momentum)?;
        let initial_energy = -state.log_density() + initial_kinetic;
        if !initial_energy.is_finite() {
            return Err(McmcError::InvalidConfig(
                "HMC initial Hamiltonian is non-finite".to_string(),
            ));
        }

        let integration = self.integrator.integrate(
            target,
            &self.metric,
            &mut self.phase_point,
            used_step_size,
            self.leapfrog_steps,
        )?;
        target_evaluations = target_evaluations.saturating_add(integration.target_evaluations);
        gradient_evaluations =
            gradient_evaluations.saturating_add(integration.gradient_evaluations);

        let mut divergent = integration.invalid_trajectory;
        let mut energy_error = None;
        let mut log_acceptance_probability = None;
        let mut acceptance_probability = 0.0;
        if !divergent {
            let final_kinetic = self.metric.kinetic_energy(&self.phase_point.momentum)?;
            let final_energy = -self.phase_point.log_density + final_kinetic;
            let difference = final_energy - initial_energy;
            if difference.is_finite() {
                energy_error = Some(difference);
                divergent = difference.abs() > self.max_energy_error;
                if !divergent {
                    let log_acceptance = (-difference).min(0.0);
                    log_acceptance_probability = Some(log_acceptance);
                    acceptance_probability = log_acceptance.exp();
                }
            } else {
                divergent = true;
            }
        }

        let accepted = !divergent
            && (acceptance_probability >= 1.0
                || rng.random::<f64>().max(f64::MIN_POSITIVE).ln()
                    < log_acceptance_probability.unwrap_or(f64::NEG_INFINITY));
        if accepted {
            state.commit_hamiltonian_proposal(
                &mut self.phase_point.position,
                self.phase_point.log_density,
                &self.phase_point.gradient,
            )?;
        } else {
            state.mark_rejected_transition();
        }

        if phase == SamplingPhase::Warmup {
            let mut metric_updated = false;
            if let Some(warmup) = &mut self.warmup {
                let observation = warmup.observe(acceptance_probability, state.position())?;
                self.step_size = observation.step_size;
                if let Some(update) = observation.metric_update {
                    self.apply_metric_update(update)?;
                    metric_updated = true;
                }
            }
            if metric_updated
                && self
                    .step_size_search
                    .is_some_and(|search| search.repeat_after_metric_update)
            {
                self.step_size_search_complete = false;
            }
        }

        Ok(TransitionReport {
            accepted: Some(accepted),
            log_acceptance_probability,
            acceptance_statistic: Some(acceptance_probability),
            proposals: 1,
            acceptances: if accepted { 1 } else { 0 },
            target_evaluations,
            gradient_evaluations,
            divergent,
            energy: Some(initial_energy),
            energy_error,
            leapfrog_steps: integration.completed_steps,
            tree_depth: None,
            max_tree_depth_reached: false,
            proposal_scale: Some(used_step_size),
            subtransitions: 1,
        })
    }

    fn on_phase_start(
        &mut self,
        _target: &mut T,
        phase: SamplingPhase,
        _state: &EuclideanState,
    ) -> Result<(), McmcError> {
        if phase == SamplingPhase::Sampling {
            self.prepare_sampling()?;
        }
        Ok(())
    }

    fn on_phase_end(
        &mut self,
        _target: &mut T,
        phase: SamplingPhase,
        _state: &EuclideanState,
    ) -> Result<(), McmcError> {
        if phase == SamplingPhase::Warmup {
            self.prepare_sampling()?;
        }
        Ok(())
    }

    fn name(&self, _target: &T) -> &'static str {
        "StaticHmc"
    }
}

const fn step_size_search_complete_default() -> bool {
    true
}

fn validate_hmc_config(
    dimension: usize,
    step_size: f64,
    leapfrog_steps: usize,
    max_energy_error: f64,
) -> Result<(), McmcError> {
    if dimension == 0 {
        return Err(McmcError::InvalidConfig(
            "HMC dimension must be positive".to_string(),
        ));
    }
    if !step_size.is_finite() || step_size <= 0.0 || leapfrog_steps == 0 {
        return Err(McmcError::InvalidConfig(
            "HMC step size and leapfrog count must be positive".to_string(),
        ));
    }
    if !max_energy_error.is_finite() || max_energy_error <= 0.0 {
        return Err(McmcError::InvalidConfig(
            "maximum HMC energy error must be finite and positive".to_string(),
        ));
    }
    Ok(())
}

fn validate_metric_adaptation(
    metric: MetricKind,
    adaptation: MetricAdaptation,
) -> Result<(), McmcError> {
    let compatible = matches!(adaptation, MetricAdaptation::None)
        || matches!(
            (metric, adaptation),
            (MetricKind::Diagonal, MetricAdaptation::Diagonal { .. })
                | (MetricKind::Dense, MetricAdaptation::Diagonal { .. })
                | (MetricKind::Dense, MetricAdaptation::Dense { .. })
        );
    if compatible {
        Ok(())
    } else {
        Err(McmcError::InvalidConfig(
            "requested metric adaptation is incompatible with the HMC metric".to_string(),
        ))
    }
}
