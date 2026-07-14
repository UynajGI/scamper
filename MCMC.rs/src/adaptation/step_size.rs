use rand::Rng;
use serde::{Deserialize, Serialize};

use crate::integrator::{LeapfrogIntegrator, PhasePoint};
use crate::{DifferentiableLogDensity, McmcError, Metric};

/// Configuration for the one-step search used to initialize an HMC step size.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub struct StepSizeSearch {
    pub acceptance_threshold: f64,
    pub minimum_step_size: f64,
    pub maximum_step_size: f64,
    pub max_iterations: u16,
    pub repeat_after_metric_update: bool,
}

impl Default for StepSizeSearch {
    fn default() -> Self {
        Self {
            acceptance_threshold: 0.5,
            minimum_step_size: 1.0e-6,
            maximum_step_size: 1.0e2,
            max_iterations: 30,
            repeat_after_metric_update: true,
        }
    }
}

impl StepSizeSearch {
    pub fn validate(self) -> Result<Self, McmcError> {
        if !self.acceptance_threshold.is_finite()
            || !(0.0..1.0).contains(&self.acceptance_threshold)
        {
            return Err(McmcError::InvalidConfig(
                "step-size search acceptance threshold must lie strictly between zero and one"
                    .to_string(),
            ));
        }
        if !self.minimum_step_size.is_finite()
            || !self.maximum_step_size.is_finite()
            || self.minimum_step_size <= 0.0
            || self.maximum_step_size <= self.minimum_step_size
        {
            return Err(McmcError::InvalidConfig(
                "step-size search bounds must be finite with 0 < minimum < maximum".to_string(),
            ));
        }
        if self.max_iterations == 0 {
            return Err(McmcError::InvalidConfig(
                "step-size search must allow at least one iteration".to_string(),
            ));
        }
        Ok(self)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct StepSizeSearchResult {
    pub step_size: f64,
    pub target_evaluations: u32,
    pub gradient_evaluations: u32,
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn find_reasonable_step_size<T, M, R>(
    target: &mut T,
    metric: &M,
    position: &[f64],
    log_density: f64,
    gradient: &[f64],
    initial_step_size: f64,
    config: StepSizeSearch,
    rng: &mut R,
    integrator: &mut LeapfrogIntegrator,
) -> Result<StepSizeSearchResult, McmcError>
where
    T: DifferentiableLogDensity + ?Sized,
    M: Metric,
    R: Rng + ?Sized,
{
    let config = config.validate()?;
    let dimension = metric.dimension();
    if position.len() != dimension || gradient.len() != dimension {
        return Err(McmcError::DimensionMismatch {
            expected: dimension,
            actual: position.len().min(gradient.len()),
        });
    }
    if !log_density.is_finite()
        || position.iter().any(|value| !value.is_finite())
        || gradient.iter().any(|value| !value.is_finite())
    {
        return Err(McmcError::InvalidConfig(
            "step-size search requires a finite accepted phase point".to_string(),
        ));
    }

    let mut base = PhasePoint::with_dimension(dimension);
    base.position.copy_from_slice(position);
    base.gradient.copy_from_slice(gradient);
    base.log_density = log_density;
    metric.sample_momentum(&mut base.momentum, rng)?;
    let initial_energy = -log_density + metric.kinetic_energy(&base.momentum)?;
    if !initial_energy.is_finite() {
        return Err(McmcError::InvalidConfig(
            "step-size search initial Hamiltonian is non-finite".to_string(),
        ));
    }

    let mut target_evaluations = 0_u32;
    let mut gradient_evaluations = 0_u32;
    let mut step_size = initial_step_size.clamp(config.minimum_step_size, config.maximum_step_size);
    let mut acceptance = one_step_acceptance(
        target,
        metric,
        &base,
        initial_energy,
        step_size,
        integrator,
        &mut target_evaluations,
        &mut gradient_evaluations,
    )?;
    let increase = acceptance > config.acceptance_threshold;

    for _ in 0..config.max_iterations {
        let next_step_size = if increase {
            (step_size * 2.0).min(config.maximum_step_size)
        } else {
            (step_size * 0.5).max(config.minimum_step_size)
        };
        if next_step_size == step_size {
            break;
        }
        step_size = next_step_size;
        acceptance = one_step_acceptance(
            target,
            metric,
            &base,
            initial_energy,
            step_size,
            integrator,
            &mut target_evaluations,
            &mut gradient_evaluations,
        )?;
        let crossed = if increase {
            acceptance <= config.acceptance_threshold
        } else {
            acceptance >= config.acceptance_threshold
        };
        if crossed {
            break;
        }
    }

    Ok(StepSizeSearchResult {
        step_size,
        target_evaluations,
        gradient_evaluations,
    })
}

#[allow(clippy::too_many_arguments)]
fn one_step_acceptance<T, M>(
    target: &mut T,
    metric: &M,
    base: &PhasePoint,
    initial_energy: f64,
    step_size: f64,
    integrator: &mut LeapfrogIntegrator,
    target_evaluations: &mut u32,
    gradient_evaluations: &mut u32,
) -> Result<f64, McmcError>
where
    T: DifferentiableLogDensity + ?Sized,
    M: Metric,
{
    let mut trial = base.clone();
    let report = integrator.integrate(target, metric, &mut trial, step_size, 1)?;
    *target_evaluations = target_evaluations.saturating_add(report.target_evaluations);
    *gradient_evaluations = gradient_evaluations.saturating_add(report.gradient_evaluations);
    if report.invalid_trajectory {
        return Ok(0.0);
    }
    let final_energy = -trial.log_density + metric.kinetic_energy(&trial.momentum)?;
    let difference = final_energy - initial_energy;
    if difference.is_finite() {
        Ok((-difference).min(0.0).exp())
    } else {
        Ok(0.0)
    }
}
