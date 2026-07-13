use serde::{Deserialize, Serialize};

use crate::{DifferentiableLogDensity, McmcError, Metric};

/// Mutable Hamiltonian phase point owned by an HMC kernel workspace.
#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct PhasePoint {
    pub position: Vec<f64>,
    pub momentum: Vec<f64>,
    pub gradient: Vec<f64>,
    pub log_density: f64,
}

impl PhasePoint {
    pub fn with_dimension(dimension: usize) -> Self {
        Self {
            position: vec![0.0; dimension],
            momentum: vec![0.0; dimension],
            gradient: vec![0.0; dimension],
            log_density: f64::NEG_INFINITY,
        }
    }

    pub fn validate_dimension(&self, dimension: usize) -> Result<(), McmcError> {
        for actual in [
            self.position.len(),
            self.momentum.len(),
            self.gradient.len(),
        ] {
            if actual != dimension {
                return Err(McmcError::DimensionMismatch {
                    expected: dimension,
                    actual,
                });
            }
        }
        Ok(())
    }
}

/// Result of a fixed-length symplectic integration.
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize, PartialEq, Eq)]
pub struct IntegrationReport {
    pub completed_steps: u32,
    pub target_evaluations: u32,
    pub gradient_evaluations: u32,
    /// The numerical trajectory left finite target support or produced a
    /// non-finite position/gradient. The accepted chain state is not affected.
    pub invalid_trajectory: bool,
}

/// Velocity-Verlet/leapfrog integrator with reusable velocity workspace.
#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct LeapfrogIntegrator {
    #[serde(default)]
    velocity: Vec<f64>,
}

impl LeapfrogIntegrator {
    pub fn with_dimension(dimension: usize) -> Self {
        Self {
            velocity: vec![0.0; dimension],
        }
    }

    pub fn integrate<T, M>(
        &mut self,
        target: &mut T,
        metric: &M,
        point: &mut PhasePoint,
        step_size: f64,
        steps: usize,
    ) -> Result<IntegrationReport, McmcError>
    where
        T: DifferentiableLogDensity + ?Sized,
        M: Metric,
    {
        let dimension = metric.dimension();
        point.validate_dimension(dimension)?;
        if self.velocity.is_empty() {
            self.velocity.resize(dimension, 0.0);
        } else if self.velocity.len() != dimension {
            return Err(McmcError::DimensionMismatch {
                expected: dimension,
                actual: self.velocity.len(),
            });
        }
        if !step_size.is_finite() || step_size <= 0.0 || steps == 0 {
            return Err(McmcError::InvalidConfig(
                "leapfrog step size and step count must be positive".to_string(),
            ));
        }
        if point.position.iter().any(|value| !value.is_finite())
            || point.momentum.iter().any(|value| !value.is_finite())
            || point.gradient.iter().any(|value| !value.is_finite())
            || !point.log_density.is_finite()
        {
            return Err(McmcError::InvalidConfig(
                "leapfrog initial phase point must be finite".to_string(),
            ));
        }

        for (momentum, gradient) in point.momentum.iter_mut().zip(point.gradient.iter()) {
            *momentum += 0.5 * step_size * gradient;
        }

        let mut report = IntegrationReport::default();
        for step in 0..steps {
            metric.velocity(&point.momentum, &mut self.velocity)?;
            for (position, velocity) in point.position.iter_mut().zip(self.velocity.iter()) {
                *position += step_size * velocity;
            }
            if point.position.iter().any(|value| !value.is_finite()) {
                report.invalid_trajectory = true;
                return Ok(report);
            }

            let log_density = target.log_density_and_gradient(&point.position, &mut point.gradient);
            report.target_evaluations = report.target_evaluations.saturating_add(1);
            report.gradient_evaluations = report.gradient_evaluations.saturating_add(1);
            if !log_density.is_finite() || point.gradient.iter().any(|value| !value.is_finite()) {
                report.invalid_trajectory = true;
                return Ok(report);
            }
            point.log_density = log_density;
            report.completed_steps = report.completed_steps.saturating_add(1);

            let momentum_scale = if step + 1 == steps { 0.5 } else { 1.0 };
            for (momentum, gradient) in point.momentum.iter_mut().zip(point.gradient.iter()) {
                *momentum += momentum_scale * step_size * gradient;
            }
        }
        Ok(report)
    }
}
