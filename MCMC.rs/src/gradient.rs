//! Finite-difference validation for analytic target gradients.

use crate::{DifferentiableLogDensity, LogDensity, McmcError};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GradientCheckConfig {
    pub step: f64,
    pub absolute_tolerance: f64,
    pub relative_tolerance: f64,
}

impl Default for GradientCheckConfig {
    fn default() -> Self {
        Self {
            step: f64::EPSILON.cbrt(),
            absolute_tolerance: 1.0e-6,
            relative_tolerance: 1.0e-5,
        }
    }
}

impl GradientCheckConfig {
    pub fn validate(self) -> Result<Self, McmcError> {
        if !self.step.is_finite() || self.step <= 0.0 {
            return Err(McmcError::InvalidConfig(
                "gradient-check step must be finite and positive".to_string(),
            ));
        }
        if !self.absolute_tolerance.is_finite() || self.absolute_tolerance < 0.0 {
            return Err(McmcError::InvalidConfig(
                "gradient-check absolute tolerance must be finite and non-negative".to_string(),
            ));
        }
        if !self.relative_tolerance.is_finite() || self.relative_tolerance < 0.0 {
            return Err(McmcError::InvalidConfig(
                "gradient-check relative tolerance must be finite and non-negative".to_string(),
            ));
        }
        Ok(self)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct GradientComponentCheck {
    pub index: usize,
    pub analytic: f64,
    pub finite_difference: f64,
    pub absolute_error: f64,
    pub relative_error: f64,
    pub passed: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct GradientCheckReport {
    pub log_density: f64,
    pub components: Vec<GradientComponentCheck>,
    pub maximum_absolute_error: f64,
    pub maximum_relative_error: f64,
    pub target_evaluations: u32,
    pub passed: bool,
}

pub fn check_gradient<T>(
    target: &mut T,
    position: &[f64],
    config: GradientCheckConfig,
) -> Result<GradientCheckReport, McmcError>
where
    T: DifferentiableLogDensity + ?Sized,
{
    let config = config.validate()?;
    if position.is_empty() {
        return Err(McmcError::InvalidConfig(
            "gradient checks require a non-empty position".to_string(),
        ));
    }
    if position.iter().any(|value| !value.is_finite()) {
        return Err(McmcError::InvalidConfig(
            "gradient checks require finite positions".to_string(),
        ));
    }

    let dimension = position.len();
    let mut analytic = vec![0.0; dimension];
    let log_density = target.log_density_and_gradient(position, &mut analytic);
    if !log_density.is_finite() || analytic.iter().any(|value| !value.is_finite()) {
        return Err(McmcError::InvalidConfig(
            "gradient checks require a finite target value and analytic gradient".to_string(),
        ));
    }

    let mut perturbed = position.to_vec();
    let mut components = Vec::with_capacity(dimension);
    let mut maximum_absolute_error: f64 = 0.0;
    let mut maximum_relative_error: f64 = 0.0;
    for index in 0..dimension {
        let scale = position[index].abs().max(1.0);
        let step = config.step * scale;
        perturbed[index] = position[index] + step;
        let plus = LogDensity::log_density(target, &perturbed);
        perturbed[index] = position[index] - step;
        let minus = LogDensity::log_density(target, &perturbed);
        perturbed[index] = position[index];
        if !plus.is_finite() || !minus.is_finite() {
            return Err(McmcError::InvalidConfig(format!(
                "gradient finite difference left the target support at component {index}"
            )));
        }

        let finite_difference = (plus - minus) / (2.0 * step);
        let absolute_error = (analytic[index] - finite_difference).abs();
        let denominator = analytic[index].abs().max(finite_difference.abs()).max(1.0);
        let relative_error = absolute_error / denominator;
        let passed =
            absolute_error <= config.absolute_tolerance + config.relative_tolerance * denominator;
        maximum_absolute_error = maximum_absolute_error.max(absolute_error);
        maximum_relative_error = maximum_relative_error.max(relative_error);
        components.push(GradientComponentCheck {
            index,
            analytic: analytic[index],
            finite_difference,
            absolute_error,
            relative_error,
            passed,
        });
    }

    Ok(GradientCheckReport {
        log_density,
        passed: components.iter().all(|component| component.passed),
        components,
        maximum_absolute_error,
        maximum_relative_error,
        target_evaluations: u32::try_from(dimension)
            .unwrap_or(u32::MAX)
            .saturating_mul(2)
            .saturating_add(1),
    })
}
