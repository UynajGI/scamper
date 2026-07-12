use rand::{Rng, RngExt};
use serde::{Deserialize, Serialize};

use crate::target::{validate_log_density, LogDensity};
use crate::{EuclideanState, McmcError, SamplingPhase, TransitionKernel, TransitionReport};

/// Coordinate-wise univariate slice sampler with stepping-out and shrinkage.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SliceSampler {
    widths: Vec<f64>,
    max_steps_out: usize,
    max_shrink_steps: usize,
    #[serde(default)]
    working_position: Vec<f64>,
}

impl SliceSampler {
    pub fn new(widths: Vec<f64>) -> Result<Self, McmcError> {
        if widths.is_empty()
            || widths
                .iter()
                .any(|width| !width.is_finite() || *width <= 0.0)
        {
            return Err(McmcError::InvalidConfig(
                "slice widths must be non-empty, finite and positive".to_string(),
            ));
        }
        let working_position = vec![0.0; widths.len()];
        Ok(Self {
            widths,
            max_steps_out: 100,
            max_shrink_steps: 1_000,
            working_position,
        })
    }

    pub const fn with_limits(mut self, max_steps_out: usize, max_shrink_steps: usize) -> Self {
        self.max_steps_out = max_steps_out;
        self.max_shrink_steps = max_shrink_steps;
        self
    }

    fn evaluate<T>(&self, target: &mut T, report: &mut TransitionReport) -> Result<f64, McmcError>
    where
        T: LogDensity<[f64]> + ?Sized,
    {
        if self.working_position.iter().all(|value| value.is_finite()) {
            report.target_evaluations = report.target_evaluations.saturating_add(1);
            validate_log_density(target.log_density(&self.working_position))
        } else {
            Ok(f64::NEG_INFINITY)
        }
    }
}

impl<T> TransitionKernel<T> for SliceSampler
where
    T: LogDensity<[f64]> + ?Sized,
{
    fn transition<R>(
        &mut self,
        target: &mut T,
        state: &mut EuclideanState,
        rng: &mut R,
        _phase: SamplingPhase,
    ) -> Result<TransitionReport, McmcError>
    where
        R: Rng + ?Sized,
    {
        state.validate()?;
        if state.dimension() != self.widths.len() {
            return Err(McmcError::DimensionMismatch {
                expected: self.widths.len(),
                actual: state.dimension(),
            });
        }
        if self.working_position.is_empty() {
            self.working_position.resize(state.dimension(), 0.0);
        } else if self.working_position.len() != state.dimension() {
            return Err(McmcError::DimensionMismatch {
                expected: self.working_position.len(),
                actual: state.dimension(),
            });
        }
        if self.max_shrink_steps == 0 {
            return Err(McmcError::InvalidConfig(
                "slice max_shrink_steps must be positive".to_string(),
            ));
        }

        self.working_position.copy_from_slice(state.position());
        let mut current_log_density = state.log_density();
        let mut report = TransitionReport::default();

        for index in 0..self.widths.len() {
            let original = self.working_position[index];
            let slice_level = current_log_density + rng.random::<f64>().max(f64::MIN_POSITIVE).ln();
            let width = self.widths[index];
            let offset = rng.random::<f64>() * width;
            let mut left = original - offset;
            let mut right = left + width;

            let mut steps_left = rng.random_range(0..=self.max_steps_out);
            let mut steps_right = self.max_steps_out - steps_left;
            while steps_left > 0 {
                self.working_position[index] = left;
                let value = self.evaluate(target, &mut report)?;
                if value <= slice_level {
                    break;
                }
                left -= width;
                steps_left -= 1;
            }
            while steps_right > 0 {
                self.working_position[index] = right;
                let value = self.evaluate(target, &mut report)?;
                if value <= slice_level {
                    break;
                }
                right += width;
                steps_right -= 1;
            }

            let mut accepted_log_density = None;
            for _ in 0..self.max_shrink_steps {
                let interval = right - left;
                if !interval.is_finite() || interval <= 0.0 {
                    return Err(McmcError::InvalidConfig(
                        "slice bracket became non-finite or empty".to_string(),
                    ));
                }
                let proposal = left + rng.random::<f64>() * interval;
                self.working_position[index] = proposal;
                let proposed_log_density = self.evaluate(target, &mut report)?;
                if proposed_log_density >= slice_level {
                    accepted_log_density = Some(proposed_log_density);
                    break;
                }
                if proposal < original {
                    left = proposal;
                } else {
                    right = proposal;
                }
            }

            let Some(proposed_log_density) = accepted_log_density else {
                return Err(McmcError::InvalidConfig(
                    "slice shrinkage exceeded configured iteration limit".to_string(),
                ));
            };
            current_log_density = proposed_log_density;
            report.proposals = report.proposals.saturating_add(1);
            report.acceptances = report.acceptances.saturating_add(1);
        }

        state.swap_position(&mut self.working_position, current_log_density);
        state.cache_mut().invalidate_gradient();
        report.accepted = None;
        report.subtransitions = 1;
        Ok(report)
    }

    fn name(&self, _target: &T) -> &'static str {
        "SliceSampler"
    }
}
