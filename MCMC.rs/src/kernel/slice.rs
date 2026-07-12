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
        Ok(Self {
            widths,
            max_steps_out: 100,
            max_shrink_steps: 1_000,
        })
    }

    pub const fn with_limits(mut self, max_steps_out: usize, max_shrink_steps: usize) -> Self {
        self.max_steps_out = max_steps_out;
        self.max_shrink_steps = max_shrink_steps;
        self
    }
}

impl TransitionKernel for SliceSampler {
    fn transition<T, R>(
        &mut self,
        target: &mut T,
        state: &mut EuclideanState,
        rng: &mut R,
        _phase: SamplingPhase,
    ) -> Result<TransitionReport, McmcError>
    where
        T: LogDensity<[f64]>,
        R: Rng + ?Sized,
    {
        if state.dimension() != self.widths.len() {
            return Err(McmcError::DimensionMismatch {
                expected: self.widths.len(),
                actual: state.dimension(),
            });
        }
        let mut report = TransitionReport::default();
        for index in 0..self.widths.len() {
            let original = state.position()[index];
            let slice_level = state.log_density() + rng.random::<f64>().max(f64::MIN_POSITIVE).ln();
            let width = self.widths[index];
            let offset = rng.random::<f64>() * width;
            let mut left = original - offset;
            let mut right = left + width;

            let mut steps_left = rng.random_range(0..=self.max_steps_out);
            let mut steps_right = self.max_steps_out - steps_left;
            while steps_left > 0 {
                state.position_mut_for_cache_rebuild()[index] = left;
                let value = validate_log_density(target.log_density(state.position()))?;
                report.target_evaluations = report.target_evaluations.saturating_add(1);
                if value <= slice_level {
                    break;
                }
                left -= width;
                steps_left -= 1;
            }
            while steps_right > 0 {
                state.position_mut_for_cache_rebuild()[index] = right;
                let value = validate_log_density(target.log_density(state.position()))?;
                report.target_evaluations = report.target_evaluations.saturating_add(1);
                if value <= slice_level {
                    break;
                }
                right += width;
                steps_right -= 1;
            }

            let mut accepted = false;
            for _ in 0..self.max_shrink_steps {
                let proposal = left + rng.random::<f64>() * (right - left);
                state.position_mut_for_cache_rebuild()[index] = proposal;
                let proposed_log_density =
                    validate_log_density(target.log_density(state.position()))?;
                report.target_evaluations = report.target_evaluations.saturating_add(1);
                if proposed_log_density >= slice_level {
                    let position = state.position().clone();
                    state.replace(position, proposed_log_density);
                    state.cache_mut().invalidate_gradient();
                    accepted = true;
                    break;
                }
                if proposal < original {
                    left = proposal;
                } else {
                    right = proposal;
                }
            }
            if !accepted {
                state.position_mut_for_cache_rebuild()[index] = original;
                return Err(McmcError::InvalidConfig(
                    "slice shrinkage exceeded configured iteration limit".to_string(),
                ));
            }
            report.proposals = report.proposals.saturating_add(1);
            report.acceptances = report.acceptances.saturating_add(1);
        }
        report.accepted = None;
        Ok(report)
    }

    fn name(&self) -> &'static str {
        "SliceSampler"
    }
}
