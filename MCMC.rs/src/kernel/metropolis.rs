use rand::{Rng, RngExt};
use serde::{Deserialize, Serialize};

use crate::adaptation::{DiagonalCovarianceAdaptation, RobbinsMonroScale};
use crate::proposal::{standard_normal, GaussianScale};
use crate::target::{validate_log_density, LogDensity};
use crate::{EuclideanState, McmcError, SamplingPhase, TransitionKernel, TransitionReport};

/// Multivariate Gaussian random-walk Metropolis kernel.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RandomWalkMetropolis {
    scale: GaussianScale,
    scale_adaptation: Option<RobbinsMonroScale>,
    covariance_adaptation: Option<DiagonalCovarianceAdaptation>,
    proposed_position: Vec<f64>,
}

impl RandomWalkMetropolis {
    pub fn isotropic(dimension: usize, scale: f64) -> Result<Self, McmcError> {
        let proposal_scale = GaussianScale::Isotropic(scale);
        proposal_scale.validate(dimension)?;
        Ok(Self {
            scale: proposal_scale,
            scale_adaptation: None,
            covariance_adaptation: None,
            proposed_position: vec![0.0; dimension],
        })
    }

    pub fn diagonal(scales: Vec<f64>) -> Result<Self, McmcError> {
        let dimension = scales.len();
        let proposal_scale = GaussianScale::Diagonal(scales);
        proposal_scale.validate(dimension)?;
        Ok(Self {
            scale: proposal_scale,
            scale_adaptation: None,
            covariance_adaptation: None,
            proposed_position: vec![0.0; dimension],
        })
    }

    pub fn with_scale_adaptation(mut self, target_acceptance: f64) -> Result<Self, McmcError> {
        self.scale_adaptation = Some(RobbinsMonroScale::new(target_acceptance)?);
        Ok(self)
    }

    pub fn with_diagonal_covariance_adaptation(
        mut self,
        minimum_scale: f64,
    ) -> Result<Self, McmcError> {
        self.covariance_adaptation = Some(DiagonalCovarianceAdaptation::new(
            self.proposed_position.len(),
            minimum_scale,
        )?);
        Ok(self)
    }

    pub fn effective_global_multiplier(&self) -> f64 {
        self.scale_adaptation
            .as_ref()
            .map_or(1.0, RobbinsMonroScale::multiplier)
    }

    pub fn scale(&self) -> &GaussianScale {
        &self.scale
    }

    pub fn adaptation_is_frozen(&self) -> bool {
        self.scale_adaptation
            .as_ref()
            .is_none_or(RobbinsMonroScale::is_frozen)
            && self
                .covariance_adaptation
                .as_ref()
                .is_none_or(DiagonalCovarianceAdaptation::is_frozen)
    }

    fn freeze_adaptation(&mut self) {
        if let Some(adaptation) = &mut self.covariance_adaptation {
            if let Some(scales) = adaptation.finalize_scales() {
                self.scale.set_diagonal(scales);
            }
        }
        if let Some(adaptation) = &mut self.scale_adaptation {
            adaptation.freeze();
        }
    }
}

impl TransitionKernel for RandomWalkMetropolis {
    fn transition<T, R>(
        &mut self,
        target: &mut T,
        state: &mut EuclideanState,
        rng: &mut R,
        phase: SamplingPhase,
    ) -> Result<TransitionReport, McmcError>
    where
        T: LogDensity<[f64]>,
        R: Rng + ?Sized,
    {
        if phase == SamplingPhase::Sampling && !self.adaptation_is_frozen() {
            self.freeze_adaptation();
        }
        let dimension = state.dimension();
        self.scale.validate(dimension)?;
        if self.proposed_position.len() != dimension {
            return Err(McmcError::DimensionMismatch {
                expected: self.proposed_position.len(),
                actual: dimension,
            });
        }
        let multiplier = self.effective_global_multiplier();
        for (index, (proposed, current)) in self
            .proposed_position
            .iter_mut()
            .zip(state.position().iter().copied())
            .enumerate()
        {
            *proposed = current + multiplier * self.scale.component(index) * standard_normal(rng);
        }

        let (proposed_log_density, target_evaluations) =
            if self.proposed_position.iter().all(|value| value.is_finite()) {
                (
                    validate_log_density(target.log_density(&self.proposed_position))?,
                    1,
                )
            } else {
                (f64::NEG_INFINITY, 0)
            };
        let log_acceptance = proposed_log_density - state.log_density();
        if log_acceptance.is_nan() {
            return Err(McmcError::InvalidLogDensity {
                value: log_acceptance,
            });
        }
        let acceptance_probability = log_acceptance.min(0.0).exp();
        let accepted = log_acceptance >= 0.0
            || rng.random::<f64>().max(f64::MIN_POSITIVE).ln() < log_acceptance;
        if accepted {
            state.swap_position(&mut self.proposed_position, proposed_log_density);
            state.cache_mut().invalidate_gradient();
        } else {
            state.mark_rejected_transition();
        }

        if phase == SamplingPhase::Warmup {
            if let Some(adaptation) = &mut self.scale_adaptation {
                adaptation.observe(acceptance_probability)?;
            }
            if let Some(adaptation) = &mut self.covariance_adaptation {
                adaptation.observe(state.position())?;
            }
        }

        Ok(TransitionReport {
            accepted: Some(accepted),
            log_acceptance_probability: if log_acceptance.is_finite() {
                Some(log_acceptance.min(0.0))
            } else {
                None
            },
            proposals: 1,
            acceptances: if accepted { 1 } else { 0 },
            target_evaluations,
            proposal_scale: Some(multiplier),
            ..TransitionReport::default()
        })
    }

    fn on_phase_start(
        &mut self,
        phase: SamplingPhase,
        _state: &EuclideanState,
    ) -> Result<(), McmcError> {
        if phase == SamplingPhase::Sampling {
            self.freeze_adaptation();
        }
        Ok(())
    }

    fn on_phase_end(
        &mut self,
        phase: SamplingPhase,
        _state: &EuclideanState,
    ) -> Result<(), McmcError> {
        if phase == SamplingPhase::Warmup {
            self.freeze_adaptation();
        }
        Ok(())
    }

    fn name(&self) -> &'static str {
        "RandomWalkMetropolis"
    }
}
