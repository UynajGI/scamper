use rand::{Rng, RngExt};
use serde::{Deserialize, Serialize};

use crate::adaptation::{
    DenseCovarianceAdaptation, DiagonalCovarianceAdaptation, RobbinsMonroScale,
};
use crate::proposal::GaussianScale;
use crate::target::{validate_log_density, LogDensity};
use crate::{EuclideanState, McmcError, SamplingPhase, TransitionKernel, TransitionReport};

/// Multivariate Gaussian random-walk Metropolis kernel.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RandomWalkMetropolis {
    scale: GaussianScale,
    scale_adaptation: Option<RobbinsMonroScale>,
    covariance_adaptation: Option<DiagonalCovarianceAdaptation>,
    #[serde(default)]
    dense_covariance_adaptation: Option<DenseCovarianceAdaptation>,
    proposed_position: Vec<f64>,
    #[serde(default)]
    normal_buffer: Vec<f64>,
}

impl RandomWalkMetropolis {
    pub fn isotropic(dimension: usize, scale: f64) -> Result<Self, McmcError> {
        let proposal_scale = GaussianScale::Isotropic(scale);
        proposal_scale.validate(dimension)?;
        Ok(Self {
            scale: proposal_scale,
            scale_adaptation: None,
            covariance_adaptation: None,
            dense_covariance_adaptation: None,
            proposed_position: vec![0.0; dimension],
            normal_buffer: vec![0.0; dimension],
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
            dense_covariance_adaptation: None,
            proposed_position: vec![0.0; dimension],
            normal_buffer: vec![0.0; dimension],
        })
    }

    pub fn dense_cholesky(dimension: usize, cholesky: Vec<f64>) -> Result<Self, McmcError> {
        let proposal_scale = GaussianScale::dense_cholesky(dimension, cholesky)?;
        Ok(Self {
            scale: proposal_scale,
            scale_adaptation: None,
            covariance_adaptation: None,
            dense_covariance_adaptation: None,
            proposed_position: vec![0.0; dimension],
            normal_buffer: vec![0.0; dimension],
        })
    }

    pub fn dense_covariance(
        dimension: usize,
        covariance: &[f64],
        jitter: f64,
    ) -> Result<Self, McmcError> {
        let proposal_scale = GaussianScale::dense_from_covariance(dimension, covariance, jitter)?;
        Ok(Self {
            scale: proposal_scale,
            scale_adaptation: None,
            covariance_adaptation: None,
            dense_covariance_adaptation: None,
            proposed_position: vec![0.0; dimension],
            normal_buffer: vec![0.0; dimension],
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
        self.dense_covariance_adaptation = None;
        self.covariance_adaptation = Some(DiagonalCovarianceAdaptation::new(
            self.proposed_position.len(),
            minimum_scale,
        )?);
        Ok(self)
    }

    pub fn with_dense_covariance_adaptation(
        mut self,
        regularization: f64,
    ) -> Result<Self, McmcError> {
        self.covariance_adaptation = None;
        self.dense_covariance_adaptation = Some(DenseCovarianceAdaptation::new(
            self.proposed_position.len(),
            regularization,
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
            && self
                .dense_covariance_adaptation
                .as_ref()
                .is_none_or(DenseCovarianceAdaptation::is_frozen)
    }

    fn freeze_adaptation(&mut self) {
        if let Some(adaptation) = &mut self.dense_covariance_adaptation {
            if let Some(cholesky) = adaptation.finalize_cholesky() {
                self.scale
                    .set_dense_cholesky(self.proposed_position.len(), cholesky);
            }
        } else if let Some(adaptation) = &mut self.covariance_adaptation {
            if let Some(scales) = adaptation.finalize_scales() {
                self.scale.set_diagonal(scales);
            }
        }
        if let Some(adaptation) = &mut self.scale_adaptation {
            adaptation.freeze();
        }
    }
}

impl<T> TransitionKernel<T> for RandomWalkMetropolis
where
    T: LogDensity<[f64]> + ?Sized,
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
            self.freeze_adaptation();
        }
        state.validate()?;
        let dimension = state.dimension();
        self.scale.validate(dimension)?;
        if self.proposed_position.len() != dimension {
            return Err(McmcError::DimensionMismatch {
                expected: self.proposed_position.len(),
                actual: dimension,
            });
        }
        if self.normal_buffer.is_empty() {
            self.normal_buffer.resize(dimension, 0.0);
        } else if self.normal_buffer.len() != dimension {
            return Err(McmcError::DimensionMismatch {
                expected: self.normal_buffer.len(),
                actual: dimension,
            });
        }

        self.scale
            .fill_displacement(rng, &mut self.normal_buffer, &mut self.proposed_position)?;
        let multiplier = self.effective_global_multiplier();
        for (proposed, current) in self
            .proposed_position
            .iter_mut()
            .zip(state.position().iter().copied())
        {
            *proposed = current + multiplier * *proposed;
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
            if let Some(adaptation) = &mut self.dense_covariance_adaptation {
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
            acceptance_statistic: Some(acceptance_probability),
            proposals: 1,
            acceptances: if accepted { 1 } else { 0 },
            target_evaluations,
            proposal_scale: Some(multiplier),
            subtransitions: 1,
            ..TransitionReport::default()
        })
    }

    fn on_phase_start(
        &mut self,
        _target: &mut T,
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
        _target: &mut T,
        phase: SamplingPhase,
        _state: &EuclideanState,
    ) -> Result<(), McmcError> {
        if phase == SamplingPhase::Warmup {
            self.freeze_adaptation();
        }
        Ok(())
    }

    fn name(&self, _target: &T) -> &'static str {
        "RandomWalkMetropolis"
    }
}
