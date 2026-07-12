use rand::{Rng, RngExt};
use serde::{Deserialize, Serialize};

use crate::adaptation::RobbinsMonroScale;
use crate::proposal::standard_normal;
use crate::target::{validate_log_density, LogDensity};
use crate::{EuclideanState, McmcError, SamplingPhase, TransitionKernel, TransitionReport};

/// Metropolis sweep that updates one coordinate at a time.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ComponentWiseMetropolis {
    scales: Vec<f64>,
    adaptations: Vec<Option<RobbinsMonroScale>>,
    #[serde(default)]
    proposed_position: Vec<f64>,
}

impl ComponentWiseMetropolis {
    pub fn new(scales: Vec<f64>) -> Result<Self, McmcError> {
        if scales.is_empty()
            || scales
                .iter()
                .any(|scale| !scale.is_finite() || *scale <= 0.0)
        {
            return Err(McmcError::InvalidConfig(
                "component scales must be non-empty, finite and positive".to_string(),
            ));
        }
        let adaptations = vec![None; scales.len()];
        let proposed_position = vec![0.0; scales.len()];
        Ok(Self {
            scales,
            adaptations,
            proposed_position,
        })
    }

    pub fn with_scale_adaptation(mut self, target_acceptance: f64) -> Result<Self, McmcError> {
        self.adaptations = (0..self.scales.len())
            .map(|_| RobbinsMonroScale::new(target_acceptance).map(Some))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(self)
    }

    pub fn adaptation_is_frozen(&self) -> bool {
        self.adaptations
            .iter()
            .flatten()
            .all(RobbinsMonroScale::is_frozen)
    }

    fn freeze(&mut self) {
        for adaptation in self.adaptations.iter_mut().flatten() {
            adaptation.freeze();
        }
    }
}

impl TransitionKernel for ComponentWiseMetropolis {
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
            self.freeze();
        }
        if state.dimension() != self.scales.len() {
            return Err(McmcError::DimensionMismatch {
                expected: self.scales.len(),
                actual: state.dimension(),
            });
        }
        if self.proposed_position.is_empty() {
            self.proposed_position.resize(state.dimension(), 0.0);
        } else if self.proposed_position.len() != state.dimension() {
            return Err(McmcError::DimensionMismatch {
                expected: self.proposed_position.len(),
                actual: state.dimension(),
            });
        }

        self.proposed_position.copy_from_slice(state.position());
        let mut current_log_density = state.log_density();
        let mut any_accepted = false;
        let mut report = TransitionReport::default();

        for index in 0..self.scales.len() {
            let multiplier = self.adaptations[index]
                .as_ref()
                .map_or(1.0, RobbinsMonroScale::multiplier);
            let old_value = self.proposed_position[index];
            let proposed_value = old_value + self.scales[index] * multiplier * standard_normal(rng);
            self.proposed_position[index] = proposed_value;

            let (proposed_log_density, target_evaluations) = if proposed_value.is_finite() {
                (
                    validate_log_density(target.log_density(&self.proposed_position))?,
                    1,
                )
            } else {
                (f64::NEG_INFINITY, 0)
            };
            let log_acceptance = proposed_log_density - current_log_density;
            if log_acceptance.is_nan() {
                return Err(McmcError::InvalidLogDensity {
                    value: log_acceptance,
                });
            }
            let acceptance_probability = log_acceptance.min(0.0).exp();
            let accepted = log_acceptance >= 0.0
                || rng.random::<f64>().max(f64::MIN_POSITIVE).ln() < log_acceptance;
            if accepted {
                current_log_density = proposed_log_density;
                any_accepted = true;
            } else {
                self.proposed_position[index] = old_value;
            }

            if phase == SamplingPhase::Warmup {
                if let Some(adaptation) = &mut self.adaptations[index] {
                    adaptation.observe(acceptance_probability)?;
                }
            }
            report.proposals = report.proposals.saturating_add(1);
            report.acceptances = report
                .acceptances
                .saturating_add(if accepted { 1 } else { 0 });
            report.target_evaluations =
                report.target_evaluations.saturating_add(target_evaluations);
        }

        if any_accepted {
            state.swap_position(&mut self.proposed_position, current_log_density);
            state.cache_mut().invalidate_gradient();
        } else {
            state.mark_rejected_transition();
        }

        report.accepted = None;
        report.log_acceptance_probability = None;
        report.proposal_scale = Some(
            self.scales
                .iter()
                .zip(&self.adaptations)
                .map(|(scale, adaptation)| {
                    scale
                        * adaptation
                            .as_ref()
                            .map_or(1.0, RobbinsMonroScale::multiplier)
                })
                .sum::<f64>()
                / self.scales.len() as f64,
        );
        Ok(report)
    }

    fn on_phase_start(
        &mut self,
        phase: SamplingPhase,
        _state: &EuclideanState,
    ) -> Result<(), McmcError> {
        if phase == SamplingPhase::Sampling {
            self.freeze();
        }
        Ok(())
    }

    fn on_phase_end(
        &mut self,
        phase: SamplingPhase,
        _state: &EuclideanState,
    ) -> Result<(), McmcError> {
        if phase == SamplingPhase::Warmup {
            self.freeze();
        }
        Ok(())
    }

    fn name(&self) -> &'static str {
        "ComponentWiseMetropolis"
    }
}
