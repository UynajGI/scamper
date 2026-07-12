use rand::Rng;
use serde::{Deserialize, Serialize};

use crate::target::{validate_log_density, LogDensity};
use crate::{EuclideanState, McmcError, SamplingPhase, TransitionKernel, TransitionReport};

/// Result returned by a target-specific Gibbs or exact block update.
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize, PartialEq)]
pub struct GibbsUpdateResult {
    /// Updated log density when the conditional update already computed it.
    /// When absent, [`GibbsKernel`] evaluates the target once.
    pub log_density: Option<f64>,
    /// Target evaluations performed inside the updater itself.
    pub target_evaluations: u32,
}

impl GibbsUpdateResult {
    pub const fn requiring_target_evaluation() -> Self {
        Self {
            log_density: None,
            target_evaluations: 0,
        }
    }

    pub const fn with_log_density(log_density: f64, target_evaluations: u32) -> Self {
        Self {
            log_density: Some(log_density),
            target_evaluations,
        }
    }
}

/// Target-specific exact conditional or blocked update.
///
/// The updater receives an immutable accepted state and a private proposal
/// workspace initialized from that state. It may change any subset of the
/// proposal. The accepted state is committed only after the updater and final
/// density validation succeed, preserving state atomicity on errors.
pub trait GibbsUpdate<T>: Send
where
    T: LogDensity<[f64]> + ?Sized,
{
    fn update<R>(
        &mut self,
        target: &mut T,
        current: &EuclideanState,
        proposed_position: &mut [f64],
        rng: &mut R,
        phase: SamplingPhase,
    ) -> Result<GibbsUpdateResult, McmcError>
    where
        R: Rng + ?Sized;

    fn on_phase_start(
        &mut self,
        _target: &mut T,
        _phase: SamplingPhase,
        _state: &EuclideanState,
    ) -> Result<(), McmcError> {
        Ok(())
    }

    fn on_phase_end(
        &mut self,
        _target: &mut T,
        _phase: SamplingPhase,
        _state: &EuclideanState,
    ) -> Result<(), McmcError> {
        Ok(())
    }

    fn name(&self, _target: &T) -> &'static str {
        "GibbsUpdate"
    }
}

/// Atomic adapter from a target-specific [`GibbsUpdate`] to a transition kernel.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GibbsKernel<U> {
    updater: U,
    #[serde(default)]
    proposed_position: Vec<f64>,
}

impl<U> GibbsKernel<U> {
    pub const fn new(updater: U) -> Self {
        Self {
            updater,
            proposed_position: Vec::new(),
        }
    }

    pub const fn updater(&self) -> &U {
        &self.updater
    }

    pub fn updater_mut(&mut self) -> &mut U {
        &mut self.updater
    }

    pub fn into_inner(self) -> U {
        self.updater
    }
}

impl<T, U> TransitionKernel<T> for GibbsKernel<U>
where
    T: LogDensity<[f64]> + ?Sized,
    U: GibbsUpdate<T>,
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
        state.validate()?;
        let dimension = state.dimension();
        if self.proposed_position.is_empty() {
            self.proposed_position.resize(dimension, 0.0);
        } else if self.proposed_position.len() != dimension {
            return Err(McmcError::DimensionMismatch {
                expected: self.proposed_position.len(),
                actual: dimension,
            });
        }
        self.proposed_position.copy_from_slice(state.position());

        let update = self
            .updater
            .update(target, state, &mut self.proposed_position, rng, phase)?;
        if self
            .proposed_position
            .iter()
            .any(|value| !value.is_finite())
        {
            return Err(McmcError::InvalidConfig(
                "Gibbs update produced a non-finite position".to_string(),
            ));
        }

        let (log_density, target_evaluations) = match update.log_density {
            Some(value) => (validate_log_density(value)?, update.target_evaluations),
            None => (
                validate_log_density(target.log_density(&self.proposed_position))?,
                update.target_evaluations.saturating_add(1),
            ),
        };
        if log_density == f64::NEG_INFINITY {
            return Err(McmcError::InvalidConfig(
                "Gibbs update produced a state outside target support".to_string(),
            ));
        }

        state.swap_position(&mut self.proposed_position, log_density);
        state.cache_mut().invalidate_gradient();
        Ok(TransitionReport {
            accepted: Some(true),
            log_acceptance_probability: Some(0.0),
            proposals: 1,
            acceptances: 1,
            target_evaluations,
            subtransitions: 1,
            ..TransitionReport::default()
        })
    }

    fn on_phase_start(
        &mut self,
        target: &mut T,
        phase: SamplingPhase,
        state: &EuclideanState,
    ) -> Result<(), McmcError> {
        self.updater.on_phase_start(target, phase, state)
    }

    fn on_phase_end(
        &mut self,
        target: &mut T,
        phase: SamplingPhase,
        state: &EuclideanState,
    ) -> Result<(), McmcError> {
        self.updater.on_phase_end(target, phase, state)
    }

    fn name(&self, target: &T) -> &'static str {
        self.updater.name(target)
    }
}
