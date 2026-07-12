mod component;
mod metropolis;
mod slice;

pub use component::ComponentWiseMetropolis;
pub use metropolis::RandomWalkMetropolis;
pub use slice::SliceSampler;

use rand::Rng;
use serde::{Deserialize, Serialize};

use crate::target::LogDensity;
use crate::{EuclideanState, McmcError, SamplingPhase};

/// Fixed-layout diagnostics for one transition.
#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct TransitionReport {
    pub accepted: Option<bool>,
    pub log_acceptance_probability: Option<f64>,
    pub proposals: u32,
    pub acceptances: u32,
    pub target_evaluations: u32,
    pub gradient_evaluations: u32,
    pub divergent: bool,
    pub energy_error: Option<f64>,
    pub leapfrog_steps: u32,
    pub tree_depth: Option<u16>,
    pub proposal_scale: Option<f64>,
}

impl TransitionReport {
    pub fn acceptance_rate(&self) -> Option<f64> {
        (self.proposals > 0).then(|| f64::from(self.acceptances) / f64::from(self.proposals))
    }

    pub fn validate(&self) -> Result<(), McmcError> {
        if self.acceptances > self.proposals
            || self
                .log_acceptance_probability
                .is_some_and(|value| !value.is_finite() || value > 0.0)
            || self.energy_error.is_some_and(|value| !value.is_finite())
            || self
                .proposal_scale
                .is_some_and(|value| !value.is_finite() || value <= 0.0)
        {
            return Err(McmcError::InvalidConfig(
                "transition report contains invalid values".to_string(),
            ));
        }
        Ok(())
    }
}

/// Allocation-conscious transition kernel over a Euclidean chain state.
pub trait TransitionKernel: Send {
    fn transition<T, R>(
        &mut self,
        target: &mut T,
        state: &mut EuclideanState,
        rng: &mut R,
        phase: SamplingPhase,
    ) -> Result<TransitionReport, McmcError>
    where
        T: LogDensity<[f64]>,
        R: Rng + ?Sized;

    fn on_phase_start(
        &mut self,
        _phase: SamplingPhase,
        _state: &EuclideanState,
    ) -> Result<(), McmcError> {
        Ok(())
    }

    fn on_phase_end(
        &mut self,
        _phase: SamplingPhase,
        _state: &EuclideanState,
    ) -> Result<(), McmcError> {
        Ok(())
    }

    fn name(&self) -> &'static str;
}
