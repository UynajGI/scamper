mod component;
mod compose;
mod gibbs;
mod hmc;
mod metropolis;
mod nuts;
mod slice;

pub use component::ComponentWiseMetropolis;
pub use compose::{Mixture, Repeat, Then};
pub use gibbs::{GibbsKernel, GibbsUpdate, GibbsUpdateResult};
pub use hmc::StaticHmc;
pub use metropolis::RandomWalkMetropolis;
pub use nuts::Nuts;
pub use slice::SliceSampler;

use rand::Rng;
use serde::{Deserialize, Serialize};

use crate::target::LogDensity;
use crate::{EuclideanState, McmcError, SamplingPhase};

/// Fixed-layout diagnostics for one transition or composed transition sweep.
#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct TransitionReport {
    pub accepted: Option<bool>,
    pub log_acceptance_probability: Option<f64>,
    /// Mean Metropolis acceptance statistic used for Hamiltonian adaptation.
    #[serde(default)]
    pub acceptance_statistic: Option<f64>,
    pub proposals: u32,
    pub acceptances: u32,
    pub target_evaluations: u32,
    pub gradient_evaluations: u32,
    pub divergent: bool,
    /// Hamiltonian energy associated with this iteration, when available.
    #[serde(default)]
    pub energy: Option<f64>,
    pub energy_error: Option<f64>,
    pub leapfrog_steps: u32,
    pub tree_depth: Option<u16>,
    /// Whether a dynamic Hamiltonian trajectory exhausted its depth limit.
    #[serde(default)]
    pub max_tree_depth_reached: bool,
    pub proposal_scale: Option<f64>,
    /// Number of elementary kernel transitions represented by this report.
    #[serde(default)]
    pub subtransitions: u32,
}

impl TransitionReport {
    pub fn acceptance_rate(&self) -> Option<f64> {
        (self.proposals > 0).then(|| f64::from(self.acceptances) / f64::from(self.proposals))
    }

    /// Merge another elementary/composed report into this one.
    ///
    /// Scalar diagnostics that are not well-defined after composition become
    /// `None`; extrema and counters are aggregated deterministically.
    pub fn merge(&mut self, other: Self) {
        if *self == Self::default() {
            *self = other;
            return;
        }

        let previous_subtransitions = self.subtransitions.max(1);
        let other_subtransitions = other.subtransitions.max(1);
        self.accepted = None;
        self.log_acceptance_probability = None;
        self.acceptance_statistic = weighted_optional_mean(
            self.acceptance_statistic,
            previous_subtransitions,
            other.acceptance_statistic,
            other_subtransitions,
        );
        self.proposals = self.proposals.saturating_add(other.proposals);
        self.acceptances = self.acceptances.saturating_add(other.acceptances);
        self.target_evaluations = self
            .target_evaluations
            .saturating_add(other.target_evaluations);
        self.gradient_evaluations = self
            .gradient_evaluations
            .saturating_add(other.gradient_evaluations);
        self.divergent |= other.divergent;
        if self.energy != other.energy {
            self.energy = None;
        }
        self.energy_error = worst_energy_error(self.energy_error, other.energy_error);
        self.leapfrog_steps = self.leapfrog_steps.saturating_add(other.leapfrog_steps);
        self.tree_depth = match (self.tree_depth, other.tree_depth) {
            (Some(left), Some(right)) => Some(left.max(right)),
            (Some(value), None) | (None, Some(value)) => Some(value),
            (None, None) => None,
        };
        self.max_tree_depth_reached |= other.max_tree_depth_reached;
        if self.proposal_scale != other.proposal_scale {
            self.proposal_scale = None;
        }
        self.subtransitions = previous_subtransitions.saturating_add(other_subtransitions);
    }

    fn normalize_subtransitions(&mut self) {
        self.subtransitions = self.subtransitions.max(1);
    }

    pub fn validate(&self) -> Result<(), McmcError> {
        if self.acceptances > self.proposals
            || self
                .log_acceptance_probability
                .is_some_and(|value| !value.is_finite() || value > 0.0)
            || self
                .acceptance_statistic
                .is_some_and(|value| !value.is_finite() || !(0.0..=1.0).contains(&value))
            || self.energy.is_some_and(|value| !value.is_finite())
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

fn weighted_optional_mean(
    left: Option<f64>,
    left_weight: u32,
    right: Option<f64>,
    right_weight: u32,
) -> Option<f64> {
    match (left, right) {
        (Some(left), Some(right)) => {
            let denominator = f64::from(left_weight.saturating_add(right_weight));
            Some((left * f64::from(left_weight) + right * f64::from(right_weight)) / denominator)
        }
        (Some(value), None) | (None, Some(value)) => Some(value),
        (None, None) => None,
    }
}

fn worst_energy_error(left: Option<f64>, right: Option<f64>) -> Option<f64> {
    match (left, right) {
        (Some(left), Some(right)) if right.abs() > left.abs() => Some(right),
        (Some(left), Some(_)) => Some(left),
        (Some(value), None) | (None, Some(value)) => Some(value),
        (None, None) => None,
    }
}

/// Allocation-conscious transition kernel over a Euclidean chain state.
///
/// The target type is a trait parameter rather than a method parameter so
/// target-specific Gibbs and blocked updates can implement this interface.
pub trait TransitionKernel<T>: Send
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

    fn name(&self, _target: &T) -> &'static str;
}
