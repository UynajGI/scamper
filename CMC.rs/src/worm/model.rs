//! Model contract for persistent classical worm updates.

use super::{WormError, WormState};
use rand::Rng;

/// Local worm-step proposal before target-weight evaluation.
#[derive(Debug, Clone, PartialEq)]
pub struct WormStepProposal<S> {
    pub step: S,
    /// `ln q(reverse) - ln q(forward)` excluding the kernel's close/step
    /// branch-probability correction.
    pub log_reverse_over_forward: f64,
}

impl<S> WormStepProposal<S> {
    pub const fn new(step: S, log_reverse_over_forward: f64) -> Self {
        Self {
            step,
            log_reverse_over_forward,
        }
    }
}

/// Evaluated local step in the extended state space.
#[derive(Debug, Clone, PartialEq)]
pub struct WormStepDelta<D> {
    pub new_head: D,
    /// `ln W(new) - ln W(old)` for the model configuration.
    pub log_weight_ratio: f64,
}

/// Model-provided local mechanics for a persistent worm chain.
///
/// The generic kernel owns open/close selection, branch Hastings factors,
/// log-domain acceptance and lifecycle statistics. A model owns available
/// local steps, configuration-weight ratios, cache patches and defect rules.
pub trait WormModel: Send {
    type Configuration;
    type Defect: Clone + PartialEq;
    type Step;
    type Patch: Default;

    /// Number of defects that may be selected when opening from a physical
    /// configuration. Every index in `0..count` must map to one defect.
    fn open_defect_count(&self, configuration: &Self::Configuration) -> usize;

    fn open_defect(
        &self,
        configuration: &Self::Configuration,
        index: usize,
    ) -> Result<Self::Defect, WormError>;

    /// Propose a local head move. `None` is a valid bounce for a defect with no
    /// available local step.
    fn propose_step(
        &self,
        state: &WormState<Self::Configuration, Self::Defect>,
        rng: &mut impl Rng,
    ) -> Result<Option<WormStepProposal<Self::Step>>, WormError>;

    /// Evaluate without modifying the accepted state.
    fn evaluate_step(
        &self,
        state: &WormState<Self::Configuration, Self::Defect>,
        step: &Self::Step,
        patch: &mut Self::Patch,
    ) -> Result<WormStepDelta<Self::Defect>, WormError>;

    /// Commit exactly one accepted step and its cache patch.
    fn commit_step(
        &self,
        state: &mut WormState<Self::Configuration, Self::Defect>,
        step: &Self::Step,
        patch: &Self::Patch,
    );

    /// Full configuration/defect/cache audit.
    fn validate_state(
        &self,
        state: &WormState<Self::Configuration, Self::Defect>,
    ) -> Result<(), WormError>;

    /// Optional dense endpoint-bin support for open-sector observables.
    fn endpoint_bin_count(&self) -> usize {
        0
    }

    fn endpoint_bin(&self, _defect: &Self::Defect) -> Option<usize> {
        None
    }
}
