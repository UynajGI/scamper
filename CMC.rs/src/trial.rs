//! Transactional trial-move execution.
//!
//! A trial is evaluated without mutating the accepted configuration.  On
//! acceptance, the backend commits the move and its cache patch exactly once;
//! rejection simply discards the patch contents.  This prevents rollback bugs
//! and gives state backends a reusable workspace for incremental caches.

use crate::ensemble::Ensemble;
use rand::{Rng, RngExt};

/// Proposal plus the Metropolis-Hastings proposal-density correction.
#[derive(Debug, Clone)]
pub struct ProposedMove<M> {
    pub movement: M,
    /// `ln q(old | new) - ln q(new | old)`.
    pub log_reverse_over_forward: f64,
}

impl<M> ProposedMove<M> {
    #[inline]
    pub const fn symmetric(movement: M) -> Self {
        Self {
            movement,
            log_reverse_over_forward: 0.0,
        }
    }

    #[inline]
    pub const fn new(movement: M, log_reverse_over_forward: f64) -> Self {
        Self {
            movement,
            log_reverse_over_forward,
        }
    }
}

/// State-side evaluator for one move type.
///
/// `Patch` is owned and reused by the update kernel.  Evaluation may fill it
/// with incremental-energy, neighbor-list or other cache updates.  The state is
/// not mutated until `commit_trial` is called.
pub trait TrialEvaluator<Model, Movement> {
    type Delta;
    type Patch;

    fn evaluate_trial(
        &self,
        model: &Model,
        movement: &Movement,
        patch: &mut Self::Patch,
    ) -> Self::Delta;

    fn commit_trial(&mut self, movement: &Movement, patch: &Self::Patch);
}

/// Result of one attempted transition.
#[derive(Debug, Clone)]
pub struct TrialOutcome<D> {
    pub accepted: bool,
    pub delta: D,
    pub log_acceptance: f64,
}

/// Evaluate and execute one generic Metropolis-Hastings transition.
#[inline]
pub fn metropolis_hastings_step<State, Model, Movement, Target, RngType>(
    state: &mut State,
    model: &Model,
    proposal: &ProposedMove<Movement>,
    target: &Target,
    patch: &mut <State as TrialEvaluator<Model, Movement>>::Patch,
    rng: &mut RngType,
) -> TrialOutcome<<State as TrialEvaluator<Model, Movement>>::Delta>
where
    State: TrialEvaluator<Model, Movement>,
    Target: Ensemble<<State as TrialEvaluator<Model, Movement>>::Delta>,
    RngType: Rng,
{
    let delta = state.evaluate_trial(model, &proposal.movement, patch);
    let target_ratio = target.log_weight_ratio(&delta);
    let log_acceptance = target_ratio + proposal.log_reverse_over_forward;
    assert!(
        !target_ratio.is_nan()
            && !proposal.log_reverse_over_forward.is_nan()
            && !log_acceptance.is_nan(),
        "trial produced NaN log weight"
    );

    let accepted =
        log_acceptance >= 0.0 || rng.random::<f64>().max(f64::MIN_POSITIVE).ln() < log_acceptance;
    if accepted {
        state.commit_trial(&proposal.movement, patch);
    }

    TrialOutcome {
        accepted,
        delta,
        log_acceptance,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    #[derive(Default)]
    struct Patch {
        value: i32,
    }

    struct State(i32);

    struct Model;

    struct DeterministicTarget;

    struct FiniteTarget;

    impl Ensemble<i32> for FiniteTarget {
        fn log_weight_ratio(&self, _delta: &i32) -> f64 {
            -2.0
        }
    }

    impl Ensemble<i32> for DeterministicTarget {
        fn log_weight_ratio(&self, delta: &i32) -> f64 {
            if *delta <= 0 {
                0.0
            } else {
                f64::NEG_INFINITY
            }
        }
    }

    impl TrialEvaluator<Model, i32> for State {
        type Delta = i32;
        type Patch = Patch;

        fn evaluate_trial(&self, _model: &Model, movement: &i32, patch: &mut Patch) -> i32 {
            patch.value = *movement - self.0;
            patch.value
        }

        fn commit_trial(&mut self, movement: &i32, _patch: &Patch) {
            self.0 = *movement;
        }
    }

    #[test]
    fn rejected_trial_does_not_mutate_state() {
        let mut state = State(1);
        let proposal = ProposedMove::symmetric(2);
        let mut patch = Patch::default();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(4);
        let outcome = metropolis_hastings_step(
            &mut state,
            &Model,
            &proposal,
            &DeterministicTarget,
            &mut patch,
            &mut rng,
        );
        assert!(!outcome.accepted);
        assert_eq!(state.0, 1);
    }

    #[test]
    fn hastings_correction_participates_in_acceptance() {
        let mut state = State(1);
        let proposal = ProposedMove::new(2, 3.0);
        let mut patch = Patch::default();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(5);
        let outcome = metropolis_hastings_step(
            &mut state,
            &Model,
            &proposal,
            &FiniteTarget,
            &mut patch,
            &mut rng,
        );
        assert!(outcome.accepted);
        assert_eq!(state.0, 2);
    }
}
