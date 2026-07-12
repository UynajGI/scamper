//! Acceptance-rule abstraction for Metropolis-Hastings transition kernels.
//!
//! Separating the acceptance formula from trial evaluation lets future kernels
//! (Barker, rejection-free, etc.) reuse the same proposal and ensemble logic
//! without modifying the trial layer.

use crate::core::ensemble::Ensemble;

/// Converts ensemble weight ratios and proposal asymmetries into a log
/// acceptance probability.
pub trait AcceptanceRule<D> {
    fn log_acceptance(
        &self,
        ensemble: &impl Ensemble<D>,
        delta: &D,
        log_proposal_ratio: f64,
    ) -> f64;
}

/// Standard Metropolis-Hastings formula.
///
/// `log P_accept = ln π(new)/π(old) + ln q(old|new)/q(new|old)`
pub struct MetropolisHastingsAcceptance;

impl<D> AcceptanceRule<D> for MetropolisHastingsAcceptance {
    #[inline]
    fn log_acceptance(
        &self,
        ensemble: &impl Ensemble<D>,
        delta: &D,
        log_proposal_ratio: f64,
    ) -> f64 {
        ensemble.log_weight_ratio(delta) + log_proposal_ratio
    }
}
