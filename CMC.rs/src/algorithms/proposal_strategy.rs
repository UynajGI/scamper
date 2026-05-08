//! Proposal strategies for Metropolis algorithm.
//!
//! A strategy defines how spins are proposed and energy changes computed.
//! This allows different proposal mechanisms (standard, OPSS, etc.) to plug
//! into the same MetropolisCore algorithm.

use crate::models::ModelMC;
use rand::Rng;

/// A proposal strategy for use with MetropolisCore.
pub trait ProposalStrategy<MC: ModelMC>: Send {
    /// Propose a new spin at the given site.
    fn propose_flip(&mut self, model: &MC, site: usize, rng: &mut impl Rng) -> (Vec<f64>, Vec<f64>);

    /// Compute energy change for the proposed flip.
    fn compute_delta_e(&self, model: &MC, site: usize, old_spin: &[f64], new_spin: &[f64]) -> f64;

    /// Called after each sweep for strategies that need per-sweep adaptation.
    /// Default: no-op. OPSSStrategy overrides to adapt sigma.
    fn adapt_after_sweep(&mut self, _model: &mut MC) {}
}
