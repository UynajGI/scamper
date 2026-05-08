//! Standard proposal strategy — delegates to model's propose_flip_spin / local_energy_change_spin.

use crate::algorithms::proposal_strategy::ProposalStrategy;
use crate::models::ModelMC;
use rand::Rng;

/// Standard proposal strategy using the model's native proposal mechanism.
pub struct StandardStrategy;

impl StandardStrategy {
    pub fn new() -> Self {
        StandardStrategy
    }
}

impl Default for StandardStrategy {
    fn default() -> Self {
        Self::new()
    }
}

impl<MC: ModelMC> ProposalStrategy<MC> for StandardStrategy {
    fn propose_flip(&mut self, model: &MC, site: usize, rng: &mut impl Rng) -> (Vec<f64>, Vec<f64>) {
        model.propose_flip_spin(site, rng)
    }

    fn compute_delta_e(&self, model: &MC, site: usize, old_spin: &[f64], new_spin: &[f64]) -> f64 {
        model.local_energy_change_spin(site, old_spin, new_spin)
    }
}
