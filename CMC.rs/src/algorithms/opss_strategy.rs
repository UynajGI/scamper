//! OPSS (Optimal Phase Space Sampling) strategy for continuous spin models.
//!
//! Gaussian move with adaptive sigma. S'_i = (S_i + sigma * F) / |S_i + sigma * F|
//! where F is a Gaussian random vector. Based on Alzate-Cardona et al. (2018).
//! Applicable to any O(N) model with spin_dim() >= 2.

use crate::algorithms::proposal_strategy::ProposalStrategy;
use crate::models::ModelMC;
use rand::Rng;
use rand::RngExt;

/// Sample from standard normal distribution using Box-Muller transform.
fn sample_gaussian(rng: &mut impl Rng) -> f64 {
    let u1: f64 = rng.random_range(0.0001..1.0);
    let u2 = rng.random::<f64>();
    (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
}

/// OPSS proposal strategy for continuous spin models (O(N), N >= 2).
/// Uses Gaussian perturbation + normalization to unit sphere, with adaptive sigma.
pub struct OPSSStrategy {
    sigma: f64,
    accepted: u64,
    total: u64,
    sigma_min: f64,
    sigma_max: f64,
}

impl OPSSStrategy {
    pub fn new(initial_sigma: f64) -> Self {
        OPSSStrategy {
            sigma: initial_sigma,
            accepted: 0,
            total: 0,
            sigma_min: 1e-6,
            sigma_max: 60.0,
        }
    }

    pub fn sigma(&self) -> f64 {
        self.sigma
    }

    fn adapt_sigma(&mut self) {
        if self.total == 0 {
            return;
        }
        let rate = self.accepted as f64 / self.total as f64;
        if rate >= 1.0 {
            self.sigma = self.sigma_max;
        } else {
            let f = 0.5 / (1.0 - rate);
            self.sigma *= f;
        }
        self.sigma = self.sigma.clamp(self.sigma_min, self.sigma_max);
        self.accepted = 0;
        self.total = 0;
    }
}

impl<MC: ModelMC> ProposalStrategy<MC> for OPSSStrategy {
    fn propose_flip(
        &mut self,
        model: &MC,
        site: usize,
        rng: &mut impl Rng,
    ) -> (Vec<f64>, Vec<f64>) {
        let dim = model.spin_dim();
        let old: Vec<f64> = (0..dim)
            .map(|d| model.spins()[site * dim + d])
            .collect();

        let new: Vec<f64> = old
            .iter()
            .map(|&s| s + self.sigma * sample_gaussian(rng))
            .collect();

        let norm = new.iter().map(|&x| x * x).sum::<f64>().sqrt();
        let normalized: Vec<f64> = new.iter().map(|&x| x / norm).collect();
        (old, normalized)
    }

    fn compute_delta_e(
        &self,
        model: &MC,
        site: usize,
        old_spin: &[f64],
        new_spin: &[f64],
    ) -> f64 {
        model.local_energy_change_spin(site, old_spin, new_spin)
    }

    fn adapt_after_sweep(&mut self, _model: &mut MC) {
        self.adapt_sigma();
    }
}
