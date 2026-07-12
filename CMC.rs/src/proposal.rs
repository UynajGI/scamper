//! Metropolis-Hastings proposal policies.

use crate::hamiltonian::{Hamiltonian, Proposable, Spin};
use crate::system::System;
use rand::{Rng, RngExt};

/// Complete proposal including the Hastings correction.
#[derive(Debug, Clone)]
pub struct ProposedSpin {
    pub spin: Spin,
    /// `ln q(old | new) - ln q(new | old)`.
    pub log_reverse_over_forward: f64,
}

impl ProposedSpin {
    pub fn symmetric(spin: Spin) -> Self {
        Self {
            spin,
            log_reverse_over_forward: 0.0,
        }
    }
}

/// Proposal strategy for [`MetropolisCore`](crate::MetropolisCore).
pub trait ProposalStrategy<H: Hamiltonian>: Send {
    fn propose(
        &mut self,
        model: &H,
        system: &System,
        site: usize,
        rng: &mut impl Rng,
    ) -> ProposedSpin;

    /// Receives every accept/reject decision, enabling correct adaptation.
    fn record_result(&mut self, _accepted: bool) {}

    /// Called after a full sweep.  `adaptation_enabled` is true only during
    /// Carlo.rs thermalization; the transition kernel is frozen afterwards.
    fn finish_sweep(&mut self, _adaptation_enabled: bool) {}
}

/// Independent model proposal.  It is assumed symmetric.
#[derive(Debug, Clone, Default)]
pub struct StandardStrategy;

impl StandardStrategy {
    pub const fn new() -> Self {
        Self
    }
}

impl<H: Hamiltonian + Proposable> ProposalStrategy<H> for StandardStrategy {
    fn propose(
        &mut self,
        model: &H,
        system: &System,
        site: usize,
        rng: &mut impl Rng,
    ) -> ProposedSpin {
        ProposedSpin::symmetric(model.propose_from(system.spin_at(site, model.spin_dim()), rng))
    }
}

/// Backward-compatible adaptive local-rotation strategy.
///
/// The old implementation used a non-involutive normalized reflection and did
/// not receive acceptance results.  This version proposes an exactly symmetric
/// random plane rotation for vector spins, records all decisions, and adapts
/// only during Carlo.rs thermalization.  For scalar models it falls back to the
/// model's standard proposal.
#[derive(Debug, Clone)]
pub struct OPSSStrategy {
    /// Maximum absolute rotation angle in radians.
    pub sigma: f64,
    pub target_acceptance: f64,
    accepted: u64,
    attempted: u64,
    adaptation_rate: f64,
}

impl OPSSStrategy {
    pub fn new() -> Self {
        Self {
            sigma: 0.5,
            target_acceptance: 0.5,
            accepted: 0,
            attempted: 0,
            adaptation_rate: 0.08,
        }
    }

    pub fn with_target(mut self, target: f64) -> Self {
        self.target_acceptance = target.clamp(0.05, 0.95);
        self
    }

    pub fn with_sigma(mut self, sigma: f64) -> Self {
        self.sigma = sigma.clamp(1e-4, std::f64::consts::PI);
        self
    }

    pub fn acceptance_counts(&self) -> (u64, u64) {
        (self.accepted, self.attempted)
    }
}

impl Default for OPSSStrategy {
    fn default() -> Self {
        Self::new()
    }
}

impl<H: Hamiltonian + Proposable> ProposalStrategy<H> for OPSSStrategy {
    fn propose(
        &mut self,
        model: &H,
        system: &System,
        site: usize,
        rng: &mut impl Rng,
    ) -> ProposedSpin {
        let dimension = model.spin_dim();
        let old = system.spin_at(site, dimension);
        if dimension < 2 {
            return ProposedSpin::symmetric(model.propose_from(old, rng));
        }

        let mut spin = Spin::from_slice(old);
        let first = rng.random_range(0..dimension);
        let mut second = rng.random_range(0..dimension - 1);
        if second >= first {
            second += 1;
        }
        let sigma = if self.sigma.is_finite() {
            self.sigma.abs().clamp(1e-4, std::f64::consts::PI)
        } else {
            0.5
        };
        let angle = rng.random_range(-sigma..sigma);
        let (sine, cosine) = angle.sin_cos();
        let left = spin[first];
        let right = spin[second];
        spin[first] = cosine * left - sine * right;
        spin[second] = sine * left + cosine * right;
        model.normalize_spin(&mut spin);
        ProposedSpin::symmetric(spin)
    }

    fn record_result(&mut self, accepted: bool) {
        self.attempted += 1;
        if accepted {
            self.accepted += 1;
        }
    }

    fn finish_sweep(&mut self, adaptation_enabled: bool) {
        if self.attempted == 0 {
            return;
        }
        if adaptation_enabled {
            let rate = self.accepted as f64 / self.attempted as f64;
            // Multiplicative Robbins-Monro style update keeps sigma positive.
            let shift = self.adaptation_rate * (rate - self.target_acceptance);
            self.sigma = (self.sigma * shift.exp()).clamp(1e-4, std::f64::consts::PI);
        }
        self.accepted = 0;
        self.attempted = 0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{build_chain, models::XYModel, System};
    use rand::SeedableRng;

    #[test]
    fn standard_ising_proposal_is_never_a_noop() {
        let model = crate::models::IsingModel::new(1.0);
        let system = System::new(build_chain(2, false), 1, 1.0, 1.0);
        let mut strategy = StandardStrategy::new();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(4);
        let proposal = strategy.propose(&model, &system, 0, &mut rng);
        assert_eq!(proposal.spin.as_slice(), &[-1.0]);
    }

    #[test]
    fn adaptive_strategy_receives_results() {
        let model = XYModel::new(1.0);
        let system = System::new(build_chain(2, false), 2, 0.0, 1.0);
        let mut strategy = OPSSStrategy::new();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(3);
        let proposal = strategy.propose(&model, &system, 0, &mut rng);
        assert_eq!(proposal.spin.len(), 2);
        <OPSSStrategy as ProposalStrategy<XYModel>>::record_result(&mut strategy, true);
        assert_eq!(strategy.acceptance_counts(), (1, 1));
    }

    #[test]
    fn adaptive_strategy_freezes_outside_thermalization() {
        let mut strategy = OPSSStrategy::new().with_sigma(0.4);
        let initial = strategy.sigma;
        <OPSSStrategy as ProposalStrategy<XYModel>>::record_result(&mut strategy, true);
        <OPSSStrategy as ProposalStrategy<XYModel>>::finish_sweep(&mut strategy, false);
        assert_eq!(strategy.sigma, initial);

        <OPSSStrategy as ProposalStrategy<XYModel>>::record_result(&mut strategy, true);
        <OPSSStrategy as ProposalStrategy<XYModel>>::finish_sweep(&mut strategy, true);
        assert!(strategy.sigma > initial);
    }
}
