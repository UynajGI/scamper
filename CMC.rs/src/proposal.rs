//! Proposal strategies for Metropolis algorithm.
//!
//! A [`ProposalStrategy`] defines how new spins are proposed and can adapt
//! its parameters based on sweep history (e.g., OPSS sigma tuning).

use crate::hamiltonian::{Hamiltonian, Proposable};
use crate::system::System;
use rand::Rng;
use smallvec::{smallvec, SmallVec};

/// A proposal strategy for use with [`MetropolisCore`](crate::MetropolisCore).
pub trait ProposalStrategy<H: Hamiltonian>: Send {
    /// Propose a new spin at the given site.
    fn propose(
        &mut self,
        model: &H,
        _system: &System,
        _site: usize,
        rng: &mut impl Rng,
    ) -> SmallVec<[f64; 3]>;

    /// Called after each sweep for strategies that need per-sweep adaptation.
    fn adapt_after_sweep(&mut self, _model: &H) {}
}

// ── Standard ────────────────────────────────────────────────

/// Standard proposal strategy — delegates entirely to `model.propose()`.
#[derive(Debug, Clone, Default)]
pub struct StandardStrategy;

impl StandardStrategy {
    pub fn new() -> Self {
        Self
    }
}

impl<H: Hamiltonian + Proposable> ProposalStrategy<H> for StandardStrategy {
    fn propose(&mut self, model: &H, _system: &System, _site: usize, rng: &mut impl Rng) -> SmallVec<[f64; 3]> {
        model.propose(rng)
    }
}

// ── OPSS (Over-relaxation) ──────────────────────────────────

/// Over-relaxed pseudo-spin-flip strategy (Yunoki & Sorella).
///
/// Mirrors the spin about the local mean-field direction. The reflection angle
/// is controlled by `sigma`, which is auto-tuned to achieve a target acceptance rate.
///
/// Works for vector spins (XY, Heisenberg) and scalar spins (Ising) — for Ising,
/// OPSS reduces to a deterministic sign flip when the local field is small.
#[derive(Debug, Clone)]
pub struct OPSSStrategy {
    /// Reflection strength (auto-adapted).
    pub sigma: f64,
    /// Target acceptance rate.
    pub target_acceptance: f64,
    /// Running acceptance counter.
    accepted: u64,
    attempted: u64,
}

impl OPSSStrategy {
    pub fn new() -> Self {
        Self {
            sigma: 1.0,
            target_acceptance: 0.5,
            accepted: 0,
            attempted: 0,
        }
    }

    pub fn with_target(mut self, target: f64) -> Self {
        self.target_acceptance = target;
        self
    }
}

impl Default for OPSSStrategy {
    fn default() -> Self {
        Self::new()
    }
}

impl<H: Hamiltonian + Proposable> ProposalStrategy<H> for OPSSStrategy {
    fn propose(&mut self, model: &H, system: &System, site: usize, rng: &mut impl Rng) -> SmallVec<[f64; 3]> {
        let sd = model.spin_dim();
        let old = system.spin_at(site, sd);

        if sd == 1 {
            // Scalar: reflect about local field
            let local_field: f64 = system.lattice.neighbors(site)
                .iter()
                .map(|&nb| system.spins[nb])
                .sum();
            let proposed = -self.sigma * old[0] + (1.0 + self.sigma) * local_field.signum();
            smallvec![proposed.signum()] // always ±1 for Ising
        } else {
            // Vector: reflect spin about local field direction
            let mut h = vec![0.0; sd];
            for &nb in system.lattice.neighbors(site) {
                let base = nb * sd;
                for (k, hk) in h.iter_mut().enumerate() {
                    *hk += system.spins[base + k];
                }
            }

            let h_norm: f64 = h.iter().map(|&x| x * x).sum::<f64>().sqrt();
            if h_norm < 1e-12 {
                // No local field → random unit vector
                let mut v = model.propose(rng);
                model.normalize_spin(&mut v);
                return v;
            }

            let h_hat: SmallVec<[f64; 3]> = h.iter().map(|&x| x / h_norm).collect();

            // s_new = 2 (s·ĥ) ĥ - σ s
            let s_dot_h: f64 = old.iter().zip(&h_hat).map(|(&s, &h)| s * h).sum();
            let mut new: SmallVec<[f64; 3]> = SmallVec::from_elem(0.0, sd);
            for k in 0..sd {
                new[k] = 2.0 * s_dot_h * h_hat[k] - self.sigma * old[k];
            }

            model.normalize_spin(&mut new);
            new
        }
    }

    fn adapt_after_sweep(&mut self, _model: &H) {
        if self.attempted == 0 {
            return;
        }
        let rate = self.accepted as f64 / self.attempted as f64;
        // Adjust sigma to approach target acceptance
        if rate > self.target_acceptance + 0.05 {
            self.sigma = (self.sigma * 1.05).min(2.0);
        } else if rate < self.target_acceptance - 0.05 {
            self.sigma = (self.sigma * 0.95).max(0.1);
        }
        self.accepted = 0;
        self.attempted = 0;
    }
}

impl OPSSStrategy {
    /// Called by MetropolisCore to track acceptance.
    pub fn record_acceptance(&mut self, accepted: bool) {
        self.attempted += 1;
        if accepted {
            self.accepted += 1;
        }
    }
}
