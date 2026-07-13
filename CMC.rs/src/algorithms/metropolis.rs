//! Metropolis-Hastings local-update kernel.

use crate::algorithms::common::{Algorithm, SimulationPhase};
use crate::audit::{audit_lattice_cache, should_audit_cache};
use crate::core::acceptance::MetropolisHastingsAcceptance;
use crate::core::cache::EnergyPatch;
use crate::core::r#move::SiteSpinMove;
use crate::core::trial::{metropolis_hastings_step, ProposedMove};
use crate::core::visit::{SiteOrder, VisitSchedule};
use crate::lattice::interaction::{Hamiltonian, Proposable};
use crate::lattice::proposal::{ProposalStrategy, StandardStrategy};
use crate::lattice::state::System;
use rand::Rng;

#[derive(Debug, Clone)]
pub struct MetropolisCore<S = StandardStrategy> {
    pub strategy: S,
    order: SiteOrder,
    visit_schedule: VisitSchedule,
    patch: EnergyPatch,
    energy_check_interval: u64,
    sweeps: u64,
}

impl MetropolisCore<StandardStrategy> {
    pub fn new() -> Self {
        Self::with_strategy(StandardStrategy::new())
    }
}

impl<S: Default> Default for MetropolisCore<S> {
    fn default() -> Self {
        Self::with_strategy(S::default())
    }
}

impl<S> MetropolisCore<S> {
    pub fn with_strategy(strategy: S) -> Self {
        Self {
            strategy,
            order: SiteOrder::new(),
            visit_schedule: VisitSchedule::RandomPermutation,
            patch: EnergyPatch::default(),
            energy_check_interval: 0,
            sweeps: 0,
        }
    }

    pub fn with_visit_schedule(mut self, schedule: VisitSchedule) -> Self {
        self.visit_schedule = schedule;
        self
    }

    /// Periodically validate the cached energy against an exact recomputation.
    /// Zero selects the build-mode automatic audit policy.
    pub fn with_energy_check_interval(mut self, interval: u64) -> Self {
        self.energy_check_interval = interval;
        self
    }
}

impl<H, S> Algorithm<H> for MetropolisCore<S>
where
    H: Hamiltonian + Proposable,
    S: ProposalStrategy<H>,
{
    fn sweep_with_phase(
        &mut self,
        system: &mut System,
        model: &H,
        rng: &mut impl Rng,
        phase: SimulationPhase,
    ) {
        let n_sites = system.n_sites();
        let spin_dim = model.spin_dim();
        let ensemble = system.canonical_ensemble();
        let sites = self.order.prepare(n_sites, self.visit_schedule, rng);

        let acceptance = MetropolisHastingsAcceptance;
        for &site in sites {
            let proposal = self.strategy.propose(model, system, site, rng);
            assert_eq!(
                proposal.spin.len(),
                spin_dim,
                "proposal dimension does not match the model"
            );
            let movement = SiteSpinMove::new(site, proposal.spin);
            let proposal = ProposedMove::new(movement, proposal.log_reverse_over_forward);
            let outcome = metropolis_hastings_step(
                system,
                model,
                &proposal,
                &ensemble,
                &acceptance,
                &mut self.patch,
                rng,
            );
            self.strategy.record_result(outcome.accepted);
        }

        self.strategy.finish_sweep(phase.allows_adaptation());
        self.sweeps = self.sweeps.wrapping_add(1);
        if should_audit_cache(self.sweeps, self.energy_check_interval) {
            audit_lattice_cache(system, model).expect("Metropolis lattice cache audit failed");
        }
    }

    fn name(&self) -> &'static str {
        "Metropolis-Hastings"
    }
}
