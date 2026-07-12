//! Heat-bath (exact conditional sampler) kernels.

use crate::algorithms::common::{Algorithm, SimulationPhase};
use crate::core::cache::EnergyPatch;
use crate::core::r#move::SiteSpinMove;
use crate::core::trial::TrialEvaluator;
use crate::core::visit::{SiteOrder, VisitSchedule};
use crate::lattice::interaction::{ContinuousHeatBathable, Hamiltonian, HeatBathable};
use crate::lattice::state::System;
use rand::Rng;

#[derive(Debug, Clone)]
pub struct HeatBathCore {
    order: SiteOrder,
    visit_schedule: VisitSchedule,
    patch: EnergyPatch,
}

impl Default for HeatBathCore {
    fn default() -> Self {
        Self::new()
    }
}

impl HeatBathCore {
    pub const fn new() -> Self {
        Self {
            order: SiteOrder::new(),
            visit_schedule: VisitSchedule::RandomPermutation,
            patch: EnergyPatch { delta_energy: 0.0 },
        }
    }

    pub fn with_visit_schedule(mut self, schedule: VisitSchedule) -> Self {
        self.visit_schedule = schedule;
        self
    }
}

impl<H: Hamiltonian + HeatBathable> Algorithm<H> for HeatBathCore {
    fn sweep_with_phase(
        &mut self,
        system: &mut System,
        model: &H,
        rng: &mut impl Rng,
        _phase: SimulationPhase,
    ) {
        let n_sites = system.n_sites();
        let spin_dim = model.spin_dim();
        let sites = self.order.prepare(n_sites, self.visit_schedule, rng);
        for &site in sites {
            let proposed =
                model.heat_bath_sample_site(&system.spins, &system.lattice, site, system.beta, rng);
            assert_eq!(
                proposed.len(),
                spin_dim,
                "heat-bath sample dimension mismatch"
            );
            let movement = SiteSpinMove::new(site, proposed);
            system.evaluate_trial(model, &movement, &mut self.patch);
            <System as TrialEvaluator<H, SiteSpinMove>>::commit_trial(
                system,
                &movement,
                &self.patch,
            );
        }
    }

    fn name(&self) -> &'static str {
        "HeatBath"
    }
}

#[derive(Debug, Clone)]
pub struct ContinuousHeatBathCore {
    order: SiteOrder,
    visit_schedule: VisitSchedule,
    patch: EnergyPatch,
}

impl Default for ContinuousHeatBathCore {
    fn default() -> Self {
        Self::new()
    }
}

impl ContinuousHeatBathCore {
    pub const fn new() -> Self {
        Self {
            order: SiteOrder::new(),
            visit_schedule: VisitSchedule::RandomPermutation,
            patch: EnergyPatch { delta_energy: 0.0 },
        }
    }

    pub fn with_visit_schedule(mut self, schedule: VisitSchedule) -> Self {
        self.visit_schedule = schedule;
        self
    }
}

impl<H: Hamiltonian + ContinuousHeatBathable> Algorithm<H> for ContinuousHeatBathCore {
    fn sweep_with_phase(
        &mut self,
        system: &mut System,
        model: &H,
        rng: &mut impl Rng,
        _phase: SimulationPhase,
    ) {
        let n_sites = system.n_sites();
        let spin_dim = model.spin_dim();
        let sites = self.order.prepare(n_sites, self.visit_schedule, rng);
        for &site in sites {
            let proposed =
                model.heat_bath_sample_site(&system.spins, &system.lattice, site, system.beta, rng);
            assert_eq!(
                proposed.len(),
                spin_dim,
                "heat-bath sample dimension mismatch"
            );
            let movement = SiteSpinMove::new(site, proposed);
            system.evaluate_trial(model, &movement, &mut self.patch);
            <System as TrialEvaluator<H, SiteSpinMove>>::commit_trial(
                system,
                &movement,
                &self.patch,
            );
        }
    }

    fn name(&self) -> &'static str {
        "ContinuousHeatBath"
    }
}
