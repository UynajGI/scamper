//! Exact microcanonical over-relaxation kernel.

use crate::algorithms::common::{Algorithm, SimulationPhase};
use crate::core::cache::EnergyPatch;
use crate::core::r#move::{SiteSpinMove, Spin};
use crate::core::trial::TrialEvaluator;
use crate::core::visit::{SiteOrder, VisitSchedule};
use crate::lattice::interaction::{Hamiltonian, LocalFieldModel};
use crate::lattice::state::System;
use rand::Rng;

#[derive(Debug, Clone)]
pub struct MicrocanonicalCore {
    order: SiteOrder,
    visit_schedule: VisitSchedule,
    field: Vec<f64>,
    patch: EnergyPatch,
}

impl Default for MicrocanonicalCore {
    fn default() -> Self {
        Self::new()
    }
}

impl MicrocanonicalCore {
    pub const fn new() -> Self {
        Self {
            order: SiteOrder::new(),
            visit_schedule: VisitSchedule::RandomPermutation,
            field: Vec::new(),
            patch: EnergyPatch { delta_energy: 0.0 },
        }
    }

    pub fn with_visit_schedule(mut self, schedule: VisitSchedule) -> Self {
        self.visit_schedule = schedule;
        self
    }
}

impl<H: Hamiltonian + LocalFieldModel> Algorithm<H> for MicrocanonicalCore {
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
        self.field.resize(spin_dim, 0.0);

        for &site in sites {
            model.local_field(&system.spins, &system.lattice, site, &mut self.field);
            let norm_squared = self.field.iter().map(|value| value * value).sum::<f64>();
            if norm_squared < 1e-28 {
                continue;
            }
            let old = Spin::from_slice(system.spin_at(site, spin_dim));
            let projection = old
                .iter()
                .zip(&self.field)
                .map(|(spin, field)| spin * field)
                .sum::<f64>()
                / norm_squared;
            let mut reflected = old.clone();
            for component in 0..spin_dim {
                reflected[component] = 2.0 * projection * self.field[component] - old[component];
            }
            let movement = SiteSpinMove::new(site, reflected);
            system.evaluate_trial(model, &movement, &mut self.patch);
            <System as TrialEvaluator<H, SiteSpinMove>>::commit_trial(
                system,
                &movement,
                &self.patch,
            );
        }
    }

    fn name(&self) -> &'static str {
        "Microcanonical"
    }
}
