//! NVT Metropolis-Hastings kernel for single-particle translations.

use crate::algorithms::SimulationPhase;
use crate::audit::{audit_particle_cache, should_audit_cache};
use crate::core::acceptance::MetropolisHastingsAcceptance;
use crate::core::trial::metropolis_hastings_step;
use crate::core::visit::{SiteOrder, VisitSchedule};
use crate::particle::{PairPotential, ParticleEnergyPatch, ParticleSystem, TranslateParticle};
use rand::Rng;

/// Marker for fixed-N, fixed-V kernels whose β-tempering weight is `-βE`.
///
/// NPT and μVT kernels deliberately do not implement this marker because their
/// replica-exchange weights require pressure-volume or activity terms.
pub trait CanonicalParticleKernel {}

/// Update policy for a continuous particle state.
pub trait ParticleAlgorithm<const D: usize, P: PairPotential>: Send {
    /// Execute one particle sweep in the supplied lifecycle phase.
    fn sweep_with_phase(
        &mut self,
        system: &mut ParticleSystem<D>,
        potential: &P,
        rng: &mut impl Rng,
        phase: SimulationPhase,
    );

    /// Direct/manual sweeps use frozen production semantics.
    fn sweep(&mut self, system: &mut ParticleSystem<D>, potential: &P, rng: &mut impl Rng) {
        self.sweep_with_phase(system, potential, rng, SimulationPhase::Measurement);
    }

    /// Human-readable kernel name.
    fn name(&self) -> &'static str {
        "Particle algorithm"
    }
}

/// One attempted translation per particle, in a configurable visit order.
#[derive(Debug, Clone)]
pub struct ParticleMetropolisCore<const D: usize> {
    translation: TranslateParticle<D>,
    order: SiteOrder,
    visit_schedule: VisitSchedule,
    patch: ParticleEnergyPatch,
    energy_check_interval: u64,
    sweeps: u64,
    last_phase: Option<SimulationPhase>,
}

impl<const D: usize> ParticleMetropolisCore<D> {
    /// Construct a translation kernel.
    pub fn new(max_displacement: f64) -> Self {
        Self {
            translation: TranslateParticle::new(max_displacement),
            order: SiteOrder::new(),
            visit_schedule: VisitSchedule::RandomPermutation,
            patch: ParticleEnergyPatch::default(),
            energy_check_interval: 0,
            sweeps: 0,
            last_phase: None,
        }
    }

    /// Replace the translation proposal, including its adaptation policy.
    pub fn with_translation(mut self, translation: TranslateParticle<D>) -> Self {
        self.translation = translation;
        self
    }

    /// Select particle visit order.
    pub fn with_visit_schedule(mut self, schedule: VisitSchedule) -> Self {
        self.visit_schedule = schedule;
        self
    }

    /// Periodically audit energy and packed cell membership; zero disables it.
    pub fn with_energy_check_interval(mut self, interval: u64) -> Self {
        self.energy_check_interval = interval;
        self
    }

    /// Translation proposal and adaptation statistics.
    #[inline]
    pub const fn translation(&self) -> &TranslateParticle<D> {
        &self.translation
    }

    /// Mutable access for explicit warmup configuration.
    #[inline]
    pub fn translation_mut(&mut self) -> &mut TranslateParticle<D> {
        &mut self.translation
    }
}

impl<const D: usize> Default for ParticleMetropolisCore<D> {
    fn default() -> Self {
        Self::new(0.1)
    }
}

impl<const D: usize> CanonicalParticleKernel for ParticleMetropolisCore<D> {}

impl<const D: usize, P: PairPotential> ParticleAlgorithm<D, P> for ParticleMetropolisCore<D> {
    fn sweep_with_phase(
        &mut self,
        system: &mut ParticleSystem<D>,
        potential: &P,
        rng: &mut impl Rng,
        phase: SimulationPhase,
    ) {
        if self.last_phase != Some(phase) {
            if phase == SimulationPhase::Measurement {
                self.translation.reset_statistics();
            }
            self.last_phase = Some(phase);
        }

        let ensemble = system.canonical_ensemble();
        let particles = self.order.prepare(system.len(), self.visit_schedule, rng);
        let acceptance = MetropolisHastingsAcceptance;
        for &particle in particles {
            let proposal = self
                .translation
                .propose(system.configuration(), particle, rng);
            let outcome = metropolis_hastings_step(
                system,
                potential,
                &proposal,
                &ensemble,
                &acceptance,
                &mut self.patch,
                rng,
            );
            self.translation.record_result(outcome.accepted);
        }
        self.translation.finish_sweep(phase.allows_adaptation());

        self.sweeps = self.sweeps.wrapping_add(1);
        if should_audit_cache(self.sweeps, self.energy_check_interval) {
            audit_particle_cache(system, potential)
                .expect("particle energy/cell-list audit failed");
        }
    }

    fn name(&self) -> &'static str {
        "Particle Metropolis-Hastings"
    }
}
