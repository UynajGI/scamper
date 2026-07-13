//! Grand-canonical insertion/deletion moves and μVT transition kernels.

use crate::algorithms::SimulationPhase;
use crate::audit::{audit_particle_cache, should_audit_cache};
use crate::core::acceptance::MetropolisHastingsAcceptance;
use crate::core::ensemble::{GrandCanonical, ThermodynamicDelta};
use crate::core::trial::{metropolis_hastings_step, ProposedMove, TrialEvaluator};
use crate::core::visit::{SiteOrder, VisitSchedule};
use crate::particle::{
    PairPotential, ParticleAlgorithm, ParticleEnergyPatch, ParticleError, ParticleSystem,
    SimulationCell, TranslateParticle,
};
use rand::{Rng, RngExt};
use std::collections::BTreeSet;

/// Insert one particle at a continuous position.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ParticleInsertion<const D: usize> {
    pub species: u16,
    pub position: [f64; D],
}

/// Delete one indexed particle using swap-remove semantics.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ParticleDeletion {
    pub particle: usize,
}

/// Particle-number-changing move.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GrandCanonicalMove<const D: usize> {
    Insert(ParticleInsertion<D>),
    Delete(ParticleDeletion),
}

/// Reusable local-energy and packed-cache patch for insertion/deletion.
#[derive(Debug, Clone, Default)]
pub struct GrandCanonicalPatch {
    /// Physical energy change.
    pub delta_energy: f64,
    cell: usize,
    candidates: Vec<usize>,
    valid: bool,
}

impl<const D: usize> ParticleSystem<D> {
    pub(crate) fn evaluate_grand_canonical<P: PairPotential>(
        &self,
        potential: &P,
        movement: &GrandCanonicalMove<D>,
        patch: &mut GrandCanonicalPatch,
    ) -> ThermodynamicDelta {
        patch.valid = false;
        match movement {
            GrandCanonicalMove::Insert(insertion) => {
                if insertion
                    .position
                    .iter()
                    .any(|coordinate| !coordinate.is_finite())
                    || !potential.supports_species(insertion.species)
                {
                    patch.delta_energy = f64::INFINITY;
                    return ThermodynamicDelta::particle(f64::INFINITY, 1);
                }
                patch.cell = self
                    .cell_list
                    .cell_of_position(&insertion.position, self.configuration.cell());
                self.cell_list
                    .fill_candidates(patch.cell, &mut patch.candidates);
                let mut local_energy = 0.0;
                for &other in &patch.candidates {
                    let distance_squared = self
                        .configuration
                        .cell()
                        .distance_squared(&insertion.position, self.configuration.position(other));
                    local_energy += potential.energy(
                        insertion.species,
                        self.configuration.species_of(other),
                        distance_squared,
                    );
                }
                if !local_energy.is_finite() {
                    patch.delta_energy = f64::INFINITY;
                    return ThermodynamicDelta::particle(f64::INFINITY, 1);
                }
                patch.delta_energy = local_energy;
                patch.valid = true;
                ThermodynamicDelta::particle(local_energy, 1)
            }
            GrandCanonicalMove::Delete(deletion) => {
                assert!(
                    deletion.particle < self.len(),
                    "deletion particle out of range"
                );
                patch.cell = self.cell_list.particle_cell(deletion.particle);
                self.cell_list
                    .fill_candidates(patch.cell, &mut patch.candidates);
                let species = self.configuration.species_of(deletion.particle);
                let mut local_energy = 0.0;
                for &other in &patch.candidates {
                    if other == deletion.particle {
                        continue;
                    }
                    let distance_squared = self.configuration.cell().distance_squared(
                        self.configuration.position(deletion.particle),
                        self.configuration.position(other),
                    );
                    local_energy += potential.energy(
                        species,
                        self.configuration.species_of(other),
                        distance_squared,
                    );
                }
                assert!(
                    local_energy.is_finite(),
                    "accepted deletion energy is non-finite"
                );
                patch.delta_energy = -local_energy;
                patch.valid = true;
                ThermodynamicDelta::particle(-local_energy, -1)
            }
        }
    }

    pub(crate) fn commit_grand_canonical(
        &mut self,
        movement: &GrandCanonicalMove<D>,
        patch: &GrandCanonicalPatch,
    ) {
        assert!(
            patch.valid,
            "invalid insertion/deletion trial must never be committed"
        );
        match movement {
            GrandCanonicalMove::Insert(insertion) => {
                let particle = self
                    .configuration
                    .push_particle(insertion.position, insertion.species);
                self.cell_list.insert_particle(particle, patch.cell);
            }
            GrandCanonicalMove::Delete(deletion) => {
                self.cell_list.remove_particle_swap(deletion.particle);
                self.configuration.swap_remove_particle(deletion.particle);
            }
        }
        self.energy += patch.delta_energy;
        assert!(
            self.energy.is_finite(),
            "μVT commit produced non-finite energy"
        );
    }
}

impl<const D: usize, P: PairPotential> TrialEvaluator<P, GrandCanonicalMove<D>>
    for ParticleSystem<D>
{
    type Delta = ThermodynamicDelta;
    type Patch = GrandCanonicalPatch;

    fn evaluate_trial(
        &self,
        model: &P,
        movement: &GrandCanonicalMove<D>,
        patch: &mut Self::Patch,
    ) -> Self::Delta {
        self.evaluate_grand_canonical(model, movement, patch)
    }

    fn commit_trial(&mut self, movement: &GrandCanonicalMove<D>, patch: &Self::Patch) {
        self.commit_grand_canonical(movement, patch);
    }
}

/// Reversible insertion/deletion proposal with explicit branch and species factors.
#[derive(Debug, Clone)]
pub struct InsertDeleteParticle {
    species: Vec<(u16, f64)>,
    total_species_weight: f64,
    insertion_probability: f64,
    minimum_particles: usize,
    maximum_particles: Option<usize>,
}

impl InsertDeleteParticle {
    /// Construct a one-species proposal with equal insertion/deletion branches.
    pub fn new(species: u16) -> Self {
        Self {
            species: vec![(species, 1.0)],
            total_species_weight: 1.0,
            insertion_probability: 0.5,
            minimum_particles: 0,
            maximum_particles: None,
        }
    }

    /// Construct a multi-species proposal. Weights affect proposals, not target activity.
    pub fn from_species(species: Vec<(u16, f64)>) -> Result<Self, ParticleError> {
        if species.is_empty() {
            return Err(ParticleError::InvalidMove(
                "insertion proposal requires at least one species".to_string(),
            ));
        }
        let mut seen = BTreeSet::new();
        let mut total = 0.0;
        for &(label, weight) in &species {
            if !seen.insert(label) {
                return Err(ParticleError::InvalidMove(format!(
                    "duplicate insertion species {label}"
                )));
            }
            if !weight.is_finite() || weight <= 0.0 {
                return Err(ParticleError::InvalidMove(
                    "species proposal weights must be finite and positive".to_string(),
                ));
            }
            total += weight;
        }
        if !total.is_finite() {
            return Err(ParticleError::InvalidMove(
                "species proposal weight sum is non-finite".to_string(),
            ));
        }
        Ok(Self {
            species,
            total_species_weight: total,
            insertion_probability: 0.5,
            minimum_particles: 0,
            maximum_particles: None,
        })
    }

    pub fn with_insertion_probability(mut self, probability: f64) -> Self {
        assert!(probability.is_finite() && (0.0..1.0).contains(&probability));
        self.insertion_probability = probability;
        self
    }

    pub fn with_particle_bounds(
        mut self,
        minimum_particles: usize,
        maximum_particles: Option<usize>,
    ) -> Result<Self, ParticleError> {
        if maximum_particles.is_some_and(|maximum| maximum <= minimum_particles) {
            return Err(ParticleError::InvalidMove(
                "maximum particle count must exceed the minimum".to_string(),
            ));
        }
        self.minimum_particles = minimum_particles;
        self.maximum_particles = maximum_particles;
        Ok(self)
    }

    /// Verify that every proposed insertion species is defined by the potential.
    pub fn validate_potential<P: PairPotential>(&self, potential: &P) -> Result<(), ParticleError> {
        for &(species, _) in &self.species {
            if !potential.supports_species(species) {
                return Err(ParticleError::InvalidMove(format!(
                    "insertion species {species} is not supported by the pair potential"
                )));
            }
        }
        Ok(())
    }

    pub fn validate_state<const D: usize>(
        &self,
        system: &ParticleSystem<D>,
    ) -> Result<(), ParticleError> {
        if system.len() < self.minimum_particles
            || self
                .maximum_particles
                .is_some_and(|maximum| system.len() > maximum)
        {
            return Err(ParticleError::InvalidMove(
                "state particle count lies outside insertion/deletion bounds".to_string(),
            ));
        }
        for &species in system.configuration().species() {
            if self.species_probability(species).is_none() {
                return Err(ParticleError::InvalidMove(format!(
                    "accepted species {species} has zero reverse insertion probability"
                )));
            }
        }
        Ok(())
    }

    pub fn propose<const D: usize>(
        &self,
        system: &ParticleSystem<D>,
        rng: &mut impl Rng,
    ) -> ProposedMove<GrandCanonicalMove<D>> {
        let particles = system.len();
        let (insert_probability, delete_probability) = self.branch_probabilities(particles);
        let insert = delete_probability == 0.0
            || (insert_probability > 0.0 && rng.random::<f64>() < insert_probability);
        if insert {
            let species = self.sample_species(rng);
            let species_probability = self
                .species_probability(species)
                .expect("sampled species has positive probability");
            let mut position = [0.0; D];
            for (axis, coordinate) in position.iter_mut().enumerate() {
                *coordinate = rng.random_range(0.0..system.configuration().cell().lengths()[axis]);
            }
            let next_particles = particles + 1;
            let (_, reverse_delete_probability) = self.branch_probabilities(next_particles);
            let log_ratio = reverse_delete_probability.ln()
                - (next_particles as f64).ln()
                - insert_probability.ln()
                - species_probability.ln()
                + system.configuration().cell().volume().ln();
            ProposedMove::new(
                GrandCanonicalMove::Insert(ParticleInsertion { species, position }),
                log_ratio,
            )
        } else {
            let particle = rng.random_range(0..particles);
            let species = system.configuration().species_of(particle);
            let species_probability = self.species_probability(species).unwrap_or(0.0);
            let previous_particles = particles - 1;
            let (reverse_insert_probability, _) = self.branch_probabilities(previous_particles);
            let log_ratio = reverse_insert_probability.ln() + species_probability.ln()
                - system.configuration().cell().volume().ln()
                - delete_probability.ln()
                + (particles as f64).ln();
            ProposedMove::new(
                GrandCanonicalMove::Delete(ParticleDeletion { particle }),
                log_ratio,
            )
        }
    }

    fn branch_probabilities(&self, particles: usize) -> (f64, f64) {
        if particles <= self.minimum_particles {
            (1.0, 0.0)
        } else if self
            .maximum_particles
            .is_some_and(|maximum| particles >= maximum)
        {
            (0.0, 1.0)
        } else {
            (self.insertion_probability, 1.0 - self.insertion_probability)
        }
    }

    fn sample_species(&self, rng: &mut impl Rng) -> u16 {
        let threshold = rng.random_range(0.0..self.total_species_weight);
        let mut cumulative = 0.0;
        for &(species, weight) in &self.species {
            cumulative += weight;
            if threshold < cumulative {
                return species;
            }
        }
        self.species.last().expect("non-empty species table").0
    }

    fn species_probability(&self, species: u16) -> Option<f64> {
        self.species
            .iter()
            .find(|(label, _)| *label == species)
            .map(|(_, weight)| *weight / self.total_species_weight)
    }
}

/// μVT kernel: canonical translations plus insertion/deletion attempts.
#[derive(Debug, Clone)]
pub struct ParticleGrandCanonicalCore<const D: usize> {
    translation: TranslateParticle<D>,
    exchange: InsertDeleteParticle,
    log_activity: f64,
    exchange_attempts_per_sweep: u64,
    order: SiteOrder,
    visit_schedule: VisitSchedule,
    translation_patch: ParticleEnergyPatch,
    exchange_patch: GrandCanonicalPatch,
    energy_check_interval: u64,
    sweeps: u64,
}

impl<const D: usize> ParticleGrandCanonicalCore<D> {
    pub fn new(max_displacement: f64, exchange: InsertDeleteParticle, log_activity: f64) -> Self {
        assert!(log_activity.is_finite(), "log activity must be finite");
        Self {
            translation: TranslateParticle::new(max_displacement),
            exchange,
            log_activity,
            exchange_attempts_per_sweep: 1,
            order: SiteOrder::new(),
            visit_schedule: VisitSchedule::RandomPermutation,
            translation_patch: ParticleEnergyPatch::default(),
            exchange_patch: GrandCanonicalPatch::default(),
            energy_check_interval: 0,
            sweeps: 0,
        }
    }

    pub fn with_exchange_attempts(mut self, attempts: u64) -> Self {
        self.exchange_attempts_per_sweep = attempts;
        self
    }

    pub fn with_translation(mut self, translation: TranslateParticle<D>) -> Self {
        self.translation = translation;
        self
    }

    pub fn with_visit_schedule(mut self, schedule: VisitSchedule) -> Self {
        self.visit_schedule = schedule;
        self
    }

    pub fn with_energy_check_interval(mut self, interval: u64) -> Self {
        self.energy_check_interval = interval;
        self
    }

    #[inline]
    pub const fn log_activity(&self) -> f64 {
        self.log_activity
    }
}

impl<const D: usize, P: PairPotential> ParticleAlgorithm<D, P> for ParticleGrandCanonicalCore<D> {
    fn sweep_with_phase(
        &mut self,
        system: &mut ParticleSystem<D>,
        potential: &P,
        rng: &mut impl Rng,
        phase: SimulationPhase,
    ) {
        self.exchange
            .validate_potential(potential)
            .expect("invalid μVT species/potential combination");
        self.exchange
            .validate_state(system)
            .expect("invalid μVT state/proposal combination");
        let target = GrandCanonical::new(system.beta, self.log_activity);
        let acceptance = MetropolisHastingsAcceptance;
        let particles = self.order.prepare(system.len(), self.visit_schedule, rng);
        for &particle in particles {
            let proposal = self
                .translation
                .propose(system.configuration(), particle, rng);
            let outcome = metropolis_hastings_step(
                system,
                potential,
                &proposal,
                &target,
                &acceptance,
                &mut self.translation_patch,
                rng,
            );
            self.translation.record_result(outcome.accepted);
        }
        for _ in 0..self.exchange_attempts_per_sweep {
            let proposal = self.exchange.propose(system, rng);
            metropolis_hastings_step(
                system,
                potential,
                &proposal,
                &target,
                &acceptance,
                &mut self.exchange_patch,
                rng,
            );
        }
        self.translation.finish_sweep(phase.allows_adaptation());

        self.sweeps = self.sweeps.wrapping_add(1);
        if should_audit_cache(self.sweeps, self.energy_check_interval) {
            audit_particle_cache(system, potential).expect("μVT particle cache audit failed");
        }
    }

    fn name(&self) -> &'static str {
        "Particle grand-canonical Metropolis-Hastings"
    }
}
