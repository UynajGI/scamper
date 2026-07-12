//! Transactional multi-particle moves used by rigid molecules and local groups.

use crate::core::ensemble::ThermodynamicDelta;
use crate::core::trial::TrialEvaluator;
use crate::particle::{PairPotential, ParticleError, ParticleSystem, SimulationCell};

/// Atomically replace the positions of a distinct set of particles.
#[derive(Debug, Clone, PartialEq)]
pub struct ParticleBatchMove<const D: usize> {
    particles: Vec<usize>,
    positions: Vec<[f64; D]>,
}

impl<const D: usize> ParticleBatchMove<D> {
    /// Build a batch move. Particle range is checked when evaluated against a state.
    pub fn new(particles: Vec<usize>, positions: Vec<[f64; D]>) -> Result<Self, ParticleError> {
        if particles.is_empty() {
            return Err(ParticleError::InvalidMove(
                "batch move must contain at least one particle".to_string(),
            ));
        }
        if particles.len() != positions.len() {
            return Err(ParticleError::InvalidMove(
                "batch particle and position buffers differ in length".to_string(),
            ));
        }
        for (index, position) in positions.iter().enumerate() {
            if position.iter().any(|coordinate| !coordinate.is_finite()) {
                return Err(ParticleError::InvalidMove(format!(
                    "batch trial position {index} contains a non-finite coordinate"
                )));
            }
        }
        let mut sorted = particles.clone();
        sorted.sort_unstable();
        if sorted.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err(ParticleError::InvalidMove(
                "batch move contains duplicate particle indices".to_string(),
            ));
        }
        Ok(Self {
            particles,
            positions,
        })
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.particles.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.particles.is_empty()
    }

    #[inline]
    pub fn particles(&self) -> &[usize] {
        &self.particles
    }

    #[inline]
    pub fn positions(&self) -> &[[f64; D]] {
        &self.positions
    }
}

/// Reusable local-energy and spatial-cache workspace for a batch move.
#[derive(Debug, Clone, Default)]
pub struct ParticleBatchPatch {
    /// Physical energy change, not multiplied by β.
    pub delta_energy: f64,
    pub(crate) new_cells: Vec<usize>,
    candidate_scratch: Vec<usize>,
    unique_candidates: Vec<usize>,
    candidate_marks: Vec<u32>,
    candidate_generation: u32,
    move_marks: Vec<u32>,
    move_slots: Vec<usize>,
    move_generation: u32,
}

impl ParticleBatchPatch {
    fn prepare<const D: usize>(&mut self, movement: &ParticleBatchMove<D>, particles: usize) {
        self.new_cells.clear();
        self.new_cells.reserve(movement.len());
        self.candidate_marks.resize(particles, 0);
        self.move_marks.resize(particles, 0);
        self.move_slots.resize(particles, 0);
        self.move_generation = next_generation(self.move_generation, &mut self.move_marks);
        for (slot, &particle) in movement.particles().iter().enumerate() {
            assert!(particle < particles, "batch trial particle out of range");
            assert_ne!(
                self.move_marks[particle], self.move_generation,
                "batch trial contains a duplicate particle"
            );
            self.move_marks[particle] = self.move_generation;
            self.move_slots[particle] = slot;
        }
    }

    #[inline]
    fn moved_slot(&self, particle: usize) -> Option<usize> {
        (self.move_marks[particle] == self.move_generation).then_some(self.move_slots[particle])
    }

    fn begin_candidates(&mut self) {
        self.candidate_generation =
            next_generation(self.candidate_generation, &mut self.candidate_marks);
        self.unique_candidates.clear();
    }

    fn mark_candidate(&mut self, particle: usize) {
        if self.candidate_marks[particle] != self.candidate_generation {
            self.candidate_marks[particle] = self.candidate_generation;
            self.unique_candidates.push(particle);
        }
    }
}

fn next_generation(current: u32, marks: &mut [u32]) -> u32 {
    if current == u32::MAX {
        marks.fill(0);
        1
    } else {
        current + 1
    }
}

impl<const D: usize> ParticleSystem<D> {
    pub(crate) fn evaluate_batch<P: PairPotential>(
        &self,
        potential: &P,
        movement: &ParticleBatchMove<D>,
        patch: &mut ParticleBatchPatch,
    ) -> ThermodynamicDelta {
        patch.prepare(movement, self.len());
        assert_eq!(
            potential.cutoff_squared().to_bits(),
            self.cell_list.cutoff_squared().to_bits(),
            "batch trial potential cutoff differs from the state cell list"
        );

        for position in movement.positions() {
            patch.new_cells.push(
                self.cell_list
                    .cell_of_position(position, self.configuration.cell()),
            );
        }

        let mut delta_energy = 0.0;

        // Changed-changed interactions are evaluated directly because both new
        // positions are absent from the accepted cell list.
        for left_slot in 0..movement.len() {
            let left = movement.particles()[left_slot];
            let species_left = self.configuration.species_of(left);
            for right_slot in left_slot + 1..movement.len() {
                let right = movement.particles()[right_slot];
                let species_right = self.configuration.species_of(right);
                let old_distance = self.configuration.cell().distance_squared(
                    self.configuration.position(left),
                    self.configuration.position(right),
                );
                let new_distance = self.configuration.cell().distance_squared(
                    &movement.positions()[left_slot],
                    &movement.positions()[right_slot],
                );
                let old_energy = potential.energy(species_left, species_right, old_distance);
                let new_energy = potential.energy(species_left, species_right, new_distance);
                if !new_energy.is_finite() {
                    patch.delta_energy = f64::INFINITY;
                    return ThermodynamicDelta::energy(f64::INFINITY);
                }
                delta_energy += new_energy - old_energy;
            }
        }

        // Changed-unchanged interactions use the union of old and new neighbor
        // cells. Each pair is owned by its changed endpoint and counted once.
        for slot in 0..movement.len() {
            let particle = movement.particles()[slot];
            patch.begin_candidates();
            self.cell_list.fill_candidates(
                self.cell_list.particle_cell(particle),
                &mut patch.candidate_scratch,
            );
            for index in 0..patch.candidate_scratch.len() {
                let candidate = patch.candidate_scratch[index];
                patch.mark_candidate(candidate);
            }
            self.cell_list
                .fill_candidates(patch.new_cells[slot], &mut patch.candidate_scratch);
            for index in 0..patch.candidate_scratch.len() {
                let candidate = patch.candidate_scratch[index];
                patch.mark_candidate(candidate);
            }

            let species = self.configuration.species_of(particle);
            for candidate_index in 0..patch.unique_candidates.len() {
                let other = patch.unique_candidates[candidate_index];
                if other == particle || patch.moved_slot(other).is_some() {
                    continue;
                }
                let other_species = self.configuration.species_of(other);
                let old_distance = self.configuration.cell().distance_squared(
                    self.configuration.position(particle),
                    self.configuration.position(other),
                );
                let new_distance = self.configuration.cell().distance_squared(
                    &movement.positions()[slot],
                    self.configuration.position(other),
                );
                let old_energy = potential.energy(species, other_species, old_distance);
                let new_energy = potential.energy(species, other_species, new_distance);
                if !new_energy.is_finite() {
                    patch.delta_energy = f64::INFINITY;
                    return ThermodynamicDelta::energy(f64::INFINITY);
                }
                delta_energy += new_energy - old_energy;
            }
        }

        assert!(!delta_energy.is_nan(), "batch trial energy change is NaN");
        patch.delta_energy = delta_energy;
        ThermodynamicDelta::energy(delta_energy)
    }

    pub(crate) fn commit_batch(
        &mut self,
        movement: &ParticleBatchMove<D>,
        patch: &ParticleBatchPatch,
    ) {
        assert!(
            patch.delta_energy.is_finite(),
            "an infinite-energy batch trial must never be committed"
        );
        assert_eq!(movement.len(), patch.new_cells.len());
        for (slot, &particle) in movement.particles().iter().enumerate() {
            let mut position = movement.positions()[slot];
            self.configuration.cell().wrap(&mut position);
            self.configuration.set_position(particle, position);
        }
        for (slot, &particle) in movement.particles().iter().enumerate() {
            self.cell_list
                .move_particle(particle, patch.new_cells[slot]);
        }
        self.energy += patch.delta_energy;
        assert!(
            self.energy.is_finite(),
            "batch commit produced non-finite energy"
        );
    }
}

impl<const D: usize, P: PairPotential> TrialEvaluator<P, ParticleBatchMove<D>>
    for ParticleSystem<D>
{
    type Delta = ThermodynamicDelta;
    type Patch = ParticleBatchPatch;

    fn evaluate_trial(
        &self,
        model: &P,
        movement: &ParticleBatchMove<D>,
        patch: &mut Self::Patch,
    ) -> Self::Delta {
        self.evaluate_batch(model, movement, patch)
    }

    fn commit_trial(&mut self, movement: &ParticleBatchMove<D>, patch: &Self::Patch) {
        self.commit_batch(movement, patch);
    }
}
