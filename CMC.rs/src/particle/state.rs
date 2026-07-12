//! Transactional particle state and local pair-energy cache patches.

use crate::core::ensemble::{CanonicalEnsemble, ThermodynamicDelta};
use crate::particle::{
    CellList, PairPotential, ParticleConfiguration, ParticleError, ParticleTranslation,
    SimulationCell,
};

/// Reusable patch for one particle translation.
#[derive(Debug, Clone, Default)]
pub struct ParticleEnergyPatch {
    /// Physical energy change, not multiplied by β.
    pub delta_energy: f64,
    pub(crate) old_cell: usize,
    pub(crate) new_cell: usize,
    pub(crate) candidates: Vec<usize>,
}

/// Accepted NVT particle state with packed cell-list acceleration.
#[derive(Debug, Clone)]
pub struct ParticleSystem<const D: usize> {
    configuration: ParticleConfiguration<D>,
    cell_list: CellList<D>,
    /// Cached physical potential energy, never multiplied by β.
    pub energy: f64,
    /// Canonical inverse temperature.
    pub beta: f64,
}

impl<const D: usize> ParticleSystem<D> {
    /// Construct an accepted state, validating species, cutoff and finite energy.
    pub fn new<P: PairPotential>(
        configuration: ParticleConfiguration<D>,
        potential: &P,
        beta: f64,
    ) -> Result<Self, ParticleError> {
        if !beta.is_finite() || beta < 0.0 {
            return Err(ParticleError::InvalidPotential(
                "beta must be finite and non-negative".to_string(),
            ));
        }
        for (particle, &species) in configuration.species().iter().enumerate() {
            if !potential.supports_species(species) {
                return Err(ParticleError::UnsupportedSpecies { particle, species });
            }
        }
        let cell_list = CellList::new(&configuration, potential.cutoff_squared())?;
        let energy = compute_total_energy(&configuration, potential);
        if !energy.is_finite() {
            return Err(ParticleError::NonFiniteAcceptedEnergy);
        }
        Ok(Self {
            configuration,
            cell_list,
            energy,
            beta,
        })
    }

    /// Accepted particle configuration.
    #[inline]
    pub const fn configuration(&self) -> &ParticleConfiguration<D> {
        &self.configuration
    }

    /// Packed cell-list cache.
    #[inline]
    pub const fn cell_list(&self) -> &CellList<D> {
        &self.cell_list
    }

    /// Number of particles.
    #[inline]
    pub fn len(&self) -> usize {
        self.configuration.len()
    }

    /// Whether the state contains no particles.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.configuration.is_empty()
    }

    /// NVT target distribution for the current β.
    #[inline]
    pub fn canonical_ensemble(&self) -> CanonicalEnsemble {
        CanonicalEnsemble::new(self.beta)
    }

    /// Update β without changing physical energy or spatial caches.
    pub fn set_beta(&mut self, beta: f64) -> Result<(), ParticleError> {
        if !beta.is_finite() || beta < 0.0 {
            return Err(ParticleError::InvalidPotential(
                "beta must be finite and non-negative".to_string(),
            ));
        }
        self.beta = beta;
        Ok(())
    }

    /// Exact O(N²) energy recomputation used for initialization and audits.
    pub fn recompute_energy<P: PairPotential>(&mut self, potential: &P) -> f64 {
        self.energy = compute_total_energy(&self.configuration, potential);
        assert!(
            self.energy.is_finite(),
            "accepted particle configuration has non-finite energy"
        );
        self.energy
    }

    /// Cached-minus-exact physical energy.
    pub fn energy_error<P: PairPotential>(&self, potential: &P) -> f64 {
        self.energy - compute_total_energy(&self.configuration, potential)
    }

    /// Audit configuration, potential compatibility, energy and cell-list invariants.
    pub fn validate<P: PairPotential>(&self, potential: &P) -> Result<(), ParticleError> {
        if potential.cutoff_squared().to_bits() != self.cell_list.cutoff_squared().to_bits() {
            return Err(ParticleError::InvalidCellList(
                "potential cutoff differs from cached cutoff".to_string(),
            ));
        }
        if !self.beta.is_finite() || self.beta < 0.0 {
            return Err(ParticleError::InvalidPotential(
                "beta must be finite and non-negative".to_string(),
            ));
        }
        if !self.energy.is_finite() {
            return Err(ParticleError::NonFiniteAcceptedEnergy);
        }
        for (particle, &species) in self.configuration.species().iter().enumerate() {
            if !potential.supports_species(species) {
                return Err(ParticleError::UnsupportedSpecies { particle, species });
            }
        }
        self.cell_list.validate(&self.configuration)?;
        let exact = compute_total_energy(&self.configuration, potential);
        if !exact.is_finite() {
            return Err(ParticleError::NonFiniteAcceptedEnergy);
        }
        let tolerance = 1e-10 * (1.0 + exact.abs());
        if (self.energy - exact).abs() > tolerance {
            return Err(ParticleError::EnergyCacheMismatch {
                cached: self.energy,
                exact,
            });
        }
        Ok(())
    }

    pub(crate) fn evaluate_translation<P: PairPotential>(
        &self,
        potential: &P,
        movement: &ParticleTranslation<D>,
        patch: &mut ParticleEnergyPatch,
    ) -> ThermodynamicDelta {
        assert!(
            movement.particle < self.len(),
            "trial particle out of range"
        );
        assert!(
            movement
                .position
                .iter()
                .all(|coordinate| coordinate.is_finite()),
            "trial position contains a non-finite coordinate"
        );
        assert_eq!(
            potential.cutoff_squared().to_bits(),
            self.cell_list.cutoff_squared().to_bits(),
            "trial potential cutoff differs from the state cell list"
        );

        let particle = movement.particle;
        let species_i = self.configuration.species_of(particle);
        patch.old_cell = self.cell_list.particle_cell(particle);
        patch.new_cell = self
            .cell_list
            .cell_of_position(&movement.position, self.configuration.cell());

        self.cell_list
            .fill_candidates(patch.old_cell, &mut patch.candidates);
        let mut old_energy = 0.0;
        for &other in &patch.candidates {
            if other == particle {
                continue;
            }
            let distance_squared = self.configuration.cell().distance_squared(
                self.configuration.position(particle),
                self.configuration.position(other),
            );
            old_energy += potential.energy(
                species_i,
                self.configuration.species_of(other),
                distance_squared,
            );
        }
        assert!(
            old_energy.is_finite(),
            "accepted local energy is non-finite"
        );

        self.cell_list
            .fill_candidates(patch.new_cell, &mut patch.candidates);
        let mut new_energy = 0.0;
        for &other in &patch.candidates {
            if other == particle {
                continue;
            }
            let distance_squared = self
                .configuration
                .cell()
                .distance_squared(&movement.position, self.configuration.position(other));
            new_energy += potential.energy(
                species_i,
                self.configuration.species_of(other),
                distance_squared,
            );
        }
        assert!(!new_energy.is_nan(), "trial local energy is NaN");

        patch.delta_energy = if new_energy.is_infinite() {
            f64::INFINITY
        } else {
            new_energy - old_energy
        };
        ThermodynamicDelta::energy(patch.delta_energy)
    }

    pub(crate) fn commit_translation(
        &mut self,
        movement: &ParticleTranslation<D>,
        patch: &ParticleEnergyPatch,
    ) {
        assert!(
            patch.delta_energy.is_finite(),
            "an infinite-energy particle trial must never be committed"
        );
        let mut wrapped_position = movement.position;
        self.configuration.cell().wrap(&mut wrapped_position);
        self.configuration
            .set_position(movement.particle, wrapped_position);
        self.cell_list
            .move_particle(movement.particle, patch.new_cell);
        self.energy += patch.delta_energy;
        assert!(
            self.energy.is_finite(),
            "accepted particle trial produced non-finite energy"
        );
    }
}

/// Exact O(N²) physical pair energy of a configuration.
pub fn compute_total_energy<const D: usize, P: PairPotential>(
    configuration: &ParticleConfiguration<D>,
    potential: &P,
) -> f64 {
    let mut energy = 0.0;
    for left in 0..configuration.len() {
        for right in left + 1..configuration.len() {
            let distance_squared = configuration
                .cell()
                .distance_squared(configuration.position(left), configuration.position(right));
            energy += potential.energy(
                configuration.species_of(left),
                configuration.species_of(right),
                distance_squared,
            );
        }
    }
    energy
}
