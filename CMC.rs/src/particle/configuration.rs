//! Accepted particle coordinates, species and periodic cell.

use crate::particle::{OrthorhombicCell, ParticleError, SimulationCell};

/// Continuous particle configuration using array-of-structures coordinates.
#[derive(Debug, Clone, PartialEq)]
pub struct ParticleConfiguration<const D: usize> {
    positions: Vec<[f64; D]>,
    species: Vec<u16>,
    cell: OrthorhombicCell<D>,
}

impl<const D: usize> ParticleConfiguration<D> {
    /// Build a configuration and wrap all coordinates into the primary cell.
    pub fn new(
        mut positions: Vec<[f64; D]>,
        species: Vec<u16>,
        cell: OrthorhombicCell<D>,
    ) -> Result<Self, ParticleError> {
        if positions.len() != species.len() {
            return Err(ParticleError::BufferLengthMismatch {
                positions: positions.len(),
                species: species.len(),
            });
        }
        for (particle, position) in positions.iter_mut().enumerate() {
            for (axis, coordinate) in position.iter().enumerate() {
                if !coordinate.is_finite() {
                    return Err(ParticleError::NonFinitePosition { particle, axis });
                }
            }
            cell.wrap(position);
        }
        Ok(Self {
            positions,
            species,
            cell,
        })
    }

    /// Number of particles.
    #[inline]
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    /// Whether the configuration contains no particles.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }

    /// Accepted positions.
    #[inline]
    pub fn positions(&self) -> &[[f64; D]] {
        &self.positions
    }

    /// Accepted position of one particle.
    #[inline]
    pub fn position(&self, particle: usize) -> &[f64; D] {
        &self.positions[particle]
    }

    /// Species labels parallel to [`Self::positions`].
    #[inline]
    pub fn species(&self) -> &[u16] {
        &self.species
    }

    /// Species label of one particle.
    #[inline]
    pub fn species_of(&self, particle: usize) -> u16 {
        self.species[particle]
    }

    /// Periodic simulation cell.
    #[inline]
    pub const fn cell(&self) -> &OrthorhombicCell<D> {
        &self.cell
    }

    pub(crate) fn set_position(&mut self, particle: usize, position: [f64; D]) {
        self.positions[particle] = position;
    }
}
