//! Packed periodic cell list with O(1) accepted membership updates.

use crate::particle::{OrthorhombicCell, ParticleConfiguration, ParticleError};

/// Cell-list acceleration structure for a fixed cutoff.
#[derive(Debug, Clone)]
pub struct CellList<const D: usize> {
    cutoff_squared: f64,
    cell_counts: [usize; D],
    cell_widths: [f64; D],
    strides: [usize; D],
    buckets: Vec<Vec<usize>>,
    particle_cell: Vec<usize>,
    particle_slot: Vec<usize>,
    neighbor_cells: Vec<Vec<usize>>,
}

impl<const D: usize> CellList<D> {
    /// Build a packed cell list from an accepted configuration.
    pub fn new(
        configuration: &ParticleConfiguration<D>,
        cutoff_squared: f64,
    ) -> Result<Self, ParticleError> {
        if !cutoff_squared.is_finite() || cutoff_squared <= 0.0 {
            return Err(ParticleError::InvalidPotential(
                "cell-list cutoff must be finite and positive".to_string(),
            ));
        }
        let cutoff = cutoff_squared.sqrt();
        let cell = configuration.cell();
        let maximum = 0.5 * cell.minimum_length();
        if cutoff > maximum * (1.0 + 16.0 * f64::EPSILON) {
            return Err(ParticleError::CutoffTooLarge { cutoff, maximum });
        }

        let mut cell_counts = [0; D];
        let mut cell_widths = [0.0; D];
        let mut strides = [0; D];
        let mut total_cells = 1usize;
        for axis in 0..D {
            let count = (cell.lengths()[axis] / cutoff).floor() as usize;
            cell_counts[axis] = count.max(1);
            cell_widths[axis] = cell.lengths()[axis] / cell_counts[axis] as f64;
            strides[axis] = total_cells;
            total_cells = total_cells
                .checked_mul(cell_counts[axis])
                .ok_or_else(|| ParticleError::InvalidCellList("cell count overflow".to_string()))?;
        }

        let mut result = Self {
            cutoff_squared,
            cell_counts,
            cell_widths,
            strides,
            buckets: vec![Vec::new(); total_cells],
            particle_cell: vec![0; configuration.len()],
            particle_slot: vec![0; configuration.len()],
            neighbor_cells: Vec::new(),
        };
        result.neighbor_cells = result.build_neighbor_cells();
        result.rebuild(configuration);
        result.validate(configuration)?;
        Ok(result)
    }

    /// Squared cutoff used to build this cache.
    #[inline]
    pub const fn cutoff_squared(&self) -> f64 {
        self.cutoff_squared
    }

    /// Number of cells along each axis.
    #[inline]
    pub const fn cell_counts(&self) -> &[usize; D] {
        &self.cell_counts
    }

    /// Physical cell width along each axis.
    #[inline]
    pub const fn cell_widths(&self) -> &[f64; D] {
        &self.cell_widths
    }

    /// Packed cell containing one accepted particle.
    #[inline]
    pub fn particle_cell(&self, particle: usize) -> usize {
        self.particle_cell[particle]
    }

    /// Map a wrapped or unwrapped trial position to a periodic cell.
    pub fn cell_of_position(&self, position: &[f64; D], cell: &OrthorhombicCell<D>) -> usize {
        let mut linear = 0usize;
        for (axis, &value) in position.iter().enumerate() {
            let wrapped = value.rem_euclid(cell.lengths()[axis]);
            let coordinate = (wrapped / self.cell_widths[axis]).floor() as usize;
            linear += coordinate.min(self.cell_counts[axis] - 1) * self.strides[axis];
        }
        linear
    }

    /// Fill reusable scratch with all particles in neighboring cells.
    pub fn fill_candidates(&self, cell_index: usize, output: &mut Vec<usize>) {
        output.clear();
        for &neighbor_cell in &self.neighbor_cells[cell_index] {
            output.extend_from_slice(&self.buckets[neighbor_cell]);
        }
    }

    /// Rebuild all packed memberships from the accepted configuration.
    fn rebuild(&mut self, configuration: &ParticleConfiguration<D>) {
        for bucket in &mut self.buckets {
            bucket.clear();
        }
        self.particle_cell.resize(configuration.len(), 0);
        self.particle_slot.resize(configuration.len(), 0);
        for particle in 0..configuration.len() {
            let cell_index =
                self.cell_of_position(configuration.position(particle), configuration.cell());
            let slot = self.buckets[cell_index].len();
            self.buckets[cell_index].push(particle);
            self.particle_cell[particle] = cell_index;
            self.particle_slot[particle] = slot;
        }
    }

    /// Apply one accepted particle's membership change without rebuilding.
    pub(crate) fn move_particle(&mut self, particle: usize, new_cell: usize) {
        let old_cell = self.particle_cell[particle];
        if old_cell == new_cell {
            return;
        }

        let old_slot = self.particle_slot[particle];
        let removed = self.buckets[old_cell].swap_remove(old_slot);
        debug_assert_eq!(removed, particle);
        if old_slot < self.buckets[old_cell].len() {
            let swapped_particle = self.buckets[old_cell][old_slot];
            self.particle_slot[swapped_particle] = old_slot;
        }

        let new_slot = self.buckets[new_cell].len();
        self.buckets[new_cell].push(particle);
        self.particle_cell[particle] = new_cell;
        self.particle_slot[particle] = new_slot;
    }

    /// Audit every packed index and position-to-cell mapping.
    pub fn validate(&self, configuration: &ParticleConfiguration<D>) -> Result<(), ParticleError> {
        if self.particle_cell.len() != configuration.len()
            || self.particle_slot.len() != configuration.len()
        {
            return Err(ParticleError::InvalidCellList(
                "particle index buffers have the wrong length".to_string(),
            ));
        }
        let mut seen = vec![false; configuration.len()];
        for (cell_index, bucket) in self.buckets.iter().enumerate() {
            for (slot, &particle) in bucket.iter().enumerate() {
                if particle >= configuration.len() {
                    return Err(ParticleError::InvalidCellList(
                        "bucket contains out-of-range particle".to_string(),
                    ));
                }
                if seen[particle] {
                    return Err(ParticleError::InvalidCellList(
                        "particle appears in multiple buckets".to_string(),
                    ));
                }
                seen[particle] = true;
                if self.particle_cell[particle] != cell_index
                    || self.particle_slot[particle] != slot
                {
                    return Err(ParticleError::InvalidCellList(
                        "reverse packed index mismatch".to_string(),
                    ));
                }
                let expected =
                    self.cell_of_position(configuration.position(particle), configuration.cell());
                if expected != cell_index {
                    return Err(ParticleError::InvalidCellList(
                        "particle is stored in the wrong spatial cell".to_string(),
                    ));
                }
            }
        }
        if seen.iter().any(|present| !present) {
            return Err(ParticleError::InvalidCellList(
                "particle is missing from all buckets".to_string(),
            ));
        }
        Ok(())
    }

    fn build_neighbor_cells(&self) -> Vec<Vec<usize>> {
        let offsets = neighbor_offsets::<D>();
        let mut result = Vec::with_capacity(self.buckets.len());
        for linear in 0..self.buckets.len() {
            let coordinates = self.coordinates(linear);
            let mut neighbors = Vec::with_capacity(offsets.len());
            for offset in &offsets {
                let mut neighbor = [0; D];
                for axis in 0..D {
                    let count = self.cell_counts[axis] as isize;
                    neighbor[axis] =
                        (coordinates[axis] as isize + offset[axis]).rem_euclid(count) as usize;
                }
                let neighbor_linear = self.linear_index(&neighbor);
                if !neighbors.contains(&neighbor_linear) {
                    neighbors.push(neighbor_linear);
                }
            }
            result.push(neighbors);
        }
        result
    }

    fn coordinates(&self, mut linear: usize) -> [usize; D] {
        let mut coordinates = [0; D];
        for axis in (0..D).rev() {
            coordinates[axis] = linear / self.strides[axis];
            linear %= self.strides[axis];
        }
        coordinates
    }

    fn linear_index(&self, coordinates: &[usize; D]) -> usize {
        coordinates
            .iter()
            .zip(self.strides)
            .map(|(&coordinate, stride)| coordinate * stride)
            .sum()
    }
}

fn neighbor_offsets<const D: usize>() -> Vec<[isize; D]> {
    fn recurse<const D: usize>(
        axis: usize,
        current: &mut [isize; D],
        output: &mut Vec<[isize; D]>,
    ) {
        if axis == D {
            output.push(*current);
            return;
        }
        for offset in -1..=1 {
            current[axis] = offset;
            recurse(axis + 1, current, output);
        }
    }

    let mut result = Vec::new();
    recurse(0, &mut [0; D], &mut result);
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::particle::OrthorhombicCell;

    #[test]
    fn small_cell_counts_do_not_duplicate_candidates() {
        let cell = OrthorhombicCell::new([5.0, 5.0]).unwrap();
        let configuration =
            ParticleConfiguration::new(vec![[0.1, 0.1], [2.6, 2.6], [4.9, 4.9]], vec![0; 3], cell)
                .unwrap();
        let list = CellList::new(&configuration, 2.5f64.powi(2)).unwrap();
        let mut candidates = Vec::new();
        list.fill_candidates(0, &mut candidates);
        candidates.sort_unstable();
        candidates.dedup();
        assert_eq!(candidates.len(), 3);
    }
}
