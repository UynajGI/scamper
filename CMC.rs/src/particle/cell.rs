//! Periodic simulation cells for continuous coordinates.

use crate::particle::ParticleError;

/// Periodic simulation-cell operations used by particle backends.
pub trait SimulationCell<const D: usize>: Send + Sync {
    /// Wrap a position into the primary cell.
    fn wrap(&self, position: &mut [f64; D]);

    /// Minimum-image displacement `to - from`.
    fn displacement(&self, from: &[f64; D], to: &[f64; D]) -> [f64; D];

    /// Cell volume (area in two dimensions).
    fn volume(&self) -> f64;
}

/// Axis-aligned periodic cell with side lengths known at compile-time dimension.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OrthorhombicCell<const D: usize> {
    lengths: [f64; D],
    inverse_lengths: [f64; D],
    volume: f64,
}

impl<const D: usize> OrthorhombicCell<D> {
    /// Construct a periodic orthorhombic cell.
    pub fn new(lengths: [f64; D]) -> Result<Self, ParticleError> {
        if D == 0 {
            return Err(ParticleError::ZeroDimension);
        }
        let mut inverse_lengths = [0.0; D];
        let mut volume = 1.0;
        for (axis, &length) in lengths.iter().enumerate() {
            if !length.is_finite() || length <= 0.0 {
                return Err(ParticleError::InvalidCellLength { axis, length });
            }
            inverse_lengths[axis] = length.recip();
            if !inverse_lengths[axis].is_finite() {
                return Err(ParticleError::InvalidCellLength { axis, length });
            }
            volume *= length;
            if !volume.is_finite() {
                return Err(ParticleError::NonFiniteCellVolume);
            }
        }
        Ok(Self {
            lengths,
            inverse_lengths,
            volume,
        })
    }

    /// Side lengths of the primary cell.
    #[inline]
    pub const fn lengths(&self) -> &[f64; D] {
        &self.lengths
    }

    /// Shortest side length.
    pub fn minimum_length(&self) -> f64 {
        self.lengths.iter().copied().fold(f64::INFINITY, f64::min)
    }

    /// Squared norm of a minimum-image displacement.
    #[inline]
    pub fn distance_squared(&self, left: &[f64; D], right: &[f64; D]) -> f64 {
        self.displacement(left, right)
            .iter()
            .map(|component| component * component)
            .sum()
    }
}

impl<const D: usize> SimulationCell<D> for OrthorhombicCell<D> {
    #[inline]
    fn wrap(&self, position: &mut [f64; D]) {
        for (coordinate, &length) in position.iter_mut().zip(&self.lengths) {
            *coordinate = coordinate.rem_euclid(length);
        }
    }

    #[inline]
    fn displacement(&self, from: &[f64; D], to: &[f64; D]) -> [f64; D] {
        let mut result = [0.0; D];
        for axis in 0..D {
            let raw = to[axis] - from[axis];
            result[axis] = raw - self.lengths[axis] * (raw * self.inverse_lengths[axis]).round();
        }
        result
    }

    #[inline]
    fn volume(&self) -> f64 {
        self.volume
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wraps_large_positive_and_negative_offsets() {
        let cell = OrthorhombicCell::new([4.0, 5.0]).unwrap();
        let mut position = [9.5, -11.0];
        cell.wrap(&mut position);
        assert_eq!(position, [1.5, 4.0]);
    }

    #[test]
    fn rejects_non_finite_derived_geometry() {
        assert!(OrthorhombicCell::new([f64::from_bits(1)]).is_err());
        assert!(OrthorhombicCell::new([f64::MAX, 2.0]).is_err());
    }

    #[test]
    fn minimum_image_is_antisymmetric() {
        let cell = OrthorhombicCell::new([10.0, 8.0, 6.0]).unwrap();
        let left = [0.2, 7.7, 0.1];
        let right = [9.8, 0.3, 5.9];
        let forward = cell.displacement(&left, &right);
        let reverse = cell.displacement(&right, &left);
        for axis in 0..3 {
            assert!((forward[axis] + reverse[axis]).abs() < 1e-14);
        }
        assert!((forward[0] + 0.4).abs() < 1e-14);
        assert!((forward[1] - 0.6).abs() < 1e-14);
        assert!((forward[2] + 0.2).abs() < 1e-14);
    }
}
