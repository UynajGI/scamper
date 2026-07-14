//! Local spin operators, basis maps, and positive vertex kinds.

use crate::impurity::core::local_hilbert::Spin;
use crate::impurity::ImpurityError;

/// Cartesian spin axis in the physical or sampled basis.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PhysicalAxis {
    X,
    Y,
    Z,
}

impl PhysicalAxis {
    /// Short uppercase label used in observable names.
    pub const fn label(self) -> &'static str {
        match self {
            Self::X => "X",
            Self::Y => "Y",
            Self::Z => "Z",
        }
    }

    const fn index(self) -> usize {
        match self {
            Self::X => 0,
            Self::Y => 1,
            Self::Z => 2,
        }
    }
}

/// A signed Cartesian axis.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SignedAxis {
    pub axis: PhysicalAxis,
    pub sign: i8,
}

impl SignedAxis {
    pub const fn positive(axis: PhysicalAxis) -> Self {
        Self { axis, sign: 1 }
    }

    pub const fn negative(axis: PhysicalAxis) -> Self {
        Self { axis, sign: -1 }
    }
}

/// Signed permutation taking sampled spin components to physical components.
///
/// `sampled_to_physical[s]` identifies the physical axis represented by the
/// sampled axis `s`.  Signs matter for one-point functions and cancel in
/// same-axis two-point functions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BasisTransform {
    sampled_to_physical: [SignedAxis; 3],
}

impl BasisTransform {
    /// Sampled and physical axes coincide.
    pub const fn identity() -> Self {
        Self {
            sampled_to_physical: [
                SignedAxis::positive(PhysicalAxis::X),
                SignedAxis::positive(PhysicalAxis::Y),
                SignedAxis::positive(PhysicalAxis::Z),
            ],
        }
    }

    /// Rotation used by the Rabi/spin-boson catalog:
    /// sampled `Z = physical X`, sampled `X = -physical Z`, and `Y` is fixed.
    pub const fn rotated_rabi() -> Self {
        Self {
            sampled_to_physical: [
                SignedAxis::negative(PhysicalAxis::Z),
                SignedAxis::positive(PhysicalAxis::Y),
                SignedAxis::positive(PhysicalAxis::X),
            ],
        }
    }

    /// Global pi/2 rotation around `Z`, used to make a negative XYZ pair-flip
    /// coefficient positive.  Same-axis `X` and `Y` correlations are swapped.
    pub const fn swap_xy_gauge() -> Self {
        Self {
            sampled_to_physical: [
                SignedAxis::positive(PhysicalAxis::Y),
                SignedAxis::negative(PhysicalAxis::X),
                SignedAxis::positive(PhysicalAxis::Z),
            ],
        }
    }

    /// Physical signed axis represented by a sampled axis.
    pub const fn physical_for_sampled(self, sampled: PhysicalAxis) -> SignedAxis {
        self.sampled_to_physical[sampled.index()]
    }

    /// Sampled signed axis corresponding to a requested physical axis.
    pub fn sampled_for_physical(self, physical: PhysicalAxis) -> SignedAxis {
        for (sampled_index, mapped) in self.sampled_to_physical.into_iter().enumerate() {
            if mapped.axis == physical {
                let sampled = match sampled_index {
                    0 => PhysicalAxis::X,
                    1 => PhysicalAxis::Y,
                    _ => PhysicalAxis::Z,
                };
                return SignedAxis {
                    axis: sampled,
                    sign: mapped.sign,
                };
            }
        }
        unreachable!("basis transform is a permutation")
    }

    /// Map sampled one-point components into physical component order.
    pub fn map_one_point(self, sampled: [f64; 3]) -> [f64; 3] {
        let mut physical = [0.0; 3];
        for sampled_axis in [PhysicalAxis::X, PhysicalAxis::Y, PhysicalAxis::Z] {
            let mapped = self.physical_for_sampled(sampled_axis);
            physical[mapped.axis.index()] = f64::from(mapped.sign) * sampled[sampled_axis.index()];
        }
        physical
    }

    /// Map same-axis sampled correlations into physical component order.
    pub fn map_same_axis_correlations(self, sampled: [f64; 3]) -> [f64; 3] {
        let mut physical = [0.0; 3];
        for sampled_axis in [PhysicalAxis::X, PhysicalAxis::Y, PhysicalAxis::Z] {
            let mapped = self.physical_for_sampled(sampled_axis);
            physical[mapped.axis.index()] = sampled[sampled_axis.index()];
        }
        physical
    }
}

impl Default for BasisTransform {
    fn default() -> Self {
        Self::identity()
    }
}

/// Local leg numbering of one retarded vertex.
pub const A_IN: usize = 0;
pub const A_OUT: usize = 1;
pub const B_IN: usize = 2;
pub const B_OUT: usize = 3;
pub const LEGS_PER_VERTEX: usize = 4;

/// Immutable positive local vertex type supplied by a spin-boson model.
#[derive(Debug, Clone, PartialEq)]
pub struct VertexKind {
    name: String,
    legs: [Spin; LEGS_PER_VERTEX],
    weight: f64,
    diagonal: bool,
}

impl VertexKind {
    pub fn new(
        name: impl Into<String>,
        legs: [Spin; LEGS_PER_VERTEX],
        weight: f64,
        diagonal: bool,
    ) -> Result<Self, ImpurityError> {
        if legs.iter().any(|spin| !matches!(spin, -1 | 1)) {
            return Err(ImpurityError::parameter(
                "vertex legs",
                "spin-1/2 legs must be encoded as -1 or +1",
            ));
        }
        if !weight.is_finite() || weight <= 0.0 {
            return Err(ImpurityError::parameter(
                "vertex weight",
                format!("must be finite and positive, got {weight}"),
            ));
        }
        let inferred_diagonal = legs[A_IN] == legs[A_OUT] && legs[B_IN] == legs[B_OUT];
        if inferred_diagonal != diagonal {
            return Err(ImpurityError::parameter(
                "diagonal",
                "diagonal flag does not match the leg pattern",
            ));
        }
        Ok(Self {
            name: name.into(),
            legs,
            weight,
            diagonal,
        })
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn legs(&self) -> &[Spin; LEGS_PER_VERTEX] {
        &self.legs
    }

    pub fn weight(&self) -> f64 {
        self.weight
    }

    pub fn is_diagonal(&self) -> bool {
        self.diagonal
    }

    pub fn spin(&self, leg: usize) -> Spin {
        self.legs[leg]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rotated_rabi_maps_sampled_z_to_physical_x() {
        let map = BasisTransform::rotated_rabi();
        assert_eq!(
            map.physical_for_sampled(PhysicalAxis::Z),
            SignedAxis::positive(PhysicalAxis::X)
        );
        assert_eq!(
            map.physical_for_sampled(PhysicalAxis::X),
            SignedAxis::negative(PhysicalAxis::Z)
        );
    }

    #[test]
    fn xyz_gauge_swaps_same_axis_correlations() {
        let mapped = BasisTransform::swap_xy_gauge().map_same_axis_correlations([1.0, 2.0, 3.0]);
        assert_eq!(mapped, [2.0, 1.0, 3.0]);
    }
}
