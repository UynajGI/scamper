//! Bond type enum for direction-dependent Hamiltonian parameters.

use std::hash::Hash;

/// Bond type enum for direction-dependent Hamiltonian parameters.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum BondType {
    /// 1D chain: horizontal bond
    ChainX,

    /// 2D square lattice
    SquareX,
    SquareY,

    /// 2D triangular lattice (0°, 60°, 120°)
    TriX,
    TriY,
    TriZ,

    /// 2D honeycomb lattice
    HoneyX,
    HoneyY,
    HoneyZ,

    /// Custom bond type for arbitrary networks
    Custom(u8),
}