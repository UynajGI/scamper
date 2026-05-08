//! Bond type for direction-dependent Hamiltonian parameters.

/// Bond type for direction-dependent Hamiltonian parameters.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum BondType {
    /// 1D chain
    ChainX,
    /// 2D square lattice
    SquareX,
    SquareY,
    SquareZ,
    /// 2D triangular lattice
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bond_type_equality() {
        let a = BondType::ChainX;
        let b = BondType::ChainX;
        assert_eq!(a, b);
    }

    #[test]
    fn test_bond_type_custom() {
        let a = BondType::Custom(42);
        let b = BondType::Custom(42);
        let c = BondType::Custom(99);
        assert_eq!(a, b);
        assert_ne!(a, c);
    }
}
