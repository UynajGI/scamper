//! Tests for lattice module.

use qmc_rs::lattice::{BondType, Lattice, Neighbor};
use qmc_rs::lattice::builders::{build_chain, build_square};

#[test]
fn test_bond_type_variants() {
    // 1D chain
    let bx = BondType::ChainX;
    assert_eq!(bx, BondType::ChainX);

    // 2D square
    let bx = BondType::SquareX;
    let by = BondType::SquareY;
    assert_ne!(bx, by);

    // 2D triangular
    assert_eq!(BondType::TriX, BondType::TriX);

    // Custom
    let custom = BondType::Custom(42);
    assert_eq!(custom, BondType::Custom(42));
}

#[test]
fn test_bond_type_hashable() {
    use std::collections::HashMap;
    let mut weights: HashMap<BondType, f64> = HashMap::new();
    weights.insert(BondType::SquareX, 1.0);
    assert_eq!(weights.get(&BondType::SquareX), Some(&1.0));
}

#[test]
fn test_neighbor_struct() {
    let n = Neighbor {
        target: 5,
        bond_type: BondType::SquareX,
    };
    assert_eq!(n.target, 5);
    assert_eq!(n.bond_type, BondType::SquareX);
}

#[test]
fn test_lattice_basic() {
    let lattice = Lattice {
        sites: vec![
            vec![Neighbor { target: 1, bond_type: BondType::ChainX }],
            vec![Neighbor { target: 0, bond_type: BondType::ChainX }],
        ],
        n_sites: 2,
        n_bonds: 2,
    };
    assert_eq!(lattice.n_sites, 2);
    assert_eq!(lattice.n_bonds, 2);
    assert_eq!(lattice.sites[0].len(), 1);
}

#[test]
fn test_lattice_clone() {
    let lattice = Lattice {
        sites: vec![vec![]],
        n_sites: 1,
        n_bonds: 0,
    };
    let cloned = lattice.clone();
    assert_eq!(cloned.n_sites, 1);
}

#[test]
fn test_build_chain_open() {
    let lattice = build_chain(4, false); // 4 sites, open boundary

    assert_eq!(lattice.n_sites, 4);
    assert_eq!(lattice.n_bonds, 3); // N-1 bonds for open chain

    // Site 0 has 1 neighbor (site 1)
    assert_eq!(lattice.sites[0].len(), 1);
    assert_eq!(lattice.sites[0][0].target, 1);

    // Site 1 has 2 neighbors (sites 0 and 2)
    assert_eq!(lattice.sites[1].len(), 2);

    // Site 3 has 1 neighbor (site 2)
    assert_eq!(lattice.sites[3].len(), 1);
}

#[test]
fn test_build_chain_periodic() {
    let lattice = build_chain(4, true); // 4 sites, periodic boundary

    assert_eq!(lattice.n_sites, 4);
    assert_eq!(lattice.n_bonds, 4); // N bonds for periodic chain

    // Every site has 2 neighbors
    for i in 0..4 {
        assert_eq!(lattice.sites[i].len(), 2);
    }
}

#[test]
fn test_build_chain_bond_types() {
    let lattice = build_chain(10, true);

    // All bonds should be ChainX
    for site in &lattice.sites {
        for neighbor in site {
            assert_eq!(neighbor.bond_type, BondType::ChainX);
        }
    }
}

#[test]
fn test_build_square_basic() {
    let lattice = build_square(4, 4, true); // 4x4 square, periodic

    assert_eq!(lattice.n_sites, 16);
    // 16 sites × 4 bonds per site / 2 (each bond counted twice) = 32
    assert_eq!(lattice.n_bonds, 32);
}

#[test]
fn test_build_square_open() {
    let lattice = build_square(3, 3, false); // 3x3, open

    assert_eq!(lattice.n_sites, 9);
    // Open: horizontal bonds = (3-1)*3 = 6, vertical bonds = 3*(3-1) = 6
    assert_eq!(lattice.n_bonds, 12);
}

#[test]
fn test_build_square_bond_types() {
    let lattice = build_square(2, 2, true);

    // Check that X and Y bond types exist
    let mut has_x = false;
    let mut has_y = false;

    for site in &lattice.sites {
        for neighbor in site {
            if neighbor.bond_type == BondType::SquareX {
                has_x = true;
            }
            if neighbor.bond_type == BondType::SquareY {
                has_y = true;
            }
        }
    }

    assert!(has_x, "Should have SquareX bonds");
    assert!(has_y, "Should have SquareY bonds");
}

#[test]
fn test_build_square_neighbors() {
    let lattice = build_square(4, 4, true);

    // Every site in periodic should have 4 neighbors
    for i in 0..16 {
        assert_eq!(lattice.sites[i].len(), 4);
    }
}