//! Lattice geometry builders.

use super::{BondType, Neighbor, Lattice};

/// Build a 1D chain lattice.
///
/// # Arguments
/// * `n_sites` - Number of sites in the chain
/// * `pbc` - Periodic boundary condition (true = ring, false = open chain)
///
/// # Panics
/// Panics if `n_sites < 2`.
pub fn build_chain(n_sites: usize, pbc: bool) -> Lattice {
    assert!(n_sites >= 2, "Chain must have at least 2 sites");

    let n_bonds = if pbc { n_sites } else { n_sites - 1 };

    let mut sites = Vec::with_capacity(n_sites);

    for i in 0..n_sites {
        let mut neighbors = Vec::new();

        // Left neighbor (i-1)
        if i > 0 || pbc {
            let left = if i > 0 { i - 1 } else { n_sites - 1 };
            neighbors.push(Neighbor {
                target: left,
                bond_type: BondType::ChainX,
            });
        }

        // Right neighbor (i+1)
        if i < n_sites - 1 || pbc {
            let right = if i < n_sites - 1 { i + 1 } else { 0 };
            neighbors.push(Neighbor {
                target: right,
                bond_type: BondType::ChainX,
            });
        }

        sites.push(neighbors);
    }

    Lattice {
        sites,
        n_sites,
        n_bonds,
    }
}

/// Build a 2D square lattice.
///
/// # Arguments
/// * `lx` - Number of sites in x direction
/// * `ly` - Number of sites in y direction
/// * `pbc` - Periodic boundary condition
///
/// # Panics
/// Panics if `lx < 2` or `ly < 2`.
pub fn build_square(lx: usize, ly: usize, pbc: bool) -> Lattice {
    assert!(lx >= 2 && ly >= 2, "Square lattice must have at least 2x2 sites");

    let n_sites = lx * ly;

    // Count bonds: horizontal bonds + vertical bonds
    let h_bonds = if pbc { lx * ly } else { (lx - 1) * ly };
    let v_bonds = if pbc { lx * ly } else { lx * (ly - 1) };
    let n_bonds = h_bonds + v_bonds;

    let mut sites = Vec::with_capacity(n_sites);

    for y in 0..ly {
        for x in 0..lx {
            let i = y * lx + x;
            let mut neighbors = Vec::new();

            // X-direction: right neighbor (x+1)
            if x < lx - 1 || pbc {
                let x_right = if x < lx - 1 { i + 1 } else { y * lx };
                neighbors.push(Neighbor {
                    target: x_right,
                    bond_type: BondType::SquareX,
                });
            }

            // X-direction: left neighbor (x-1)
            if x > 0 || pbc {
                let x_left = if x > 0 { i - 1 } else { y * lx + lx - 1 };
                neighbors.push(Neighbor {
                    target: x_left,
                    bond_type: BondType::SquareX,
                });
            }

            // Y-direction: down neighbor (y+1)
            if y < ly - 1 || pbc {
                let y_down = if y < ly - 1 { i + lx } else { x };
                neighbors.push(Neighbor {
                    target: y_down,
                    bond_type: BondType::SquareY,
                });
            }

            // Y-direction: up neighbor (y-1)
            if y > 0 || pbc {
                let y_up = if y > 0 { i - lx } else { (ly - 1) * lx + x };
                neighbors.push(Neighbor {
                    target: y_up,
                    bond_type: BondType::SquareY,
                });
            }

            sites.push(neighbors);
        }
    }

    Lattice {
        sites,
        n_sites,
        n_bonds,
    }
}