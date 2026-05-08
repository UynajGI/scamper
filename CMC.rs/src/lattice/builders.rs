//! Geometry builders for constructing lattices.

use super::{BondType, Lattice, Neighbor};

/// Build N-dimensional hypercubic lattice with **bidirectional** bonds.
///
/// For each bond in the positive direction, a reciprocal bond in the negative
/// direction is also added. This means each site has `2 * n_dims` neighbors
/// (or fewer at boundaries for non-PBC).
///
/// `dimensions` specifies the size along each axis (e.g., `[8]` for 1D chain,
/// `[4, 4]` for 2D square). `bond_types` provides the `BondType` for each
/// dimension direction — must have the same length as `dimensions`.
/// `pbc` applies uniformly to all directions.
pub fn build_hypercubic(dimensions: &[usize], bond_types: &[BondType], pbc: bool) -> Lattice {
    assert_eq!(dimensions.len(), bond_types.len(), "bond_types must match dimensions.len()");
    assert!(!dimensions.is_empty(), "dimensions must not be empty");
    assert!(dimensions.iter().all(|&d| d > 0), "all dimensions must be > 0");

    let n_sites: usize = dimensions.iter().product();
    let n_dims = dimensions.len();

    // Precompute strides: stride[k] = product of dimensions[0..k]
    let mut strides = vec![1usize; n_dims];
    for k in 1..n_dims {
        strides[k] = strides[k - 1] * dimensions[k - 1];
    }

    let mut sites = vec![vec![]; n_sites];
    let mut n_bonds = 0;

    for (idx, neighbors) in sites.iter_mut().enumerate() {
        // Decode coordinates from linear index
        let mut coords = vec![0usize; n_dims];
        let mut remaining = idx;
        for k in (0..n_dims).rev() {
            coords[k] = remaining / strides[k];
            remaining %= strides[k];
        }

        // Connect to neighbors in both positive and negative directions
        for dim in 0..n_dims {
            let coord = coords[dim];

            // Positive direction
            let neighbor_coord_pos = if pbc {
                (coord + 1) % dimensions[dim]
            } else {
                coord.saturating_add(1)
            };

            if neighbor_coord_pos < dimensions[dim] {
                let offset = (neighbor_coord_pos as isize - coord as isize) * strides[dim] as isize;
                let target = (idx as isize + offset) as usize;
                neighbors.push(Neighbor {
                    target,
                    bond_type: bond_types[dim],
                });
                n_bonds += 1;
            }

            // Negative direction (reciprocal bond)
            let neighbor_coord_neg = if pbc {
                Some((coord + dimensions[dim] - 1) % dimensions[dim])
            } else {
                coord.checked_sub(1)
            };

            if let Some(neg) = neighbor_coord_neg {
                if neg < dimensions[dim] {
                    let offset = (neg as isize - coord as isize) * strides[dim] as isize;
                    let target = (idx as isize + offset) as usize;
                    // Use the same BondType for the reciprocal direction
                    neighbors.push(Neighbor {
                        target,
                        bond_type: bond_types[dim],
                    });
                    n_bonds += 1;
                }
            }
        }
    }

    Lattice {
        sites,
        n_sites,
        n_bonds,
    }
}

/// Build 1D chain lattice.
pub fn build_chain(n_sites: usize, pbc: bool) -> Lattice {
    build_hypercubic(&[n_sites], &[BondType::ChainX], pbc)
}

/// Build 2D square lattice.
pub fn build_square(lx: usize, ly: usize, pbc: bool) -> Lattice {
    build_hypercubic(
        &[lx, ly],
        &[BondType::SquareX, BondType::SquareY],
        pbc,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chain_4_open() {
        let lattice = build_chain(4, false);
        assert_eq!(lattice.n_sites, 4);
        // 3 physical bonds × 2 (bidirectional) = 6 directed bonds
        assert_eq!(lattice.n_bonds, 6);
    }

    #[test]
    fn test_chain_4_pbc() {
        let lattice = build_chain(4, true);
        assert_eq!(lattice.n_sites, 4);
        // 4 physical bonds × 2 (bidirectional) = 8 directed bonds
        assert_eq!(lattice.n_bonds, 8);
    }

    #[test]
    fn test_square_2x2_open() {
        let lattice = build_square(2, 2, false);
        assert_eq!(lattice.n_sites, 4);
        // 4 physical bonds × 2 = 8 directed bonds
        assert_eq!(lattice.n_bonds, 8);
    }

    #[test]
    fn test_square_2x2_pbc() {
        let lattice = build_square(2, 2, true);
        assert_eq!(lattice.n_sites, 4);
        // 8 physical bonds × 2 = 16 directed bonds
        assert_eq!(lattice.n_bonds, 16);
    }

    #[test]
    fn test_square_4x4_pbc() {
        let lattice = build_square(4, 4, true);
        assert_eq!(lattice.n_sites, 16);
        // 32 physical bonds × 2 = 64 directed bonds
        assert_eq!(lattice.n_bonds, 64);
    }

    #[test]
    fn test_hypercubic_1d_is_chain() {
        let lattice = build_hypercubic(&[8], &[BondType::ChainX], true);
        assert_eq!(lattice.n_sites, 8);
        // 8 physical bonds × 2 = 16
        assert_eq!(lattice.n_bonds, 16);
    }

    #[test]
    fn test_hypercubic_2d_is_square() {
        let lattice = build_hypercubic(&[3, 4], &[BondType::SquareX, BondType::SquareY], false);
        assert_eq!(lattice.n_sites, 12);
        // open 3x4: physical x-bonds = (3-1)*4=8, y-bonds = 3*(4-1)=9 → 17 physical × 2 = 34
        assert_eq!(lattice.n_bonds, 34);
    }

    #[test]
    fn test_hypercubic_3d_pbc() {
        let lattice = build_hypercubic(
            &[2, 3, 4],
            &[BondType::SquareX, BondType::SquareY, BondType::SquareZ],
            true,
        );
        assert_eq!(lattice.n_sites, 24);
        // 3D PBC: physical bonds = 24*3=72, directed = 144
        assert_eq!(lattice.n_bonds, 144);
    }

    #[test]
    fn test_hypercubic_3d_open() {
        let lattice = build_hypercubic(
            &[2, 3, 4],
            &[BondType::SquareX, BondType::SquareY, BondType::SquareZ],
            false,
        );
        assert_eq!(lattice.n_sites, 24);
        // open physical: x:(2-1)*3*4=12, y:2*(3-1)*4=16, z:2*3*(4-1)=18 → 46 × 2 = 92
        assert_eq!(lattice.n_bonds, 92);
    }
}
