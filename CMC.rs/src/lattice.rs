//! Lattice topology — adjacency-list representation.

/// Bond direction labels for hypercubic lattices.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BondType {
    ChainX,
    SquareX,
    SquareY,
    SquareZ,
    CubicX,
    CubicY,
    CubicZ,
}

/// Neighbor entry: target site index + bond type.
#[derive(Debug, Clone)]
pub struct Neighbor {
    pub target: usize,
    pub bond_type: BondType,
}

/// Lattice as adjacency list. `sites[i]` lists neighbors of site `i`.
#[derive(Debug, Clone)]
pub struct Lattice {
    pub sites: Vec<Vec<Neighbor>>,
    pub n_sites: usize,
    pub n_bonds: usize,
}

// ── Builders ────────────────────────────────────────────────

/// 1D chain with optional PBC.
pub fn build_chain(n: usize, pbc: bool) -> Lattice {
    let mut sites: Vec<Vec<Neighbor>> = vec![Vec::new(); n];
    let mut n_bonds = 0usize;

    for i in 0..n {
        // right neighbor
        if i + 1 < n {
            sites[i].push(Neighbor {
                target: i + 1,
                bond_type: BondType::ChainX,
            });
            sites[i + 1].push(Neighbor {
                target: i,
                bond_type: BondType::ChainX,
            });
            n_bonds += 2;
        } else if pbc && n > 1 {
            sites[i].push(Neighbor {
                target: 0,
                bond_type: BondType::ChainX,
            });
            sites[0].push(Neighbor {
                target: i,
                bond_type: BondType::ChainX,
            });
            n_bonds += 2;
        }
    }

    Lattice {
        sites,
        n_sites: n,
        n_bonds,
    }
}

/// 2D square lattice (width × height) with optional PBC.
pub fn build_square(w: usize, h: usize, pbc: bool) -> Lattice {
    build_hypercubic(&[w, h], &[BondType::SquareX, BondType::SquareY], pbc)
}

/// N-dimensional hypercubic lattice.
///
/// `dims` gives size along each axis. `bond_types[i]` labels bonds in direction `i`.
/// Must have `dims.len() == bond_types.len()`.
pub fn build_hypercubic(dims: &[usize], bond_types: &[BondType], pbc: bool) -> Lattice {
    assert_eq!(dims.len(), bond_types.len());
    let n_dims = dims.len();
    let n_sites: usize = dims.iter().product();
    let mut sites: Vec<Vec<Neighbor>> = vec![Vec::new(); n_sites];

    // Precompute stride for each dimension
    let mut strides = vec![1usize; n_dims];
    for k in 1..n_dims {
        strides[k] = strides[k - 1] * dims[k - 1];
    }

    let mut n_bonds = 0usize;
    for idx in 0..n_sites {
        // Decode linear index → coordinates
        let mut remaining = idx;
        let mut coords = vec![0usize; n_dims];
        for k in (0..n_dims).rev() {
            coords[k] = remaining / strides[k];
            remaining %= strides[k];
        }

        for dim in 0..n_dims {
            let c = coords[dim];
            // positive direction
            let neighbor_coord = if pbc {
                (c + 1) % dims[dim]
            } else if c + 1 < dims[dim] {
                c + 1
            } else {
                continue; // boundary, skip
            };

            let mut neighbor_coords = coords.clone();
            neighbor_coords[dim] = neighbor_coord;
            let neighbor_idx: usize = neighbor_coords
                .iter()
                .enumerate()
                .map(|(k, &v)| v * strides[k])
                .sum();

            sites[idx].push(Neighbor {
                target: neighbor_idx,
                bond_type: bond_types[dim],
            });
            sites[neighbor_idx].push(Neighbor {
                target: idx,
                bond_type: bond_types[dim],
            });
            n_bonds += 2;
        }
    }

    Lattice {
        sites,
        n_sites,
        n_bonds,
    }
}

// ── Tests ───────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chain_open() {
        let l = build_chain(4, false);
        assert_eq!(l.n_sites, 4);
        // 3 physical bonds × 2 directed = 6
        assert_eq!(l.n_bonds, 6);
    }

    #[test]
    fn test_chain_pbc() {
        let l = build_chain(4, true);
        assert_eq!(l.n_sites, 4);
        assert_eq!(l.n_bonds, 8);
    }

    #[test]
    fn test_chain_single_site() {
        let l = build_chain(1, true);
        assert_eq!(l.n_sites, 1);
        assert_eq!(l.n_bonds, 0);
        assert!(l.sites[0].is_empty());
    }

    #[test]
    fn test_square_2x2_open() {
        let l = build_square(2, 2, false);
        assert_eq!(l.n_sites, 4);
        assert_eq!(l.n_bonds, 8); // 4 physical × 2
    }

    #[test]
    fn test_square_2x2_pbc() {
        let l = build_square(2, 2, true);
        assert_eq!(l.n_sites, 4);
        assert_eq!(l.n_bonds, 16); // 8 physical × 2
    }

    #[test]
    fn test_square_4x4_pbc() {
        let l = build_square(4, 4, true);
        assert_eq!(l.n_sites, 16);
        assert_eq!(l.n_bonds, 64);
    }

    #[test]
    fn test_hypercubic_3d_pbc() {
        let l = build_hypercubic(
            &[2, 3, 4],
            &[BondType::SquareX, BondType::SquareY, BondType::SquareZ],
            true,
        );
        assert_eq!(l.n_sites, 24);
        // 3D PBC: every site has 6 directed neighbors → 24*6/2 physical = 72, directed = 144
        assert_eq!(l.n_bonds, 144);
    }
}
