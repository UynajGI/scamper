//! Lattice topology — adjacency-list representation.

/// Bond direction labels.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BondType {
    ChainX,
    SquareX,
    SquareY,
    SquareZ,
    CubicX,
    CubicY,
    CubicZ,
    /// Triangular lattice horizontal bond (x direction).
    TriX,
    /// Triangular lattice vertical bond (y direction).
    TriY,
    /// Triangular lattice diagonal bond (x+y direction).
    TriDiag,
    /// Honeycomb lattice horizontal bond.
    HoneyX,
    /// Honeycomb lattice vertical bond.
    HoneyY,
    /// Kagome lattice bond (all equivalent).
    Kagome,
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

/// 2D triangular lattice (Lx × Ly) with PBC.
///
/// Each site has 6 neighbors: ±x, ±y, ±(x+y) diagonal.
/// Total directed bonds = 6 * Lx * Ly.
pub fn build_triangular(lx: usize, ly: usize) -> Lattice {
    assert!(lx >= 2 && ly >= 2, "triangular lattice needs Lx,Ly >= 2");
    let n_sites = lx * ly;
    let mut sites: Vec<Vec<Neighbor>> = vec![Vec::with_capacity(6); n_sites];

    let idx = |x: usize, y: usize| -> usize { y * lx + x };
    let mut n_bonds = 0usize;

    for y in 0..ly {
        for x in 0..lx {
            let i = idx(x, y);
            let neighbors = [
                (idx((x + 1) % lx, y), BondType::TriX),
                (idx((x + lx - 1) % lx, y), BondType::TriX),
                (idx(x, (y + 1) % ly), BondType::TriY),
                (idx(x, (y + ly - 1) % ly), BondType::TriY),
                (idx((x + 1) % lx, (y + 1) % ly), BondType::TriDiag),
                (idx((x + lx - 1) % lx, (y + ly - 1) % ly), BondType::TriDiag),
            ];
            for &(target, bt) in &neighbors {
                sites[i].push(Neighbor {
                    target,
                    bond_type: bt,
                });
                n_bonds += 1;
            }
        }
    }

    Lattice {
        sites,
        n_sites,
        n_bonds,
    }
}

/// 2D honeycomb lattice (Lx × Ly) with PBC.
///
/// Brick-wall representation: each site has 3 neighbors — 2 horizontal,
/// 1 vertical (direction alternates with column parity).
/// Total directed bonds = 3 * Lx * Ly.
pub fn build_honeycomb(lx: usize, ly: usize) -> Lattice {
    assert!(lx >= 2, "honeycomb lattice needs Lx >= 2");
    assert!(lx.is_multiple_of(2), "honeycomb lattice needs even Lx");
    let n_sites = lx * ly;
    let mut sites: Vec<Vec<Neighbor>> = vec![Vec::with_capacity(3); n_sites];

    let idx = |x: usize, y: usize| -> usize { y * lx + x };
    let mut n_bonds = 0usize;

    for y in 0..ly {
        for x in 0..lx {
            let i = idx(x, y);
            // horizontal neighbors (both directions = 2 neighbors)
            let right = idx((x + 1) % lx, y);
            sites[i].push(Neighbor {
                target: right,
                bond_type: BondType::HoneyX,
            });
            n_bonds += 1;

            let left = idx((x + lx - 1) % lx, y);
            sites[i].push(Neighbor {
                target: left,
                bond_type: BondType::HoneyX,
            });
            n_bonds += 1;

            // vertical neighbor: direction alternates with column parity
            let vert = if x % 2 == 0 {
                idx(x, (y + 1) % ly)
            } else {
                idx(x, (y + ly - 1) % ly)
            };
            sites[i].push(Neighbor {
                target: vert,
                bond_type: BondType::HoneyY,
            });
            n_bonds += 1;
        }
    }

    Lattice {
        sites,
        n_sites,
        n_bonds,
    }
}

/// 2D kagome lattice (Lx × Ly unit cells) with PBC.
///
/// 3 sites per unit cell (A, B, C). Total sites = 3 * Lx * Ly.
/// Each site has 4 neighbors via corner-sharing triangles.
/// Total directed bonds = 12 * Lx * Ly.
pub fn build_kagome(lx: usize, ly: usize) -> Lattice {
    assert!(lx >= 2 && ly >= 2, "kagome lattice needs Lx,Ly >= 2");
    let n_sites = 3 * lx * ly;
    let mut sites: Vec<Vec<Neighbor>> = vec![Vec::with_capacity(4); n_sites];

    let idx = |sublat: usize, ux: usize, uy: usize| -> usize { sublat + 3 * (ux + uy * lx) };
    let bt = BondType::Kagome;
    let mut n_bonds = 0usize;

    for uy in 0..ly {
        for ux in 0..lx {
            // Sublattice A (0): B(same), C(same), B(left), C(down)
            {
                let i = idx(0, ux, uy);
                let neighbors = [
                    idx(1, ux, uy),
                    idx(2, ux, uy),
                    idx(1, (ux + lx - 1) % lx, uy),
                    idx(2, ux, (uy + ly - 1) % ly),
                ];
                for &t in &neighbors {
                    sites[i].push(Neighbor {
                        target: t,
                        bond_type: bt,
                    });
                    n_bonds += 1;
                }
            }
            // Sublattice B (1): A(same), C(same), A(right), C(left)
            {
                let i = idx(1, ux, uy);
                let neighbors = [
                    idx(0, ux, uy),
                    idx(2, ux, uy),
                    idx(0, (ux + 1) % lx, uy),
                    idx(2, (ux + lx - 1) % lx, uy),
                ];
                for &t in &neighbors {
                    sites[i].push(Neighbor {
                        target: t,
                        bond_type: bt,
                    });
                    n_bonds += 1;
                }
            }
            // Sublattice C (2): A(same), B(same), A(up), B(right)
            {
                let i = idx(2, ux, uy);
                let neighbors = [
                    idx(0, ux, uy),
                    idx(1, ux, uy),
                    idx(0, ux, (uy + 1) % ly),
                    idx(1, (ux + 1) % lx, uy),
                ];
                for &t in &neighbors {
                    sites[i].push(Neighbor {
                        target: t,
                        bond_type: bt,
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

    // ── Triangular tests ──────────────────────────────────

    #[test]
    fn test_triangular_2x2() {
        let l = build_triangular(2, 2);
        assert_eq!(l.n_sites, 4);
        assert_eq!(l.n_bonds, 24); // 4 sites × 6 neighbors
        for site in &l.sites {
            assert_eq!(site.len(), 6);
        }
    }

    #[test]
    fn test_triangular_3x3_degree() {
        let l = build_triangular(3, 3);
        assert_eq!(l.n_sites, 9);
        assert_eq!(l.n_bonds, 54); // 9 × 6
        for site in &l.sites {
            assert_eq!(site.len(), 6);
        }
    }

    #[test]
    fn test_triangular_bond_types() {
        let l = build_triangular(4, 4);
        let mut counts = [0usize; 3]; // TriX, TriY, TriDiag
        for site in &l.sites {
            for nb in site {
                match nb.bond_type {
                    BondType::TriX => counts[0] += 1,
                    BondType::TriY => counts[1] += 1,
                    BondType::TriDiag => counts[2] += 1,
                    _ => panic!("unexpected bond type"),
                }
            }
        }
        // Each site has 2 TriX, 2 TriY, 2 TriDiag = 16*2 each = 32 each
        assert_eq!(counts[0], 32);
        assert_eq!(counts[1], 32);
        assert_eq!(counts[2], 32);
    }

    // ── Honeycomb tests ───────────────────────────────────

    #[test]
    fn test_honeycomb_2x2() {
        let l = build_honeycomb(2, 2);
        assert_eq!(l.n_sites, 4);
        assert_eq!(l.n_bonds, 12); // 4 sites × 3 neighbors
        for site in &l.sites {
            assert_eq!(site.len(), 3);
        }
    }

    #[test]
    fn test_honeycomb_4x4() {
        let l = build_honeycomb(4, 4);
        assert_eq!(l.n_sites, 16);
        assert_eq!(l.n_bonds, 48); // 16 × 3
        for site in &l.sites {
            assert_eq!(site.len(), 3);
        }
    }

    #[test]
    fn test_honeycomb_bond_types() {
        let l = build_honeycomb(4, 4);
        let mut nx = 0usize;
        let mut ny = 0usize;
        for site in &l.sites {
            for nb in site {
                match nb.bond_type {
                    BondType::HoneyX => nx += 1,
                    BondType::HoneyY => ny += 1,
                    _ => panic!("unexpected bond type"),
                }
            }
        }
        // 16 sites × 2 HoneyX = 32, 16 sites × 1 HoneyY = 16
        assert_eq!(nx, 32);
        assert_eq!(ny, 16);
    }

    // ── Kagome tests ──────────────────────────────────────

    #[test]
    fn test_kagome_2x2() {
        let l = build_kagome(2, 2);
        assert_eq!(l.n_sites, 12); // 3 × 2 × 2
        assert_eq!(l.n_bonds, 48); // 12 sites × 4 neighbors
        for site in &l.sites {
            assert_eq!(site.len(), 4);
        }
    }

    #[test]
    fn test_kagome_3x3() {
        let l = build_kagome(3, 3);
        assert_eq!(l.n_sites, 27); // 3 × 3 × 3
        assert_eq!(l.n_bonds, 108); // 27 × 4
        for site in &l.sites {
            assert_eq!(site.len(), 4);
        }
    }

    #[test]
    fn test_kagome_bond_types() {
        let l = build_kagome(3, 3);
        for site in &l.sites {
            for nb in site {
                assert_eq!(nb.bond_type, BondType::Kagome);
            }
        }
    }
}
