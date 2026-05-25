//! Lattice topology — CSR (Compressed Sparse Row) representation.
//!
//! The CSR format stores all neighbor indices in a single flat array,
//! giving cache-friendly iteration and eliminating per-site Vec overhead.

/// Bond direction labels (construction-time metadata).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BondType {
    ChainX,
    SquareX,
    SquareY,
    SquareZ,
    CubicX,
    CubicY,
    CubicZ,
    TriX,
    TriY,
    TriDiag,
    HoneyX,
    HoneyY,
    Kagome,
}

/// CSR-format lattice topology.
///
/// Neighbor indices for site `i` are stored contiguously in
/// `neighbors[offsets[i]..offsets[i+1]]`. This gives cache-friendly
/// iteration and O(1) access to any site's adjacency list.
#[derive(Debug, Clone)]
pub struct CsrLattice {
    /// `offsets[i]..offsets[i+1]` is the range in `neighbors` for site `i`.
    /// Length = `n_sites + 1`.
    pub offsets: Vec<usize>,
    /// Flat array of neighbor site indices.
    pub neighbors: Vec<usize>,
    /// Number of lattice sites.
    pub n_sites: usize,
    /// Number of directed bonds (= neighbors.len()).
    pub n_bonds: usize,
}

impl CsrLattice {
    /// Slice of neighbor site indices for the given site.
    #[inline]
    pub fn neighbors(&self, site: usize) -> &[usize] {
        &self.neighbors[self.offsets[site]..self.offsets[site + 1]]
    }

    /// Degree (number of neighbors) of a site.
    #[inline]
    pub fn degree(&self, site: usize) -> usize {
        self.offsets[site + 1] - self.offsets[site]
    }
}

/// Convert adjacency list to CSR format.
fn from_adjacency(sites: &[Vec<usize>]) -> CsrLattice {
    let n_sites = sites.len();
    let n_bonds: usize = sites.iter().map(|s| s.len()).sum();
    let mut offsets = Vec::with_capacity(n_sites + 1);
    let mut neighbors = Vec::with_capacity(n_bonds);
    offsets.push(0);
    for site in sites {
        neighbors.extend_from_slice(site);
        offsets.push(neighbors.len());
    }
    CsrLattice {
        offsets,
        neighbors,
        n_sites,
        n_bonds,
    }
}

// ── Builders ────────────────────────────────────────────────

/// 1D chain with optional PBC.
pub fn build_chain(n: usize, pbc: bool) -> CsrLattice {
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for i in 0..n {
        if i + 1 < n {
            adj[i].push(i + 1);
            adj[i + 1].push(i);
        } else if pbc && n > 1 {
            adj[i].push(0);
            adj[0].push(i);
        }
    }
    from_adjacency(&adj)
}

/// 2D square lattice (width × height) with optional PBC.
pub fn build_square(w: usize, h: usize, pbc: bool) -> CsrLattice {
    build_hypercubic(&[w, h], &[BondType::SquareX, BondType::SquareY], pbc)
}

/// N-dimensional hypercubic lattice.
///
/// `dims` gives size along each axis. `bond_types` labels bonds per direction
/// (kept for API compatibility). Must have `dims.len() == bond_types.len()`.
pub fn build_hypercubic(dims: &[usize], bond_types: &[BondType], pbc: bool) -> CsrLattice {
    assert_eq!(dims.len(), bond_types.len());
    let n_dims = dims.len();
    let n_sites: usize = dims.iter().product();
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n_sites];

    let mut strides = vec![1usize; n_dims];
    for k in 1..n_dims {
        strides[k] = strides[k - 1] * dims[k - 1];
    }

    for idx in 0..n_sites {
        let mut remaining = idx;
        let mut coords = vec![0usize; n_dims];
        for k in (0..n_dims).rev() {
            coords[k] = remaining / strides[k];
            remaining %= strides[k];
        }

        for dim in 0..n_dims {
            let c = coords[dim];
            let neighbor_coord = if pbc {
                (c + 1) % dims[dim]
            } else if c + 1 < dims[dim] {
                c + 1
            } else {
                continue;
            };

            let mut neighbor_coords = coords.clone();
            neighbor_coords[dim] = neighbor_coord;
            let neighbor_idx: usize = neighbor_coords
                .iter()
                .enumerate()
                .map(|(k, &v)| v * strides[k])
                .sum();

            adj[idx].push(neighbor_idx);
            adj[neighbor_idx].push(idx);
        }
    }

    from_adjacency(&adj)
}

/// 2D triangular lattice (Lx × Ly) with PBC.
///
/// Each site has 6 neighbors: ±x, ±y, ±(x+y) diagonal.
pub fn build_triangular(lx: usize, ly: usize) -> CsrLattice {
    assert!(lx >= 2 && ly >= 2, "triangular lattice needs Lx,Ly >= 2");
    let n_sites = lx * ly;
    let mut adj: Vec<Vec<usize>> = vec![Vec::with_capacity(6); n_sites];

    let idx = |x: usize, y: usize| -> usize { y * lx + x };

    for y in 0..ly {
        for x in 0..lx {
            let i = idx(x, y);
            let neighbors = [
                idx((x + 1) % lx, y),
                idx((x + lx - 1) % lx, y),
                idx(x, (y + 1) % ly),
                idx(x, (y + ly - 1) % ly),
                idx((x + 1) % lx, (y + 1) % ly),
                idx((x + lx - 1) % lx, (y + ly - 1) % ly),
            ];
            adj[i].extend_from_slice(&neighbors);
        }
    }

    from_adjacency(&adj)
}

/// 2D honeycomb lattice (Lx × Ly) with PBC.
///
/// Brick-wall representation: each site has 3 neighbors.
pub fn build_honeycomb(lx: usize, ly: usize) -> CsrLattice {
    assert!(lx >= 2, "honeycomb lattice needs Lx >= 2");
    assert!(lx.is_multiple_of(2), "honeycomb lattice needs even Lx");
    let n_sites = lx * ly;
    let mut adj: Vec<Vec<usize>> = vec![Vec::with_capacity(3); n_sites];

    let idx = |x: usize, y: usize| -> usize { y * lx + x };

    for y in 0..ly {
        for x in 0..lx {
            let i = idx(x, y);
            let right = idx((x + 1) % lx, y);
            let left = idx((x + lx - 1) % lx, y);
            let vert = if x % 2 == 0 {
                idx(x, (y + 1) % ly)
            } else {
                idx(x, (y + ly - 1) % ly)
            };
            adj[i].extend_from_slice(&[right, left, vert]);
        }
    }

    from_adjacency(&adj)
}

/// 2D kagome lattice (Lx × Ly unit cells) with PBC.
///
/// 3 sites per unit cell (A, B, C). Each site has 4 neighbors.
pub fn build_kagome(lx: usize, ly: usize) -> CsrLattice {
    assert!(lx >= 2 && ly >= 2, "kagome lattice needs Lx,Ly >= 2");
    let n_sites = 3 * lx * ly;
    let mut adj: Vec<Vec<usize>> = vec![Vec::with_capacity(4); n_sites];

    let idx = |sublat: usize, ux: usize, uy: usize| -> usize { sublat + 3 * (ux + uy * lx) };

    for uy in 0..ly {
        for ux in 0..lx {
            // Sublattice A
            {
                let i = idx(0, ux, uy);
                adj[i].extend_from_slice(&[
                    idx(1, ux, uy),
                    idx(2, ux, uy),
                    idx(1, (ux + lx - 1) % lx, uy),
                    idx(2, ux, (uy + ly - 1) % ly),
                ]);
            }
            // Sublattice B
            {
                let i = idx(1, ux, uy);
                adj[i].extend_from_slice(&[
                    idx(0, ux, uy),
                    idx(2, ux, uy),
                    idx(0, (ux + 1) % lx, uy),
                    idx(2, (ux + lx - 1) % lx, uy),
                ]);
            }
            // Sublattice C
            {
                let i = idx(2, ux, uy);
                adj[i].extend_from_slice(&[
                    idx(0, ux, uy),
                    idx(1, ux, uy),
                    idx(0, ux, (uy + 1) % ly),
                    idx(1, (ux + 1) % lx, uy),
                ]);
            }
        }
    }

    from_adjacency(&adj)
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
        assert!(l.neighbors(0).is_empty());
    }

    #[test]
    fn test_square_2x2_open() {
        let l = build_square(2, 2, false);
        assert_eq!(l.n_sites, 4);
        assert_eq!(l.n_bonds, 8);
    }

    #[test]
    fn test_square_2x2_pbc() {
        let l = build_square(2, 2, true);
        assert_eq!(l.n_sites, 4);
        assert_eq!(l.n_bonds, 16);
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
        assert_eq!(l.n_bonds, 144);
    }

    #[test]
    fn test_triangular_2x2() {
        let l = build_triangular(2, 2);
        assert_eq!(l.n_sites, 4);
        assert_eq!(l.n_bonds, 24);
        for i in 0..l.n_sites {
            assert_eq!(l.degree(i), 6);
        }
    }

    #[test]
    fn test_triangular_3x3_degree() {
        let l = build_triangular(3, 3);
        assert_eq!(l.n_sites, 9);
        assert_eq!(l.n_bonds, 54);
        for i in 0..l.n_sites {
            assert_eq!(l.degree(i), 6);
        }
    }

    #[test]
    fn test_honeycomb_2x2() {
        let l = build_honeycomb(2, 2);
        assert_eq!(l.n_sites, 4);
        assert_eq!(l.n_bonds, 12);
        for i in 0..l.n_sites {
            assert_eq!(l.degree(i), 3);
        }
    }

    #[test]
    fn test_honeycomb_4x4() {
        let l = build_honeycomb(4, 4);
        assert_eq!(l.n_sites, 16);
        assert_eq!(l.n_bonds, 48);
        for i in 0..l.n_sites {
            assert_eq!(l.degree(i), 3);
        }
    }

    #[test]
    fn test_kagome_2x2() {
        let l = build_kagome(2, 2);
        assert_eq!(l.n_sites, 12);
        assert_eq!(l.n_bonds, 48);
        for i in 0..l.n_sites {
            assert_eq!(l.degree(i), 4);
        }
    }

    #[test]
    fn test_kagome_3x3() {
        let l = build_kagome(3, 3);
        assert_eq!(l.n_sites, 27);
        assert_eq!(l.n_bonds, 108);
        for i in 0..l.n_sites {
            assert_eq!(l.degree(i), 4);
        }
    }

    #[test]
    fn test_csr_cache_friendly() {
        // Verify the CSR layout is contiguous
        let l = build_square(4, 4, true);
        assert_eq!(l.offsets.len(), l.n_sites + 1);
        assert_eq!(l.offsets[0], 0);
        assert_eq!(l.offsets[l.n_sites], l.n_bonds);
        // Every offset range is valid and non-decreasing
        for i in 0..l.n_sites {
            assert!(l.offsets[i] <= l.offsets[i + 1]);
        }
    }
}
