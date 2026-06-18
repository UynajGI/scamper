//! Minimal lattice topology for QMC.rs.
//!
//! Deliberately small and decoupled from CMC.rs's [`CsrLattice`](../../cmc_rs/lattice/struct.CsrLattice.html)
//! to keep the two toolbox crates independent. If a shared lattice ever
//! becomes desirable, lift this into `carlo-rs`.
//!
//! For now only the 1D chain with periodic boundaries is provided — enough
//! for the Heisenberg-chain worm validation. Higher dimensions can be added
//! here without touching the worm layer.

/// 1D chain lattice with periodic boundary conditions.
///
/// Site `i` connects to `(i + 1) % n` and `(i + n - 1) % n`. Each site has
/// degree 2 (for `n >= 2`). Stored as a CSR-like flat neighbor list so the
/// worm's hot loop iterates a contiguous slice.
#[derive(Debug, Clone)]
pub struct ChainLattice {
    /// `offsets[i]..offsets[i+1]` indexes into `neighbors` for site `i`.
    pub offsets: Vec<usize>,
    /// Flat neighbor indices.
    pub neighbors: Vec<usize>,
    /// Number of sites.
    pub n_sites: usize,
    /// Number of directed bonds (= neighbors.len()).
    pub n_bonds: usize,
}

impl ChainLattice {
    /// Build a chain of `n` sites with periodic boundaries.
    ///
    /// `n == 1` is degenerate (the single site is its own neighbor twice);
    /// we forbid it since it's never a useful QMC lattice. Use `n >= 2`.
    pub fn new(n: usize) -> Self {
        assert!(n >= 2, "ChainLattice needs n >= 2, got {n}");
        // Each site has exactly 2 neighbors.
        let mut offsets = Vec::with_capacity(n + 1);
        let mut neighbors = Vec::with_capacity(2 * n);
        offsets.push(0);
        for i in 0..n {
            neighbors.push((i + 1) % n);
            neighbors.push((i + n - 1) % n);
            offsets.push(neighbors.len());
        }
        Self {
            offsets,
            neighbors,
            n_sites: n,
            n_bonds: 2 * n,
        }
    }

    /// Neighbor site indices for site `i`.
    #[inline]
    pub fn neighbors(&self, i: usize) -> &[usize] {
        &self.neighbors[self.offsets[i]..self.offsets[i + 1]]
    }

    /// Coordination number (2 for the chain).
    #[inline]
    pub fn degree(&self, _i: usize) -> usize {
        2
    }

    /// Iterate over undirected bonds `(i, j)` with `i < j`, once each.
    /// Useful for diagonal measurements (energy, staggered magnetization).
    pub fn bonds(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        (0..self.n_sites).map(move |i| (i, (i + 1) % self.n_sites))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chain_n4_structure() {
        let l = ChainLattice::new(4);
        assert_eq!(l.n_sites, 4);
        assert_eq!(l.n_bonds, 8); // 4 sites × degree 2
        assert_eq!(l.neighbors(0), &[1, 3]);
        assert_eq!(l.neighbors(1), &[2, 0]);
        assert_eq!(l.neighbors(3), &[0, 2]); // PBC wrap
    }

    #[test]
    fn chain_uniform_degree() {
        let l = ChainLattice::new(8);
        for i in 0..8 {
            assert_eq!(l.degree(i), 2);
        }
    }

    #[test]
    fn chain_bonds_undirected_unique() {
        let l = ChainLattice::new(4);
        let bonds: Vec<_> = l.bonds().collect();
        assert_eq!(bonds.len(), 4); // 4 undirected bonds for a 4-ring
        for (i, j) in &bonds {
            assert!(i < j || *i == 3); // (3,0) is the PBC wrap, i>j
        }
    }

    #[test]
    #[should_panic(expected = "n >= 2")]
    fn chain_rejects_single_site() {
        let _ = ChainLattice::new(1);
    }
}
