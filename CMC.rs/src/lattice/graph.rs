//! Graph/lattice topology in compressed sparse-row form.
//!
//! `CsrLattice` deliberately separates **physical undirected bonds** from
//! adjacency incidences.  A physical bond is stored exactly once in `edges`,
//! while every endpoint receives one CSR incidence.  This removes the old
//! "sum directed neighbours and divide by two" convention and makes weighted,
//! typed and parallel bonds unambiguous.

use std::collections::BTreeMap;

/// Construction-time bond labels.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum BondType {
    #[default]
    Generic,
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

impl BondType {
    /// Stable serialization label (NOT derived from Debug).
    pub fn as_label(self) -> &'static str {
        match self {
            BondType::Generic => "generic",
            BondType::ChainX => "chain_x",
            BondType::SquareX => "square_x",
            BondType::SquareY => "square_y",
            BondType::SquareZ => "square_z",
            BondType::CubicX => "cubic_x",
            BondType::CubicY => "cubic_y",
            BondType::CubicZ => "cubic_z",
            BondType::TriX => "tri_x",
            BondType::TriY => "tri_y",
            BondType::TriDiag => "tri_diag",
            BondType::HoneyX => "honey_x",
            BondType::HoneyY => "honey_y",
            BondType::Kagome => "kagome",
        }
    }

    /// Inverse of `as_label`.
    pub fn from_label(label: &str) -> Option<Self> {
        match label {
            "generic" => Some(Self::Generic),
            "chain_x" => Some(Self::ChainX),
            "square_x" => Some(Self::SquareX),
            "square_y" => Some(Self::SquareY),
            "square_z" => Some(Self::SquareZ),
            "cubic_x" => Some(Self::CubicX),
            "cubic_y" => Some(Self::CubicY),
            "cubic_z" => Some(Self::CubicZ),
            "tri_x" => Some(Self::TriX),
            "tri_y" => Some(Self::TriY),
            "tri_diag" => Some(Self::TriDiag),
            "honey_x" => Some(Self::HoneyX),
            "honey_y" => Some(Self::HoneyY),
            "kagome" => Some(Self::Kagome),
            _ => None,
        }
    }
}

/// One physical undirected bond.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bond {
    pub source: usize,
    pub target: usize,
    pub kind: BondType,
    /// Multiplicative coupling weight.  Built-in lattices use `1.0`.
    pub weight: f64,
}

impl Bond {
    pub const fn new(source: usize, target: usize, kind: BondType, weight: f64) -> Self {
        Self {
            source,
            target,
            kind,
            weight,
        }
    }

    #[inline]
    pub fn other(self, site: usize) -> Option<usize> {
        if site == self.source {
            Some(self.target)
        } else if site == self.target {
            Some(self.source)
        } else {
            None
        }
    }
}

/// CSR-format arbitrary undirected multigraph.
///
/// Compatibility fields `offsets`, `neighbors`, `n_sites`, and `n_bonds` are
/// retained.  `n_bonds` means the number of directed incidences, exactly as in
/// the original crate.  Use [`CsrLattice::n_edges`] for physical bond count.
#[derive(Debug, Clone)]
pub struct CsrLattice {
    /// `offsets[i]..offsets[i+1]` is site `i`'s incidence range.
    pub offsets: Vec<usize>,
    /// Neighbor endpoint for each incidence.
    pub neighbors: Vec<usize>,
    /// Physical edge id corresponding to each entry of `neighbors`.
    pub edge_ids: Vec<usize>,
    /// Physical undirected edges, each stored exactly once.
    pub edges: Vec<Bond>,
    pub n_sites: usize,
    /// Number of directed incidences (`neighbors.len()`).
    pub n_bonds: usize,
}

impl CsrLattice {
    /// Build an arbitrary graph from physical bonds.
    ///
    /// Parallel bonds and self-loops are supported.  A self-loop contributes
    /// two incidences to preserve the usual graph-theoretic degree convention.
    pub fn try_from_edges(n_sites: usize, edges: Vec<Bond>) -> Result<Self, String> {
        if n_sites == 0 {
            return Err("lattice must contain at least one site".to_string());
        }

        let mut degree = vec![0usize; n_sites];
        for (edge_id, edge) in edges.iter().enumerate() {
            if edge.source >= n_sites || edge.target >= n_sites {
                return Err(format!(
                    "edge {edge_id} endpoint out of range: ({}, {}) for {n_sites} sites",
                    edge.source, edge.target
                ));
            }
            if !edge.weight.is_finite() {
                return Err(format!("edge {edge_id} has non-finite weight"));
            }
            degree[edge.source] += 1;
            degree[edge.target] += 1;
        }

        let mut offsets = Vec::with_capacity(n_sites + 1);
        offsets.push(0);
        for d in degree {
            offsets.push(offsets.last().copied().unwrap_or(0) + d);
        }

        let n_bonds = offsets[n_sites];
        let mut neighbors = vec![0usize; n_bonds];
        let mut edge_ids = vec![0usize; n_bonds];
        let mut cursor = offsets[..n_sites].to_vec();

        for (edge_id, edge) in edges.iter().enumerate() {
            let left = cursor[edge.source];
            neighbors[left] = edge.target;
            edge_ids[left] = edge_id;
            cursor[edge.source] += 1;

            let right = cursor[edge.target];
            neighbors[right] = edge.source;
            edge_ids[right] = edge_id;
            cursor[edge.target] += 1;
        }

        let lattice = Self {
            offsets,
            neighbors,
            edge_ids,
            edges,
            n_sites,
            n_bonds,
        };
        lattice.validate()?;
        Ok(lattice)
    }

    /// Infallible convenience constructor for programmatically valid edges.
    pub fn from_edges(n_sites: usize, edges: Vec<Bond>) -> Self {
        Self::try_from_edges(n_sites, edges).expect("invalid lattice edge list")
    }

    /// Convert a symmetric adjacency list into a physical multigraph.
    ///
    /// For each pair `(i,j)`, the number of `i→j` and `j→i` entries must match.
    /// This is useful for importing existing neighbor-list based geometries.
    pub fn try_from_adjacency(sites: &[Vec<usize>]) -> Result<Self, String> {
        if sites.is_empty() {
            return Err("lattice must contain at least one site".to_string());
        }

        let n_sites = sites.len();
        let mut counts: BTreeMap<(usize, usize), (usize, usize)> = BTreeMap::new();
        let mut self_counts = vec![0usize; n_sites];

        for (source, row) in sites.iter().enumerate() {
            for &target in row {
                if target >= n_sites {
                    return Err(format!(
                        "adjacency endpoint out of range: {source} -> {target} for {n_sites} sites"
                    ));
                }
                if source == target {
                    self_counts[source] += 1;
                } else {
                    let key = if source < target {
                        (source, target)
                    } else {
                        (target, source)
                    };
                    let entry = counts.entry(key).or_insert((0, 0));
                    if source < target {
                        entry.0 += 1;
                    } else {
                        entry.1 += 1;
                    }
                }
            }
        }

        let mut edges = Vec::new();
        for ((source, target), (forward, reverse)) in counts {
            if forward != reverse {
                return Err(format!(
                    "adjacency is not symmetric for ({source}, {target}): {forward} vs {reverse}"
                ));
            }
            for _ in 0..forward {
                edges.push(Bond::new(source, target, BondType::Generic, 1.0));
            }
        }
        for (site, count) in self_counts.into_iter().enumerate() {
            if count % 2 != 0 {
                return Err(format!(
                    "self-loop adjacency at site {site} must contain an even number of incidences"
                ));
            }
            for _ in 0..count / 2 {
                edges.push(Bond::new(site, site, BondType::Generic, 1.0));
            }
        }

        Self::try_from_edges(n_sites, edges)
    }

    pub fn from_adjacency(sites: &[Vec<usize>]) -> Self {
        Self::try_from_adjacency(sites).expect("invalid symmetric adjacency list")
    }

    /// Validate all CSR and physical-edge invariants.
    pub fn validate(&self) -> Result<(), String> {
        if self.n_sites == 0 {
            return Err("n_sites must be positive".to_string());
        }
        if self.offsets.len() != self.n_sites + 1 {
            return Err(format!(
                "offset length mismatch: expected {}, got {}",
                self.n_sites + 1,
                self.offsets.len()
            ));
        }
        if self.offsets.first().copied() != Some(0) {
            return Err("offsets must start at zero".to_string());
        }
        if self.offsets.windows(2).any(|w| w[0] > w[1]) {
            return Err("offsets must be non-decreasing".to_string());
        }
        if self.neighbors.len() != self.edge_ids.len() || self.neighbors.len() != self.n_bonds {
            return Err("incidence array length mismatch".to_string());
        }
        if self.offsets[self.n_sites] != self.n_bonds {
            return Err("last offset must equal n_bonds".to_string());
        }

        let mut incidence_counts = vec![0usize; self.edges.len()];
        let mut endpoint_counts = vec![[0usize; 2]; self.edges.len()];
        for site in 0..self.n_sites {
            for incidence in self.offsets[site]..self.offsets[site + 1] {
                let neighbor = self.neighbors[incidence];
                let edge_id = self.edge_ids[incidence];
                if neighbor >= self.n_sites || edge_id >= self.edges.len() {
                    return Err(format!("invalid incidence {incidence} at site {site}"));
                }
                let edge = self.edges[edge_id];
                if edge.other(site) != Some(neighbor) {
                    return Err(format!(
                        "incidence {incidence} does not match physical edge {edge_id}"
                    ));
                }
                incidence_counts[edge_id] += 1;
                if edge.source == edge.target || site == edge.source {
                    endpoint_counts[edge_id][0] += 1;
                } else {
                    endpoint_counts[edge_id][1] += 1;
                }
            }
        }
        for (edge_id, count) in incidence_counts.into_iter().enumerate() {
            if count != 2 {
                return Err(format!(
                    "physical edge {edge_id} must have exactly two incidences, got {count}"
                ));
            }
            let edge = self.edges[edge_id];
            let expected = if edge.source == edge.target {
                [2, 0]
            } else {
                [1, 1]
            };
            if endpoint_counts[edge_id] != expected {
                return Err(format!(
                    "physical edge {edge_id} has invalid endpoint incidence counts {:?}",
                    endpoint_counts[edge_id]
                ));
            }
        }
        Ok(())
    }

    #[inline]
    pub fn neighbors(&self, site: usize) -> &[usize] {
        &self.neighbors[self.offsets[site]..self.offsets[site + 1]]
    }

    #[inline]
    pub fn edge_ids(&self, site: usize) -> &[usize] {
        &self.edge_ids[self.offsets[site]..self.offsets[site + 1]]
    }

    #[inline]
    pub fn incidences(&self, site: usize) -> impl Iterator<Item = (usize, usize)> + '_ {
        let range = self.offsets[site]..self.offsets[site + 1];
        self.neighbors[range.clone()]
            .iter()
            .copied()
            .zip(self.edge_ids[range].iter().copied())
    }

    #[inline]
    pub fn degree(&self, site: usize) -> usize {
        self.offsets[site + 1] - self.offsets[site]
    }

    #[inline]
    pub fn n_edges(&self) -> usize {
        self.edges.len()
    }
}

fn validate_dims(dims: &[usize], bond_types: &[BondType]) {
    assert!(
        !dims.is_empty(),
        "at least one lattice dimension is required"
    );
    assert_eq!(dims.len(), bond_types.len());
    assert!(
        dims.iter().all(|&d| d > 0),
        "lattice dimensions must be positive"
    );
}

/// 1D chain with optional periodic boundaries.
pub fn build_chain(n: usize, pbc: bool) -> CsrLattice {
    assert!(n > 0, "chain length must be positive");
    build_hypercubic(&[n], &[BondType::ChainX], pbc)
}

/// 2D square lattice (`width × height`).
pub fn build_square(width: usize, height: usize, pbc: bool) -> CsrLattice {
    build_hypercubic(
        &[width, height],
        &[BondType::SquareX, BondType::SquareY],
        pbc,
    )
}

/// N-dimensional hypercubic lattice represented as an arbitrary graph.
pub fn build_hypercubic(dims: &[usize], bond_types: &[BondType], pbc: bool) -> CsrLattice {
    validate_dims(dims, bond_types);
    let n_dims = dims.len();
    let n_sites: usize = dims.iter().product();
    let mut strides = vec![1usize; n_dims];
    for axis in 1..n_dims {
        strides[axis] = strides[axis - 1] * dims[axis - 1];
    }

    let mut edges = Vec::new();
    for site in 0..n_sites {
        let mut coords = vec![0usize; n_dims];
        let mut remaining = site;
        for axis in (0..n_dims).rev() {
            coords[axis] = remaining / strides[axis];
            remaining %= strides[axis];
        }

        for axis in 0..n_dims {
            let coordinate = coords[axis];
            let next = if coordinate + 1 < dims[axis] {
                coordinate + 1
            } else if pbc && dims[axis] > 1 {
                0
            } else {
                continue;
            };

            let target = site - coordinate * strides[axis] + next * strides[axis];
            edges.push(Bond::new(site, target, bond_types[axis], 1.0));
        }
    }

    CsrLattice::from_edges(n_sites, edges)
}

/// 2D triangular lattice with periodic boundaries.
pub fn build_triangular(lx: usize, ly: usize) -> CsrLattice {
    assert!(lx >= 2 && ly >= 2, "triangular lattice needs Lx,Ly >= 2");
    let index = |x: usize, y: usize| y * lx + x;
    let mut edges = Vec::with_capacity(3 * lx * ly);
    for y in 0..ly {
        for x in 0..lx {
            let site = index(x, y);
            edges.push(Bond::new(site, index((x + 1) % lx, y), BondType::TriX, 1.0));
            edges.push(Bond::new(site, index(x, (y + 1) % ly), BondType::TriY, 1.0));
            edges.push(Bond::new(
                site,
                index((x + 1) % lx, (y + 1) % ly),
                BondType::TriDiag,
                1.0,
            ));
        }
    }
    CsrLattice::from_edges(lx * ly, edges)
}

/// 2D honeycomb lattice in brick-wall representation with periodic boundaries.
pub fn build_honeycomb(lx: usize, ly: usize) -> CsrLattice {
    assert!(lx >= 2, "honeycomb lattice needs Lx >= 2");
    assert!(ly >= 2, "honeycomb lattice needs Ly >= 2");
    assert!(lx.is_multiple_of(2), "honeycomb lattice needs even Lx");
    let index = |x: usize, y: usize| y * lx + x;
    let mut adjacency = vec![Vec::with_capacity(3); lx * ly];

    for y in 0..ly {
        for x in 0..lx {
            let matched = if x % 2 == 0 {
                index((x + 1) % lx, (y + 1) % ly)
            } else {
                index((x + lx - 1) % lx, (y + ly - 1) % ly)
            };
            adjacency[index(x, y)].extend_from_slice(&[
                index((x + 1) % lx, y),
                index((x + lx - 1) % lx, y),
                matched,
            ]);
        }
    }

    let mut lattice = CsrLattice::from_adjacency(&adjacency);
    for edge in &mut lattice.edges {
        edge.kind = if edge.source / lx == edge.target / lx {
            BondType::HoneyX
        } else {
            BondType::HoneyY
        };
    }
    lattice
}

/// 2D kagome lattice (`3 × Lx × Ly` sites) with periodic boundaries.
pub fn build_kagome(lx: usize, ly: usize) -> CsrLattice {
    assert!(lx >= 2 && ly >= 2, "kagome lattice needs Lx,Ly >= 2");
    let n_sites = 3 * lx * ly;
    let mut adjacency = vec![Vec::with_capacity(4); n_sites];
    let index = |sublattice: usize, x: usize, y: usize| sublattice + 3 * (x + y * lx);

    for y in 0..ly {
        for x in 0..lx {
            adjacency[index(0, x, y)].extend_from_slice(&[
                index(1, x, y),
                index(2, x, y),
                index(1, (x + lx - 1) % lx, y),
                index(2, x, (y + ly - 1) % ly),
            ]);
            adjacency[index(1, x, y)].extend_from_slice(&[
                index(0, x, y),
                index(2, x, y),
                index(0, (x + 1) % lx, y),
                index(2, (x + lx - 1) % lx, y),
            ]);
            adjacency[index(2, x, y)].extend_from_slice(&[
                index(0, x, y),
                index(1, x, y),
                index(0, x, (y + 1) % ly),
                index(1, (x + 1) % lx, y),
            ]);
        }
    }

    let mut lattice = CsrLattice::from_adjacency(&adjacency);
    for edge in &mut lattice.edges {
        edge.kind = BondType::Kagome;
    }
    lattice
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chain_counts() {
        assert_eq!(build_chain(4, false).n_bonds, 6);
        assert_eq!(build_chain(4, true).n_bonds, 8);
        assert_eq!(build_chain(1, true).n_bonds, 0);
    }

    #[test]
    fn square_counts() {
        assert_eq!(build_square(2, 2, false).n_bonds, 8);
        assert_eq!(build_square(2, 2, true).n_bonds, 16);
        assert_eq!(build_square(4, 4, true).n_edges(), 32);
    }

    #[test]
    fn non_bravais_degrees() {
        for site in 0..build_triangular(3, 3).n_sites {
            assert_eq!(build_triangular(3, 3).degree(site), 6);
        }
        for (lx, ly) in [(2, 3), (4, 4), (6, 5)] {
            let honey = build_honeycomb(lx, ly);
            assert!(honey.validate().is_ok());
            assert!((0..honey.n_sites).all(|site| honey.degree(site) == 3));
        }
        let kagome = build_kagome(3, 3);
        assert!((0..kagome.n_sites).all(|site| kagome.degree(site) == 4));
    }

    #[test]
    fn weighted_parallel_edges_are_preserved() {
        let graph = CsrLattice::from_edges(
            2,
            vec![
                Bond::new(0, 1, BondType::Generic, 1.0),
                Bond::new(0, 1, BondType::Generic, 2.0),
            ],
        );
        assert_eq!(graph.n_edges(), 2);
        assert_eq!(graph.degree(0), 2);
        assert_eq!(graph.degree(1), 2);
    }
}
