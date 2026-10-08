//! Graph/lattice topology in compressed sparse-row form.
//!
//! `CsrLattice` deliberately separates **physical undirected bonds** from
//! adjacency incidences.  A physical bond is stored exactly once in `edges`,
//! while every endpoint receives one CSR incidence.  This removes the old
//! "sum directed neighbours and divide by two" convention and makes weighted,
//! typed and parallel bonds unambiguous.

use std::collections::BTreeMap;

use crate::topology::{EmbeddingError, LatticeEmbedding};

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

    /// Connected components of the bond graph, each as an ascending site list.
    ///
    /// Components are ordered by their smallest site; isolated sites form
    /// their own single-site component. Every physical edge lies fully inside
    /// exactly one component, so ensembles that factorize over components can
    /// be sampled by decomposing along this partition.
    pub fn connected_components(&self) -> Vec<Vec<usize>> {
        let mut components = Vec::new();
        let mut seen = vec![false; self.n_sites];
        for seed in 0..self.n_sites {
            if seen[seed] {
                continue;
            }
            seen[seed] = true;
            let mut stack = vec![seed];
            let mut members = Vec::new();
            while let Some(site) = stack.pop() {
                members.push(site);
                for &neighbor in self.neighbors(site) {
                    if !seen[neighbor] {
                        seen[neighbor] = true;
                        stack.push(neighbor);
                    }
                }
            }
            members.sort_unstable();
            components.push(members);
        }
        components
    }
}

impl crate::topology::GraphView for CsrLattice {
    #[inline(always)]
    fn vertex_count(&self) -> usize {
        self.n_sites
    }
}

impl crate::topology::UndirectedGraphView for CsrLattice {
    #[inline(always)]
    fn edge_count(&self) -> usize {
        self.edges.len()
    }

    #[inline(always)]
    fn edge_endpoints(&self, edge: crate::topology::EdgeId) -> [crate::topology::VertexId; 2] {
        let edge = self.edges[edge.index()];
        [
            crate::topology::VertexId::from_valid_index(edge.source),
            crate::topology::VertexId::from_valid_index(edge.target),
        ]
    }

    #[inline(always)]
    fn incidences(
        &self,
        vertex: crate::topology::VertexId,
    ) -> impl Iterator<Item = crate::topology::Incidence> + '_ {
        let range = self.offsets[vertex.index()]..self.offsets[vertex.index() + 1];
        self.neighbors[range.clone()]
            .iter()
            .copied()
            .zip(self.edge_ids[range].iter().copied())
            .map(|(neighbor, edge)| crate::topology::Incidence {
                neighbor: crate::topology::VertexId::from_valid_index(neighbor),
                edge: crate::topology::EdgeId::from_valid_index(edge),
            })
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

/// 1D ring plus its periodic embedding; see
/// [`build_hypercubic_with_embedding`] for the boundary semantics.
pub fn build_chain_with_embedding(
    n: usize,
    pbc: bool,
) -> Result<(CsrLattice, LatticeEmbedding), EmbeddingError> {
    assert!(n > 0, "chain length must be positive");
    build_hypercubic_with_embedding(&[n], &[BondType::ChainX], pbc)
}

/// 2D square lattice (`width × height`).
pub fn build_square(width: usize, height: usize, pbc: bool) -> CsrLattice {
    build_hypercubic(
        &[width, height],
        &[BondType::SquareX, BondType::SquareY],
        pbc,
    )
}

/// 2D torus plus its periodic embedding; see
/// [`build_hypercubic_with_embedding`] for the boundary semantics.
pub fn build_square_with_embedding(
    width: usize,
    height: usize,
    pbc: bool,
) -> Result<(CsrLattice, LatticeEmbedding), EmbeddingError> {
    build_hypercubic_with_embedding(
        &[width, height],
        &[BondType::SquareX, BondType::SquareY],
        pbc,
    )
}

/// N-dimensional hypercubic lattice represented as an arbitrary graph.
pub fn build_hypercubic(dims: &[usize], bond_types: &[BondType], pbc: bool) -> CsrLattice {
    let (n_sites, edges, _) = hypercubic_bonds(dims, bond_types, pbc, false);
    CsrLattice::from_edges(n_sites, edges)
}

/// Hypercubic lattice plus its periodic embedding (W1).
///
/// With `pbc = true` this returns the same lattice as [`build_hypercubic`]
/// together with a [`LatticeEmbedding`] whose per-edge cell displacements were
/// recorded while the bonds were generated — the only stage where the wrap
/// information exists unambiguously. Displacements are never reconstructed
/// afterwards from `BondType` or vertex IDs: on small periodic cells the same
/// endpoint pair is joined by parallel edges with different displacements, so
/// endpoint geometry alone cannot identify them.
///
/// `pbc = false` is rejected with [`EmbeddingError::OpenBoundaries`]: an
/// open-boundary graph is legitimately non-periodic and has no wrapping
/// queries, so it has no embedding at all rather than a vacuous one.
pub fn build_hypercubic_with_embedding(
    dims: &[usize],
    bond_types: &[BondType],
    pbc: bool,
) -> Result<(CsrLattice, LatticeEmbedding), EmbeddingError> {
    if !pbc {
        return Err(EmbeddingError::OpenBoundaries);
    }
    let (n_sites, edges, displacements) = hypercubic_bonds(dims, bond_types, true, true);
    let lattice = CsrLattice::from_edges(n_sites, edges);
    let embedding =
        LatticeEmbedding::try_from_edge_displacements(&lattice, dims.len(), &displacements)?;
    Ok((lattice, embedding))
}

/// Shared hypercubic bond loop.
///
/// Every bond moves one fundamental cell forward along one axis from source to
/// target. When `track_displacements` is set, the source-to-target cell
/// displacement is appended per bond (flat, edge-major): `+1` on the bond's
/// axis exactly when the bond wraps the periodic boundary, `0` otherwise.
/// [`BondType`] is carried through untouched and never consulted for geometry.
fn hypercubic_bonds(
    dims: &[usize],
    bond_types: &[BondType],
    pbc: bool,
    track_displacements: bool,
) -> (usize, Vec<Bond>, Vec<i32>) {
    validate_dims(dims, bond_types);
    let n_dims = dims.len();
    let n_sites: usize = dims.iter().product();
    let mut strides = vec![1usize; n_dims];
    for axis in 1..n_dims {
        strides[axis] = strides[axis - 1] * dims[axis - 1];
    }

    let mut edges = Vec::new();
    let mut displacements = Vec::new();
    for site in 0..n_sites {
        let mut coords = vec![0usize; n_dims];
        let mut remaining = site;
        for axis in (0..n_dims).rev() {
            coords[axis] = remaining / strides[axis];
            remaining %= strides[axis];
        }

        for axis in 0..n_dims {
            let coordinate = coords[axis];
            let (next, wraps) = if coordinate + 1 < dims[axis] {
                (coordinate + 1, false)
            } else if pbc && dims[axis] > 1 {
                (0, true)
            } else {
                continue;
            };

            let target = site - coordinate * strides[axis] + next * strides[axis];
            edges.push(Bond::new(site, target, bond_types[axis], 1.0));
            if track_displacements {
                displacements.extend((0..n_dims).map(|a| i32::from(a == axis && wraps)));
            }
        }
    }

    (n_sites, edges, displacements)
}

/// 2D triangular lattice with periodic boundaries.
pub fn build_triangular(lx: usize, ly: usize) -> CsrLattice {
    let (n_sites, edges, _) = triangular_bonds(lx, ly, false);
    CsrLattice::from_edges(n_sites, edges)
}

/// 2D triangular torus plus its periodic embedding; see
/// [`build_hypercubic_with_embedding`] for the boundary semantics.
///
/// Diagonal bonds lift `(x, y)` to `(x + 1, y + 1)` in the covering lattice,
/// so they wrap each axis independently exactly when the emitting site sits
/// on the last coordinate of that axis.
pub fn build_triangular_with_embedding(
    lx: usize,
    ly: usize,
) -> Result<(CsrLattice, LatticeEmbedding), EmbeddingError> {
    let (n_sites, edges, displacements) = triangular_bonds(lx, ly, true);
    let lattice = CsrLattice::from_edges(n_sites, edges);
    let embedding = LatticeEmbedding::try_from_edge_displacements(&lattice, 2, &displacements)?;
    Ok((lattice, embedding))
}

/// Shared triangular bond loop; see [`hypercubic_bonds`] for the
/// displacement-tracking convention.
fn triangular_bonds(
    lx: usize,
    ly: usize,
    track_displacements: bool,
) -> (usize, Vec<Bond>, Vec<i32>) {
    assert!(lx >= 2 && ly >= 2, "triangular lattice needs Lx,Ly >= 2");
    let index = |x: usize, y: usize| y * lx + x;
    let mut edges = Vec::with_capacity(3 * lx * ly);
    let mut displacements = Vec::new();
    for y in 0..ly {
        for x in 0..lx {
            let site = index(x, y);
            let (x_next, x_wraps) = ((x + 1) % lx, x + 1 == lx);
            let (y_next, y_wraps) = ((y + 1) % ly, y + 1 == ly);
            edges.push(Bond::new(site, index(x_next, y), BondType::TriX, 1.0));
            edges.push(Bond::new(site, index(x, y_next), BondType::TriY, 1.0));
            edges.push(Bond::new(
                site,
                index(x_next, y_next),
                BondType::TriDiag,
                1.0,
            ));
            if track_displacements {
                displacements.extend([
                    i32::from(x_wraps),
                    0,
                    0,
                    i32::from(y_wraps),
                    i32::from(x_wraps),
                    i32::from(y_wraps),
                ]);
            }
        }
    }
    (lx * ly, edges, displacements)
}

/// 2D honeycomb lattice in brick-wall representation with periodic boundaries.
///
/// No `_with_embedding` variant exists: this builder assembles a symmetric
/// adjacency list, and `CsrLattice::from_adjacency` drops the per-bond
/// construction direction while deduplicating to physical edges. On `lx = 2`
/// the two horizontal bonds between a site pair carry different cell
/// displacements that the emitted bond list can no longer distinguish, so
/// recovering them would be guessing (W1 forbids that). A honeycomb embedding
/// requires a builder that records displacements while emitting bonds.
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
///
/// No `_with_embedding` variant exists, for the same reason as
/// [`build_honeycomb`]: the symmetric-adjacency construction loses the
/// per-bond direction, and on small cells distinct covering bonds collapse
/// onto the same endpoint pair, so per-edge cell displacements cannot be
/// attached without guessing.
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
    use crate::topology::{GraphView, PeriodicEmbedding, UndirectedGraphView};

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

    #[test]
    fn connected_components_partition_disconnected_graphs() {
        // Two disjoint bonds plus an isolated site.
        let graph = CsrLattice::from_edges(
            5,
            vec![
                Bond::new(3, 4, BondType::Generic, 1.0),
                Bond::new(0, 1, BondType::Generic, 1.0),
            ],
        );
        let components = graph.connected_components();
        assert_eq!(
            components,
            vec![vec![0, 1], vec![2], vec![3, 4]],
            "components must be ordered by smallest site, ascending within"
        );
        // Every edge lies fully inside exactly one component.
        let covered: usize = components
            .iter()
            .map(|sites| {
                graph
                    .edges
                    .iter()
                    .filter(|edge| sites.contains(&edge.source))
                    .count()
            })
            .sum();
        assert_eq!(covered, graph.n_edges());

        // Connected graphs decompose into a single component.
        assert_eq!(build_chain(6, true).connected_components().len(), 1);
        assert_eq!(
            build_square(3, 3, true).connected_components(),
            vec![(0..9).collect::<Vec<_>>()]
        );
    }

    // ---- W1 periodic embedding ----

    /// Displacement of one directed incidence, as an owned vector.
    fn displacement(
        lattice: &CsrLattice,
        embedding: &LatticeEmbedding,
        edge_index: usize,
        from_index: usize,
    ) -> Vec<i32> {
        let edge = lattice.edge_id(edge_index).expect("dense edge index");
        let from = lattice.vertex_id(from_index).expect("dense vertex index");
        embedding
            .edge_displacement(edge, from)
            .expect("endpoint query")
            .iter()
            .collect()
    }

    /// Physical bond index joining two sites, if any (unique at L >= 3).
    fn bond_between(lattice: &CsrLattice, a: usize, b: usize) -> Option<usize> {
        lattice
            .incidences(a)
            .find_map(|(neighbor, edge)| (neighbor == b).then_some(edge))
    }

    /// Accumulated displacement around a closed site walk.
    fn cycle_sum(lattice: &CsrLattice, embedding: &LatticeEmbedding, walk: &[usize]) -> Vec<i32> {
        let mut sum = vec![0i32; embedding.periodic_dimensions()];
        for index in 0..walk.len() {
            let (from, to) = (walk[index], walk[(index + 1) % walk.len()]);
            let edge_index = bond_between(lattice, from, to).expect("cycle bond");
            let edge = lattice.edge_id(edge_index).expect("dense edge index");
            let from_id = lattice.vertex_id(from).expect("dense vertex index");
            let displacement = embedding
                .edge_displacement(edge, from_id)
                .expect("endpoint query");
            for (axis, total) in sum.iter_mut().enumerate() {
                *total += displacement.axis(axis).expect("periodic axis");
            }
        }
        sum
    }

    /// Hypercubic strides: `strides[0] = 1`, `strides[a] = strides[a-1] * dims[a-1]`.
    fn hypercubic_strides(dims: &[usize]) -> Vec<usize> {
        let mut strides = vec![1usize; dims.len()];
        for axis in 1..dims.len() {
            strides[axis] = strides[axis - 1] * dims[axis - 1];
        }
        strides
    }

    /// Bond-type labels matching a hypercubic dimension count.
    fn bond_types_for(n_dims: usize) -> Vec<BondType> {
        match n_dims {
            1 => vec![BondType::ChainX],
            2 => vec![BondType::SquareX, BondType::SquareY],
            _ => vec![BondType::CubicX, BondType::CubicY, BondType::CubicZ],
        }
    }

    /// Mixed-radix coordinates under the hypercubic site convention.
    fn site_coords(dims: &[usize], site: usize) -> Vec<usize> {
        let strides = hypercubic_strides(dims);
        let mut coords = vec![0usize; dims.len()];
        let mut remaining = site;
        for axis in (0..dims.len()).rev() {
            coords[axis] = remaining / strides[axis];
            remaining %= strides[axis];
        }
        coords
    }

    #[test]
    fn square_2x2_torus_displacements_match_hand_derivation() {
        let (lattice, embedding) =
            build_square_with_embedding(2, 2, true).expect("torus embedding");
        assert_eq!(embedding.periodic_dimensions(), 2);
        // Sites are indexed y*2+x. Each of the four vertex pairs carries two
        // parallel bonds: an in-cell one and a wrapping one, in emission order.
        let expected: [(usize, usize, [i32; 2]); 8] = [
            (0, 0, [0, 0]), // 0 -> 1, in-cell x
            (1, 0, [0, 0]), // 0 -> 2, in-cell y
            (2, 1, [1, 0]), // 1 -> 0, wraps +x
            (3, 1, [0, 0]), // 1 -> 3, in-cell y
            (4, 2, [0, 0]), // 2 -> 3, in-cell x
            (5, 2, [0, 1]), // 2 -> 0, wraps +y
            (6, 3, [1, 0]), // 3 -> 2, wraps +x
            (7, 3, [0, 1]), // 3 -> 1, wraps +y
        ];
        for (edge_index, from, vector) in expected {
            assert_eq!(
                displacement(&lattice, &embedding, edge_index, from),
                vector.to_vec(),
                "edge {edge_index} from site {from}"
            );
            let [source, target] =
                lattice.edge_endpoints(lattice.edge_id(edge_index).expect("dense edge index"));
            let other = if from == source.index() {
                target.index()
            } else {
                source.index()
            };
            assert_eq!(
                displacement(&lattice, &embedding, edge_index, other),
                vector.map(|component| -component).to_vec(),
                "edge {edge_index} from the other endpoint must negate exactly"
            );
        }
    }

    #[test]
    fn parallel_edges_carry_distinct_displacements() {
        // L = 2 ring: both physical bonds join sites 0 and 1, but one stays
        // in-cell while the other wraps. Endpoint geometry cannot tell them
        // apart; only the per-edge table can.
        let (lattice, embedding) = build_chain_with_embedding(2, true).expect("ring embedding");
        assert_eq!(lattice.n_edges(), 2);
        assert_eq!(displacement(&lattice, &embedding, 0, 0), vec![0]);
        assert_eq!(displacement(&lattice, &embedding, 1, 1), vec![1]);
        assert_eq!(displacement(&lattice, &embedding, 0, 1), vec![0]);
        assert_eq!(displacement(&lattice, &embedding, 1, 0), vec![-1]);

        // 2x2 square: pair {0, 1} carries d(0->1) = (0,0) and (-1,0).
        let (lattice, embedding) =
            build_square_with_embedding(2, 2, true).expect("torus embedding");
        assert_eq!(displacement(&lattice, &embedding, 0, 0), vec![0, 0]);
        assert_eq!(displacement(&lattice, &embedding, 2, 0), vec![-1, 0]);

        // 2x2 triangular: the two diagonals between sites 1 and 2 wrap
        // different axes — (1,0) versus (0,-1) from site 1.
        let (lattice, embedding) = build_triangular_with_embedding(2, 2).expect("torus embedding");
        assert_eq!(displacement(&lattice, &embedding, 5, 1), vec![1, 0]);
        assert_eq!(displacement(&lattice, &embedding, 8, 2), vec![0, 1]);
        assert_eq!(displacement(&lattice, &embedding, 8, 1), vec![0, -1]);
    }

    #[test]
    fn torus_displacements_follow_builder_geometry() {
        let cases: [Vec<usize>; 8] = [
            vec![2],
            vec![3],
            vec![5],
            vec![2, 2],
            vec![3, 3],
            vec![4, 3],
            vec![2, 3, 4],
            vec![3, 1, 2],
        ];
        for dims in cases {
            let n_dims = dims.len();
            let bond_types = bond_types_for(n_dims);
            let (lattice, embedding) =
                build_hypercubic_with_embedding(&dims, &bond_types, true).expect("torus embedding");
            assert_eq!(embedding.periodic_dimensions(), n_dims);
            let strides = hypercubic_strides(&dims);
            for (edge_index, bond) in lattice.edges.iter().enumerate() {
                let coords = site_coords(&dims, bond.source);
                // Independently reconstruct which axis this bond advances and
                // whether it wraps, from the builder's coordinate convention.
                let mut matches = Vec::new();
                for axis in 0..n_dims {
                    if dims[axis] == 1 {
                        continue; // no bonds along a length-1 axis
                    }
                    let next = (coords[axis] + 1) % dims[axis];
                    let target = bond.source - coords[axis] * strides[axis] + next * strides[axis];
                    if target == bond.target {
                        matches.push(axis);
                    }
                }
                assert_eq!(matches.len(), 1, "one advancing axis per hypercubic bond");
                let axis = matches[0];
                let mut expected = vec![0i32; n_dims];
                expected[axis] = i32::from(coords[axis] + 1 == dims[axis]);

                let edge = lattice.edge_id(edge_index).expect("dense edge index");
                let source = lattice.vertex_id(bond.source).expect("dense vertex index");
                let target = lattice.vertex_id(bond.target).expect("dense vertex index");
                let forward = embedding
                    .edge_displacement(edge, source)
                    .expect("source endpoint");
                let backward = embedding
                    .edge_displacement(edge, target)
                    .expect("target endpoint");
                for (periodic_axis, &component) in expected.iter().enumerate() {
                    assert_eq!(
                        forward.axis(periodic_axis).expect("axis"),
                        component,
                        "dims {dims:?} edge {edge_index}"
                    );
                    assert_eq!(
                        backward.axis(periodic_axis).expect("axis"),
                        -component,
                        "reverse direction negates, dims {dims:?} edge {edge_index}"
                    );
                }
                assert_eq!(
                    forward.is_zero(),
                    expected.iter().all(|&component| component == 0)
                );
            }
        }
    }

    #[test]
    fn open_boundaries_have_no_embedding() {
        assert_eq!(
            build_chain_with_embedding(4, false).unwrap_err(),
            EmbeddingError::OpenBoundaries
        );
        assert_eq!(
            build_square_with_embedding(2, 3, false).unwrap_err(),
            EmbeddingError::OpenBoundaries
        );
        assert_eq!(
            build_hypercubic_with_embedding(&[2, 2, 2], &bond_types_for(3), false).unwrap_err(),
            EmbeddingError::OpenBoundaries
        );
    }

    #[test]
    fn square_plaquette_displacements_sum_to_zero_on_3x3() {
        let (lattice, embedding) =
            build_square_with_embedding(3, 3, true).expect("torus embedding");
        let index = |x: usize, y: usize| 3 * y + x;
        for y in 0..3 {
            for x in 0..3 {
                let plaquette = [
                    index(x, y),
                    index((x + 1) % 3, y),
                    index((x + 1) % 3, (y + 1) % 3),
                    index(x, (y + 1) % 3),
                ];
                assert_eq!(
                    cycle_sum(&lattice, &embedding, &plaquette),
                    vec![0, 0],
                    "plaquette at ({x},{y}) is contractible"
                );
            }
        }
    }

    #[test]
    fn triangular_elementary_triangles_sum_to_zero_on_3x3() {
        let (lattice, embedding) = build_triangular_with_embedding(3, 3).expect("torus embedding");
        let index = |x: usize, y: usize| 3 * y + x;
        for y in 0..3 {
            for x in 0..3 {
                let a = index(x, y);
                let b = index((x + 1) % 3, y);
                let c = index((x + 1) % 3, (y + 1) % 3);
                let d = index(x, (y + 1) % 3);
                let up = [a, b, c]; // TriX forward, TriY forward, TriDiag reversed
                let down = [a, d, c]; // TriY forward, TriDiag forward, TriX reversed
                assert_eq!(
                    cycle_sum(&lattice, &embedding, &up),
                    vec![0, 0],
                    "up triangle at ({x},{y})"
                );
                assert_eq!(
                    cycle_sum(&lattice, &embedding, &down),
                    vec![0, 0],
                    "down triangle at ({x},{y})"
                );
            }
        }
    }

    #[test]
    fn axis_rings_wind_exactly_once() {
        // Walking one full ring along each axis of a torus must accumulate
        // exactly one fundamental-cell crossing on that axis and none on the
        // others: the winding number of the fundamental cycle is one.
        fn walk_axis(
            lattice: &CsrLattice,
            embedding: &LatticeEmbedding,
            dims: &[usize],
            axis: usize,
        ) -> Vec<i32> {
            let strides = hypercubic_strides(dims);
            let mut total = vec![0i32; dims.len()];
            let mut site = 0usize;
            for _ in 0..dims[axis] {
                let coords = site_coords(dims, site);
                let next = (coords[axis] + 1) % dims[axis];
                let target = site - coords[axis] * strides[axis] + next * strides[axis];
                // The ring bond emitted FROM this site; at L = 2 a parallel
                // bond joins the same pair in the opposite direction.
                let edge_index = lattice
                    .edges
                    .iter()
                    .position(|bond| bond.source == site && bond.target == target)
                    .expect("ring bond emitted from this site");
                for (sum, component) in total
                    .iter_mut()
                    .zip(displacement(lattice, embedding, edge_index, site))
                {
                    *sum += component;
                }
                site = target;
            }
            assert_eq!(site, 0, "ring walk returns to the origin");
            total
        }

        let (chain, chain_embedding) = build_chain_with_embedding(5, true).expect("ring embedding");
        assert_eq!(walk_axis(&chain, &chain_embedding, &[5], 0), vec![1]);

        for n in [2usize, 3, 4] {
            let (square, square_embedding) =
                build_square_with_embedding(n, n, true).expect("torus embedding");
            assert_eq!(
                walk_axis(&square, &square_embedding, &[n, n], 0),
                vec![1, 0]
            );
            assert_eq!(
                walk_axis(&square, &square_embedding, &[n, n], 1),
                vec![0, 1]
            );
        }

        let dims = [2usize, 3, 2];
        let (cube, cube_embedding) =
            build_hypercubic_with_embedding(&dims, &bond_types_for(3), true)
                .expect("torus embedding");
        assert_eq!(walk_axis(&cube, &cube_embedding, &dims, 0), vec![1, 0, 0]);
        assert_eq!(walk_axis(&cube, &cube_embedding, &dims, 1), vec![0, 1, 0]);
        assert_eq!(walk_axis(&cube, &cube_embedding, &dims, 2), vec![0, 0, 1]);
    }

    #[test]
    fn embedding_builders_reproduce_the_plain_lattices_full_structure() {
        let (chain, _) = build_chain_with_embedding(4, true).expect("ring embedding");
        let plain_chain = build_chain(4, true);
        assert_eq!(chain.edges, plain_chain.edges);
        assert_eq!(chain.n_sites, plain_chain.n_sites);
        assert_eq!(chain.n_bonds, plain_chain.n_bonds);
        assert_eq!(chain.offsets, plain_chain.offsets);
        assert_eq!(chain.neighbors, plain_chain.neighbors);
        assert_eq!(chain.edge_ids, plain_chain.edge_ids);

        let (square, _) = build_square_with_embedding(3, 4, true).expect("torus embedding");
        let plain_square = build_square(3, 4, true);
        assert_eq!(square.edges, plain_square.edges);
        assert_eq!(square.n_sites, plain_square.n_sites);
        assert_eq!(square.offsets, plain_square.offsets);
        assert_eq!(square.neighbors, plain_square.neighbors);
        assert_eq!(square.edge_ids, plain_square.edge_ids);

        let (triangular, _) = build_triangular_with_embedding(3, 3).expect("torus embedding");
        let plain_triangular = build_triangular(3, 3);
        assert_eq!(triangular.edges, plain_triangular.edges);
        assert_eq!(triangular.n_sites, plain_triangular.n_sites);
        assert_eq!(triangular.offsets, plain_triangular.offsets);
        assert_eq!(triangular.neighbors, plain_triangular.neighbors);
        assert_eq!(triangular.edge_ids, plain_triangular.edge_ids);

        let dims = [2usize, 3, 2];
        let bond_types = bond_types_for(3);
        let (cube, _) =
            build_hypercubic_with_embedding(&dims, &bond_types, true).expect("torus embedding");
        let plain_cube = build_hypercubic(&dims, &bond_types, true);
        assert_eq!(cube.edges, plain_cube.edges);
        assert_eq!(cube.n_sites, plain_cube.n_sites);
        assert_eq!(cube.offsets, plain_cube.offsets);
        assert_eq!(cube.neighbors, plain_cube.neighbors);
        assert_eq!(cube.edge_ids, plain_cube.edge_ids);
    }
}
