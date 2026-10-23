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
/// Site numbering is a public contract (the embedding builders and their
/// tests rely on it): site `(x, y)` carries index `y * lx + x`, with `x`
/// even/odd selecting the brick-wall sublattice. Every site has three
/// neighbours, listed here in family order:
///
/// - right `(x + 1 mod lx, y)` — horizontal, [`BondType::HoneyX`];
/// - left `(x + lx - 1 mod lx, y)` — horizontal, [`BondType::HoneyX`];
/// - matched: `(x + 1 mod lx, y + 1 mod ly)` for even `x`,
///   `(x + lx - 1 mod lx, y + ly - 1 mod ly)` for odd `x` — the inter-row
///   bond, [`BondType::HoneyY`].
///
/// On `lx = 2` the right and left steps of one site reach the same
/// neighbour, so every row carries a doubled horizontal bond pair; the two
/// physical bonds keep distinct cell displacements
/// (see [`build_honeycomb_with_embedding`]).
pub fn build_honeycomb(lx: usize, ly: usize) -> CsrLattice {
    let (n_sites, edges, _) = honeycomb_bonds(lx, ly, false);
    CsrLattice::from_edges(n_sites, edges)
}

/// Honeycomb torus plus its periodic embedding; see
/// [`build_hypercubic_with_embedding`] for the boundary semantics.
///
/// Returns exactly the lattice of [`build_honeycomb`] (pinned bit-identical
/// by test) with a [`LatticeEmbedding`] whose per-bond fundamental-cell
/// displacements were recorded at emission. A step moving `+1` along an axis
/// from coordinate `c` crosses one fundamental-cell boundary iff `c + 1` is
/// the axis length (displacement `+1`), a `-1` step crosses one iff `c == 0`
/// (displacement `-1`), and otherwise the component is `0`. Within a
/// doubled horizontal pair (`lx = 2`) the non-wrapping bond `(0, 0)` is
/// emitted before the wrapping one `(-1, 0)`, which fixes the
/// displacement-table row order inside parallel groups.
pub fn build_honeycomb_with_embedding(
    lx: usize,
    ly: usize,
) -> Result<(CsrLattice, LatticeEmbedding), EmbeddingError> {
    let (n_sites, edges, displacements) = honeycomb_bonds(lx, ly, true);
    let lattice = CsrLattice::from_edges(n_sites, edges);
    let embedding = LatticeEmbedding::try_from_edge_displacements(&lattice, 2, &displacements)?;
    Ok((lattice, embedding))
}

/// Shared honeycomb bond loop; see [`hypercubic_bonds`] for the
/// displacement-tracking convention and [`build_honeycomb`] for the site
/// numbering and neighbour families.
///
/// Each physical bond is emitted exactly once, from its smaller-index
/// endpoint, visiting pairs in ascending `(source, target)` order — the order
/// the historical `from_adjacency` BTreeMap produced, so the refactored plain
/// builder stayed bit-identical and that order doubles as the
/// displacement-table row order. Within one site the bonds sort by target; a
/// stable sort keeps the family order for parallel pairs (right before left
/// at `lx = 2`).
fn honeycomb_bonds(
    lx: usize,
    ly: usize,
    track_displacements: bool,
) -> (usize, Vec<Bond>, Vec<i32>) {
    assert!(lx >= 2, "honeycomb lattice needs Lx >= 2");
    assert!(ly >= 2, "honeycomb lattice needs Ly >= 2");
    assert!(lx.is_multiple_of(2), "honeycomb lattice needs even Lx");
    let index = |x: usize, y: usize| y * lx + x;

    let mut edges = Vec::with_capacity(3 * lx * ly / 2);
    let mut displacements = Vec::new();
    for y in 0..ly {
        for x in 0..lx {
            let site = index(x, y);
            let mut candidates = [
                (
                    index((x + 1) % lx, y),
                    [i32::from(x + 1 == lx), 0],
                    BondType::HoneyX,
                ),
                (
                    index((x + lx - 1) % lx, y),
                    [-i32::from(x == 0), 0],
                    BondType::HoneyX,
                ),
                if x.is_multiple_of(2) {
                    (
                        index((x + 1) % lx, (y + 1) % ly),
                        [i32::from(x + 1 == lx), i32::from(y + 1 == ly)],
                        BondType::HoneyY,
                    )
                } else {
                    (
                        index((x + lx - 1) % lx, (y + ly - 1) % ly),
                        [-i32::from(x == 0), -i32::from(y == 0)],
                        BondType::HoneyY,
                    )
                },
            ];
            candidates.sort_by_key(|&(target, _, _)| target);
            for (target, [dx, dy], kind) in candidates {
                if target > site {
                    edges.push(Bond::new(site, target, kind, 1.0));
                    if track_displacements {
                        displacements.extend([dx, dy]);
                    }
                }
            }
        }
    }
    (lx * ly, edges, displacements)
}

/// 2D kagome lattice (`3 × Lx × Ly` sites) with periodic boundaries.
///
/// Site numbering is a public contract (the embedding builders and their
/// tests rely on it): unit cell `(x, y)` hosts three sublattice sites, with
/// sublattice `S` at index `S + 3 * (x + y * lx)`. Writing `S(x, y)` for
/// sublattice `S` in cell `(x, y)`, the neighbours are:
///
/// - `A(x, y)`: `B(x, y)`, `C(x, y)`, `B(x - 1, y)`, `C(x, y - 1)`;
/// - `B(x, y)`: `A(x, y)`, `C(x, y)`, `A(x + 1, y)`, `C(x - 1, y)`;
/// - `C(x, y)`: `A(x, y)`, `B(x, y)`, `A(x, y + 1)`, `B(x + 1, y)`;
///
/// with cell coordinates taken mod `lx`/`ly`; every bond is
/// [`BondType::Kagome`]. No parallel bonds arise for `Lx, Ly >= 2` (each
/// site's four neighbours are distinct); the small-cell obstruction recorded
/// in W1 was reconstructability — at `Lx = 2` the left/right cross families
/// join the same sublattice pairs with swapped cells, so endpoint data alone
/// cannot attribute which pair wraps. [`build_kagome_with_embedding`]
/// records the attribution at emission instead.
pub fn build_kagome(lx: usize, ly: usize) -> CsrLattice {
    let (n_sites, edges, _) = kagome_bonds(lx, ly, false);
    CsrLattice::from_edges(n_sites, edges)
}

/// Kagome torus plus its periodic embedding; see
/// [`build_hypercubic_with_embedding`] for the boundary semantics and
/// [`build_honeycomb_with_embedding`] for the displacement convention
/// (identical here: ±1 exactly on the axis whose boundary the step crosses).
pub fn build_kagome_with_embedding(
    lx: usize,
    ly: usize,
) -> Result<(CsrLattice, LatticeEmbedding), EmbeddingError> {
    let (n_sites, edges, displacements) = kagome_bonds(lx, ly, true);
    let lattice = CsrLattice::from_edges(n_sites, edges);
    let embedding = LatticeEmbedding::try_from_edge_displacements(&lattice, 2, &displacements)?;
    Ok((lattice, embedding))
}

/// Shared kagome bond loop; see [`hypercubic_bonds`] for the
/// displacement-tracking convention, [`build_kagome`] for the site numbering
/// and neighbour families, and [`honeycomb_bonds`] for the emission-order
/// contract (identical: each bond once, from its smaller endpoint, ascending
/// `(source, target)`).
fn kagome_bonds(lx: usize, ly: usize, track_displacements: bool) -> (usize, Vec<Bond>, Vec<i32>) {
    assert!(lx >= 2 && ly >= 2, "kagome lattice needs Lx,Ly >= 2");
    let index = |sublattice: usize, x: usize, y: usize| sublattice + 3 * (x + y * lx);

    let mut edges = Vec::with_capacity(6 * lx * ly);
    let mut displacements = Vec::new();
    for y in 0..ly {
        for x in 0..lx {
            for sublattice in 0..3 {
                let site = index(sublattice, x, y);
                let mut candidates = match sublattice {
                    0 => [
                        (index(1, x, y), [0, 0]),
                        (index(2, x, y), [0, 0]),
                        (index(1, (x + lx - 1) % lx, y), [-i32::from(x == 0), 0]),
                        (index(2, x, (y + ly - 1) % ly), [0, -i32::from(y == 0)]),
                    ],
                    1 => [
                        (index(0, x, y), [0, 0]),
                        (index(2, x, y), [0, 0]),
                        (index(0, (x + 1) % lx, y), [i32::from(x + 1 == lx), 0]),
                        (index(2, (x + lx - 1) % lx, y), [-i32::from(x == 0), 0]),
                    ],
                    _ => [
                        (index(0, x, y), [0, 0]),
                        (index(1, x, y), [0, 0]),
                        (index(0, x, (y + 1) % ly), [0, i32::from(y + 1 == ly)]),
                        (index(1, (x + 1) % lx, y), [i32::from(x + 1 == lx), 0]),
                    ],
                };
                candidates.sort_by_key(|&(target, _)| target);
                for (target, [dx, dy]) in candidates {
                    if target > site {
                        edges.push(Bond::new(site, target, BondType::Kagome, 1.0));
                        if track_displacements {
                            displacements.extend([dx, dy]);
                        }
                    }
                }
            }
        }
    }
    (3 * lx * ly, edges, displacements)
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

    // ---- T3 honeycomb/kagome embeddings ----

    #[test]
    fn honeycomb_and_kagome_embedding_builders_reproduce_the_plain_lattices() {
        let honeycomb_cases: [(usize, usize); 10] = [
            (2, 2),
            (2, 3),
            (2, 4),
            (4, 2),
            (4, 3),
            (4, 4),
            (6, 2),
            (6, 3),
            (6, 5),
            (8, 4),
        ];
        for (lx, ly) in honeycomb_cases {
            let (lattice, _) = build_honeycomb_with_embedding(lx, ly).expect("honeycomb torus");
            let plain = build_honeycomb(lx, ly);
            assert_eq!(lattice.edges, plain.edges, "honeycomb {lx}x{ly} edges");
            assert_eq!(lattice.n_sites, plain.n_sites, "honeycomb {lx}x{ly} sites");
            assert_eq!(lattice.n_bonds, plain.n_bonds, "honeycomb {lx}x{ly} bonds");
            assert_eq!(
                lattice.offsets, plain.offsets,
                "honeycomb {lx}x{ly} offsets"
            );
            assert_eq!(
                lattice.neighbors, plain.neighbors,
                "honeycomb {lx}x{ly} neighbors"
            );
            assert_eq!(
                lattice.edge_ids, plain.edge_ids,
                "honeycomb {lx}x{ly} edge ids"
            );
        }

        let kagome_cases: [(usize, usize); 10] = [
            (2, 2),
            (2, 3),
            (2, 4),
            (3, 2),
            (3, 3),
            (3, 4),
            (4, 2),
            (4, 3),
            (5, 3),
            (6, 4),
        ];
        for (lx, ly) in kagome_cases {
            let (lattice, _) = build_kagome_with_embedding(lx, ly).expect("kagome torus");
            let plain = build_kagome(lx, ly);
            assert_eq!(lattice.edges, plain.edges, "kagome {lx}x{ly} edges");
            assert_eq!(lattice.n_sites, plain.n_sites, "kagome {lx}x{ly} sites");
            assert_eq!(lattice.n_bonds, plain.n_bonds, "kagome {lx}x{ly} bonds");
            assert_eq!(lattice.offsets, plain.offsets, "kagome {lx}x{ly} offsets");
            assert_eq!(
                lattice.neighbors, plain.neighbors,
                "kagome {lx}x{ly} neighbors"
            );
            assert_eq!(
                lattice.edge_ids, plain.edge_ids,
                "kagome {lx}x{ly} edge ids"
            );
        }
    }

    #[test]
    fn honeycomb_2x2_displacements_match_hand_derivation() {
        let (lattice, embedding) = build_honeycomb_with_embedding(2, 2).expect("honeycomb torus");
        assert_eq!(embedding.periodic_dimensions(), 2);
        assert_eq!(lattice.n_edges(), 6);
        // Sites are y*2+x. Both rows carry one doubled horizontal pair; the
        // matched bonds are the two inter-row diagonals. Derived on paper
        // from the documented families: right/left steps wrap only from
        // x = 0 (left) or x = lx - 1 (right, never emitted at lx = 2); the
        // odd matched step from y = 0 wraps y.
        let expected: [(usize, usize, [i32; 2]); 6] = [
            (0, 0, [0, 0]),  // 0 -> 1, right, in-cell
            (1, 0, [-1, 0]), // 0 -> 1, left, wraps x (doubled pair)
            (2, 0, [0, 0]),  // 0 -> 3, matched up-right
            (3, 1, [0, -1]), // 1 -> 2, matched down-left, wraps y
            (4, 2, [0, 0]),  // 2 -> 3, right, in-cell
            (5, 2, [-1, 0]), // 2 -> 3, left, wraps x (doubled pair)
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
    fn kagome_2x2_displacements_match_hand_derivation() {
        let (lattice, embedding) = build_kagome_with_embedding(2, 2).expect("kagome torus");
        assert_eq!(lattice.n_edges(), 24);
        // Cells (x, y) host A/B/C at s + 3(x + 2y): 0-2, 3-5, 6-8, 9-11.
        // Within-cell bonds never wrap; the wrapping bonds are exactly six:
        // the x = 0 left steps of A (A(0,y)-B(1,y)) and of B
        // (B(0,y)-C(1,y)) for y in 0..2, plus the y = 0 down steps of A
        // (A(0,0)-C(0,1), A(1,0)-C(1,1)) — matching the 2*ly + lx closed
        // form at lx = ly = 2. Hand-derived on paper from the documented
        // neighbour lists before any assertion.
        let expected: [(usize, usize, [i32; 2]); 24] = [
            (0, 0, [0, 0]),   // A(0,0)-B(0,0)
            (1, 0, [0, 0]),   // A(0,0)-C(0,0)
            (2, 0, [-1, 0]),  // A(0,0)-B(1,0), left, wraps x
            (3, 0, [0, -1]),  // A(0,0)-C(0,1), down, wraps y
            (4, 1, [0, 0]),   // B(0,0)-C(0,0)
            (5, 1, [0, 0]),   // B(0,0)-A(1,0), right, in-cell at Lx = 2
            (6, 1, [-1, 0]),  // B(0,0)-C(1,0), left, wraps x
            (7, 2, [0, 0]),   // C(0,0)-B(1,0), right, in-cell at Lx = 2
            (8, 2, [0, 0]),   // C(0,0)-A(0,1), up, in-cell at Ly = 2
            (9, 3, [0, 0]),   // A(1,0)-B(1,0)
            (10, 3, [0, 0]),  // A(1,0)-C(1,0)
            (11, 3, [0, -1]), // A(1,0)-C(1,1), down, wraps y
            (12, 4, [0, 0]),  // B(1,0)-C(1,0)
            (13, 5, [0, 0]),  // C(1,0)-A(1,1), up, in-cell
            (14, 6, [0, 0]),  // A(0,1)-B(0,1)
            (15, 6, [0, 0]),  // A(0,1)-C(0,1)
            (16, 6, [-1, 0]), // A(0,1)-B(1,1), left, wraps x
            (17, 7, [0, 0]),  // B(0,1)-C(0,1)
            (18, 7, [0, 0]),  // B(0,1)-A(1,1), right, in-cell
            (19, 7, [-1, 0]), // B(0,1)-C(1,1), left, wraps x
            (20, 8, [0, 0]),  // C(0,1)-B(1,1), right, in-cell
            (21, 9, [0, 0]),  // A(1,1)-B(1,1)
            (22, 9, [0, 0]),  // A(1,1)-C(1,1)
            (23, 10, [0, 0]), // B(1,1)-C(1,1)
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
    fn honeycomb_4x2_and_kagome_3x2_4x3_displacement_spot_checks() {
        // Honeycomb 4x2, full table (12 bonds, lx >= 4 so nothing doubles):
        // the only wrapping bonds are the left steps from x = 0 (one per row)
        // and the odd matched steps from y = 0 (one per odd column).
        let (lattice, embedding) = build_honeycomb_with_embedding(4, 2).expect("honeycomb torus");
        assert_eq!(lattice.n_edges(), 12);
        let expected: [(usize, usize, [i32; 2]); 12] = [
            (0, 0, [0, 0]),
            (1, 0, [-1, 0]), // 0 -> 3, left from x = 0
            (2, 0, [0, 0]),  // 0 -> 5, matched up-right
            (3, 1, [0, 0]),
            (4, 1, [0, -1]), // 1 -> 4, matched down-left from y = 0
            (5, 2, [0, 0]),
            (6, 2, [0, 0]),  // 2 -> 7, matched up-right
            (7, 3, [0, -1]), // 3 -> 6, matched down-left from y = 0
            (8, 4, [0, 0]),
            (9, 4, [-1, 0]), // 4 -> 7, left from x = 0
            (10, 5, [0, 0]),
            (11, 6, [0, 0]),
        ];
        for (edge_index, from, vector) in expected {
            assert_eq!(
                displacement(&lattice, &embedding, edge_index, from),
                vector.to_vec(),
                "honeycomb 4x2 edge {edge_index}"
            );
        }

        // Kagome 3x2: every wrapping bond asserted by endpoint pair (all
        // pairs are unique). Left steps from x = 0 wrap x; the y = 0 down
        // steps of A wrap y (three of them: one per column).
        let (lattice, embedding) = build_kagome_with_embedding(3, 2).expect("kagome torus");
        assert_eq!(lattice.n_edges(), 36);
        for (a, b, vector) in [
            (0usize, 7usize, [-1, 0]), // A(0,0)-B(2,0)
            (9, 16, [-1, 0]),          // A(0,1)-B(2,1)
            (1, 8, [-1, 0]),           // B(0,0)-C(2,0)
            (10, 17, [-1, 0]),         // B(0,1)-C(2,1)
            (0, 11, [0, -1]),          // A(0,0)-C(0,1)
            (3, 14, [0, -1]),          // A(1,0)-C(1,1)
            (6, 17, [0, -1]),          // A(2,0)-C(2,1)
        ] {
            let edge_index = bond_between(&lattice, a, b).expect("wrapping bond");
            assert_eq!(
                displacement(&lattice, &embedding, edge_index, a),
                vector.to_vec(),
                "kagome 3x2 bond ({a},{b})"
            );
            assert_eq!(
                displacement(&lattice, &embedding, edge_index, b),
                vector.map(|component| -component).to_vec(),
                "kagome 3x2 bond ({a},{b}) reverse"
            );
        }

        // 4x3 spot checks for both families.
        let (lattice, embedding) = build_honeycomb_with_embedding(4, 3).expect("honeycomb torus");
        for (a, b, vector) in [
            (0usize, 3usize, [-1, 0]), // left from x = 0, row 0
            (1, 8, [0, -1]),           // matched down-left from (1,0)
            (3, 10, [0, -1]),          // matched down-left from (3,0)
        ] {
            let edge_index = bond_between(&lattice, a, b).expect("honeycomb 4x3 bond");
            assert_eq!(
                displacement(&lattice, &embedding, edge_index, a),
                vector.to_vec()
            );
        }
        let (lattice, embedding) = build_kagome_with_embedding(4, 3).expect("kagome torus");
        for (a, b, vector) in [
            (0usize, 10usize, [-1, 0]), // A(0,0)-B(3,0)
            (1, 11, [-1, 0]),           // B(0,0)-C(3,0)
            (0, 26, [0, -1]),           // A(0,0)-C(0,2)
        ] {
            let edge_index = bond_between(&lattice, a, b).expect("kagome 4x3 bond");
            assert_eq!(
                displacement(&lattice, &embedding, edge_index, a),
                vector.to_vec()
            );
        }
    }

    /// One directed incidence per documented neighbour-family step of every
    /// site: `(site, target, displacement site -> target)`, with the target
    /// recovered by mixed-radix decomposition of the endpoint IDs and the
    /// displacement from wrapping cell arithmetic. This enumerates every
    /// direction independently of the builder's once-per-bond emission.
    fn honeycomb_directed_incidences(lx: usize, ly: usize) -> Vec<(usize, usize, [i32; 2])> {
        let mut incidences = Vec::with_capacity(3 * lx * ly);
        for site in 0..lx * ly {
            let (x, y) = (site % lx, site / lx);
            let families: [(isize, isize); 3] = [
                (1, 0),
                (-1, 0),
                if x.is_multiple_of(2) {
                    (1, 1)
                } else {
                    (-1, -1)
                },
            ];
            for (step_x, step_y) in families {
                let target = (y as isize + step_y).rem_euclid(ly as isize) as usize * lx
                    + (x as isize + step_x).rem_euclid(lx as isize) as usize;
                let dx = match step_x {
                    1 => i32::from(x + 1 == lx),
                    -1 => -i32::from(x == 0),
                    _ => 0,
                };
                let dy = match step_y {
                    1 => i32::from(y + 1 == ly),
                    -1 => -i32::from(y == 0),
                    _ => 0,
                };
                incidences.push((site, target, [dx, dy]));
            }
        }
        incidences
    }

    /// Same independent enumeration for kagome, from the documented
    /// neighbour lists of the three sublattices.
    fn kagome_directed_incidences(lx: usize, ly: usize) -> Vec<(usize, usize, [i32; 2])> {
        let decompose = |site: usize| (site % 3, (site / 3) % lx, site / (3 * lx));
        let mut incidences = Vec::with_capacity(4 * 3 * lx * ly);
        for site in 0..3 * lx * ly {
            let (sublattice, x, y) = decompose(site);
            // (target sublattice, step) pairs per sublattice.
            let families: [(usize, isize, isize); 4] = match sublattice {
                0 => [(1, 0, 0), (2, 0, 0), (1, -1, 0), (2, 0, -1)],
                1 => [(0, 0, 0), (2, 0, 0), (0, 1, 0), (2, -1, 0)],
                _ => [(0, 0, 0), (1, 0, 0), (0, 0, 1), (1, 1, 0)],
            };
            for (target_sublattice, step_x, step_y) in families {
                let target = target_sublattice
                    + 3 * ((y as isize + step_y).rem_euclid(ly as isize) as usize * lx
                        + (x as isize + step_x).rem_euclid(lx as isize) as usize);
                let dx = match step_x {
                    1 => i32::from(x + 1 == lx),
                    -1 => -i32::from(x == 0),
                    _ => 0,
                };
                let dy = match step_y {
                    1 => i32::from(y + 1 == ly),
                    -1 => -i32::from(y == 0),
                    _ => 0,
                };
                incidences.push((site, target, [dx, dy]));
            }
        }
        incidences
    }

    #[test]
    fn honeycomb_and_kagome_displacements_follow_builder_geometry() {
        // Per (unordered pair, displacement normalized to the smaller
        // endpoint), the independent enumeration counts one directed
        // incidence per direction: two exact negations per physical bond,
        // so every count is even and halves to the emission's table.
        fn check(
            lattice: &CsrLattice,
            embedding: &LatticeEmbedding,
            incidences: Vec<(usize, usize, [i32; 2])>,
            label: &str,
        ) {
            let mut reference: BTreeMap<((usize, usize), [i32; 2]), usize> = BTreeMap::new();
            for (site, target, vector) in incidences {
                let key = (site.min(target), site.max(target));
                let normalized = if site == key.0 {
                    vector
                } else {
                    vector.map(|component| -component)
                };
                *reference.entry((key, normalized)).or_insert(0) += 1;
            }

            let mut produced: BTreeMap<((usize, usize), [i32; 2]), usize> = BTreeMap::new();
            for edge_index in 0..lattice.n_edges() {
                let [source, target] =
                    lattice.edge_endpoints(lattice.edge_id(edge_index).expect("dense edge index"));
                assert!(
                    source.index() < target.index(),
                    "{label}: emission must key bonds by ascending endpoints"
                );
                let forward = embedding
                    .edge_displacement(
                        lattice.edge_id(edge_index).expect("dense edge index"),
                        source,
                    )
                    .expect("source endpoint");
                let backward = embedding
                    .edge_displacement(
                        lattice.edge_id(edge_index).expect("dense edge index"),
                        target,
                    )
                    .expect("target endpoint");
                let vector: Vec<i32> = forward.iter().collect();
                for (forward_component, backward_component) in forward.iter().zip(backward.iter()) {
                    assert_eq!(
                        backward_component, -forward_component,
                        "{label} edge {edge_index}: reverse direction must negate exactly"
                    );
                }
                assert_eq!(
                    forward.is_zero(),
                    vector.iter().all(|&component| component == 0),
                    "{label} edge {edge_index}: is_zero convention"
                );
                *produced
                    .entry(((source.index(), target.index()), [vector[0], vector[1]]))
                    .or_insert(0) += 1;
            }

            assert_eq!(
                produced.len(),
                reference.len(),
                "{label}: distinct (pair, displacement) entries"
            );
            for (key, &count) in &reference {
                assert_eq!(
                    count % 2,
                    0,
                    "{label} {key:?}: the two directions of a bond must negate exactly"
                );
                assert_eq!(
                    produced.get(key).copied().unwrap_or(0),
                    count / 2,
                    "{label} {key:?}: displacement multiplicity"
                );
            }
        }

        for (lx, ly) in [
            (2, 2),
            (2, 3),
            (2, 4),
            (4, 2),
            (4, 3),
            (6, 3),
            (6, 5),
            (8, 4),
        ] {
            let (lattice, embedding) =
                build_honeycomb_with_embedding(lx, ly).expect("honeycomb torus");
            assert_eq!(
                lattice.n_edges(),
                3 * lx * ly / 2,
                "honeycomb {lx}x{ly} edge count"
            );
            check(
                &lattice,
                &embedding,
                honeycomb_directed_incidences(lx, ly),
                "honeycomb",
            );
        }
        for (lx, ly) in [
            (2, 2),
            (2, 3),
            (2, 4),
            (3, 2),
            (3, 3),
            (4, 3),
            (5, 3),
            (5, 4),
            (6, 4),
        ] {
            let (lattice, embedding) = build_kagome_with_embedding(lx, ly).expect("kagome torus");
            assert_eq!(
                lattice.n_edges(),
                6 * lx * ly,
                "kagome {lx}x{ly} edge count"
            );
            check(
                &lattice,
                &embedding,
                kagome_directed_incidences(lx, ly),
                "kagome",
            );
        }
    }

    #[test]
    fn honeycomb_and_kagome_small_cell_bond_content() {
        // Honeycomb: exactly the lx = 2 rows double their horizontal pair
        // (displacements (0,0) and (-1,0)); wrapping bonds are the ly left
        // steps from x = 0 plus the lx/2 odd matched steps from y = 0.
        for (lx, ly) in [(2, 2), (2, 3), (2, 4), (4, 2), (4, 3), (6, 5)] {
            let (lattice, embedding) =
                build_honeycomb_with_embedding(lx, ly).expect("honeycomb torus");
            let mut pairs = std::collections::BTreeSet::new();
            let mut doubled = 0usize;
            let mut wrapping = 0usize;
            for edge_index in 0..lattice.n_edges() {
                let [source, target] =
                    lattice.edge_endpoints(lattice.edge_id(edge_index).expect("dense edge index"));
                let vector: Vec<i32> = embedding
                    .edge_displacement(
                        lattice.edge_id(edge_index).expect("dense edge index"),
                        source,
                    )
                    .expect("source endpoint")
                    .iter()
                    .collect();
                wrapping += usize::from(vector != vec![0, 0]);
                if !pairs.insert((source.index(), target.index())) {
                    doubled += 1;
                }
            }
            let expected_doubled = if lx == 2 { ly } else { 0 };
            assert_eq!(
                lattice.n_edges(),
                3 * lx * ly / 2,
                "honeycomb {lx}x{ly} edge count"
            );
            assert_eq!(
                lattice.n_edges() - pairs.len(),
                expected_doubled,
                "honeycomb {lx}x{ly} doubled horizontal pairs"
            );
            assert_eq!(doubled, expected_doubled);
            assert_eq!(
                wrapping,
                ly + lx / 2,
                "honeycomb {lx}x{ly} wrapping bonds (left steps + odd matched steps)"
            );
        }

        // Kagome: no parallel bonds at any size >= 2 (each site's four
        // neighbours are distinct); wrapping bonds are the A/B left steps
        // from x = 0 (2*ly) plus the A down steps from y = 0 (lx) — at
        // Lx = 2 these are exactly the collapsed cross-family bonds whose
        // wrap could not be attributed from endpoints alone.
        for (lx, ly) in [(2, 2), (2, 3), (3, 2), (3, 3), (4, 3), (5, 4)] {
            let (lattice, embedding) = build_kagome_with_embedding(lx, ly).expect("kagome torus");
            let mut pairs = std::collections::BTreeSet::new();
            let mut wrapping = 0usize;
            for edge_index in 0..lattice.n_edges() {
                let [source, target] =
                    lattice.edge_endpoints(lattice.edge_id(edge_index).expect("dense edge index"));
                let vector: Vec<i32> = embedding
                    .edge_displacement(
                        lattice.edge_id(edge_index).expect("dense edge index"),
                        source,
                    )
                    .expect("source endpoint")
                    .iter()
                    .collect();
                wrapping += usize::from(vector != vec![0, 0]);
                assert!(
                    pairs.insert((source.index(), target.index())),
                    "kagome {lx}x{ly} must not carry parallel bonds"
                );
            }
            assert_eq!(
                lattice.n_edges(),
                6 * lx * ly,
                "kagome {lx}x{ly} edge count"
            );
            assert_eq!(
                wrapping,
                2 * ly + lx,
                "kagome {lx}x{ly} wrapping bonds (A/B left steps + A down steps)"
            );
        }
    }
}
