//! Periodic embedding capability for winding-capable topology (W1).
//!
//! A periodic graph (for example a hypercubic torus) embeds into a covering
//! lattice `Z^D`: every directed incidence carries an explicit integer cell
//! displacement, and a closed walk whose accumulated displacement is nonzero
//! winds around the torus. This module stores that displacement as data. It is
//! never inferred from [`crate::BondType`] or from vertex IDs: small periodic
//! cells, parallel edges between the same endpoint pair, and non-Bravais
//! lattices all make such inference ambiguous (on a 2-wide torus the two
//! parallel edges between a pair of sites carry different displacements,
//! `0` and `-1`).
//!
//! The capability is the data foundation for the W2 winding analyzer: one pass
//! over physical edges can call
//! `embedding.edge_displacement(edge, from)` next to
//! `UndirectedGraphView::edge_endpoints` and run a union-find with integer
//! potentials. W1 deliberately provides no analysis; accumulation semantics
//! for arbitrary user-supplied magnitudes (overflow/saturation policy) are the
//! W2 analyzer's contract, while W1 stores vectors verbatim.

use core::fmt;

use super::{EdgeId, VertexId};
use crate::CsrLattice;

/// Signed view of one displacement vector, from one endpoint toward the other.
///
/// The backing storage holds the vector in one fixed direction; this view
/// negates it lazily when the displacement was queried from the other
/// endpoint. Antisymmetry is therefore structural, not a validated copy.
#[derive(Clone, Copy, Debug)]
pub struct DirectedDisplacement<'a> {
    vector: &'a [i32],
    negated: bool,
}

impl<'a> DirectedDisplacement<'a> {
    #[inline]
    pub(crate) fn new(vector: &'a [i32], negated: bool) -> DirectedDisplacement<'a> {
        DirectedDisplacement { vector, negated }
    }

    /// Iterates the signed component per periodic axis.
    #[inline]
    pub fn iter(self) -> impl Iterator<Item = i32> + 'a {
        self.vector
            .iter()
            .map(move |&component| if self.negated { -component } else { component })
    }

    /// Returns the signed component for `axis`, or `None` outside `[0, D)`.
    #[inline]
    pub fn axis(self, axis: usize) -> Option<i32> {
        self.vector
            .get(axis)
            .map(|&component| if self.negated { -component } else { component })
    }

    /// Number of periodic axes (the embedding dimension `D`).
    #[inline]
    pub fn dimensions(self) -> usize {
        self.vector.len()
    }

    /// True when no cell boundary is crossed in either direction.
    #[inline]
    pub fn is_zero(self) -> bool {
        self.vector.iter().all(|&component| component == 0)
    }
}

impl<const N: usize> PartialEq<[i32; N]> for DirectedDisplacement<'_> {
    fn eq(&self, other: &[i32; N]) -> bool {
        self.vector.len() == N && self.iter().zip(other.iter().copied()).all(|(a, b)| a == b)
    }
}

impl PartialEq<&[i32]> for DirectedDisplacement<'_> {
    fn eq(&self, other: &&[i32]) -> bool {
        self.vector.len() == other.len()
            && self.iter().zip(other.iter().copied()).all(|(a, b)| a == b)
    }
}

/// Read-only periodic-embedding capability.
///
/// Displacements are integer fundamental-cell crossings. Moving from `from`
/// across `edge` toward the other endpoint changes the covering-space cell by
/// the returned vector: `+1` on an axis means one cell boundary crossed in the
/// positive axis direction, `0` means the edge stays inside one cell.
///
/// The lookup key is the pair (physical edge, endpoint) — one directed
/// incidence — and never the endpoint pair or a bond label. The same physical
/// edge queried from either endpoint returns exact negations, and parallel
/// edges between the same endpoint pair may legitimately carry different
/// displacements. `edge_displacement` returns `None` when `from` is not an
/// endpoint of `edge`. IDs must come from the graph this embedding was
/// constructed for (same provenance contract as [`super::GraphView`]);
/// out-of-range edge IDs may panic at this infallible accessor.
///
/// `periodic_dimensions()` is at least 1 for every implementation. A graph
/// without periodicity (open boundaries) is represented by *not* having this
/// capability, not by a zero-dimensional embedding: a winding query on a
/// non-embedded graph is an explicit error in consumers, not a vacuous
/// "no winding" answer.
///
/// Self-loop edges cannot carry a directed displacement under the
/// (edge, endpoint) key — their two directed incidences would need opposite
/// displacements but share one key — so embedding construction rejects them.
///
/// [`super::BorrowedUndirectedCsr`] intentionally does not implement this
/// capability: a CSR topology carries no displacement payload, and
/// displacements cannot be reconstructed from endpoints after the fact
/// without the construction-time wrap information (the same ambiguity that
/// forbids `BondType`/vertex-ID inference). An embedding is separate data
/// that must be produced by whoever knows the geometry, beside the borrowed
/// topology, in the same way [`LatticeEmbedding`] accompanies a
/// [`CsrLattice`].
pub trait PeriodicEmbedding {
    /// Number of independent periodic directions `D` of the embedding.
    fn periodic_dimensions(&self) -> usize;

    /// Displacement from `from` toward the other endpoint of `edge`, in cell
    /// units per periodic axis. `None` when `from` is not an endpoint.
    fn edge_displacement(&self, edge: EdgeId, from: VertexId) -> Option<DirectedDisplacement<'_>>;
}

/// Owned displacement table for one [`CsrLattice`].
///
/// Built alongside a lattice by the periodic hypercubic builders
/// (`build_chain_with_embedding`, `build_square_with_embedding`,
/// `build_hypercubic_with_embedding`) or from explicitly supplied per-edge
/// vectors via [`LatticeEmbedding::try_from_edge_displacements`]. It stores one
/// `source -> target` vector per physical edge (flat, edge-major, `D` components
/// per edge) plus the endpoint pair it was validated against; the reverse
/// direction is the exact negation produced by [`DirectedDisplacement`].
///
/// The table is detached from its lattice: both must stay unchanged together,
/// like the ID provenance contract of the view capabilities.
///
/// Storage is `16 + 4 * D` bytes per physical edge; lookups are O(1) with no
/// allocation, so W1 needs no dedicated benchmark.
#[derive(Clone, Debug)]
pub struct LatticeEmbedding {
    dimensions: usize,
    endpoints: Vec<[usize; 2]>,
    displacements: Vec<i32>,
}

impl LatticeEmbedding {
    /// Builds an embedding table from explicit per-edge displacements.
    ///
    /// `displacements` is flat and edge-major: the vector of lattice edge `e`
    /// (the `Bond` order in `lattice.edges`) occupies
    /// `displacements[e * periodic_dimensions .. (e + 1) * periodic_dimensions]`
    /// and is measured from that bond's `source` toward its `target`. The
    /// lattice must satisfy its documented invariants (as produced by its
    /// constructors or revalidated with `CsrLattice::validate`).
    ///
    /// Rejections, all before any state is kept:
    /// - [`EmbeddingError::DimensionsZero`] for `periodic_dimensions == 0`
    ///   (non-periodic graphs have no embedding at all);
    /// - [`EmbeddingError::DisplacementLength`] when the table length is not
    ///   `n_edges * periodic_dimensions`;
    /// - [`EmbeddingError::SelfLoopEdge`] — self-loops have no unambiguous
    ///   directed displacement;
    /// - [`EmbeddingError::NegationOverflow`] — a component of `i32::MIN`
    ///   cannot be negated for the reverse direction;
    /// - [`EmbeddingError::ValidationCapacity`] when the table cannot reserve
    ///   its storage.
    ///
    /// This constructor does **not** verify that the supplied vectors form a
    /// consistent `Z^D` cocycle (zero sum around contractible cycles); that is
    /// a property of whoever computes the vectors, checked for the built-in
    /// builders by plaquette-enumeration tests, and exploited by the W2
    /// analyzer rather than re-validated here. It also cannot check that
    /// `periodic_dimensions` matches the lattice's intrinsic geometry — the
    /// lattice carries no dimension information, so `D` is caller data.
    pub fn try_from_edge_displacements(
        lattice: &CsrLattice,
        periodic_dimensions: usize,
        displacements: &[i32],
    ) -> Result<Self, EmbeddingError> {
        let edge_count = lattice.n_edges();
        if periodic_dimensions == 0 {
            return Err(EmbeddingError::DimensionsZero);
        }
        let expected = edge_count.saturating_mul(periodic_dimensions);
        if displacements.len() != expected {
            return Err(EmbeddingError::DisplacementLength {
                expected,
                actual: displacements.len(),
            });
        }
        for (edge, bond) in lattice.edges.iter().enumerate() {
            if bond.source == bond.target {
                return Err(EmbeddingError::SelfLoopEdge { edge });
            }
        }
        if let Some(component) = displacements
            .iter()
            .position(|&component| component == i32::MIN)
        {
            return Err(EmbeddingError::NegationOverflow {
                edge: component / periodic_dimensions,
            });
        }

        let mut endpoints = Vec::new();
        endpoints
            .try_reserve_exact(edge_count)
            .map_err(|_| EmbeddingError::ValidationCapacity { edge_count })?;
        let mut table = Vec::new();
        table
            .try_reserve_exact(displacements.len())
            .map_err(|_| EmbeddingError::ValidationCapacity { edge_count })?;
        endpoints.extend(lattice.edges.iter().map(|bond| [bond.source, bond.target]));
        table.extend_from_slice(displacements);

        Ok(LatticeEmbedding {
            dimensions: periodic_dimensions,
            endpoints,
            displacements: table,
        })
    }
}

impl PeriodicEmbedding for LatticeEmbedding {
    #[inline(always)]
    fn periodic_dimensions(&self) -> usize {
        self.dimensions
    }

    #[inline]
    fn edge_displacement(&self, edge: EdgeId, from: VertexId) -> Option<DirectedDisplacement<'_>> {
        let index = edge.index();
        let [source, target] = self.endpoints[index];
        let vector = &self.displacements[index * self.dimensions..][..self.dimensions];
        if from.index() == source {
            Some(DirectedDisplacement::new(vector, false))
        } else if from.index() == target {
            Some(DirectedDisplacement::new(vector, true))
        } else {
            None
        }
    }
}

/// Validation failure while constructing a periodic embedding.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum EmbeddingError {
    /// An embedding was requested for an open-boundary lattice. Open graphs
    /// are legitimately non-periodic; construct them with the plain builders.
    OpenBoundaries,
    /// The periodic dimension count was zero.
    DimensionsZero,
    /// The flat displacement table length was not `n_edges * dimensions`.
    DisplacementLength { expected: usize, actual: usize },
    /// A self-loop edge cannot carry a directed displacement keyed by
    /// endpoint: its two directed incidences would need opposite values.
    SelfLoopEdge { edge: usize },
    /// A displacement component of `i32::MIN` cannot be negated for the
    /// reverse direction.
    NegationOverflow { edge: usize },
    /// The embedding table could not reserve its per-edge storage.
    ValidationCapacity { edge_count: usize },
}

impl fmt::Display for EmbeddingError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::OpenBoundaries => write!(
                formatter,
                "open-boundary lattices have no periodic embedding; use the plain builders"
            ),
            Self::DimensionsZero => {
                write!(formatter, "periodic dimension count must be at least one")
            }
            Self::DisplacementLength { expected, actual } => write!(
                formatter,
                "displacement table has {actual} components, expected {expected} \
                 (n_edges * periodic_dimensions)"
            ),
            Self::SelfLoopEdge { edge } => write!(
                formatter,
                "self-loop edge {edge} cannot carry a directed displacement"
            ),
            Self::NegationOverflow { edge } => write!(
                formatter,
                "edge {edge} has a displacement component of i32::MIN, which cannot be negated"
            ),
            Self::ValidationCapacity { edge_count } => write!(
                formatter,
                "could not reserve embedding storage for {edge_count} physical edges"
            ),
        }
    }
}

impl std::error::Error for EmbeddingError {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Bond, BondType, GraphView, UndirectedGraphView};

    /// Hand-built 4-site ring with one wrap edge, bypassing the built-in
    /// builders: bulk edges 0-1, 1-2, 2-3 stay in-cell, edge 3-0 wraps +1.
    fn ring() -> (CsrLattice, LatticeEmbedding) {
        let lattice = CsrLattice::from_edges(
            4,
            vec![
                Bond::new(0, 1, BondType::Generic, 1.0),
                Bond::new(1, 2, BondType::Generic, 1.0),
                Bond::new(2, 3, BondType::Generic, 1.0),
                Bond::new(3, 0, BondType::Generic, 1.0),
            ],
        );
        let embedding = LatticeEmbedding::try_from_edge_displacements(&lattice, 1, &[0, 0, 0, 1])
            .expect("valid ring embedding");
        (lattice, embedding)
    }

    fn vertex(lattice: &CsrLattice, index: usize) -> VertexId {
        lattice.vertex_id(index).expect("dense vertex index")
    }

    #[test]
    fn ring_displacements_are_explicit_data() {
        let (lattice, embedding) = ring();
        assert_eq!(embedding.periodic_dimensions(), 1);

        for edge_index in 0..lattice.edge_count() {
            let edge = lattice.edge_id(edge_index).expect("dense edge index");
            let [source, target] = lattice.edge_endpoints(edge);
            let forward = embedding
                .edge_displacement(edge, source)
                .expect("source is an endpoint");
            let backward = embedding
                .edge_displacement(edge, target)
                .expect("target is an endpoint");
            let expected = [0, 0, 0, 1][edge_index];
            assert_eq!(forward, [expected], "edge {edge_index} source direction");
            assert_eq!(
                backward,
                [-expected],
                "edge {edge_index} target direction must negate exactly"
            );
        }

        // The full ring is one winding cycle: one net +1 crossing.
        let total: i32 = (0..lattice.edge_count())
            .map(|index| {
                let edge = lattice.edge_id(index).expect("dense edge index");
                let [source, _] = lattice.edge_endpoints(edge);
                embedding
                    .edge_displacement(edge, source)
                    .expect("endpoint")
                    .axis(0)
                    .expect("one periodic axis")
            })
            .sum();
        assert_eq!(total, 1, "the ring winds exactly once around the torus");
    }

    #[test]
    fn non_endpoint_query_returns_none() {
        let (lattice, embedding) = ring();
        let edge = lattice.edge_id(0).expect("dense edge index"); // 0-1
        assert!(
            embedding
                .edge_displacement(edge, vertex(&lattice, 2))
                .is_none(),
            "vertex 2 is not an endpoint of edge 0"
        );
    }

    #[test]
    fn generic_code_accepts_the_capability() {
        fn winding_norm<E: PeriodicEmbedding>(embedding: &E) -> usize {
            embedding.periodic_dimensions()
        }
        let (_, embedding) = ring();
        assert_eq!(winding_norm(&embedding), 1);
    }

    #[test]
    fn directed_displacement_view_helpers() {
        let (lattice, embedding) = ring();
        let edge = lattice.edge_id(3).expect("dense edge index"); // 3-0 wrap
        let forward = embedding
            .edge_displacement(edge, vertex(&lattice, 3))
            .expect("endpoint");
        assert_eq!(forward.iter().collect::<Vec<i32>>(), vec![1]);
        assert_eq!(forward.axis(0), Some(1));
        assert_eq!(forward.axis(1), None);
        assert_eq!(forward.dimensions(), 1);
        assert!(!forward.is_zero());

        let bulk = lattice.edge_id(0).expect("dense edge index");
        let zero = embedding
            .edge_displacement(bulk, vertex(&lattice, 0))
            .expect("endpoint");
        assert!(zero.is_zero());
        assert_eq!(zero, &[0][..]);
    }

    #[test]
    fn dimensions_zero_is_rejected() {
        let lattice = CsrLattice::from_edges(2, vec![Bond::new(0, 1, BondType::Generic, 1.0)]);
        assert_eq!(
            LatticeEmbedding::try_from_edge_displacements(&lattice, 0, &[]).unwrap_err(),
            EmbeddingError::DimensionsZero
        );
    }

    #[test]
    fn displacement_table_length_is_rejected() {
        let lattice = CsrLattice::from_edges(2, vec![Bond::new(0, 1, BondType::Generic, 1.0)]);
        assert_eq!(
            LatticeEmbedding::try_from_edge_displacements(&lattice, 2, &[0, 0, 0]).unwrap_err(),
            EmbeddingError::DisplacementLength {
                expected: 2,
                actual: 3
            }
        );
    }

    #[test]
    fn self_loop_edges_are_rejected() {
        let lattice = CsrLattice::from_edges(
            2,
            vec![
                Bond::new(0, 1, BondType::Generic, 1.0),
                Bond::new(1, 1, BondType::Generic, 1.0),
            ],
        );
        assert_eq!(
            LatticeEmbedding::try_from_edge_displacements(&lattice, 1, &[0, 0]).unwrap_err(),
            EmbeddingError::SelfLoopEdge { edge: 1 }
        );
    }

    #[test]
    fn non_negatable_component_is_rejected() {
        let lattice = CsrLattice::from_edges(2, vec![Bond::new(0, 1, BondType::Generic, 1.0)]);
        assert_eq!(
            LatticeEmbedding::try_from_edge_displacements(&lattice, 1, &[i32::MIN]).unwrap_err(),
            EmbeddingError::NegationOverflow { edge: 0 }
        );
    }
}
