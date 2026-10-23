//! External bond-configuration analysis (roadmap §15 R1).
//!
//! This module is the **decoupled analysis contract** of the FK/random-cluster
//! track: it analyzes a *caller-supplied* bond configuration through the
//! shared F3 analyzer capabilities ([`crate::analyze`],
//! [`crate::analyze_with_labels`], [`crate::analyze_winding`]) and is
//! deliberately ignorant of the sampler that produced the bonds. Independent
//! bond Bernoulli, a future random-cluster q = 1 reduction (R2), or the actual
//! FK bond draw of an Edwards-Sokal/Swendsen-Wang adapter (R3) all arrive as
//! the same thing: one `bool` per physical edge of an
//! [`UndirectedGraphView`]. Nothing here imports, constrains, or couples to
//! any sampler; the worm high-temperature even subgraph stays out by design
//! (it is not an FK configuration).
//!
//! Semantics follow F3 exactly: every vertex is active (bond percolation —
//! the configuration is an edge *subset*), an edge is active iff its mask
//! entry is open, and an open self-loop is a union no-op while parallel edges
//! stay distinct physical edges with independent mask entries.

use core::fmt;

use super::{
    analyze, AnalysisError, AnalysisResult, BoundaryQuery, ComponentWorkspace, EdgeActivity,
    VertexActivity,
};
use crate::{EdgeId, UndirectedGraphView};

/// A caller-supplied bond mask, ready for the shared analyzers.
///
/// Construction pins the dimensions of the reference graph; the analyzer
/// re-validates those dimensions against the graph it is actually given, so a
/// configuration built for one graph and analyzed against another is a typed
/// rejection ([`AnalysisError::EdgeActivityLength`] or
/// [`AnalysisError::VertexActivityLength`]), never silent misuse. The mask is
/// borrowed, not copied: this is the zero-copy, zero-allocation path for
/// researchers who already hold the bond draw.
#[derive(Clone, Copy, Debug)]
pub struct ExternalBondConfiguration<'a> {
    vertex_count: usize,
    edge_count: usize,
    open: &'a [bool],
}

impl<'a> ExternalBondConfiguration<'a> {
    /// Validate `open_bonds` against `graph` and borrow it.
    ///
    /// One entry per physical edge in dense `EdgeId` order; `true` means the
    /// bond is present in the configuration. An empty mask on a graph with at
    /// least one edge is the typed [`BondConfigurationError::MaskLength`]
    /// rejection — "no open bond" is an all-`false` mask of the *correct*
    /// length, which is legal and leaves every vertex an isolate.
    pub fn new<G: UndirectedGraphView>(
        graph: &G,
        open_bonds: &'a [bool],
    ) -> Result<Self, BondConfigurationError> {
        if open_bonds.len() != graph.edge_count() {
            return Err(BondConfigurationError::MaskLength {
                expected: graph.edge_count(),
                actual: open_bonds.len(),
            });
        }
        Ok(Self {
            vertex_count: graph.vertex_count(),
            edge_count: graph.edge_count(),
            open: open_bonds,
        })
    }

    /// Bind a mask whose length already matched the graph.
    fn from_length_checked(vertex_count: usize, edge_count: usize, open: &'a [bool]) -> Self {
        Self {
            vertex_count,
            edge_count,
            open,
        }
    }

    /// Vertex count of the construction-time graph (all vertices active).
    #[inline]
    pub const fn vertex_count(&self) -> usize {
        self.vertex_count
    }

    /// Physical edge count of the construction-time graph (mask length).
    #[inline]
    pub const fn edge_count(&self) -> usize {
        self.edge_count
    }

    /// The borrowed bond mask in dense `EdgeId` order.
    #[inline]
    pub const fn mask(&self) -> &'a [bool] {
        self.open
    }

    /// Number of open bonds. `O(E)` over the mask.
    pub fn open_bond_count(&self) -> usize {
        self.open.iter().filter(|&&open| open).count()
    }
}

impl VertexActivity for ExternalBondConfiguration<'_> {
    #[inline]
    fn vertex_count(&self) -> usize {
        self.vertex_count
    }

    #[inline]
    fn vertex_active(&self, _vertex: crate::VertexId) -> bool {
        true
    }
}

impl EdgeActivity for ExternalBondConfiguration<'_> {
    #[inline]
    fn edge_count(&self) -> usize {
        self.edge_count
    }

    #[inline]
    fn edge_active(&self, edge: EdgeId) -> bool {
        *self
            .open
            .get(edge.index())
            .expect("edge ID outside bond-configuration domain")
    }
}

/// Analyze a raw bond mask in one call: `mask -> AnalysisResult`.
///
/// Sampler and analysis are deliberately decoupled: this convenience only
/// validates the mask against `graph`, wraps it in
/// [`ExternalBondConfiguration`], and runs the F3 [`analyze`] — it does not
/// know or care how the mask was drawn. Pass the same configuration as both
/// activity sources (every vertex active, edges per mask). Labels are not
/// materialized; call [`crate::analyze_with_labels`] with the
/// [`ExternalBondConfiguration`] directly when canonical labels are needed.
///
/// A mask of the wrong length surfaces as
/// [`AnalysisError::EdgeActivityLength`] — the same typed rejection the
/// analyzer produces itself.
pub fn analyze_bonds<'a, G>(
    graph: &G,
    open_bonds: &[bool],
    queries: &[BoundaryQuery],
    workspace: &'a mut ComponentWorkspace,
) -> Result<AnalysisResult<'a>, AnalysisError>
where
    G: UndirectedGraphView,
{
    if open_bonds.len() != graph.edge_count() {
        return Err(AnalysisError::EdgeActivityLength {
            expected: graph.edge_count(),
            actual: open_bonds.len(),
        });
    }
    let configuration = ExternalBondConfiguration::from_length_checked(
        graph.vertex_count(),
        graph.edge_count(),
        open_bonds,
    );
    analyze(graph, &configuration, &configuration, queries, workspace)
}

/// Owned record of *which* physical edges are open in one bond draw.
///
/// Storage is a **sorted, deduplicated `Vec<EdgeId>`**, chosen over a bitset
/// because the producer this type serves — a future R3 Edwards-Sokal/SW
/// adapter exposing its actual FK draw — naturally emits open edges as it
/// draws them, membership stays `O(log k)` via [`Self::contains`], and sparse
/// draws at small occupation cost `O(k)` instead of `O(E)`. Dense-mask views
/// for the analyzers are produced by [`Self::write_mask`] (reusable buffer)
/// or [`Self::to_mask`] (fresh allocation), then wrapped in
/// [`ExternalBondConfiguration`].
///
/// This is a record, not a sampler: it stores no probability, no RNG state,
/// and no coupling to any update kernel.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct BondSample {
    edge_count: usize,
    open: Vec<EdgeId>,
}

impl BondSample {
    /// Record the open edges of a draw on `graph`.
    ///
    /// Every ID must come from `graph` (a dense index outside its edge count
    /// is the typed [`BondConfigurationError::EdgeIndex`] rejection). Input
    /// order is irrelevant; duplicates collapse — the record is the *set* of
    /// open bonds, and a duplicated ID carries no extra information for
    /// connectivity analysis.
    pub fn from_open_edges<G: UndirectedGraphView>(
        graph: &G,
        open_edges: impl IntoIterator<Item = EdgeId>,
    ) -> Result<Self, BondConfigurationError> {
        let edge_count = graph.edge_count();
        let mut open: Vec<EdgeId> = open_edges.into_iter().collect();
        for edge in &open {
            let index = edge.index();
            if index >= edge_count {
                return Err(BondConfigurationError::EdgeIndex { index, edge_count });
            }
        }
        open.sort_unstable();
        open.dedup();
        Ok(Self { edge_count, open })
    }

    /// Record the open edges of a dense bond mask on `graph`.
    ///
    /// The mask length must equal `graph.edge_count()`
    /// ([`BondConfigurationError::MaskLength`] otherwise). Together with
    /// [`Self::to_mask`]/[`Self::write_mask`] this round-trips a
    /// configuration between both representations.
    pub fn from_mask<G: UndirectedGraphView>(
        graph: &G,
        open_bonds: &[bool],
    ) -> Result<Self, BondConfigurationError> {
        let configuration = ExternalBondConfiguration::new(graph, open_bonds)?;
        let open = (0..configuration.edge_count())
            .filter_map(|index| {
                let edge = graph.edge_id(index).expect("dense physical edge index");
                configuration.edge_active(edge).then_some(edge)
            })
            .collect();
        Ok(Self {
            edge_count: configuration.edge_count(),
            open,
        })
    }

    /// Physical edge count of the reference graph (the mask-view length).
    #[inline]
    pub const fn edge_count(&self) -> usize {
        self.edge_count
    }

    /// Number of open bonds.
    #[inline]
    pub fn open_count(&self) -> usize {
        self.open.len()
    }

    /// The open edges, sorted ascending by dense index, without duplicates.
    #[inline]
    pub fn open_edge_ids(&self) -> &[EdgeId] {
        &self.open
    }

    /// Whether `edge` is open. `O(log k)` binary search.
    pub fn contains(&self, edge: EdgeId) -> bool {
        self.open.binary_search(&edge).is_ok()
    }

    /// Write the dense mask view into a caller-owned buffer of exactly
    /// [`Self::edge_count`] entries (closed edges are `false`; entries are
    /// overwritten, so the buffer is safely reusable across draws).
    pub fn write_mask(&self, open_bonds: &mut [bool]) -> Result<(), BondConfigurationError> {
        if open_bonds.len() != self.edge_count {
            return Err(BondConfigurationError::MaskLength {
                expected: self.edge_count,
                actual: open_bonds.len(),
            });
        }
        open_bonds.fill(false);
        for edge in &self.open {
            open_bonds[edge.index()] = true;
        }
        Ok(())
    }

    /// Freshly allocated dense mask view of this sample.
    pub fn to_mask(&self) -> Vec<bool> {
        let mut open_bonds = vec![false; self.edge_count];
        for edge in &self.open {
            open_bonds[edge.index()] = true;
        }
        open_bonds
    }
}

/// Validation failure while binding a bond configuration to a graph.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BondConfigurationError {
    /// A bond mask (or mask-view buffer) length differs from the graph's
    /// physical edge count.
    MaskLength { expected: usize, actual: usize },
    /// An open-edge record named a dense edge index outside the graph.
    EdgeIndex { index: usize, edge_count: usize },
}

impl fmt::Display for BondConfigurationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MaskLength { expected, actual } => write!(
                formatter,
                "bond mask has {actual} entries, expected {expected} physical edges"
            ),
            Self::EdgeIndex { index, edge_count } => write!(
                formatter,
                "open-edge index {index} outside graph edge count {edge_count}"
            ),
        }
    }
}

impl std::error::Error for BondConfigurationError {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{build_chain, GraphView};

    #[test]
    fn all_active_vertices_and_masked_edges_feed_the_capabilities() {
        let graph = build_chain(4, false);
        let mask = [true, false, true];
        let configuration = ExternalBondConfiguration::new(&graph, &mask).unwrap();
        assert_eq!(configuration.vertex_count(), 4);
        assert_eq!(configuration.edge_count(), 3);
        assert_eq!(configuration.mask(), &mask);
        assert_eq!(configuration.open_bond_count(), 2);
        for vertex in graph.vertex_ids() {
            assert!(configuration.vertex_active(vertex));
        }
        for (index, &open) in mask.iter().enumerate() {
            let edge = graph.edge_id(index).unwrap();
            assert_eq!(configuration.edge_active(edge), open);
        }
    }

    #[test]
    fn mask_length_mismatch_is_a_typed_rejection_from_both_paths() {
        let graph = build_chain(4, false);
        assert_eq!(
            ExternalBondConfiguration::new(&graph, &[true; 2]).unwrap_err(),
            BondConfigurationError::MaskLength {
                expected: 3,
                actual: 2
            }
        );
        // Empty mask on a non-empty graph is the same rejection; the legal
        // "nothing open" configuration is an all-false mask of full length.
        assert_eq!(
            ExternalBondConfiguration::new(&graph, &[]).unwrap_err(),
            BondConfigurationError::MaskLength {
                expected: 3,
                actual: 0
            }
        );
        assert!(matches!(
            analyze_bonds(&graph, &[true; 4], &[], &mut ComponentWorkspace::new()).unwrap_err(),
            AnalysisError::EdgeActivityLength {
                expected: 3,
                actual: 4
            }
        ));
        // The rejection is the analyzer's own typed domain error.
        assert!(
            analyze_bonds(&graph, &[true; 4], &[], &mut ComponentWorkspace::new())
                .unwrap_err()
                .to_string()
                .contains("edge activity length")
        );
    }

    #[test]
    fn all_false_mask_leaves_every_vertex_an_isolate() {
        let graph = build_chain(4, false);
        let mut workspace = ComponentWorkspace::new();
        let result = analyze_bonds(&graph, &[false, false, false], &[], &mut workspace).unwrap();
        let summary = result.summary();
        assert_eq!(summary.active_vertex_count, 4);
        assert_eq!(summary.active_edge_count, 0);
        assert_eq!(summary.component_count, 4);
        assert_eq!(summary.largest_component_size, 1);
        assert_eq!(summary.raw_second_moment, 4);
    }

    #[test]
    fn open_self_loop_is_a_counted_union_no_op() {
        // Self-loop at 0, parallel edges 1-2, plus one longer path piece.
        let graph = crate::CsrLattice::from_edges(
            4,
            vec![
                crate::Bond::new(0, 0, crate::BondType::Generic, 1.0),
                crate::Bond::new(1, 2, crate::BondType::Generic, 1.0),
                crate::Bond::new(1, 2, crate::BondType::Generic, 1.0),
                crate::Bond::new(2, 3, crate::BondType::Generic, 1.0),
            ],
        );
        let mut workspace = ComponentWorkspace::new();
        let result = analyze_bonds(&graph, &[true, true, true, true], &[], &mut workspace).unwrap();
        let summary = result.summary();
        // All four edges count as active (F3), the self-loop unions nothing,
        // and the parallel pair collapses to one connection.
        assert_eq!(summary.active_edge_count, 4);
        assert_eq!(summary.component_count, 2);
        assert_eq!(summary.largest_component_size, 3);
        assert_eq!(summary.raw_second_moment, 1 + 9);
    }

    #[test]
    fn bond_sample_round_trips_both_representations() {
        let graph = build_chain(5, false);
        let mask = [true, false, false, true];
        let sample = BondSample::from_mask(&graph, &mask).unwrap();
        assert_eq!(sample.edge_count(), 4);
        assert_eq!(sample.open_count(), 2);
        assert_eq!(
            sample.open_edge_ids(),
            &[graph.edge_id(0).unwrap(), graph.edge_id(3).unwrap()]
        );
        assert_eq!(sample.to_mask(), mask);

        let mut buffer = vec![true; 4];
        sample.write_mask(&mut buffer).unwrap();
        assert_eq!(buffer, mask);
        assert_eq!(
            sample.write_mask(&mut [false; 3]).unwrap_err(),
            BondConfigurationError::MaskLength {
                expected: 4,
                actual: 3
            }
        );
    }

    #[test]
    fn bond_sample_sorts_dedups_and_rejects_foreign_indices() {
        let graph = build_chain(4, false);
        let edges = [
            graph.edge_id(2).unwrap(),
            graph.edge_id(0).unwrap(),
            graph.edge_id(2).unwrap(),
        ];
        let sample = BondSample::from_open_edges(&graph, edges).unwrap();
        assert_eq!(
            sample.open_edge_ids(),
            &[graph.edge_id(0).unwrap(), graph.edge_id(2).unwrap()],
            "duplicates collapse and IDs sort ascending"
        );
        assert!(sample.contains(graph.edge_id(2).unwrap()));
        assert!(!sample.contains(graph.edge_id(1).unwrap()));

        let foreign = EdgeId::from_valid_index(3);
        assert_eq!(
            BondSample::from_open_edges(&graph, [foreign]).unwrap_err(),
            BondConfigurationError::EdgeIndex {
                index: 3,
                edge_count: 3
            }
        );
    }
}
