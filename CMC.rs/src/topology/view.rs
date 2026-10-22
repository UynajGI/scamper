use super::{ArcId, EdgeId, VertexId};

/// One undirected incidence from the queried vertex.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Incidence {
    pub neighbor: VertexId,
    pub edge: EdgeId,
}

/// Minimal read-only graph capability with dense vertex IDs.
///
/// IDs produced by a view are valid only for that view while its topology is
/// unchanged. IDs carry no graph provenance. Implementations with externally
/// mutable topology must document and uphold validity before infallible access;
/// violating this contract may panic but cannot cause memory unsafety.
pub trait GraphView {
    fn vertex_count(&self) -> usize;

    /// Returns the ID at `index`, or `None` outside this view's dense range.
    #[inline]
    fn vertex_id(&self, index: usize) -> Option<VertexId> {
        VertexId::new(index, self.vertex_count())
    }

    /// Convenience iterator over this view's dense vertex IDs in ascending order.
    ///
    /// This default is not a performance obligation for other implementations.
    #[inline]
    fn vertex_ids(&self) -> impl Iterator<Item = VertexId> {
        (0..self.vertex_count()).map(VertexId::from_valid_index)
    }
}

/// Read-only undirected multigraph capability.
///
/// Physical edge IDs are stable dense indices. Every physical edge has exactly
/// two incidences, including a self-loop; parallel edges retain distinct IDs.
/// IDs must come from this unchanged, validated view. Passing an ID from another
/// graph or using a structurally invalid implementation may panic at infallible
/// accessors. `CsrLattice` has public fields, so callers that mutate them must
/// restore its documented invariants and call `CsrLattice::validate` before
/// using this trait implementation again.
pub trait UndirectedGraphView: GraphView {
    fn edge_count(&self) -> usize;

    /// Returns the ID at `index`, or `None` outside this view's dense range.
    #[inline]
    fn edge_id(&self, index: usize) -> Option<EdgeId> {
        EdgeId::new(index, self.edge_count())
    }

    fn edge_endpoints(&self, edge: EdgeId) -> [VertexId; 2];

    fn incidences(&self, vertex: VertexId) -> impl Iterator<Item = Incidence> + '_;
}

/// One directed incidence (arc) leaving or entering the queried vertex.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ArcIncidence {
    /// The vertex reached by following the arc in the queried direction.
    pub neighbor: VertexId,
    /// The arc traversed.
    pub arc: ArcId,
}

/// Read-only directed multigraph capability (D1).
///
/// Arc IDs are stable dense indices; arc `a` runs
/// `arc_endpoints(a)[0] -> arc_endpoints(a)[1]`. Parallel arcs (same
/// source and target) keep distinct IDs; self-loops are legal and count as
/// one out-incidence and one in-incidence at their vertex. Directedness is
/// a capability, never a storage flag: out- and in-adjacency are separate
/// traits so a caller cannot silently treat an arc list as symmetric.
pub trait DirectedGraphView: GraphView {
    fn arc_count(&self) -> usize;

    /// Returns the ID at `index`, or `None` outside this view's dense range.
    #[inline]
    fn arc_id(&self, index: usize) -> Option<ArcId> {
        ArcId::new(index, self.arc_count())
    }

    /// `[source, target]` of the arc; direction is source -> target.
    fn arc_endpoints(&self, arc: ArcId) -> [VertexId; 2];

    /// Convenience iterator over this view's dense arc IDs in ascending
    /// order. Like [`GraphView::vertex_ids`], this default is a convenience,
    /// not a performance obligation.
    #[inline]
    fn arc_ids(&self) -> impl Iterator<Item = ArcId> {
        (0..self.arc_count()).map(ArcId::from_valid_index)
    }
}

/// Outgoing-adjacency capability: iterate the arcs leaving a vertex.
pub trait DirectedOutGraphView: DirectedGraphView {
    fn out_arcs(&self, vertex: VertexId) -> impl Iterator<Item = ArcIncidence> + '_;
}

/// Incoming-adjacency capability: iterate the arcs entering a vertex,
/// backed by a transpose representation.
pub trait DirectedInGraphView: DirectedGraphView {
    fn in_arcs(&self, vertex: VertexId) -> impl Iterator<Item = ArcIncidence> + '_;
}

#[cfg(test)]
mod tests {
    use super::*;

    struct OutFixture([[VertexId; 1]; 2]);

    impl GraphView for OutFixture {
        fn vertex_count(&self) -> usize {
            2
        }
    }

    impl DirectedGraphView for OutFixture {
        fn arc_count(&self) -> usize {
            1
        }

        fn arc_endpoints(&self, _arc: ArcId) -> [VertexId; 2] {
            [VertexId::from_valid_index(0), VertexId::from_valid_index(1)]
        }
    }

    impl DirectedOutGraphView for OutFixture {
        fn out_arcs(&self, vertex: VertexId) -> impl Iterator<Item = ArcIncidence> + '_ {
            let neighbor = self.0[vertex.index()][0];
            [ArcIncidence {
                neighbor,
                arc: ArcId::from_valid_index(0),
            }]
            .into_iter()
        }
    }

    #[test]
    fn directed_out_is_a_capability_with_arc_identity() {
        let graph = OutFixture([
            [VertexId::from_valid_index(1)],
            [VertexId::from_valid_index(1)],
        ]);
        let first = graph
            .out_arcs(VertexId::from_valid_index(0))
            .next()
            .expect("fixture arc");
        assert_eq!(first.arc.index(), 0);
        assert_eq!(first.neighbor.index(), 1);
        assert_eq!(graph.arc_endpoints(first.arc)[1], first.neighbor);
    }
}
