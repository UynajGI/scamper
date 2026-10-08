use super::{EdgeId, VertexId};

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

// D1 will decide the public directed-incidence shape. Keeping these capabilities
// crate-private prevents the F1 sketch from becoming a compatibility promise.
#[allow(dead_code)]
pub(crate) trait DirectedOutGraphView: GraphView {
    fn out_neighbors(&self, vertex: VertexId) -> impl Iterator<Item = VertexId> + '_;
}

#[allow(dead_code)]
pub(crate) trait DirectedInGraphView: GraphView {
    fn in_neighbors(&self, vertex: VertexId) -> impl Iterator<Item = VertexId> + '_;
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

    impl DirectedOutGraphView for OutFixture {
        fn out_neighbors(&self, vertex: VertexId) -> impl Iterator<Item = VertexId> + '_ {
            self.0[vertex.index()].iter().copied()
        }
    }

    struct InFixture([[VertexId; 1]; 2]);

    impl GraphView for InFixture {
        fn vertex_count(&self) -> usize {
            2
        }
    }

    impl DirectedInGraphView for InFixture {
        fn in_neighbors(&self, vertex: VertexId) -> impl Iterator<Item = VertexId> + '_ {
            self.0[vertex.index()].iter().copied()
        }
    }

    fn first_out<G: DirectedOutGraphView>(graph: &G, vertex: VertexId) -> VertexId {
        graph.out_neighbors(vertex).next().expect("fixture arc")
    }

    fn first_in<G: DirectedInGraphView>(graph: &G, vertex: VertexId) -> VertexId {
        graph.in_neighbors(vertex).next().expect("fixture arc")
    }

    #[test]
    fn directed_in_and_out_sketches_are_separate_capabilities() {
        let zero = VertexId::new(0, 2).expect("fixture ID");
        let one = VertexId::new(1, 2).expect("fixture ID");
        assert_eq!(first_out(&OutFixture([[one], [one]]), zero), one);
        assert_eq!(first_in(&InFixture([[zero], [zero]]), one), zero);
    }
}
