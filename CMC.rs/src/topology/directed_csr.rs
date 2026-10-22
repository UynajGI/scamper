//! Owned directed CSR multigraph (D1).
//!
//! Arcs are stored once in a dense physical table; adjacency lives in two
//! CSR incidence structures, one outgoing (row-major by source) and one
//! transposed (row-major by target), both built by [`DirectedCsr::from_arcs`].
//! Directedness is structural: there is no symmetric fallback, and the
//! undirected [`crate::CsrLattice`] never masquerades as directed.

use super::{
    ArcId, ArcIncidence, DirectedGraphView, DirectedInGraphView, DirectedOutGraphView, GraphView,
    TopologyError, VertexId,
};

/// One physical directed arc.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Arc {
    pub source: usize,
    pub target: usize,
}

impl Arc {
    pub const fn new(source: usize, target: usize) -> Self {
        Self { source, target }
    }
}

/// Owned directed CSR multigraph with dense arc IDs.
///
/// Fields are public for construction parity with [`crate::CsrLattice`];
/// callers that mutate them must restore the documented invariants and call
/// [`DirectedCsr::validate`] before using the trait implementations again.
/// Invariants: `out_offsets` and `in_offsets` have length `vertex_count + 1`,
/// start at zero, are non-decreasing, and end at the respective incidence
/// array lengths; every arc appears exactly once in its source's out-range
/// (with matching neighbor) and exactly once in its target's in-range; arc
/// IDs inside the incidence arrays are a permutation of `0..arcs.len()`.
#[derive(Clone, Debug, Default)]
pub struct DirectedCsr {
    /// Physical arcs in dense order; arc ID = position.
    pub arcs: Vec<Arc>,
    /// `out_offsets[v]..out_offsets[v+1]` is vertex `v`'s out-incidence range.
    pub out_offsets: Vec<usize>,
    pub out_neighbors: Vec<usize>,
    pub out_arc_ids: Vec<usize>,
    /// Transpose (incoming) incidence structure, same conventions.
    pub in_offsets: Vec<usize>,
    pub in_neighbors: Vec<usize>,
    pub in_arc_ids: Vec<usize>,
    pub vertex_count: usize,
}

impl DirectedCsr {
    /// Build from arcs with dense IDs in input order.
    ///
    /// Parallel arcs and self-loops are legal (a self-loop contributes one
    /// out- and one in-incidence at its vertex). Endpoint values are
    /// validated; the construction itself cannot fail otherwise.
    pub fn from_arcs(vertex_count: usize, arcs: Vec<Arc>) -> Result<Self, TopologyError> {
        for (arc_index, arc) in arcs.iter().enumerate() {
            for &endpoint in &[arc.source, arc.target] {
                if endpoint >= vertex_count {
                    return Err(TopologyError::EndpointOutOfRange {
                        edge: arc_index,
                        endpoint,
                        vertex_count,
                    });
                }
            }
        }

        let mut out_offsets = vec![0usize; vertex_count + 1];
        let mut in_offsets = vec![0usize; vertex_count + 1];
        for arc in &arcs {
            out_offsets[arc.source + 1] += 1;
            in_offsets[arc.target + 1] += 1;
        }
        for offsets in [&mut out_offsets, &mut in_offsets] {
            for index in 1..offsets.len() {
                offsets[index] += offsets[index - 1];
            }
        }

        let mut out_neighbors = vec![0usize; arcs.len()];
        let mut out_arc_ids = vec![0usize; arcs.len()];
        let mut in_neighbors = vec![0usize; arcs.len()];
        let mut in_arc_ids = vec![0usize; arcs.len()];
        let mut out_cursor = out_offsets.clone();
        let mut in_cursor = in_offsets.clone();
        for (arc_index, arc) in arcs.iter().enumerate() {
            let out_slot = out_cursor[arc.source];
            out_cursor[arc.source] += 1;
            out_neighbors[out_slot] = arc.target;
            out_arc_ids[out_slot] = arc_index;
            let in_slot = in_cursor[arc.target];
            in_cursor[arc.target] += 1;
            in_neighbors[in_slot] = arc.source;
            in_arc_ids[in_slot] = arc_index;
        }

        Ok(Self {
            arcs,
            out_offsets,
            out_neighbors,
            out_arc_ids,
            in_offsets,
            in_neighbors,
            in_arc_ids,
            vertex_count,
        })
    }

    /// Validate every documented invariant (public-field contract).
    pub fn validate(&self) -> Result<(), TopologyError> {
        super::borrowed_directed_csr::validate_out_structure(
            self.vertex_count,
            &self.out_offsets,
            &self.out_neighbors,
            &self.out_arc_ids,
            &self
                .arcs
                .iter()
                .map(|arc| [arc.source, arc.target])
                .collect::<Vec<_>>(),
        )?;
        super::borrowed_directed_csr::validate_transpose_structure(
            self.vertex_count,
            &self.in_offsets,
            &self.in_neighbors,
            &self.in_arc_ids,
            &self
                .arcs
                .iter()
                .map(|arc| [arc.source, arc.target])
                .collect::<Vec<_>>(),
        )
    }
}

impl GraphView for DirectedCsr {
    #[inline]
    fn vertex_count(&self) -> usize {
        self.vertex_count
    }
}

impl DirectedGraphView for DirectedCsr {
    #[inline]
    fn arc_count(&self) -> usize {
        self.arcs.len()
    }

    #[inline]
    fn arc_endpoints(&self, arc: ArcId) -> [VertexId; 2] {
        let physical = self.arcs[arc.index()];
        [
            VertexId::from_valid_index(physical.source),
            VertexId::from_valid_index(physical.target),
        ]
    }
}

impl DirectedOutGraphView for DirectedCsr {
    #[inline]
    fn out_arcs(&self, vertex: VertexId) -> impl Iterator<Item = ArcIncidence> + '_ {
        let index = vertex.index();
        (self.out_offsets[index]..self.out_offsets[index + 1]).map(|incidence| ArcIncidence {
            neighbor: VertexId::from_valid_index(self.out_neighbors[incidence]),
            arc: ArcId::from_valid_index(self.out_arc_ids[incidence]),
        })
    }
}

impl DirectedInGraphView for DirectedCsr {
    #[inline]
    fn in_arcs(&self, vertex: VertexId) -> impl Iterator<Item = ArcIncidence> + '_ {
        let index = vertex.index();
        (self.in_offsets[index]..self.in_offsets[index + 1]).map(|incidence| ArcIncidence {
            neighbor: VertexId::from_valid_index(self.in_neighbors[incidence]),
            arc: ArcId::from_valid_index(self.in_arc_ids[incidence]),
        })
    }
}
