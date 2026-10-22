//! Borrowed zero-copy directed CSR views (D1).
//!
//! [`BorrowedDirectedOutCsr`] borrows a complete outgoing CSR plus the arc
//! endpoint table — the minimum that carries stable arc identity — and
//! implements [`DirectedOutGraphView`]. [`BorrowedDirectedCsr`] additionally
//! borrows the transpose arrays and implements [`DirectedInGraphView`].
//! Construction validates each representation once: offset shape, incidence
//! lengths, neighbor and arc ranges, per-arc source/target consistency, and
//! that every arc appears exactly once on each side (multiplicity).

use std::collections::TryReserveError;

use super::{
    ArcId, ArcIncidence, DirectedGraphView, DirectedInGraphView, DirectedOutGraphView, GraphView,
    TopologyError, VertexId,
};

/// Zero-copy view over a complete outgoing directed CSR.
///
/// `arc_endpoints[a] = [source, target]` is the array-of-structures arc
/// table; the out CSR need not be sorted by neighbor. Empty graphs are
/// valid (`out_offsets = [0]`, other slices empty).
#[derive(Clone, Copy, Debug)]
pub struct BorrowedDirectedOutCsr<'a> {
    vertex_count: usize,
    out_offsets: &'a [usize],
    out_neighbors: &'a [usize],
    out_arc_ids: &'a [usize],
    arc_endpoints: &'a [[usize; 2]],
}

impl<'a> BorrowedDirectedOutCsr<'a> {
    pub fn new(
        out_offsets: &'a [usize],
        out_neighbors: &'a [usize],
        out_arc_ids: &'a [usize],
        arc_endpoints: &'a [[usize; 2]],
    ) -> Result<Self, TopologyError> {
        validate_out_structure(
            usize::MAX,
            out_offsets,
            out_neighbors,
            out_arc_ids,
            arc_endpoints,
        )?;
        Ok(Self {
            vertex_count: out_offsets.len() - 1,
            out_offsets,
            out_neighbors,
            out_arc_ids,
            arc_endpoints,
        })
    }
}

/// Zero-copy view over a complete directed CSR including the transpose.
#[derive(Clone, Copy, Debug)]
pub struct BorrowedDirectedCsr<'a> {
    out: BorrowedDirectedOutCsr<'a>,
    in_offsets: &'a [usize],
    in_neighbors: &'a [usize],
    in_arc_ids: &'a [usize],
}

impl<'a> BorrowedDirectedCsr<'a> {
    pub fn new(
        out_offsets: &'a [usize],
        out_neighbors: &'a [usize],
        out_arc_ids: &'a [usize],
        arc_endpoints: &'a [[usize; 2]],
        in_offsets: &'a [usize],
        in_neighbors: &'a [usize],
        in_arc_ids: &'a [usize],
    ) -> Result<Self, TopologyError> {
        let out =
            BorrowedDirectedOutCsr::new(out_offsets, out_neighbors, out_arc_ids, arc_endpoints)?;
        let vertex_count = out.vertex_count;
        validate_transpose_structure(
            vertex_count,
            in_offsets,
            in_neighbors,
            in_arc_ids,
            arc_endpoints,
        )?;
        Ok(Self {
            out,
            in_offsets,
            in_neighbors,
            in_arc_ids,
        })
    }
}

impl GraphView for BorrowedDirectedOutCsr<'_> {
    #[inline]
    fn vertex_count(&self) -> usize {
        self.vertex_count
    }
}

impl DirectedGraphView for BorrowedDirectedOutCsr<'_> {
    #[inline]
    fn arc_count(&self) -> usize {
        self.arc_endpoints.len()
    }

    #[inline]
    fn arc_endpoints(&self, arc: ArcId) -> [VertexId; 2] {
        let [source, target] = self.arc_endpoints[arc.index()];
        [
            VertexId::from_valid_index(source),
            VertexId::from_valid_index(target),
        ]
    }
}

impl DirectedOutGraphView for BorrowedDirectedOutCsr<'_> {
    #[inline]
    fn out_arcs(&self, vertex: VertexId) -> impl Iterator<Item = ArcIncidence> + '_ {
        let index = vertex.index();
        (self.out_offsets[index]..self.out_offsets[index + 1]).map(|incidence| ArcIncidence {
            neighbor: VertexId::from_valid_index(self.out_neighbors[incidence]),
            arc: ArcId::from_valid_index(self.out_arc_ids[incidence]),
        })
    }
}

impl GraphView for BorrowedDirectedCsr<'_> {
    #[inline]
    fn vertex_count(&self) -> usize {
        self.out.vertex_count()
    }
}

impl DirectedGraphView for BorrowedDirectedCsr<'_> {
    #[inline]
    fn arc_count(&self) -> usize {
        self.out.arc_count()
    }

    #[inline]
    fn arc_endpoints(&self, arc: ArcId) -> [VertexId; 2] {
        self.out.arc_endpoints(arc)
    }
}

impl DirectedOutGraphView for BorrowedDirectedCsr<'_> {
    #[inline]
    fn out_arcs(&self, vertex: VertexId) -> impl Iterator<Item = ArcIncidence> + '_ {
        self.out.out_arcs(vertex)
    }
}

impl DirectedInGraphView for BorrowedDirectedCsr<'_> {
    #[inline]
    fn in_arcs(&self, vertex: VertexId) -> impl Iterator<Item = ArcIncidence> + '_ {
        let index = vertex.index();
        (self.in_offsets[index]..self.in_offsets[index + 1]).map(|incidence| ArcIncidence {
            neighbor: VertexId::from_valid_index(self.in_neighbors[incidence]),
            arc: ArcId::from_valid_index(self.in_arc_ids[incidence]),
        })
    }
}

/// Validate one side (out or transpose) of a directed CSR against the arc
/// table. Returns the vertex count implied by the offsets.
///
/// `declared_vertex_count` is `usize::MAX` when the offsets define the
/// count (borrowed construction) and the owner's count otherwise.
pub(super) fn validate_out_structure(
    declared_vertex_count: usize,
    out_offsets: &[usize],
    out_neighbors: &[usize],
    out_arc_ids: &[usize],
    arc_endpoints: &[[usize; 2]],
) -> Result<usize, TopologyError> {
    if out_neighbors.len() != out_arc_ids.len() {
        return Err(TopologyError::IncidenceLength {
            neighbors: out_neighbors.len(),
            edge_ids: out_arc_ids.len(),
        });
    }
    let vertex_count = validate_offsets(out_offsets, out_neighbors.len())?;
    if declared_vertex_count != usize::MAX && declared_vertex_count != vertex_count {
        return Err(TopologyError::OffsetLength {
            expected: declared_vertex_count + 1,
            actual: out_offsets.len(),
        });
    }
    for (arc, endpoints) in arc_endpoints.iter().enumerate() {
        for &endpoint in endpoints {
            if endpoint >= vertex_count {
                return Err(TopologyError::EndpointOutOfRange {
                    edge: arc,
                    endpoint,
                    vertex_count,
                });
            }
        }
    }
    let mut seen =
        validation_state(arc_endpoints.len()).map_err(|_| TopologyError::ValidationCapacity {
            edge_count: arc_endpoints.len(),
        })?;
    for vertex in 0..vertex_count {
        for incidence in out_offsets[vertex]..out_offsets[vertex + 1] {
            let neighbor = out_neighbors[incidence];
            if neighbor >= vertex_count {
                return Err(TopologyError::NeighborOutOfRange {
                    incidence,
                    neighbor,
                    vertex_count,
                });
            }
            let arc = out_arc_ids[incidence];
            if arc >= arc_endpoints.len() {
                return Err(TopologyError::EdgeIdOutOfRange {
                    incidence,
                    edge: arc,
                    edge_count: arc_endpoints.len(),
                });
            }
            let [source, target] = arc_endpoints[arc];
            if vertex != source || neighbor != target {
                return Err(TopologyError::OutArcMismatch {
                    incidence,
                    vertex,
                    neighbor,
                    arc,
                });
            }
            seen[arc] = seen[arc].saturating_add(1);
            if seen[arc] > 1 {
                return Err(TopologyError::ArcOutMultiplicity {
                    arc,
                    actual: seen[arc] as usize,
                });
            }
        }
    }
    for (arc, &count) in seen.iter().enumerate() {
        if count != 1 {
            return Err(TopologyError::ArcOutMultiplicity {
                arc,
                actual: count as usize,
            });
        }
    }
    Ok(vertex_count)
}

/// Validate the transpose side: every arc enters its target exactly once
/// with the source as neighbor.
pub(super) fn validate_transpose_structure(
    vertex_count: usize,
    in_offsets: &[usize],
    in_neighbors: &[usize],
    in_arc_ids: &[usize],
    arc_endpoints: &[[usize; 2]],
) -> Result<(), TopologyError> {
    if in_neighbors.len() != in_arc_ids.len() {
        return Err(TopologyError::IncidenceLength {
            neighbors: in_neighbors.len(),
            edge_ids: in_arc_ids.len(),
        });
    }
    let transpose_vertex_count = validate_offsets(in_offsets, in_neighbors.len())?;
    if transpose_vertex_count != vertex_count {
        return Err(TopologyError::OffsetLength {
            expected: vertex_count + 1,
            actual: in_offsets.len(),
        });
    }
    let mut seen =
        validation_state(arc_endpoints.len()).map_err(|_| TopologyError::ValidationCapacity {
            edge_count: arc_endpoints.len(),
        })?;
    for vertex in 0..vertex_count {
        for incidence in in_offsets[vertex]..in_offsets[vertex + 1] {
            let neighbor = in_neighbors[incidence];
            if neighbor >= vertex_count {
                return Err(TopologyError::NeighborOutOfRange {
                    incidence,
                    neighbor,
                    vertex_count,
                });
            }
            let arc = in_arc_ids[incidence];
            if arc >= arc_endpoints.len() {
                return Err(TopologyError::EdgeIdOutOfRange {
                    incidence,
                    edge: arc,
                    edge_count: arc_endpoints.len(),
                });
            }
            let [source, target] = arc_endpoints[arc];
            if vertex != target || neighbor != source {
                return Err(TopologyError::InArcMismatch {
                    incidence,
                    vertex,
                    neighbor,
                    arc,
                });
            }
            seen[arc] = seen[arc].saturating_add(1);
            if seen[arc] > 1 {
                return Err(TopologyError::ArcInMultiplicity {
                    arc,
                    actual: seen[arc] as usize,
                });
            }
        }
    }
    for (arc, &count) in seen.iter().enumerate() {
        if count != 1 {
            return Err(TopologyError::ArcInMultiplicity {
                arc,
                actual: count as usize,
            });
        }
    }
    Ok(())
}

fn validate_offsets(offsets: &[usize], incidence_count: usize) -> Result<usize, TopologyError> {
    let vertex_count = offsets
        .len()
        .checked_sub(1)
        .ok_or(TopologyError::OffsetLength {
            expected: 1,
            actual: 0,
        })?;
    if offsets[0] != 0 {
        return Err(TopologyError::OffsetStart { actual: offsets[0] });
    }
    if let Some(index) = offsets.windows(2).position(|window| window[0] > window[1]) {
        return Err(TopologyError::OffsetsNotMonotone { index });
    }
    if offsets[vertex_count] != incidence_count {
        return Err(TopologyError::OffsetEnd {
            expected: incidence_count,
            actual: offsets[vertex_count],
        });
    }
    Ok(vertex_count)
}

fn validation_state(arc_count: usize) -> Result<Vec<u8>, TryReserveError> {
    let mut state = Vec::new();
    state.try_reserve_exact(arc_count)?;
    state.resize(arc_count, 0);
    Ok(state)
}
