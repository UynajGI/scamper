use super::{EdgeId, GraphView, Incidence, TopologyError, UndirectedGraphView, VertexId};

const COUNT_MASK: u8 = 0b0000_0011;
const SOURCE_SEEN: u8 = 0b0000_0100;
const TARGET_SEEN: u8 = 0b0000_1000;

/// Zero-copy read-only view over a complete undirected CSR representation.
///
/// `edge_endpoints` uses an array-of-structures layout: entry `e` is the two
/// endpoints of physical edge `e`. Empty graphs are valid and use `offsets =
/// [0]` with all other slices empty. Adjacency need not be sorted.
///
/// Construction validates once using temporary state of one byte per physical
/// edge. That state is dropped before the view is returned; the returned view
/// owns no topology payload. Capacity and allocation failures reported by
/// `try_reserve_exact` become [`TopologyError::ValidationCapacity`], although an
/// allocator configured to abort on OOM may still terminate the process.
///
/// An offsets/neighbors-only CSR is intentionally insufficient: stable physical
/// edge identity cannot be inferred for parallel edges and self-loops.
#[derive(Clone, Copy, Debug)]
pub struct BorrowedUndirectedCsr<'a> {
    offsets: &'a [usize],
    neighbors: &'a [usize],
    edge_ids: &'a [usize],
    edge_endpoints: &'a [[usize; 2]],
}

impl<'a> BorrowedUndirectedCsr<'a> {
    pub fn new(
        offsets: &'a [usize],
        neighbors: &'a [usize],
        edge_ids: &'a [usize],
        edge_endpoints: &'a [[usize; 2]],
    ) -> Result<Self, TopologyError> {
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
        if neighbors.len() != edge_ids.len() {
            return Err(TopologyError::IncidenceLength {
                neighbors: neighbors.len(),
                edge_ids: edge_ids.len(),
            });
        }
        let end = offsets[vertex_count];
        if end != neighbors.len() {
            return Err(TopologyError::OffsetEnd {
                expected: neighbors.len(),
                actual: end,
            });
        }

        for (edge, endpoints) in edge_endpoints.iter().enumerate() {
            for &endpoint in endpoints {
                if endpoint >= vertex_count {
                    return Err(TopologyError::EndpointOutOfRange {
                        edge,
                        endpoint,
                        vertex_count,
                    });
                }
            }
        }

        let mut incidence_state = validation_state(edge_endpoints.len())?;
        for vertex in 0..vertex_count {
            for incidence in offsets[vertex]..offsets[vertex + 1] {
                let neighbor = neighbors[incidence];
                if neighbor >= vertex_count {
                    return Err(TopologyError::NeighborOutOfRange {
                        incidence,
                        neighbor,
                        vertex_count,
                    });
                }
                let edge = edge_ids[incidence];
                if edge >= edge_endpoints.len() {
                    return Err(TopologyError::EdgeIdOutOfRange {
                        incidence,
                        edge,
                        edge_count: edge_endpoints.len(),
                    });
                }
                let [source, target] = edge_endpoints[edge];
                if !((vertex == source && neighbor == target)
                    || (vertex == target && neighbor == source))
                {
                    return Err(TopologyError::IncidenceMismatch {
                        incidence,
                        vertex,
                        neighbor,
                        edge,
                    });
                }

                let state = &mut incidence_state[edge];
                let count = *state & COUNT_MASK;
                if count == 2 {
                    return Err(TopologyError::IncidenceMultiplicity { edge, actual: 3 });
                }
                *state = (*state & !COUNT_MASK) | (count + 1);
                if source != target {
                    *state |= if vertex == source {
                        SOURCE_SEEN
                    } else {
                        TARGET_SEEN
                    };
                }
            }
        }
        for (edge, state) in incidence_state.into_iter().enumerate() {
            let actual = usize::from(state & COUNT_MASK);
            if actual != 2 {
                return Err(TopologyError::IncidenceMultiplicity { edge, actual });
            }
            if edge_endpoints[edge][0] != edge_endpoints[edge][1]
                && state & (SOURCE_SEEN | TARGET_SEEN) != SOURCE_SEEN | TARGET_SEEN
            {
                return Err(TopologyError::EndpointIncidenceMultiplicity {
                    edge,
                    source_seen: state & SOURCE_SEEN != 0,
                    target_seen: state & TARGET_SEEN != 0,
                });
            }
        }

        Ok(Self {
            offsets,
            neighbors,
            edge_ids,
            edge_endpoints,
        })
    }
}

fn validation_state(edge_count: usize) -> Result<Vec<u8>, TopologyError> {
    let mut state = Vec::new();
    state
        .try_reserve_exact(edge_count)
        .map_err(|_| TopologyError::ValidationCapacity { edge_count })?;
    state.resize(edge_count, 0);
    Ok(state)
}

impl GraphView for BorrowedUndirectedCsr<'_> {
    #[inline(always)]
    fn vertex_count(&self) -> usize {
        self.offsets.len() - 1
    }
}

impl UndirectedGraphView for BorrowedUndirectedCsr<'_> {
    #[inline(always)]
    fn edge_count(&self) -> usize {
        self.edge_endpoints.len()
    }

    #[inline(always)]
    fn edge_endpoints(&self, edge: EdgeId) -> [VertexId; 2] {
        self.edge_endpoints[edge.index()].map(VertexId::from_valid_index)
    }

    #[inline(always)]
    fn incidences(&self, vertex: VertexId) -> impl Iterator<Item = Incidence> + '_ {
        let range = self.offsets[vertex.index()]..self.offsets[vertex.index() + 1];
        self.neighbors[range.clone()]
            .iter()
            .copied()
            .zip(self.edge_ids[range].iter().copied())
            .map(|(neighbor, edge)| Incidence {
                neighbor: VertexId::from_valid_index(neighbor),
                edge: EdgeId::from_valid_index(edge),
            })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn validation_capacity_overflow_is_reported() {
        assert_eq!(
            validation_state(usize::MAX),
            Err(TopologyError::ValidationCapacity {
                edge_count: usize::MAX
            })
        );
    }
}
