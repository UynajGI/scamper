use core::fmt;
use std::collections::TryReserveError;

use crate::{GraphView, VertexId};

/// A validated boundary-crossing query over two vertex sets.
///
/// Duplicate indices are removed. Empty sets are valid and can never cross.
/// If the sets overlap, an active vertex in the intersection is sufficient.
/// IDs remain subject to the graph-view provenance contract: equal dimensions
/// do not make a query constructed for another graph valid.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BoundaryQuery {
    from: Vec<VertexId>,
    to: Vec<VertexId>,
    vertex_count: usize,
}

impl BoundaryQuery {
    pub fn new<G: GraphView>(
        graph: &G,
        from: &[usize],
        to: &[usize],
    ) -> Result<Self, BoundaryQueryError> {
        let vertex_count = graph.vertex_count();
        validate_indices("from", from, vertex_count)?;
        validate_indices("to", to, vertex_count)?;

        let mut from_ids = Vec::new();
        from_ids
            .try_reserve_exact(from.len())
            .map_err(|source| BoundaryQueryError::Capacity {
                set: "from",
                source,
            })?;
        from_ids.extend(
            from.iter()
                .map(|&index| graph.vertex_id(index).expect("validated vertex index")),
        );
        from_ids.sort_unstable();
        from_ids.dedup();

        let mut to_ids = Vec::new();
        to_ids
            .try_reserve_exact(to.len())
            .map_err(|source| BoundaryQueryError::Capacity { set: "to", source })?;
        to_ids.extend(
            to.iter()
                .map(|&index| graph.vertex_id(index).expect("validated vertex index")),
        );
        to_ids.sort_unstable();
        to_ids.dedup();

        Ok(Self {
            from: from_ids,
            to: to_ids,
            vertex_count,
        })
    }

    #[inline]
    pub fn from(&self) -> &[VertexId] {
        &self.from
    }

    #[inline]
    pub fn to(&self) -> &[VertexId] {
        &self.to
    }

    #[inline]
    pub const fn vertex_count(&self) -> usize {
        self.vertex_count
    }
}

fn validate_indices(
    set: &'static str,
    indices: &[usize],
    vertex_count: usize,
) -> Result<(), BoundaryQueryError> {
    if let Some((position, &index)) = indices
        .iter()
        .enumerate()
        .find(|&(_, &index)| index >= vertex_count)
    {
        return Err(BoundaryQueryError::VertexOutOfRange {
            set,
            position,
            index,
            vertex_count,
        });
    }
    Ok(())
}

#[derive(Debug)]
pub enum BoundaryQueryError {
    VertexOutOfRange {
        set: &'static str,
        position: usize,
        index: usize,
        vertex_count: usize,
    },
    Capacity {
        set: &'static str,
        source: TryReserveError,
    },
}

impl fmt::Display for BoundaryQueryError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::VertexOutOfRange {
                set,
                position,
                index,
                vertex_count,
            } => write!(
                formatter,
                "boundary {set} vertex at position {position} is {index}, but vertex count is {vertex_count}"
            ),
            Self::Capacity { set, .. } => {
                write!(formatter, "could not reserve boundary {set} storage")
            }
        }
    }
}

impl std::error::Error for BoundaryQueryError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Capacity { source, .. } => Some(source),
            Self::VertexOutOfRange { .. } => None,
        }
    }
}
