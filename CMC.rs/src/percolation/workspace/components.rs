use std::collections::TryReserveError;
#[cfg(feature = "allocation-probe")]
use std::mem::size_of;

use crate::VertexId;

/// Per-vertex component label; inactive vertices have no component.
pub type ComponentLabel = Option<VertexId>;

/// Reusable storage for undirected component analysis.
///
/// Call [`Self::prepare`] before a repeated workload to make each subsequent
/// analysis allocation-free. Internal union-find and query scratch are private.
#[derive(Clone, Debug, Default)]
pub struct ComponentWorkspace {
    pub(crate) parent: Vec<usize>,
    pub(crate) size: Vec<usize>,
    pub(crate) canonical: Vec<usize>,
    pub(crate) query_root_stamps: Vec<u64>,
    pub(crate) query_outcomes: Vec<bool>,
    pub(crate) labels: Vec<ComponentLabel>,
    pub(crate) stamp: u64,
}

impl ComponentWorkspace {
    pub const fn new() -> Self {
        Self {
            parent: Vec::new(),
            size: Vec::new(),
            canonical: Vec::new(),
            query_root_stamps: Vec::new(),
            query_outcomes: Vec::new(),
            labels: Vec::new(),
            stamp: 0,
        }
    }

    /// Reserve and initialize storage for a repeated workload.
    pub fn prepare(
        &mut self,
        vertex_capacity: usize,
        query_capacity: usize,
        with_labels: bool,
    ) -> Result<(), TryReserveError> {
        reserve_len(&mut self.parent, vertex_capacity)?;
        reserve_len(&mut self.size, vertex_capacity)?;
        reserve_len(&mut self.canonical, vertex_capacity)?;
        reserve_len(&mut self.query_root_stamps, vertex_capacity)?;
        reserve_len(&mut self.query_outcomes, query_capacity)?;
        if with_labels {
            reserve_len(&mut self.labels, vertex_capacity)?;
        }
        Ok(())
    }

    /// Reserved vertex capacity shared by all mandatory analyzer arrays.
    pub fn vertex_capacity(&self) -> usize {
        self.parent
            .capacity()
            .min(self.size.capacity())
            .min(self.canonical.capacity())
            .min(self.query_root_stamps.capacity())
    }

    pub fn query_capacity(&self) -> usize {
        self.query_outcomes.capacity()
    }

    pub fn label_capacity(&self) -> usize {
        self.labels.capacity()
    }

    /// Audit-only retained-capacity accounting for the allocation probe.
    #[cfg(feature = "allocation-probe")]
    #[doc(hidden)]
    pub fn capacity_audit(&self) -> WorkspaceCapacityAudit {
        WorkspaceCapacityAudit {
            parent: vector_capacity_bytes(&self.parent),
            size: vector_capacity_bytes(&self.size),
            canonical: vector_capacity_bytes(&self.canonical),
            query_root_stamps: vector_capacity_bytes(&self.query_root_stamps),
            query_outcomes: vector_capacity_bytes(&self.query_outcomes),
            labels: vector_capacity_bytes(&self.labels),
        }
    }

    pub(crate) fn resize_for(
        &mut self,
        vertex_count: usize,
        query_count: usize,
        with_labels: bool,
    ) -> Result<(), TryReserveError> {
        self.prepare(vertex_count, query_count, with_labels)?;
        self.parent.resize(vertex_count, usize::MAX);
        self.size.resize(vertex_count, 0);
        self.canonical.resize(vertex_count, usize::MAX);
        self.query_root_stamps.resize(vertex_count, 0);
        self.query_outcomes.resize(query_count, false);
        if with_labels {
            self.labels.resize(vertex_count, None);
        }
        Ok(())
    }

    pub(crate) fn next_stamp(&mut self) -> u64 {
        if self.stamp == u64::MAX {
            self.query_root_stamps.fill(0);
            self.stamp = 1;
        } else {
            self.stamp += 1;
        }
        self.stamp
    }
}

fn reserve_len<T>(storage: &mut Vec<T>, len: usize) -> Result<(), TryReserveError> {
    if storage.capacity() < len {
        storage.try_reserve_exact(len - storage.len())?;
    }
    Ok(())
}

#[cfg(feature = "allocation-probe")]
fn vector_capacity_bytes<T>(storage: &Vec<T>) -> usize {
    storage.capacity().saturating_mul(size_of::<T>())
}

/// Audit-only byte accounting for every owned workspace vector.
#[cfg(feature = "allocation-probe")]
#[doc(hidden)]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct WorkspaceCapacityAudit {
    pub parent: usize,
    pub size: usize,
    pub canonical: usize,
    pub query_root_stamps: usize,
    pub query_outcomes: usize,
    pub labels: usize,
}

#[cfg(feature = "allocation-probe")]
impl WorkspaceCapacityAudit {
    pub const fn total_bytes(self) -> usize {
        self.parent
            .saturating_add(self.size)
            .saturating_add(self.canonical)
            .saturating_add(self.query_root_stamps)
            .saturating_add(self.query_outcomes)
            .saturating_add(self.labels)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::percolation::{analyze, BoundaryQuery, EdgeActivity, VertexActivity};
    use crate::{build_chain, EdgeId, GraphView, UndirectedGraphView, VertexId};

    struct AllVertices(usize);

    impl VertexActivity for AllVertices {
        fn vertex_count(&self) -> usize {
            self.0
        }

        fn vertex_active(&self, _vertex: VertexId) -> bool {
            true
        }
    }

    struct AllEdges(usize);

    impl EdgeActivity for AllEdges {
        fn edge_count(&self) -> usize {
            self.0
        }

        fn edge_active(&self, _edge: EdgeId) -> bool {
            true
        }
    }

    #[test]
    fn generation_wrap_clears_stamps_on_real_analysis() {
        let graph = build_chain(4, false);
        let query = BoundaryQuery::new(&graph, &[0], &[3]).unwrap();
        let queries = [query];
        let vertices = AllVertices(graph.vertex_count());
        let edges = AllEdges(graph.edge_count());

        let mut wrapped = ComponentWorkspace::new();
        wrapped.prepare(4, queries.len(), false).unwrap();
        wrapped.query_root_stamps.resize(4, 0);
        wrapped.query_root_stamps.copy_from_slice(&[7, 11, 13, 17]);
        wrapped.stamp = u64::MAX;
        let wrapped_result = analyze(&graph, &vertices, &edges, &queries, &mut wrapped).unwrap();
        let wrapped_summary = wrapped_result.summary();
        let wrapped_outcomes = wrapped_result.query_outcomes().to_vec();

        let mut clean = ComponentWorkspace::new();
        let clean_result = analyze(&graph, &vertices, &edges, &queries, &mut clean).unwrap();
        assert_eq!(wrapped_summary, clean_result.summary());
        assert_eq!(wrapped_outcomes, clean_result.query_outcomes());
        assert_eq!(wrapped_outcomes, [true]);
        assert_eq!(wrapped.stamp, 1);
    }
}
