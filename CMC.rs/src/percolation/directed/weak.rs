//! Weak components: connectivity with arc directions ignored.

use super::{find_root, DirectedAnalysisError};
use crate::{DirectedGraphView, VertexId};

/// Raw weak-component facts for one directed graph.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct WeakComponentSummary {
    /// Weak components, including isolated vertices.
    pub component_count: usize,
    /// Canonical identity (lowest member vertex ID) of the largest weak
    /// component; equal sizes resolve to the lowest identity.
    pub largest_component: Option<VertexId>,
    /// Size of [`WeakComponentSummary::largest_component`], `0` when empty.
    pub largest_component_size: usize,
    /// Raw `sum(s_i^2)` over weak-component sizes.
    pub raw_second_moment: u128,
}

/// Borrowed weak-component output: the summary plus per-vertex canonical
/// component labels (borrowed from the workspace).
#[derive(Clone, Copy, Debug)]
pub struct WeakComponentResult<'a> {
    summary: WeakComponentSummary,
    labels: &'a [VertexId],
}

impl WeakComponentResult<'_> {
    #[inline]
    pub const fn summary(&self) -> WeakComponentSummary {
        self.summary
    }

    /// Canonical component label per dense vertex index.
    #[inline]
    pub fn labels(&self) -> &[VertexId] {
        self.labels
    }
}

/// Reusable union-find workspace for weak-component analysis.
#[derive(Debug, Default)]
pub struct WeakComponentWorkspace {
    parent: Vec<usize>,
    size: Vec<usize>,
    canonical: Vec<usize>,
    labels: Vec<VertexId>,
}

impl WeakComponentWorkspace {
    /// Prepare for graphs with `vertex_count` vertices.
    pub fn prepare(&mut self, vertex_count: usize) -> Result<(), DirectedAnalysisError> {
        self.parent.clear();
        self.parent.try_reserve_exact(vertex_count)?;
        self.parent.extend(0..vertex_count);
        self.size.clear();
        self.size.try_reserve_exact(vertex_count)?;
        self.size.resize(vertex_count, 1);
        self.canonical.clear();
        self.canonical.try_reserve_exact(vertex_count)?;
        self.canonical.extend(0..vertex_count);
        self.labels.clear();
        self.labels.try_reserve_exact(vertex_count)?;
        self.labels
            .resize(vertex_count, VertexId::from_valid_index(0));
        Ok(())
    }
}

/// Analyze weak components: one pass over the physical arc table, arcs
/// unioned as undirected connections.
pub fn analyze_weak_components<'a, G: DirectedGraphView>(
    graph: &G,
    workspace: &'a mut WeakComponentWorkspace,
) -> Result<WeakComponentResult<'a>, DirectedAnalysisError> {
    let vertex_count = graph.vertex_count();
    workspace.prepare(vertex_count)?;

    for arc in graph.arc_ids() {
        let [source, target] = graph.arc_endpoints(arc);
        // Self-loops union a vertex with itself: a no-op by convention.
        super::union_by_size(
            &mut workspace.parent,
            &mut workspace.size,
            &mut workspace.canonical,
            source.index(),
            target.index(),
        );
    }

    let mut summary = WeakComponentSummary {
        component_count: 0,
        largest_component: None,
        largest_component_size: 0,
        raw_second_moment: 0,
    };
    for index in 0..vertex_count {
        let root = find_root(&mut workspace.parent, index);
        let canonical = workspace.canonical[root];
        let identity = graph.vertex_id(canonical).expect("dense canonical index");
        workspace.labels[index] = identity;
        if root == index {
            let size = workspace.size[index];
            summary.component_count += 1;
            summary.raw_second_moment += (size as u128) * (size as u128);
            if size > summary.largest_component_size
                || (size == summary.largest_component_size
                    && summary
                        .largest_component
                        .is_none_or(|current| identity < current))
            {
                summary.largest_component_size = size;
                summary.largest_component = Some(identity);
            }
        }
    }

    Ok(WeakComponentResult {
        summary,
        labels: &workspace.labels,
    })
}
