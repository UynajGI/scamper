//! Strongly connected components: iterative Tarjan SCC.

use super::DirectedAnalysisError;
use crate::{DirectedGraphView, DirectedOutGraphView, VertexId};

/// Raw SCC facts for one directed graph.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct StrongComponentSummary {
    /// Strongly connected components, including singleton components of
    /// isolated vertices.
    pub component_count: usize,
    /// Canonical identity (lowest member vertex ID) of the largest SCC;
    /// equal sizes resolve to the lowest identity.
    pub largest_scc: Option<VertexId>,
    /// Size of [`StrongComponentSummary::largest_scc`], `0` when empty.
    pub largest_scc_size: usize,
    /// Raw `sum(s_i^2)` over SCC sizes.
    pub raw_second_moment: u128,
    /// True iff the graph has no directed cycle at all: every SCC is a
    /// singleton and no self-loop exists. (A self-loop is a cycle but
    /// leaves its SCC a singleton, so it is tracked separately.)
    pub acyclic: bool,
}

/// Borrowed SCC output: the summary plus per-vertex canonical component
/// labels (borrowed from the workspace).
#[derive(Clone, Copy, Debug)]
pub struct StrongComponentResult<'a> {
    summary: StrongComponentSummary,
    labels: &'a [VertexId],
}

impl StrongComponentResult<'_> {
    #[inline]
    pub const fn summary(&self) -> StrongComponentSummary {
        self.summary
    }

    /// Canonical SCC label per dense vertex index; two vertices are in the
    /// same SCC iff their labels are equal.
    #[inline]
    pub fn labels(&self) -> &[VertexId] {
        self.labels
    }
}

const UNSET: u32 = u32::MAX;

/// Reusable Tarjan workspace (adjacency mirror, indices, lowlinks, stacks).
#[derive(Debug, Default)]
pub struct StrongComponentWorkspace {
    /// Materialized out-adjacency: `adjacency[v]` is the vertex's
    /// successor list, flat with `adj_start` offsets. Built once per
    /// analysis so the DFS walks plain slices.
    adjacency: Vec<usize>,
    adj_start: Vec<usize>,
    index: Vec<u32>,
    lowlink: Vec<u32>,
    on_stack: Vec<bool>,
    /// Tarjan's component stack of vertex indices.
    stack: Vec<usize>,
    /// Explicit DFS stack of (vertex, successor cursor).
    dfs: Vec<(usize, usize)>,
    labels: Vec<VertexId>,
    next_index: u32,
}

impl StrongComponentWorkspace {
    /// Prepare for graphs with `vertex_count` vertices and `arc_count`
    /// arcs.
    pub fn prepare(
        &mut self,
        vertex_count: usize,
        arc_count: usize,
    ) -> Result<(), DirectedAnalysisError> {
        self.adjacency.clear();
        self.adjacency.try_reserve_exact(arc_count)?;
        self.adj_start.clear();
        self.adj_start.try_reserve_exact(vertex_count + 1)?;
        self.index.clear();
        self.index.try_reserve_exact(vertex_count)?;
        self.index.resize(vertex_count, UNSET);
        self.lowlink.clear();
        self.lowlink.try_reserve_exact(vertex_count)?;
        self.lowlink.resize(vertex_count, UNSET);
        self.on_stack.clear();
        self.on_stack.try_reserve_exact(vertex_count)?;
        self.on_stack.resize(vertex_count, false);
        self.stack.clear();
        self.dfs.clear();
        self.labels.clear();
        self.labels.try_reserve_exact(vertex_count)?;
        self.labels
            .resize(vertex_count, VertexId::from_valid_index(0));
        self.next_index = 0;
        Ok(())
    }
}

/// Analyze strongly connected components with an iterative Tarjan pass.
///
/// Requires out-adjacency ([`DirectedOutGraphView`]); the workspace
/// materializes the successor lists once (O(V + A)) and the DFS then walks
/// plain slices. Non-recursive by construction: deep graphs cannot
/// overflow the call stack.
pub fn analyze_strong_components<'a, G>(
    graph: &G,
    workspace: &'a mut StrongComponentWorkspace,
) -> Result<StrongComponentResult<'a>, DirectedAnalysisError>
where
    G: DirectedGraphView + DirectedOutGraphView,
{
    let vertex_count = graph.vertex_count();
    // Tarjan indices live in u32 with UNSET = u32::MAX: the discovery
    // index space covers at most 2^32 - 1 vertices. A graph that large
    // cannot be reserved on real hardware (the typed capacity error fires
    // first), but the bound is recorded explicitly.
    debug_assert!(vertex_count < u32::MAX as usize);
    workspace.prepare(vertex_count, graph.arc_count())?;

    let mut has_self_loop = false;
    for vertex in 0..vertex_count {
        workspace.adj_start.push(workspace.adjacency.len());
        let identity = graph.vertex_id(vertex).expect("dense index");
        for incidence in graph.out_arcs(identity) {
            let successor = incidence.neighbor.index();
            has_self_loop |= successor == vertex;
            workspace.adjacency.push(successor);
        }
    }
    workspace.adj_start.push(workspace.adjacency.len());

    let mut summary = StrongComponentSummary {
        component_count: 0,
        largest_scc: None,
        largest_scc_size: 0,
        raw_second_moment: 0,
        acyclic: true,
    };

    for root in 0..vertex_count {
        if workspace.index[root] != UNSET {
            continue;
        }
        workspace.dfs.push((root, workspace.adj_start[root]));
        while let Some(&mut (vertex, ref mut cursor)) = workspace.dfs.last_mut() {
            if workspace.index[vertex] == UNSET {
                let index = workspace.next_index;
                workspace.next_index = index + 1;
                workspace.index[vertex] = index;
                workspace.lowlink[vertex] = index;
                workspace.stack.push(vertex);
                workspace.on_stack[vertex] = true;
            }

            let mut pushed = false;
            let end = workspace.adj_start[vertex + 1];
            while *cursor < end {
                let successor = workspace.adjacency[*cursor];
                *cursor += 1;
                if workspace.index[successor] == UNSET {
                    workspace
                        .dfs
                        .push((successor, workspace.adj_start[successor]));
                    pushed = true;
                    break;
                }
                if workspace.on_stack[successor]
                    && workspace.index[successor] < workspace.lowlink[vertex]
                {
                    workspace.lowlink[vertex] = workspace.index[successor];
                }
            }
            if pushed {
                continue;
            }

            workspace.dfs.pop();
            if let Some(&(parent, _)) = workspace.dfs.last() {
                let child_low = workspace.lowlink[vertex];
                if child_low < workspace.lowlink[parent] {
                    workspace.lowlink[parent] = child_low;
                }
            }
            if workspace.lowlink[vertex] == workspace.index[vertex] {
                let mut size = 0usize;
                let mut canonical = usize::MAX;
                // Locate the component extent: members run from this
                // index to the stack top, ending at `vertex`.
                let mut members_start = workspace.stack.len();
                loop {
                    members_start -= 1;
                    let member = workspace.stack[members_start];
                    workspace.on_stack[member] = false;
                    canonical = canonical.min(member);
                    size += 1;
                    if member == vertex {
                        break;
                    }
                }
                let identity = graph.vertex_id(canonical).expect("dense canonical index");
                // Label every member before dropping them from the stack.
                for &member in &workspace.stack[members_start..] {
                    workspace.labels[member] = identity;
                }
                workspace.stack.truncate(members_start);
                summary.component_count += 1;
                summary.raw_second_moment += (size as u128) * (size as u128);
                if size > 1 {
                    summary.acyclic = false;
                }
                if size > summary.largest_scc_size
                    || (size == summary.largest_scc_size
                        && summary.largest_scc.is_none_or(|current| identity < current))
                {
                    summary.largest_scc_size = size;
                    summary.largest_scc = Some(identity);
                }
            }
        }
    }
    summary.acyclic &= !has_self_loop;

    Ok(StrongComponentResult {
        summary,
        labels: &workspace.labels,
    })
}
