//! Source-to-target reachability: early-exit BFS with stamped visitation.

use super::DirectedAnalysisError;
use crate::{DirectedOutGraphView, VertexId};

/// Reusable BFS workspace.
///
/// `visited` carries per-vertex stamps so repeated queries need no clearing
/// pass: a query touches only its own stamp.
#[derive(Debug, Default)]
pub struct ReachabilityWorkspace {
    visited: Vec<u32>,
    queue: Vec<usize>,
    stamp: u32,
}

impl ReachabilityWorkspace {
    /// Prepare for graphs with `vertex_count` vertices.
    pub fn prepare(&mut self, vertex_count: usize) -> Result<(), DirectedAnalysisError> {
        if self.visited.len() < vertex_count {
            self.visited.clear();
            self.visited.try_reserve_exact(vertex_count)?;
            self.visited.resize(vertex_count, 0);
        }
        if self.queue.capacity() < vertex_count {
            self.queue = Vec::with_capacity(vertex_count);
        }
        Ok(())
    }
}

/// Whether a directed path (possibly the empty path) leads `source ->
/// target`, following out-arcs.
///
/// `reachable(v, v)` is true by definition. Early exit as soon as the
/// target is dequeued. The workspace stamp mechanism makes repeated
/// queries allocation-free and O(visited subgraph) each.
pub fn reachable<G: DirectedOutGraphView>(
    graph: &G,
    source: VertexId,
    target: VertexId,
    workspace: &mut ReachabilityWorkspace,
) -> Result<bool, DirectedAnalysisError> {
    workspace.prepare(graph.vertex_count())?;
    // Stamp rollover after 2^32 queries on one workspace: clear every slot
    // and restart from 1 so stale stamps can never alias a fresh one.
    if workspace.stamp == u32::MAX {
        workspace.visited.iter_mut().for_each(|slot| *slot = 0);
        workspace.stamp = 0;
    }
    workspace.stamp += 1;
    let stamp = workspace.stamp;

    let start = source.index();
    let goal = target.index();
    workspace.queue.clear();
    workspace.queue.push(start);
    workspace.visited[start] = stamp;
    let mut head = 0;
    while head < workspace.queue.len() {
        let vertex = workspace.queue[head];
        head += 1;
        if vertex == goal {
            return Ok(true);
        }
        let identity = graph.vertex_id(vertex).expect("dense index");
        for incidence in graph.out_arcs(identity) {
            let successor = incidence.neighbor.index();
            if workspace.visited[successor] != stamp {
                workspace.visited[successor] = stamp;
                workspace.queue.push(successor);
            }
        }
    }
    Ok(false)
}
