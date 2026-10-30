//! G1: incremental component process — generic union substrate (roadmap §13).
//!
//! Achlioptas (G2) and invasion (G3) share infrastructure, not a fabricated
//! common law trait: what they share is *incremental union with online
//! observables*. This module is that substrate and nothing more. It knows
//! nothing about candidate edges, frontiers, thresholds, or rules — the Z1
//! step (one uniformly random vertex activation), the G2 step (one accepted
//! candidate edge), and the G3 step (one frontier growth) are all just
//! sequences of `activate_vertex`/`activate_edge` calls.
//!
//! The workspace is generic over [`UndirectedGraphView`] (edge endpoints and
//! incidences come from the caller's graph) and tracks, after EVERY accepted
//! operation, the online observables of the active subgraph:
//!
//! - component count (active isolates included),
//! - largest component size with its canonical identity (the lowest dense
//!   vertex index in the component — the F3 `ComponentSummary` convention,
//!   equal-size ties resolved by the lowest identity),
//! - raw second moment `sum(s_i^2)`,
//! - active vertex and active physical-edge counts.
//!
//! ## Process coordinates are caller-defined
//!
//! Every accepted operation appends `(step, summary)` to the workspace's
//! [`ProcessTrajectory`]. The `step` value is a caller-defined coordinate that
//! this module never interprets: an activation index, an accepted-edge count,
//! an invasion time, or anything else the driving law means by "step". This is
//! the roadmap's "record process coordinates without forcing a specific law".
//! The trajectory budget is exact: each vertex and each physical edge is
//! accepted at most once, so `prepare(vertices, max_edges)` reserves
//! `vertices + max_edges` entries and an accepted operation can never overflow
//! it.
//!
//! ## Conventions (F3 parity)
//!
//! An edge is *active* iff it was activated AND both endpoints are active;
//! self-loops count once and parallel edges count separately. Activation is
//! order-independent: the summary after any sequence depends only on the sets
//! of activated vertices and edges (verified by exhaustive-ordering gates in
//! `tests/physics/incremental_process.rs` against a from-scratch reference).
//!
//! ## Reuse
//!
//! `prepare` once per topology size; `reset` between realizations
//! (deactivates everything — union-find cannot split components, so
//! per-element deactivation is out of scope — and clears the trajectory);
//! zero steady-state allocation thereafter.

use core::fmt;
use std::collections::TryReserveError;

use crate::UndirectedGraphView;

/// Edge activation state.
///
/// An activated edge is *live* (counted as active, endpoints unioned) exactly
/// when both endpoints are active; the transition happens once, at the later
/// of the two events, so self-loops (two incidences at one vertex) are counted
/// once.
const EDGE_INACTIVE: u8 = 0;
const EDGE_ACTIVATED: u8 = 1;
const EDGE_LIVE: u8 = 2;

/// Online component facts after one accepted incremental operation (or the
/// current state via [`IncrementalUnionWorkspace::summary`]).
///
/// Mirrors the F3 [`crate::ComponentSummary`] conventions, including the
/// canonical largest-component identity: the lowest dense vertex index in the
/// selected component, with equal-size ties resolved by the lowest identity.
/// Map it through `graph.vertex_id(index)` for the typed ID of the driving
/// graph.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct IncrementalSummary {
    /// Vertices activated so far.
    pub active_vertex_count: usize,
    /// Activated physical edges whose two endpoints are both active.
    /// Self-loops count once; parallel edges count separately.
    pub active_edge_count: usize,
    /// Connected components among active vertices, including active isolates.
    pub component_count: usize,
    /// Canonical identity (lowest dense vertex index) of the selected largest
    /// component; `None` when no vertex is active.
    pub largest_component: Option<usize>,
    /// Active vertices in the largest component; zero when nothing is active.
    pub largest_component_size: usize,
    /// Raw `sum(s_i^2)` over all active components.
    pub raw_second_moment: u128,
}

/// Recorded process coordinates of one realization: parallel reusable buffers
/// of `(step, summary)` pairs, one entry per accepted operation.
///
/// The `step` value is caller-defined and never interpreted here (see the
/// module documentation). Buffers are prepared with a fixed capacity budget
/// and reused across realizations; exposure is by slice.
#[derive(Clone, Debug, Default)]
pub struct ProcessTrajectory {
    steps: Vec<usize>,
    summaries: Vec<IncrementalSummary>,
}

impl ProcessTrajectory {
    /// Reserve capacity for at most `capacity` recorded operations.
    pub fn prepare(&mut self, capacity: usize) -> Result<(), TryReserveError> {
        reserve_len(&mut self.steps, capacity)?;
        reserve_len(&mut self.summaries, capacity)?;
        Ok(())
    }

    /// Drop every recorded entry, retaining capacity (O(1)).
    #[inline]
    pub fn clear(&mut self) {
        self.steps.clear();
        self.summaries.clear();
    }

    /// Number of recorded operations.
    #[inline]
    pub fn len(&self) -> usize {
        self.steps.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.steps.is_empty()
    }

    /// Caller-defined step coordinate per accepted operation.
    #[inline]
    pub fn steps(&self) -> &[usize] {
        &self.steps
    }

    /// Online summary per accepted operation.
    #[inline]
    pub fn summaries(&self) -> &[IncrementalSummary] {
        &self.summaries
    }

    /// Append one entry. The workspace guarantees the prepared budget (at
    /// most one accepted operation per vertex and per physical edge), so this
    /// never allocates in a correctly prepared workspace.
    #[inline]
    fn push(&mut self, step: usize, summary: IncrementalSummary) {
        debug_assert!(
            self.steps.len() < self.steps.capacity()
                && self.summaries.len() < self.summaries.capacity(),
            "trajectory exceeded its prepared budget"
        );
        self.steps.push(step);
        self.summaries.push(summary);
    }
}

/// Reusable workspace for incremental union with online observables.
///
/// Call [`Self::prepare`] once for a topology size (this also resets all
/// state); drive any number of realizations with `activate_vertex` /
/// `activate_edge`; call [`Self::reset`] between realizations. Zero
/// steady-state allocation: every buffer is reserved in `prepare` and the
/// trajectory budget `vertices + max_edges` cannot be exceeded by accepted
/// operations.
#[derive(Clone, Debug, Default)]
pub struct IncrementalUnionWorkspace {
    parent: Vec<usize>,
    size: Vec<usize>,
    canonical: Vec<usize>,
    vertex_active: Vec<u8>,
    edge_state: Vec<u8>,
    trajectory: ProcessTrajectory,
    // Online aggregates over the current active subgraph.
    components: usize,
    largest_size: usize,
    /// `usize::MAX` when nothing is active; otherwise the canonical identity
    /// of the selected largest component.
    largest_identity: usize,
    second_moment: u128,
    active_vertices: usize,
    active_edges: usize,
    prepared_vertices: usize,
    prepared_edges: usize,
}

impl IncrementalUnionWorkspace {
    pub const fn new() -> Self {
        Self {
            parent: Vec::new(),
            size: Vec::new(),
            canonical: Vec::new(),
            vertex_active: Vec::new(),
            edge_state: Vec::new(),
            trajectory: ProcessTrajectory {
                steps: Vec::new(),
                summaries: Vec::new(),
            },
            components: 0,
            largest_size: 0,
            largest_identity: usize::MAX,
            second_moment: 0,
            active_vertices: 0,
            active_edges: 0,
            prepared_vertices: 0,
            prepared_edges: 0,
        }
    }

    /// Prepare for topologies with exactly `vertices` vertices and at most
    /// `max_edges` physical edges, resetting every previous state.
    ///
    /// Re-preparing with a different size is supported: logical arrays are
    /// rebuilt at the new size, so no stale state survives a shrink. The
    /// trajectory budget is `vertices + max_edges` (each vertex and each edge
    /// is accepted at most once).
    pub fn prepare(&mut self, vertices: usize, max_edges: usize) -> Result<(), IncrementalError> {
        reserve_len(&mut self.parent, vertices)?;
        reserve_len(&mut self.size, vertices)?;
        reserve_len(&mut self.canonical, vertices)?;
        reserve_len(&mut self.vertex_active, vertices)?;
        reserve_len(&mut self.edge_state, max_edges)?;
        self.trajectory
            .prepare(vertices.saturating_add(max_edges))?;
        // Rebuild the logical arrays at the new size so a shrink leaves no
        // stale entries (clear + resize rewrites every slot; capacity is
        // retained).
        self.parent.clear();
        self.parent.resize(vertices, usize::MAX);
        self.size.clear();
        self.size.resize(vertices, 0);
        self.canonical.clear();
        self.canonical.resize(vertices, usize::MAX);
        self.vertex_active.clear();
        self.vertex_active.resize(vertices, 0);
        self.edge_state.clear();
        self.edge_state.resize(max_edges, EDGE_INACTIVE);
        self.prepared_vertices = vertices;
        self.prepared_edges = max_edges;
        self.reset();
        Ok(())
    }

    /// Prepared topology sizes.
    #[inline]
    pub const fn prepared_sizes(&self) -> (usize, usize) {
        (self.prepared_vertices, self.prepared_edges)
    }

    /// Deactivate everything and clear the trajectory for the next
    /// realization, retaining all capacity.
    ///
    /// O(V + E) (two byte-array clears; the `usize` union-find arrays need no
    /// rewrite because every activation writes `parent`/`size`/`canonical`
    /// for its vertex before any read). Per-element deactivation is out of
    /// scope: union-find cannot split a merged component.
    pub fn reset(&mut self) {
        self.vertex_active.fill(0);
        self.edge_state.fill(EDGE_INACTIVE);
        self.components = 0;
        self.largest_size = 0;
        self.largest_identity = usize::MAX;
        self.second_moment = 0;
        self.active_vertices = 0;
        self.active_edges = 0;
        self.trajectory.clear();
    }

    /// Current online summary (the state before/after any number of
    /// operations; not appended to the trajectory).
    pub fn summary(&self) -> IncrementalSummary {
        IncrementalSummary {
            active_vertex_count: self.active_vertices,
            active_edge_count: self.active_edges,
            component_count: self.components,
            largest_component: (self.largest_size != 0).then_some(self.largest_identity),
            largest_component_size: self.largest_size,
            raw_second_moment: self.second_moment,
        }
    }

    /// The trajectory recorded since the last `prepare`/`reset`: one
    /// `(step, summary)` entry per accepted operation, in acceptance order.
    #[inline]
    pub const fn trajectory(&self) -> &ProcessTrajectory {
        &self.trajectory
    }

    /// Size of the component containing `vertex`, or `None` when the vertex
    /// is inactive or outside the prepared domain.
    ///
    /// The scoring primitive for future rule-driven processes (G2 candidate
    /// selection): a read-only root walk, no path compression (chains are
    /// already compressed by the unions that built them).
    pub fn component_size(&self, vertex: crate::VertexId) -> Option<usize> {
        let index = vertex.index();
        if index >= self.prepared_vertices || self.vertex_active[index] == 0 {
            return None;
        }
        let mut root = index;
        while self.parent[root] != root {
            root = self.parent[root];
        }
        Some(self.size[root])
    }

    /// Activate one vertex and union it with every live-capable activated
    /// incident edge. Returns the online summary after the activation.
    ///
    /// `step` is the caller-defined process coordinate recorded with the
    /// summary (see the module documentation).
    pub fn activate_vertex<G: UndirectedGraphView>(
        &mut self,
        graph: &G,
        vertex: crate::VertexId,
        step: usize,
    ) -> Result<IncrementalSummary, IncrementalError> {
        self.check_graph_sizes(graph)?;
        let index = vertex.index();
        if index >= self.prepared_vertices {
            return Err(IncrementalError::VertexIndexOutOfBounds {
                requested: index,
                vertices: self.prepared_vertices,
            });
        }
        if self.vertex_active[index] != 0 {
            return Err(IncrementalError::VertexAlreadyActive { vertex: index });
        }

        self.vertex_active[index] = 1;
        self.parent[index] = index;
        self.size[index] = 1;
        self.canonical[index] = index;
        self.components += 1;
        self.second_moment += 1;
        self.active_vertices += 1;
        self.note_largest(1, index);

        for incidence in graph.incidences(vertex) {
            let edge = incidence.edge.index();
            // EDGE_INACTIVE: not activated by the process; EDGE_LIVE: already
            // counted (the self-loop second incidence lands here).
            if self.edge_state[edge] != EDGE_ACTIVATED {
                continue;
            }
            let neighbor = incidence.neighbor.index();
            if self.vertex_active[neighbor] == 0 {
                continue;
            }
            self.edge_state[edge] = EDGE_LIVE;
            self.active_edges += 1;
            if let Some((left_size, right_size, merged_identity)) =
                self.union_indices(index, neighbor)
            {
                self.components -= 1;
                self.second_moment += 2 * (left_size as u128) * (right_size as u128);
                self.note_largest(left_size + right_size, merged_identity);
            }
        }

        let summary = self.summary();
        self.trajectory.push(step, summary);
        Ok(summary)
    }

    /// Activate one physical edge; union its endpoints iff both are already
    /// active (otherwise the edge waits and goes live at the second endpoint's
    /// activation). Returns the online summary after the activation.
    ///
    /// `step` is the caller-defined process coordinate recorded with the
    /// summary (see the module documentation).
    pub fn activate_edge<G: UndirectedGraphView>(
        &mut self,
        graph: &G,
        edge: crate::EdgeId,
        step: usize,
    ) -> Result<IncrementalSummary, IncrementalError> {
        self.check_graph_sizes(graph)?;
        let index = edge.index();
        if index >= self.prepared_edges {
            return Err(IncrementalError::EdgeIndexOutOfBounds {
                requested: index,
                edges: self.prepared_edges,
            });
        }
        if self.edge_state[index] != EDGE_INACTIVE {
            return Err(IncrementalError::EdgeAlreadyActive { edge: index });
        }

        let [left, right] = graph.edge_endpoints(edge);
        self.edge_state[index] = EDGE_ACTIVATED;
        if self.vertex_active[left.index()] != 0 && self.vertex_active[right.index()] != 0 {
            self.edge_state[index] = EDGE_LIVE;
            self.active_edges += 1;
            if let Some((left_size, right_size, merged_identity)) =
                self.union_indices(left.index(), right.index())
            {
                self.components -= 1;
                self.second_moment += 2 * (left_size as u128) * (right_size as u128);
                self.note_largest(left_size + right_size, merged_identity);
            }
        }

        let summary = self.summary();
        self.trajectory.push(step, summary);
        Ok(summary)
    }

    fn check_graph_sizes<G: UndirectedGraphView>(&self, graph: &G) -> Result<(), IncrementalError> {
        if graph.vertex_count() != self.prepared_vertices
            || graph.edge_count() != self.prepared_edges
        {
            return Err(IncrementalError::GraphSizeMismatch {
                vertices: self.prepared_vertices,
                edges: self.prepared_edges,
                graph_vertices: graph.vertex_count(),
                graph_edges: graph.edge_count(),
            });
        }
        Ok(())
    }

    /// Update the selected-largest bookkeeping with one newly created
    /// component of `merged_size` and canonical identity `merged_identity`.
    ///
    /// Correct online because a merge only removes the old best when the
    /// merged component is strictly larger (every component has size >= 1),
    /// so the equal-size case is always a tie-break against an untouched
    /// incumbent — exactly the F3 selection rule.
    fn note_largest(&mut self, merged_size: usize, merged_identity: usize) {
        if merged_size > self.largest_size {
            self.largest_size = merged_size;
            self.largest_identity = merged_identity;
        } else if merged_size == self.largest_size {
            self.largest_identity = self.largest_identity.min(merged_identity);
        }
    }

    /// Union by size with path compression; returns the two pre-merge sizes
    /// and the merged canonical identity when the roots were distinct.
    fn union_indices(&mut self, left: usize, right: usize) -> Option<(usize, usize, usize)> {
        let left_root = self.find_root(left);
        let right_root = self.find_root(right);
        if left_root == right_root {
            return None; // includes self-loops: no component change
        }
        let (big, small) = if self.size[left_root] >= self.size[right_root] {
            (left_root, right_root)
        } else {
            (right_root, left_root)
        };
        let (big_size, small_size) = (self.size[big], self.size[small]);
        self.parent[small] = big;
        self.size[big] = big_size + small_size;
        self.canonical[big] = self.canonical[big].min(self.canonical[small]);
        Some((big_size, small_size, self.canonical[big]))
    }

    fn find_root(&mut self, node: usize) -> usize {
        let mut root = node;
        while self.parent[root] != root {
            root = self.parent[root];
        }
        let mut current = node;
        while self.parent[current] != root {
            let next = self.parent[current];
            self.parent[current] = root;
            current = next;
        }
        root
    }
}

fn reserve_len<T>(storage: &mut Vec<T>, len: usize) -> Result<(), TryReserveError> {
    if storage.capacity() < len {
        storage.try_reserve_exact(len - storage.len())?;
    }
    Ok(())
}

/// Failure of an incremental-process operation.
#[derive(Debug)]
pub enum IncrementalError {
    /// A workspace buffer could not be reserved.
    Capacity(TryReserveError),
    /// The driving graph's sizes differ from the prepared sizes.
    GraphSizeMismatch {
        vertices: usize,
        edges: usize,
        graph_vertices: usize,
        graph_edges: usize,
    },
    /// The vertex index lies outside the prepared domain.
    VertexIndexOutOfBounds { requested: usize, vertices: usize },
    /// The edge index lies outside the prepared domain.
    EdgeIndexOutOfBounds { requested: usize, edges: usize },
    /// The vertex was already active; the operation was refused without any
    /// state change.
    VertexAlreadyActive { vertex: usize },
    /// The edge was already activated; the operation was refused without any
    /// state change.
    EdgeAlreadyActive { edge: usize },
}

impl From<TryReserveError> for IncrementalError {
    fn from(source: TryReserveError) -> Self {
        Self::Capacity(source)
    }
}

impl fmt::Display for IncrementalError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Capacity(_) => {
                write!(formatter, "could not reserve incremental-process buffers")
            }
            Self::GraphSizeMismatch {
                vertices,
                edges,
                graph_vertices,
                graph_edges,
            } => write!(
                formatter,
                "graph has {graph_vertices} vertices and {graph_edges} edges, \
                 workspace was prepared for {vertices} and {edges}"
            ),
            Self::VertexIndexOutOfBounds {
                requested,
                vertices,
            } => write!(
                formatter,
                "vertex index {requested} is outside the prepared 0..{vertices} domain"
            ),
            Self::EdgeIndexOutOfBounds { requested, edges } => write!(
                formatter,
                "edge index {requested} is outside the prepared 0..{edges} domain"
            ),
            Self::VertexAlreadyActive { vertex } => {
                write!(formatter, "vertex {vertex} is already active")
            }
            Self::EdgeAlreadyActive { edge } => {
                write!(formatter, "edge {edge} is already activated")
            }
        }
    }
}

impl std::error::Error for IncrementalError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Capacity(source) => Some(source),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{build_chain, Bond, BondType, CsrLattice, GraphView, VertexId};
    use rand::{RngExt, SeedableRng};
    use rand_xoshiro::Xoshiro256PlusPlus;

    /// 3 vertices, 4 edges: self-loop (0,0), parallel pair (0,1) x2, (1,2).
    fn multigraph() -> CsrLattice {
        CsrLattice::try_from_edges(
            3,
            vec![
                Bond::new(0, 0, BondType::Generic, 1.0),
                Bond::new(0, 1, BondType::Generic, 1.0),
                Bond::new(0, 1, BondType::Generic, 1.0),
                Bond::new(1, 2, BondType::Generic, 1.0),
            ],
        )
        .expect("multigraph is valid")
    }

    #[test]
    fn hand_worked_chain_with_edges_first() {
        let graph = build_chain(4, false);
        let mut workspace = IncrementalUnionWorkspace::new();
        workspace.prepare(4, 3).unwrap();

        let empty = workspace.summary();
        assert_eq!(empty.component_count, 0);
        assert_eq!(empty.largest_component, None);
        assert_eq!(empty.largest_component_size, 0);

        // Edge 0 (0-1) waits: no endpoint is active.
        let summary = workspace
            .activate_edge(&graph, graph.edge_id(0).unwrap(), 0)
            .unwrap();
        assert_eq!(
            summary,
            IncrementalSummary {
                active_vertex_count: 0,
                active_edge_count: 0,
                component_count: 0,
                largest_component: None,
                largest_component_size: 0,
                raw_second_moment: 0,
            }
        );

        // Vertex 0: singleton; the waiting edge still has an inactive peer.
        let summary = workspace
            .activate_vertex(&graph, graph.vertex_id(0).unwrap(), 1)
            .unwrap();
        assert_eq!(summary.component_count, 1);
        assert_eq!(summary.largest_component, Some(0));
        assert_eq!(summary.active_edge_count, 0);
        assert_eq!(summary.raw_second_moment, 1);

        // Vertex 1: the edge goes live and merges {0} + {1}.
        let summary = workspace
            .activate_vertex(&graph, graph.vertex_id(1).unwrap(), 2)
            .unwrap();
        assert_eq!(summary.active_vertex_count, 2);
        assert_eq!(summary.active_edge_count, 1);
        assert_eq!(summary.component_count, 1);
        assert_eq!(summary.largest_component, Some(0));
        assert_eq!(summary.largest_component_size, 2);
        assert_eq!(summary.raw_second_moment, 4);

        assert_eq!(workspace.trajectory().steps(), &[0, 1, 2]);
        assert_eq!(workspace.trajectory().summaries().len(), 3);
        assert_eq!(workspace.trajectory().summaries()[2], summary);
    }

    #[test]
    fn self_loop_counts_once_and_parallel_edges_count_separately() {
        let graph = multigraph();
        let mut workspace = IncrementalUnionWorkspace::new();
        workspace.prepare(3, 4).unwrap();

        // Activate the self-loop and both parallel edges before any vertex.
        for edge in 0..3 {
            let summary = workspace
                .activate_edge(&graph, graph.edge_id(edge).unwrap(), edge)
                .unwrap();
            assert_eq!(summary.active_edge_count, 0, "edges wait for endpoints");
        }
        let s0 = workspace
            .activate_vertex(&graph, graph.vertex_id(0).unwrap(), 10)
            .unwrap();
        // The self-loop goes live once at its only endpoint (its second
        // incidence sees EDGE_LIVE); the parallel edges wait for vertex 1.
        assert_eq!(s0.active_vertex_count, 1);
        assert_eq!(s0.active_edge_count, 1);
        assert_eq!(s0.component_count, 1);
        assert_eq!(s0.raw_second_moment, 1);

        let s1 = workspace
            .activate_vertex(&graph, graph.vertex_id(1).unwrap(), 11)
            .unwrap();
        // Both parallel edges go live; only the first of them merges.
        assert_eq!(s1.active_vertex_count, 2);
        assert_eq!(s1.active_edge_count, 3);
        assert_eq!(s1.component_count, 1);
        assert_eq!(s1.largest_component_size, 2);
        assert_eq!(s1.raw_second_moment, 4);

        let s2 = workspace
            .activate_vertex(&graph, graph.vertex_id(2).unwrap(), 12)
            .unwrap();
        // Edge 3 (1-2) was never activated: {0,1} and {2} remain.
        assert_eq!(s2.active_vertex_count, 3);
        assert_eq!(s2.active_edge_count, 3);
        assert_eq!(s2.component_count, 2);
        assert_eq!(s2.largest_component, Some(0));
        assert_eq!(s2.raw_second_moment, 4 + 1);
    }

    #[test]
    fn equal_size_ties_resolve_by_lowest_identity() {
        let graph = build_chain(4, false); // edges (0,1), (1,2), (2,3)
        let mut workspace = IncrementalUnionWorkspace::new();
        workspace.prepare(4, 3).unwrap();

        // Build the identity-2 pair first: {2,3} becomes the incumbent
        // largest with identity 2.
        workspace
            .activate_vertex(&graph, graph.vertex_id(2).unwrap(), 0)
            .unwrap();
        workspace
            .activate_vertex(&graph, graph.vertex_id(3).unwrap(), 1)
            .unwrap();
        let first_pair = workspace
            .activate_edge(&graph, graph.edge_id(2).unwrap(), 2)
            .unwrap();
        assert_eq!(first_pair.largest_component_size, 2);
        assert_eq!(first_pair.largest_component, Some(2));

        // Now the equal-size identity-0 pair: the tie must REPLACE the
        // selected identity with the lower one.
        workspace
            .activate_vertex(&graph, graph.vertex_id(0).unwrap(), 3)
            .unwrap();
        workspace
            .activate_vertex(&graph, graph.vertex_id(1).unwrap(), 4)
            .unwrap();
        let retied = workspace
            .activate_edge(&graph, graph.edge_id(0).unwrap(), 5)
            .unwrap();
        assert_eq!(retied.largest_component_size, 2);
        assert_eq!(retied.largest_component, Some(0));
        assert_eq!(retied.component_count, 2);
        assert_eq!(retied.raw_second_moment, 8);
    }

    #[test]
    fn refusals_leave_state_untouched() {
        let graph = build_chain(3, false);
        let mut workspace = IncrementalUnionWorkspace::new();
        workspace.prepare(3, 2).unwrap();

        workspace
            .activate_vertex(&graph, graph.vertex_id(0).unwrap(), 0)
            .unwrap();
        let after_first = workspace.summary();
        assert!(matches!(
            workspace.activate_vertex(&graph, graph.vertex_id(0).unwrap(), 1),
            Err(IncrementalError::VertexAlreadyActive { vertex: 0 })
        ));
        assert_eq!(workspace.summary(), after_first);

        workspace
            .activate_edge(&graph, graph.edge_id(0).unwrap(), 2)
            .unwrap();
        let after_edge = workspace.summary();
        assert!(matches!(
            workspace.activate_edge(&graph, graph.edge_id(0).unwrap(), 3),
            Err(IncrementalError::EdgeAlreadyActive { edge: 0 })
        ));
        assert_eq!(workspace.summary(), after_edge);
        assert_eq!(workspace.trajectory().len(), 2, "only accepted ops record");

        // A forged out-of-domain vertex (an ID from a larger graph) is a
        // typed refusal, not a panic and not a state change.
        assert!(matches!(
            workspace.activate_vertex(&graph, VertexId::from_valid_index(3), 0),
            Err(IncrementalError::VertexIndexOutOfBounds { .. })
        ));
        assert!(matches!(
            workspace.activate_edge(&graph, crate::EdgeId::from_valid_index(2), 0),
            Err(IncrementalError::EdgeIndexOutOfBounds { .. })
        ));
        assert_eq!(workspace.summary(), after_edge);

        // A graph whose sizes differ from the prepared sizes is refused.
        let small = build_chain(2, false);
        assert!(matches!(
            workspace.activate_vertex(&small, small.vertex_id(0).unwrap(), 0),
            Err(IncrementalError::GraphSizeMismatch { .. })
        ));
        assert_eq!(workspace.summary(), after_edge);
    }

    #[test]
    fn reset_and_reprepare_leave_no_stale_state() {
        let graph = multigraph();
        let mut workspace = IncrementalUnionWorkspace::new();
        workspace.prepare(3, 4).unwrap();
        for index in 0..3 {
            workspace
                .activate_vertex(&graph, graph.vertex_id(index).unwrap(), index)
                .unwrap();
        }
        assert_eq!(workspace.summary().active_vertex_count, 3);

        workspace.reset();
        assert_eq!(
            workspace.summary(),
            IncrementalSummary {
                active_vertex_count: 0,
                active_edge_count: 0,
                component_count: 0,
                largest_component: None,
                largest_component_size: 0,
                raw_second_moment: 0,
            }
        );
        assert!(workspace.trajectory().is_empty());

        // Shrink to a smaller topology: stale slots from the old size must be
        // gone (index 2 is out of domain, vertex 0 inactive again).
        let small = build_chain(2, false);
        workspace.prepare(2, 1).unwrap();
        assert!(matches!(
            workspace.activate_vertex(&small, VertexId::from_valid_index(2), 0),
            Err(IncrementalError::VertexIndexOutOfBounds { .. })
        ));
        assert_eq!(workspace.summary().component_count, 0);
        workspace
            .activate_vertex(&small, small.vertex_id(0).unwrap(), 0)
            .unwrap();
        assert_eq!(workspace.summary().largest_component, Some(0));
    }

    #[test]
    fn component_size_reports_live_components() {
        let graph = build_chain(4, false);
        let mut workspace = IncrementalUnionWorkspace::new();
        workspace.prepare(4, 3).unwrap();
        assert_eq!(
            workspace.component_size(graph.vertex_id(0).unwrap()),
            None,
            "inactive vertices have no component"
        );
        for index in 0..3 {
            workspace
                .activate_vertex(&graph, graph.vertex_id(index).unwrap(), index)
                .unwrap();
        }
        assert_eq!(
            workspace.component_size(graph.vertex_id(0).unwrap()),
            Some(1),
            "no edges activated: three isolates"
        );
        workspace
            .activate_edge(&graph, graph.edge_id(0).unwrap(), 3)
            .unwrap();
        workspace
            .activate_edge(&graph, graph.edge_id(1).unwrap(), 4)
            .unwrap();
        for index in 0..3 {
            assert_eq!(
                workspace.component_size(graph.vertex_id(index).unwrap()),
                Some(3)
            );
        }
        assert_eq!(workspace.component_size(graph.vertex_id(3).unwrap()), None);
    }

    #[test]
    fn random_interleavings_converge_to_the_same_terminal_summary() {
        let graph = multigraph();
        let mut workspace = IncrementalUnionWorkspace::new();
        workspace.prepare(3, 4).unwrap();
        let mut reference = None;
        for seed in 0..16_u64 {
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
            workspace.reset();
            // Random interleaving of all vertex and all edge activations.
            let mut ops: Vec<Option<usize>> =
                (0..3).map(Some).chain((0..4).map(|_| None)).collect();
            for index in (1..ops.len()).rev() {
                let swap_with = rng.random_range(0..=index);
                ops.swap(index, swap_with);
            }
            let mut edges = 0_usize;
            for op in ops {
                match op {
                    Some(vertex) => workspace
                        .activate_vertex(&graph, graph.vertex_id(vertex).unwrap(), vertex)
                        .unwrap(),
                    None => {
                        let edge = edges;
                        edges += 1;
                        workspace
                            .activate_edge(&graph, graph.edge_id(edge).unwrap(), edge)
                            .unwrap()
                    }
                };
            }
            let terminal = workspace.summary();
            if let Some(reference) = reference {
                assert_eq!(terminal, reference, "seed {seed} diverged");
            } else {
                reference = Some(terminal);
            }
        }
        let terminal = reference.expect("at least one seed");
        // Fully activated: one component of all 3 vertices, all 4 edges live.
        assert_eq!(terminal.active_vertex_count, 3);
        assert_eq!(terminal.active_edge_count, 4);
        assert_eq!(terminal.component_count, 1);
        assert_eq!(terminal.largest_component, Some(0));
        assert_eq!(terminal.largest_component_size, 3);
        assert_eq!(terminal.raw_second_moment, 9);
    }
}
