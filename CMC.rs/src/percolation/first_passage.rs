//! F1: weighted shortest-path layer over the topology substrate (roadmap
//! v2 §3 F1).
//!
//! [`first_passage_run`] is a single-source, single-pass Dijkstra with a
//! binary heap, generic over [`UndirectedGraphView`]: topology comes from
//! the view (parallel bonds are distinct dense edges, self-loops are legal),
//! while the edge weights enter as a caller-owned slice of one `f64` per
//! physical edge in dense edge order — for a [`CsrLattice`](crate::CsrLattice)
//! that slice is exactly `lattice.edges[e].weight`. Weights are data, not
//! topology: first-passage realizations redraw them per call (F2's passage
//! times are frozen fields in exactly the [`invasion_run`](crate::invasion_run)
//! threshold-field sense), so the workspace is prepared per (V, E) sizes and
//! every call may carry a different weight field.
//!
//! ## Weight contract
//!
//! - **Non-negative and finite**, validated in full before the pass starts
//!   (a rejected run never touches the workspace): negative weights are the
//!   typed [`FirstPassageError::NegativeWeight`] naming the edge, NaN and
//!   infinities are [`FirstPassageError::WeightNotFinite`]. Zero weights are
//!   ALLOWED — non-negativity is Dijkstra's actual requirement — so
//!   zero-weight edges and zero-weight cycles are correct inputs.
//! - **Parallel edges are independent relaxations**; the minimum-weight bond
//!   between two vertices wins naturally (strict-improvement updates only),
//!   pinned by test. Self-loops can never strictly improve and are ignored.
//!
//! ## Distances, unreachable vertices, and the tie policy
//!
//! Distances are one `f64` per dense vertex. The unreachable sentinel is
//! `f64::INFINITY` (the array's reset value), so reachability is the derived
//! `distance.is_finite()` — no separate flags; a reachable vertex always has
//! a finite distance because validated weights are finite and any path sum
//! that would overflow to infinity reads as unreachable (a saturation
//! convention for weights near `f64::MAX`, documented rather than guarded).
//!
//! Heap entries are keyed `(distance, dense vertex id)` with `distance`
//! under `f64::total_cmp` (`-0.0 < +0.0`: distinct keys, not a tie) and the
//! vertex id as the tie key. The id exists because Rust's `BinaryHeap`
//! specifies NO order among equal keys: without a total tie key, which of
//! two equal-distance frontier entries settles first — and therefore which
//! of several equal-cost predecessors a vertex keeps — would depend on heap
//! internals. Lowest-id-first makes the pop order a total order (the same
//! tie discipline as the invasion frontier, and a 16-byte heap entry — no
//! separate counter field), so the relaxation order and the predecessor
//! output are deterministic functions of the input, pinned by hand-case
//! test. Keys are unique anyway: pushes happen only on strict improvement,
//! so one vertex never has two entries at the same distance.
//!
//! ## Workspace and F2 hand-off
//!
//! [`FirstPassageWorkspace`] reserves the distance array, the optional
//! predecessor arrays, and the heap (pushes are bounded by E + 1: the source
//! plus at most one per physical edge — the non-negative weights make
//! relaxing an edge from its later-settling endpoint never strictly improve)
//! once per topology size; repeated single-source calls on a correctly sized
//! workspace perform zero steady-state allocation. Predecessors (vertex and
//! physical edge per reached vertex, for geodesic reconstruction in F2/F3)
//! are opt-in through [`GeodesicTracking`] and live in the workspace beside
//! the distances. One call is one deterministic computation, like
//! [`analyze_winding`](crate::analyze_winding): no MonteCarlo adapter and no
//! scheduler integration at this layer — F2 builds its observables on top.

use std::cmp::Ordering;
use std::cmp::Reverse;
use std::collections::BinaryHeap;

use super::process::ProcessError;
use crate::{UndirectedGraphView, VertexId};

/// Predecessor-array sentinel: the source, unreached vertices, and (on a
/// tracking-off run) everything.
pub const NO_PREDECESSOR: usize = usize::MAX;

/// Whether a run also records the shortest-path tree.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum GeodesicTracking {
    /// Distances only: `prepare` rewrites every predecessor entry to
    /// `NO_PREDECESSOR` before the pass, so the arrays read as empty
    /// regardless of earlier runs (review follow-up: the previous
    /// "stale from the last tracking-on run" wording was wrong —
    /// preparation clears unconditionally).
    #[default]
    Off,
    /// Record predecessor vertex and physical edge per reached vertex
    /// ([`NO_PREDECESSOR`] elsewhere) for geodesic reconstruction.
    Predecessors,
}

/// Scalar outcome of one single-source pass.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct FirstPassageOutcome {
    /// The dense source vertex index of the pass.
    pub source: usize,
    /// Reached vertices, the source included (`distance.is_finite()`).
    pub reached: usize,
}

/// One heap entry: `(distance, dense vertex id)`.
///
/// Ordering is the documented tie policy — `total_cmp` on the distance,
/// then the vertex id ascending — so the pop order is a total order even
/// under massive ties (zero-weight edges). Rejections guarantee finite
/// weights, keeping the ordering total.
#[derive(Clone, Copy, Debug)]
struct PassageEntry {
    distance: f64,
    vertex: usize,
}

impl PartialEq for PassageEntry {
    fn eq(&self, other: &Self) -> bool {
        matches!(self.partial_cmp(other), Some(Ordering::Equal))
    }
}

impl Eq for PassageEntry {}

impl PartialOrd for PassageEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for PassageEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        self.distance
            .total_cmp(&other.distance)
            .then_with(|| self.vertex.cmp(&other.vertex))
    }
}

/// Reusable weighted-shortest-path workspace: the distance array, the
/// optional predecessor arrays, and the binary heap.
///
/// `prepare` once per topology size ([`first_passage_run`] re-prepares
/// idempotently); every buffer is reserved up front, so a correctly sized
/// workspace performs zero steady-state allocation per single-source call.
#[derive(Debug, Default)]
pub struct FirstPassageWorkspace {
    /// Min-distance frontier ([`Reverse`] pops the smallest key first).
    heap: BinaryHeap<Reverse<PassageEntry>>,
    /// One `f64` per dense vertex; `f64::INFINITY` = unreachable.
    dist: Vec<f64>,
    /// Predecessor vertex per dense vertex ([`NO_PREDECESSOR`] = none);
    /// written only under [`GeodesicTracking::Predecessors`].
    predecessor_vertex: Vec<usize>,
    /// Predecessor physical edge per dense vertex; written with the vertex.
    predecessor_edge: Vec<usize>,
    prepared_vertices: usize,
    prepared_edges: usize,
}

impl FirstPassageWorkspace {
    /// An unprepared workspace.
    pub fn new() -> Self {
        Self::default()
    }

    /// Reserve every buffer for a topology of exactly `vertices` vertices
    /// and `edges` physical edges. Re-preparing with different sizes is
    /// supported and leaves no stale state.
    pub fn prepare(&mut self, vertices: usize, edges: usize) -> Result<(), FirstPassageError> {
        // Clear BEFORE reserving: `try_reserve_exact` budgets for
        // `len + additional`, so reserving with the previous run's entries
        // still in the heap would demand nearly twice the capacity once and
        // allocate. Drained first, `edges + 1` (the source push plus at
        // most one push per physical edge) is exactly the budget every time.
        self.heap.clear();
        self.heap
            .try_reserve_exact(edges + 1)
            .map_err(|source| FirstPassageError::Capacity(ProcessError::from(source)))?;
        reserve_exact(&mut self.dist, vertices, f64::INFINITY)?;
        reserve_exact(&mut self.predecessor_vertex, vertices, NO_PREDECESSOR)?;
        reserve_exact(&mut self.predecessor_edge, vertices, NO_PREDECESSOR)?;
        self.prepared_vertices = vertices;
        self.prepared_edges = edges;
        Ok(())
    }

    /// Prepared topology sizes.
    #[inline]
    pub const fn prepared_sizes(&self) -> (usize, usize) {
        (self.prepared_vertices, self.prepared_edges)
    }

    /// Distances of the last run, one `f64` per dense vertex:
    /// `0.0` at the source, `f64::INFINITY` where unreachable.
    #[inline]
    pub fn distances(&self) -> &[f64] {
        &self.dist
    }

    /// Whether the vertex was reached by the last run.
    #[inline]
    pub fn is_reachable(&self, index: usize) -> bool {
        self.dist[index].is_finite()
    }

    /// Predecessor vertex per dense vertex of the last run
    /// ([`NO_PREDECESSOR`] at the source and unreached vertices). Meaningful
    /// only when the last run used [`GeodesicTracking::Predecessors`].
    #[inline]
    pub fn predecessor_vertices(&self) -> &[usize] {
        &self.predecessor_vertex
    }

    /// Predecessor physical edge per dense vertex, written with
    /// [`Self::predecessor_vertices`] (same tracking requirement).
    #[inline]
    pub fn predecessor_edges(&self) -> &[usize] {
        &self.predecessor_edge
    }

    /// Reset for the next run, retaining all capacity.
    fn reset(&mut self, tracking: GeodesicTracking) {
        self.heap.clear();
        self.dist.fill(f64::INFINITY);
        if tracking == GeodesicTracking::Predecessors {
            self.predecessor_vertex.fill(NO_PREDECESSOR);
            self.predecessor_edge.fill(NO_PREDECESSOR);
        }
    }
}

/// Reserve `len` logical entries (rewriting the default value so a shrink
/// leaves no stale state, capacity retained).
fn reserve_exact<T: Clone>(
    storage: &mut Vec<T>,
    len: usize,
    default: T,
) -> Result<(), FirstPassageError> {
    if storage.capacity() < len {
        storage.clear();
        storage
            .try_reserve_exact(len)
            .map_err(|source| FirstPassageError::Capacity(ProcessError::from(source)))?;
    }
    storage.clear();
    storage.resize(len, default);
    Ok(())
}

/// One single-source Dijkstra pass (the reusable, RNG-free core).
///
/// `weights` holds one `f64` per physical edge in dense edge order — for a
/// `CsrLattice`, `lattice.edges[e].weight`. Every entry must be finite and
/// non-negative (zero allowed); violations are typed rejections raised
/// during validation, before the workspace is touched, so a rejected call
/// leaves the previous run's results intact. `source` is a dense vertex
/// index. With [`GeodesicTracking::Predecessors`] the workspace's
/// predecessor arrays hold the deterministic shortest-path tree (one vertex
/// and one physical edge per reached vertex). Complexity
/// O((V + E) log E); a correctly sized workspace allocates nothing.
pub fn first_passage_run<G: UndirectedGraphView>(
    graph: &G,
    weights: &[f64],
    source: usize,
    tracking: GeodesicTracking,
    workspace: &mut FirstPassageWorkspace,
) -> Result<FirstPassageOutcome, FirstPassageError> {
    let vertices = graph.vertex_count();
    let edges = graph.edge_count();
    if vertices == 0 {
        return Err(FirstPassageError::EmptyGraph);
    }
    if weights.len() != edges {
        return Err(FirstPassageError::WeightLengthMismatch {
            expected: edges,
            actual: weights.len(),
        });
    }
    for (edge, &value) in weights.iter().enumerate() {
        if !value.is_finite() {
            return Err(FirstPassageError::WeightNotFinite { edge, value });
        }
        if value < 0.0 {
            return Err(FirstPassageError::NegativeWeight { edge, value });
        }
    }
    if source >= vertices {
        return Err(FirstPassageError::SourceDomain { source, vertices });
    }
    workspace.prepare(vertices, edges)?;
    workspace.reset(tracking);

    workspace.dist[source] = 0.0;
    workspace.heap.push(Reverse(PassageEntry {
        distance: 0.0,
        vertex: source,
    }));

    while let Some(Reverse(entry)) = workspace.heap.pop() {
        // Stale-entry check: pushes happen only on strict improvement, so a
        // vertex's pushed distances are strictly decreasing and its final
        // (smallest) entry is the only one that can pop with an equal
        // distance — each vertex settles exactly once and no settled flags
        // are needed. Entries with a larger distance were superseded.
        if entry.distance > workspace.dist[entry.vertex] {
            continue;
        }
        let u = VertexId::from_valid_index(entry.vertex);
        for incidence in graph.incidences(u) {
            let v = incidence.neighbor.index();
            let candidate = entry.distance + weights[incidence.edge.index()];
            if candidate < workspace.dist[v] {
                workspace.dist[v] = candidate;
                if tracking == GeodesicTracking::Predecessors {
                    workspace.predecessor_vertex[v] = entry.vertex;
                    workspace.predecessor_edge[v] = incidence.edge.index();
                }
                workspace.heap.push(Reverse(PassageEntry {
                    distance: candidate,
                    vertex: v,
                }));
            }
        }
    }

    let reached = workspace
        .dist
        .iter()
        .filter(|&&distance| distance.is_finite())
        .count();
    Ok(FirstPassageOutcome { source, reached })
}

/// Typed rejections for the weighted shortest-path layer.
#[derive(Debug)]
pub enum FirstPassageError {
    /// A graph with zero vertices.
    EmptyGraph,
    /// A weight field with the wrong length (one entry per physical edge).
    WeightLengthMismatch {
        /// The graph's physical edge count.
        expected: usize,
        actual: usize,
    },
    /// A NaN or infinite weight (`f64::INFINITY` included) at `edge`.
    WeightNotFinite {
        /// The dense physical-edge index of the offending weight.
        edge: usize,
        value: f64,
    },
    /// A negative weight at `edge` — Dijkstra requires non-negative.
    NegativeWeight {
        /// The dense physical-edge index of the offending weight.
        edge: usize,
        value: f64,
    },
    /// A source outside the graph's dense index range.
    SourceDomain { source: usize, vertices: usize },
    /// A buffer reservation failure.
    Capacity(ProcessError),
}

impl core::fmt::Display for FirstPassageError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::EmptyGraph => {
                write!(
                    formatter,
                    "shortest-path analysis needs at least one vertex"
                )
            }
            Self::WeightLengthMismatch { expected, actual } => write!(
                formatter,
                "weight field has {actual} entries, expected {expected} (one per physical edge)"
            ),
            Self::WeightNotFinite { edge, value } => write!(
                formatter,
                "weight {value} at edge {edge} is not finite (non-finite weights break Dijkstra)"
            ),
            Self::NegativeWeight { edge, value } => write!(
                formatter,
                "weight {value} at edge {edge} is negative (Dijkstra requires non-negative \
                 weights)"
            ),
            Self::SourceDomain { source, vertices } => {
                write!(
                    formatter,
                    "source {source} outside the graph (0..{vertices})"
                )
            }
            Self::Capacity(source) => {
                write!(
                    formatter,
                    "could not prepare shortest-path buffers: {source}"
                )
            }
        }
    }
}

impl std::error::Error for FirstPassageError {}
