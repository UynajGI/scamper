//! G3: invasion percolation — frozen-threshold minimum-frontier growth
//! (roadmap §13).
//!
//! One realization freezes a threshold field — one `f64` per vertex (site
//! mode) or per physical edge (bond mode) — and grows an invasion cluster
//! from explicit seed vertices through the **minimum-threshold frontier**:
//! at every step the frontier's smallest-threshold entry is accepted, its
//! un-invaded endpoint joins the cluster, and the entry's physical edge
//! becomes the invasion-tree edge. Non-trapping growth continues until the
//! frontier empties (the whole component of the seeds is invaded); the
//! optional Wilkinson-Willemsen trapping rule additionally requires an
//! entry's site to be connected to an explicit exit set through un-invaded
//! sites.
//!
//! ## Explicit policies (roadmap §13: "make tie and stop policies explicit")
//!
//! - **Tie policy**: equal thresholds resolve by the **lowest dense ID** —
//!   the vertex index in site mode, the physical-edge index in bond mode.
//!   The heap orders entries by `f64::total_cmp` on the threshold, then by
//!   dense ID ascending; NaN and infinite thresholds are typed rejections,
//!   so the heap ordering is total. `total_cmp` orders `-0.0 < +0.0`
//!   (distinct keys, no tie) — documented behavior, not an accident.
//! - **Stop policy** ([`InvasionStop`]): grow until the frontier empties
//!   (`GrowAll`); stop at the first moment both spanning vertices are
//!   invaded (`FirstSpanning { from, to }` — checked after the seeds are
//!   planted and after every accepted growth step, so two seeds that
//!   already satisfy it stop at step zero); stop after exactly `n`
//!   accepted growth steps (`Steps(n)`, zero is legal and accepts
//!   nothing). A frontier exhausted before the policy fires stops with
//!   [`InvasionStopReason::FrontierExhausted`].
//! - **Trapping policy** ([`InvasionTrapping`]): non-trapping by default;
//!   `Trapping { exit }` implements the Wilkinson-Willemsen rule — a
//!   candidate site is invadable iff it is connected to some **un-invaded
//!   exit vertex** through un-invaded sites (the defending phase; an
//!   invaded exit vertex no longer provides escape, and exit vertices are
//!   themselves invadable through the normal minimum rule). Because the
//!   defending phase only shrinks, connectivity to the exit is monotone
//!   decreasing: once a site is found trapped it is retired permanently,
//!   and each trapped verdict costs one DFS of the candidate's defending
//!   region with early exit at the first exit vertex hit — O(region) per
//!   check, worst-case O(V (V + E)) total. A split-maintaining union-find
//!   has the same asymptotic worst case (invading a site can split its
//!   class, forcing an O(class) rebuild), so the DFS-plus-monotone
//!   retirement cost model is the honest engineering choice, stated here
//!   rather than hidden.
//!
//! ## Semantics per mode
//!
//! - **Site mode** ([`InvasionMode::Site`]): `thresholds` has one entry
//!   per dense vertex. A vertex enters the frontier at most once — when
//!   the first invaded neighbor reaches it — carrying that neighbor's
//!   connecting edge; when accepted, that recorded edge is the tree edge.
//! - **Bond mode** ([`InvasionMode::Bond`]): `thresholds` has one entry
//!   per physical edge. An edge enters the frontier when exactly one
//!   endpoint is invaded and is accepted only if that still holds at pop
//!   time (entries whose far endpoint was invaded through another edge
//!   are stale and discarded). Parallel edges are independent frontier
//!   entries; self-loops never enter the frontier.
//!
//! The invasion graph is always a forest rooted at the seeds: growth only
//! ever adds an un-invaded vertex through one accepted edge, so the
//! accepted-edge count is exactly `invaded - seeds` and no cycle can form
//! (verified by gate).
//!
//! ## MST / Prim prefix relation (bond mode, non-trapping, unique
//! ## thresholds)
//!
//! The non-trapping bond rule accepts the minimum-threshold edge crossing
//! the cut (invaded, un-invaded) — exactly Prim's greedy step. With all
//! thresholds distinct the accepted-edge sequence of a `GrowAll`
//! realization is therefore **exactly the Prim growth order** from the
//! seed, and the terminal invasion tree is the unique MST of the seed's
//! connected component; a stopped run is the corresponding Prim prefix.
//! Verified against an independent test-side Prim implementation (gate in
//! `tests/physics/invasion_percolation.rs`). Literature anchor:
//! Sweeney & Middleton, "Minimal spanning trees at the percolation
//! threshold: a numerical calculation", Phys. Rev. E 88, 032129 (2013)
//! (arXiv:1307.0043, archived in `docs/ref/`), which states the
//! Prim/invasion-percolation correspondence (its reference [27]:
//! Barabási, Phys. Rev. Lett. 76, 3750 (1996)); the founding
//! Wilkinson-Willemsen paper (J. Phys. A 16, 3365 (1983)) predates arXiv
//! and is not archived. Evidence in `CMC.rs/PERCOLATION.md` (G3).
//!
//! ## RNG draw contract (adapter, bitwise reproducibility)
//!
//! [`InvasionPercolationMC`] draws the frozen field per realization: one
//! `f64` per element in **ascending dense order** — vertices `0..V` in
//! site mode, edges `0..E` in bond mode — uniform on `[0, 1)`. The core
//! [`invasion_run`] accepts caller-sampled thresholds (any finite values;
//! the law does not require uniformity) and is RNG-free.
//!
//! ## G1 substrate drive
//!
//! Every accepted operation drives the G1 [`IncrementalUnionWorkspace`]:
//! seeds call `activate_vertex(seed, 0)` and the `i`-th accepted growth
//! calls `activate_vertex(v, i)` then `activate_edge(e, i)` — the
//! roadmap's "trajectory steps = invasion times". The substrate's online
//! component summaries are thereby recorded along the invasion, exactly
//! as the G1 contract anticipated, with no allocation beyond its prepared
//! `V + E` trajectory budget. The invasion cluster is a forest by
//! construction, so the substrate's component bookkeeping is exact.
//!
//! ## C0 recording (adapter)
//!
//! One realization is ONE element-preserving vector sample plus scalars,
//! never per-step scalar samples:
//!
//! - `InvasionThresholdOrder[i]` (fixed length `steps`): the threshold
//!   accepted at growth step `i`; entries past the realization's growth
//!   length repeat the final accepted threshold (all `0.0` when the
//!   realization accepts nothing), so every realization keeps the same
//!   shape while the scalars carry the true lengths.
//! - `InvasionAcceptedTotal`: accepted growth steps.
//! - `InvasionInvadedVertices`: seeds + accepted growth steps.
//! - `InvasionMaxAcceptedThreshold`: the largest accepted threshold (the
//!   realization's invasion barrier); `0.0` when nothing is accepted.
//!
//! Configure `thermalization_sweeps = 0`; realizations are i.i.d.

use std::cmp::Ordering;
use std::cmp::Reverse;
use std::collections::BinaryHeap;

use carlo_rs::{Context, MonteCarlo};
use rand::{Rng, RngExt};
use rand_xoshiro::Xoshiro256PlusPlus;

use super::{IncrementalError, IncrementalUnionWorkspace, ProcessError, StampFlags};
use crate::{EdgeId, UndirectedGraphView, VertexId};

/// Which element carries the frozen per-realization threshold field.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum InvasionMode {
    /// One threshold per dense vertex; a vertex enters the frontier when
    /// the first invaded neighbor reaches it.
    Site,
    /// One threshold per physical edge; boundary edges (exactly one
    /// invaded endpoint) are the frontier entries.
    Bond,
}

/// Explicit stop policy for one invasion realization.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum InvasionStop {
    /// Grow until the frontier empties — with non-trapping growth this
    /// invades every site of the seeds' connected components.
    GrowAll,
    /// Stop at the first moment (after planting the seeds, and after every
    /// accepted growth step) where both `from` and `to` are invaded. With
    /// `from` among the seeds this is the invasion time of `to`. A
    /// frontier exhausted before both are invaded stops with
    /// [`InvasionStopReason::FrontierExhausted`].
    FirstSpanning {
        /// One spanning endpoint (dense vertex index).
        from: usize,
        /// The other spanning endpoint (dense vertex index).
        to: usize,
    },
    /// Stop after exactly `n` accepted growth steps (zero accepts
    /// nothing). A frontier exhausted earlier stops with
    /// [`InvasionStopReason::FrontierExhausted`].
    Steps(usize),
}

/// Explicit trapping policy for one invasion realization.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub enum InvasionTrapping {
    /// Every frontier entry is invadable; growth stops only at frontier
    /// exhaustion (the default).
    #[default]
    NonTrapping,
    /// Wilkinson-Willemsen trapping: `exit` is a non-empty set of dense
    /// exit (outlet) vertices; a frontier entry's site is invadable iff
    /// it is connected to some un-invaded exit vertex through un-invaded
    /// sites. Trapped sites are retired permanently (the defending phase
    /// only shrinks, so a trapped verdict is final).
    Trapping {
        /// Dense exit (outlet) vertex indices; non-empty, in-domain,
        /// duplicate-free.
        exit: Vec<usize>,
    },
}

/// Why an invasion realization stopped.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum InvasionStopReason {
    /// The frontier emptied (every reachable site invaded or retired).
    FrontierExhausted,
    /// The spanning policy fired: both endpoints invaded.
    SpanningStep,
    /// The step budget fired: exactly `n` growth steps accepted.
    StepBudget,
}

/// Scalar outcome of one invasion realization; the per-site record lives in
/// the workspace (`invaded_order`, `parent`, `accepted_edges`,
/// `accepted_thresholds`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InvasionOutcome {
    /// Planted seed vertices.
    pub seed_count: usize,
    /// Accepted growth steps (one vertex and one tree edge each).
    pub accepted_total: usize,
    /// Invaded vertices (`seed_count + accepted_total`).
    pub invaded_vertices: usize,
    /// Distinct frontier sites retired as trapped (trapping variant only).
    pub trapped_sites: usize,
    /// Largest accepted threshold; `0.0` when nothing was accepted.
    pub max_accepted_threshold: f64,
    /// Why the realization stopped.
    pub stop_reason: InvasionStopReason,
    /// Growth steps accepted when the stop fired (`0` for a step-zero
    /// spanning stop; the full accepted count under frontier exhaustion).
    pub stop_step: usize,
}

/// One frontier entry in the binary heap.
///
/// Ordering is the documented tie policy: threshold under `total_cmp`,
/// then dense ID ascending (vertex ID in site mode, edge ID in bond mode).
/// Non-finite thresholds are rejected before any entry is built.
#[derive(Clone, Copy, Debug)]
struct FrontierEntry {
    /// The frozen threshold this entry competes with.
    threshold: f64,
    /// Tie-break key: vertex index (site mode) or edge index (bond mode).
    id: usize,
    /// The un-invaded vertex acceptance would add.
    vertex: usize,
    /// The physical edge acceptance would add to the invasion tree.
    edge: usize,
}

impl FrontierEntry {
    fn key_cmp(&self, other: &Self) -> Ordering {
        self.threshold
            .total_cmp(&other.threshold)
            .then_with(|| self.id.cmp(&other.id))
    }
}

impl PartialEq for FrontierEntry {
    fn eq(&self, other: &Self) -> bool {
        matches!(self.partial_cmp(other), Some(Ordering::Equal))
    }
}

impl Eq for FrontierEntry {}

impl PartialOrd for FrontierEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for FrontierEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        self.key_cmp(other)
    }
}

/// Reusable invasion workspace: the minimum-threshold binary heap, invasion
/// state, trapping scratch, and the driven G1 substrate.
///
/// `prepare` once per topology size ([`invasion_run`] re-prepares
/// idempotently); every buffer is reserved up front, so a correctly sized
/// workspace performs **zero steady-state allocation** per realization.
/// Site-mode frontier pushes are bounded by `V`, bond-mode pushes by `E`,
/// so the heap reserves `max(V, E)`.
#[derive(Debug, Default)]
pub struct InvasionWorkspace {
    /// Min-threshold frontier ([`Reverse`] makes the binary heap pop the
    /// smallest key first).
    heap: BinaryHeap<Reverse<FrontierEntry>>,
    /// `usize::MAX` = un-invaded; a seed points to itself; every other
    /// invaded vertex points to its invasion parent.
    parent: Vec<usize>,
    /// Site mode: the vertex has entered the frontier.
    vertex_pushed: Vec<u8>,
    /// Bond mode: the edge has entered the frontier.
    edge_pushed: Vec<u8>,
    /// Trapping: the site was found trapped (permanent retirement).
    retired: Vec<u8>,
    /// Dense exit flags (trapping variant).
    exit_marks: Vec<u8>,
    /// Dense vertex indices in invasion order (seeds first, caller order).
    order: Vec<usize>,
    /// Accepted tree edge per growth step.
    accepted_edges: Vec<usize>,
    /// Accepted threshold per growth step.
    accepted_thresholds: Vec<f64>,
    /// Trapped-check DFS stack.
    dfs_stack: Vec<usize>,
    /// Trapped-check visitation stamps (one generation per check).
    visited: StampFlags,
    /// The G1 substrate driven by every accepted operation.
    substrate: IncrementalUnionWorkspace,
    prepared_vertices: usize,
    prepared_edges: usize,
}

impl InvasionWorkspace {
    /// An unprepared workspace.
    pub fn new() -> Self {
        Self::default()
    }

    /// Reserve every buffer for a topology of exactly `vertices` vertices
    /// and `edges` physical edges, resetting all state. Re-preparing with
    /// different sizes is supported and leaves no stale state.
    pub fn prepare(&mut self, vertices: usize, edges: usize) -> Result<(), InvasionError> {
        // Clear BEFORE reserving: `try_reserve_exact` budgets for
        // `len + additional`, so reserving while the previous realization's
        // entries are still in the heap would demand up to twice the
        // capacity on the first reuse and allocate once. With the heap
        // drained first, `max(vertices, edges)` is exactly the budget every
        // time (site pushes are bounded by V, bond pushes by E).
        self.heap.clear();
        self.heap
            .try_reserve_exact(vertices.max(edges))
            .map_err(|source| InvasionError::Capacity(ProcessError::from(source)))?;
        reserve_exact(&mut self.parent, vertices, usize::MAX)?;
        reserve_exact(&mut self.vertex_pushed, vertices, 0)?;
        reserve_exact(&mut self.edge_pushed, edges, 0)?;
        reserve_exact(&mut self.retired, vertices, 0)?;
        reserve_exact(&mut self.exit_marks, vertices, 0)?;
        reserve_exact(&mut self.order, vertices, 0)?;
        reserve_exact(&mut self.accepted_edges, vertices, 0)?;
        reserve_exact(&mut self.accepted_thresholds, vertices, 0.0)?;
        reserve_exact(&mut self.dfs_stack, vertices, 0)?;
        self.visited
            .prepare(vertices)
            .map_err(InvasionError::Capacity)?;
        self.substrate
            .prepare(vertices, edges)
            .map_err(|error| match error {
                IncrementalError::Capacity(source) => {
                    InvasionError::Capacity(ProcessError::from(source))
                }
                other => InvasionError::Substrate(other.to_string()),
            })?;
        self.prepared_vertices = vertices;
        self.prepared_edges = edges;
        self.reset();
        Ok(())
    }

    /// Prepared topology sizes.
    #[inline]
    pub const fn prepared_sizes(&self) -> (usize, usize) {
        (self.prepared_vertices, self.prepared_edges)
    }

    /// Reset for the next realization, retaining all capacity: the heap is
    /// drained, logical arrays cleared, the substrate reset.
    pub fn reset(&mut self) {
        self.heap.clear();
        self.parent.fill(usize::MAX);
        self.vertex_pushed.fill(0);
        self.edge_pushed.fill(0);
        self.retired.fill(0);
        self.exit_marks.fill(0);
        self.order.clear();
        self.accepted_edges.clear();
        self.accepted_thresholds.clear();
        self.dfs_stack.clear();
        self.substrate.reset();
    }

    /// Dense vertex indices in invasion order: seeds first (caller order),
    /// then one entry per accepted growth step.
    #[inline]
    pub fn invaded_order(&self) -> &[usize] {
        &self.order
    }

    /// Invasion parent per dense vertex: `usize::MAX` when un-invaded, the
    /// vertex itself for a seed, the invading neighbor otherwise.
    #[inline]
    pub fn parent(&self) -> &[usize] {
        &self.parent
    }

    /// Whether the vertex is invaded in the last run.
    #[inline]
    pub fn is_invaded(&self, index: usize) -> bool {
        self.parent
            .get(index)
            .is_some_and(|&parent| parent != usize::MAX)
    }

    /// Accepted tree edge (dense edge index) per growth step.
    #[inline]
    pub fn accepted_edges(&self) -> &[usize] {
        &self.accepted_edges
    }

    /// Accepted threshold per growth step.
    #[inline]
    pub fn accepted_thresholds(&self) -> &[f64] {
        &self.accepted_thresholds
    }

    /// The driven G1 substrate (online component summaries along the
    /// invasion; trajectory steps are the invasion times).
    #[inline]
    pub const fn substrate(&self) -> &IncrementalUnionWorkspace {
        &self.substrate
    }

    /// Whether an un-invaded candidate is connected to an un-invaded exit
    /// vertex through un-invaded sites (the trapping rule), by DFS with
    /// early exit at the first exit vertex hit.
    fn escapes_to_exit<G: UndirectedGraphView>(&mut self, graph: &G, start: usize) -> bool {
        if self.exit_marks[start] != 0 {
            return true;
        }
        self.visited.new_generation();
        self.visited.mark(start);
        self.dfs_stack.clear();
        self.dfs_stack.push(start);
        while let Some(vertex) = self.dfs_stack.pop() {
            for incidence in graph.incidences(VertexId::from_valid_index(vertex)) {
                let neighbor = incidence.neighbor.index();
                if self.parent[neighbor] != usize::MAX {
                    continue; // invaded: not part of the defending phase
                }
                if self.visited.is_marked(neighbor) {
                    continue;
                }
                if self.exit_marks[neighbor] != 0 {
                    return true; // early exit: escape found
                }
                self.visited.mark(neighbor);
                self.dfs_stack.push(neighbor);
            }
        }
        false
    }

    /// Push the frontier entries created by invading `vertex` (its
    /// incidences into un-invaded territory), under the mode's entry rule.
    fn push_growth<G: UndirectedGraphView>(
        &mut self,
        graph: &G,
        thresholds: &[f64],
        mode: InvasionMode,
        vertex: usize,
    ) {
        for incidence in graph.incidences(VertexId::from_valid_index(vertex)) {
            let neighbor = incidence.neighbor.index();
            if self.parent[neighbor] != usize::MAX {
                continue; // already invaded (self-loops land here too)
            }
            let edge = incidence.edge.index();
            match mode {
                InvasionMode::Site => {
                    if self.vertex_pushed[neighbor] != 0 {
                        continue; // already queued by an earlier invader
                    }
                    self.vertex_pushed[neighbor] = 1;
                    self.heap.push(Reverse(FrontierEntry {
                        threshold: thresholds[neighbor],
                        id: neighbor,
                        vertex: neighbor,
                        edge,
                    }));
                }
                InvasionMode::Bond => {
                    if self.edge_pushed[edge] != 0 {
                        continue;
                    }
                    self.edge_pushed[edge] = 1;
                    self.heap.push(Reverse(FrontierEntry {
                        threshold: thresholds[edge],
                        id: edge,
                        vertex: neighbor,
                        edge,
                    }));
                }
            }
        }
    }
}

/// Reserve `len` logical entries (rewriting the default value so a shrink
/// leaves no stale state, capacity retained).
fn reserve_exact<T: Clone>(
    storage: &mut Vec<T>,
    len: usize,
    default: T,
) -> Result<(), InvasionError> {
    if storage.capacity() < len {
        storage.clear();
        storage
            .try_reserve_exact(len)
            .map_err(|source| InvasionError::Capacity(ProcessError::from(source)))?;
    }
    storage.clear();
    storage.resize(len, default);
    Ok(())
}

/// Run one complete invasion realization (the reusable, RNG-free core).
///
/// `thresholds` is the frozen field — one entry per dense vertex (site
/// mode) or per physical edge (bond mode); every entry must be finite (NaN
/// and infinities are typed rejections). `seeds` must be non-empty,
/// in-domain, and duplicate-free. The workspace is re-prepared
/// idempotently for the graph's sizes, so a correctly sized workspace
/// never allocates in steady state.
pub fn invasion_run<G: UndirectedGraphView>(
    graph: &G,
    mode: InvasionMode,
    thresholds: &[f64],
    seeds: &[usize],
    stop: InvasionStop,
    trapping: &InvasionTrapping,
    workspace: &mut InvasionWorkspace,
) -> Result<InvasionOutcome, InvasionError> {
    let vertices = graph.vertex_count();
    let edges = graph.edge_count();
    if vertices == 0 {
        return Err(InvasionError::EmptyGraph);
    }
    let expected = match mode {
        InvasionMode::Site => vertices,
        InvasionMode::Bond => edges,
    };
    if thresholds.len() != expected {
        return Err(InvasionError::ThresholdLengthMismatch {
            expected,
            actual: thresholds.len(),
        });
    }
    for (index, &value) in thresholds.iter().enumerate() {
        if !value.is_finite() {
            return Err(InvasionError::ThresholdNotFinite { index, value });
        }
    }
    validate_seeds_and_policies(seeds, stop, trapping, vertices)?;
    workspace.prepare(vertices, edges)?;
    if let InvasionTrapping::Trapping { exit } = trapping {
        for &vertex in exit {
            workspace.exit_marks[vertex] = 1;
        }
    }

    // Plant the seeds (substrate step 0).
    for &seed in seeds {
        workspace.parent[seed] = seed;
        workspace.order.push(seed);
        workspace
            .substrate
            .activate_vertex(graph, VertexId::from_valid_index(seed), 0)
            .map_err(|error| InvasionError::Substrate(error.to_string()))?;
    }
    if let InvasionStop::FirstSpanning { from, to } = stop {
        if workspace.is_invaded(from) && workspace.is_invaded(to) {
            return Ok(finished(workspace, 0, InvasionStopReason::SpanningStep, 0));
        }
    }
    if let InvasionStop::Steps(0) = stop {
        return Ok(finished(workspace, 0, InvasionStopReason::StepBudget, 0));
    }

    for &seed in seeds {
        workspace.push_growth(graph, thresholds, mode, seed);
    }

    let mut trapped_sites = 0_usize;
    while let Some(Reverse(entry)) = workspace.heap.pop() {
        if workspace.parent[entry.vertex] != usize::MAX {
            continue; // stale bond entry: far endpoint invaded meanwhile
        }
        if matches!(trapping, InvasionTrapping::Trapping { .. }) {
            if workspace.retired[entry.vertex] != 0 {
                continue; // retired earlier: entry consumed without growth
            }
            if !workspace.escapes_to_exit(graph, entry.vertex) {
                // Trapped forever (the defending phase only shrinks):
                // retire the site permanently.
                workspace.retired[entry.vertex] = 1;
                trapped_sites += 1;
                continue;
            }
        }

        // Accept: derive the parent from the entry's edge (its other
        // endpoint is the already-invaded invader), then drive the
        // substrate with the invasion time (1-based growth step).
        let [left, right] = graph.edge_endpoints(EdgeId::from_valid_index(entry.edge));
        let parent = if left.index() == entry.vertex {
            right.index()
        } else {
            left.index()
        };
        let step = workspace.accepted_edges.len() + 1;
        workspace.parent[entry.vertex] = parent;
        workspace.order.push(entry.vertex);
        workspace.accepted_edges.push(entry.edge);
        workspace.accepted_thresholds.push(entry.threshold);
        workspace
            .substrate
            .activate_vertex(graph, VertexId::from_valid_index(entry.vertex), step)
            .map_err(|error| InvasionError::Substrate(error.to_string()))?;
        workspace
            .substrate
            .activate_edge(graph, EdgeId::from_valid_index(entry.edge), step)
            .map_err(|error| InvasionError::Substrate(error.to_string()))?;
        workspace.push_growth(graph, thresholds, mode, entry.vertex);

        match stop {
            InvasionStop::GrowAll => {}
            InvasionStop::Steps(budget) => {
                if budget == step {
                    return Ok(finished(
                        workspace,
                        trapped_sites,
                        InvasionStopReason::StepBudget,
                        step,
                    ));
                }
            }
            InvasionStop::FirstSpanning { from, to } => {
                if workspace.is_invaded(from) && workspace.is_invaded(to) {
                    return Ok(finished(
                        workspace,
                        trapped_sites,
                        InvasionStopReason::SpanningStep,
                        step,
                    ));
                }
            }
        }
    }

    let accepted = workspace.accepted_edges.len();
    Ok(finished(
        workspace,
        trapped_sites,
        InvasionStopReason::FrontierExhausted,
        accepted,
    ))
}

/// Shared seed/stop/trapping domain validation (core and adapter).
fn validate_seeds_and_policies(
    seeds: &[usize],
    stop: InvasionStop,
    trapping: &InvasionTrapping,
    vertices: usize,
) -> Result<(), InvasionError> {
    if seeds.is_empty() {
        return Err(InvasionError::EmptySeeds);
    }
    for &seed in seeds {
        if seed >= vertices {
            return Err(InvasionError::SeedDomain { seed, vertices });
        }
    }
    for (position, &seed) in seeds.iter().enumerate() {
        if seeds[..position].contains(&seed) {
            return Err(InvasionError::DuplicateSeed { seed });
        }
    }
    if let InvasionStop::FirstSpanning { from, to } = stop {
        for vertex in [from, to] {
            if vertex >= vertices {
                return Err(InvasionError::SpanningDomain { vertex, vertices });
            }
        }
    }
    match trapping {
        InvasionTrapping::NonTrapping => {}
        InvasionTrapping::Trapping { exit } => {
            if exit.is_empty() {
                return Err(InvasionError::EmptyExitSet);
            }
            for &vertex in exit {
                if vertex >= vertices {
                    return Err(InvasionError::ExitDomain { vertex, vertices });
                }
            }
            for (position, &vertex) in exit.iter().enumerate() {
                if exit[..position].contains(&vertex) {
                    return Err(InvasionError::DuplicateExit { vertex });
                }
            }
        }
    }
    Ok(())
}

/// Assemble the outcome from the workspace's terminal state.
fn finished(
    workspace: &InvasionWorkspace,
    trapped_sites: usize,
    reason: InvasionStopReason,
    stop_step: usize,
) -> InvasionOutcome {
    let accepted_total = workspace.accepted_edges.len();
    InvasionOutcome {
        seed_count: workspace.order.len() - accepted_total,
        accepted_total,
        invaded_vertices: workspace.order.len(),
        trapped_sites,
        // `reduce` (not `fold(0.0, max)`): core callers may freeze negative
        // thresholds, and folding from 0.0 would clamp an all-negative
        // maximum up to 0.0.
        max_accepted_threshold: workspace
            .accepted_thresholds
            .iter()
            .copied()
            .reduce(f64::max)
            .unwrap_or(0.0),
        stop_reason: reason,
        stop_step,
    }
}

/// MonteCarlo adapter for invasion percolation on an arbitrary owned
/// undirected graph view.
///
/// One `sweep()` freezes the threshold field (the documented dense-order
/// draw contract) and grows one complete realization under the configured
/// mode, stop, and trapping policies; `measure()` records it as ONE C0
/// element-preserving vector sample plus scalars. Configure
/// `thermalization_sweeps = 0`; realizations are i.i.d.
pub struct InvasionPercolationMC<G: UndirectedGraphView> {
    graph: G,
    mode: InvasionMode,
    seeds: Vec<usize>,
    stop: InvasionStop,
    trapping: InvasionTrapping,
    steps: usize,
    workspace: InvasionWorkspace,
    thresholds: Vec<f64>,
    history: Vec<f64>,
    outcome: InvasionOutcome,
}

impl<G: UndirectedGraphView> InvasionPercolationMC<G> {
    /// Build the adapter. `steps` is the recorded-array length (`>= 1`);
    /// every domain is validated here so `sweep` cannot fail.
    pub fn new(
        graph: G,
        mode: InvasionMode,
        seeds: Vec<usize>,
        stop: InvasionStop,
        trapping: InvasionTrapping,
        steps: usize,
    ) -> Result<Self, InvasionError> {
        let vertices = graph.vertex_count();
        if vertices == 0 {
            return Err(InvasionError::EmptyGraph);
        }
        if steps == 0 {
            return Err(InvasionError::NoSteps);
        }
        validate_seeds_and_policies(&seeds, stop, &trapping, vertices)?;
        let mut workspace = InvasionWorkspace::new();
        workspace.prepare(vertices, graph.edge_count())?;
        let field_len = match mode {
            InvasionMode::Site => vertices,
            InvasionMode::Bond => graph.edge_count(),
        };
        let mut thresholds = Vec::new();
        thresholds
            .try_reserve_exact(field_len)
            .map_err(|source| InvasionError::Capacity(ProcessError::from(source)))?;
        thresholds.resize(field_len, 0.0);
        let mut history = Vec::new();
        history
            .try_reserve_exact(steps)
            .map_err(|source| InvasionError::Capacity(ProcessError::from(source)))?;
        Ok(Self {
            graph,
            mode,
            seeds,
            stop,
            trapping,
            steps,
            workspace,
            thresholds,
            history,
            outcome: InvasionOutcome {
                seed_count: 0,
                accepted_total: 0,
                invaded_vertices: 0,
                trapped_sites: 0,
                max_accepted_threshold: 0.0,
                stop_reason: InvasionStopReason::FrontierExhausted,
                stop_step: 0,
            },
        })
    }

    /// The owned graph the process runs on.
    pub const fn graph(&self) -> &G {
        &self.graph
    }

    /// Which element carries the frozen threshold field.
    pub const fn mode(&self) -> InvasionMode {
        self.mode
    }

    /// The explicit seed set.
    pub fn seeds(&self) -> &[usize] {
        &self.seeds
    }

    /// The explicit stop policy.
    pub const fn stop(&self) -> InvasionStop {
        self.stop
    }

    /// The explicit trapping policy.
    pub fn trapping(&self) -> &InvasionTrapping {
        &self.trapping
    }

    /// Recorded-array length (the C0 shape of `InvasionThresholdOrder`).
    pub const fn steps(&self) -> usize {
        self.steps
    }

    /// Outcome of the last realized invasion.
    pub const fn last_outcome(&self) -> InvasionOutcome {
        self.outcome
    }

    /// The fixed-length recorded array of the last realization.
    pub fn last_history(&self) -> &[f64] {
        &self.history
    }

    /// Dense vertex indices in invasion order of the last realization.
    pub fn last_order(&self) -> &[usize] {
        self.workspace.invaded_order()
    }

    /// Accepted tree edge per growth step of the last realization.
    pub fn last_accepted_edges(&self) -> &[usize] {
        self.workspace.accepted_edges()
    }

    /// Accepted threshold per growth step of the last realization.
    pub fn last_accepted_thresholds(&self) -> &[f64] {
        self.workspace.accepted_thresholds()
    }

    /// Draw the frozen field and run one complete realization.
    ///
    /// Draw contract (bitwise reproducibility): one `f64` per element in
    /// ascending dense order — vertices `0..V` (site mode) or edges
    /// `0..E` (bond mode). Deterministic given the seed law and RNG
    /// stream; exposed for tests and RNG-free analysis paths.
    pub fn realize<R: Rng + ?Sized>(&mut self, rng: &mut R) {
        for slot in self.thresholds.iter_mut() {
            *slot = rng.random::<f64>();
        }
        let outcome = {
            let Self {
                graph,
                mode,
                seeds,
                stop,
                trapping,
                workspace,
                thresholds,
                ..
            } = self;
            invasion_run(graph, *mode, thresholds, seeds, *stop, trapping, workspace)
                .expect("constructor validated every domain")
        };
        self.outcome = outcome;
        // Compress the accepted-threshold sequence into the fixed C0 shape:
        // entries past the growth length repeat the final accepted
        // threshold (zero acceptances fill 0.0).
        let accepted = self.workspace.accepted_thresholds();
        let final_threshold = accepted.last().copied().unwrap_or(0.0);
        self.history.clear();
        self.history.extend((0..self.steps).map(|step| {
            if step < accepted.len() {
                accepted[step]
            } else {
                final_threshold
            }
        }));
    }
}

impl<G: UndirectedGraphView> MonteCarlo for InvasionPercolationMC<G> {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        self.realize(&mut ctx.rng);
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        ctx.measure_array("InvasionThresholdOrder", &self.history)
            .expect("array shape is fixed at construction");
        ctx.measure("InvasionAcceptedTotal", self.outcome.accepted_total as f64);
        ctx.measure(
            "InvasionInvadedVertices",
            self.outcome.invaded_vertices as f64,
        );
        ctx.measure(
            "InvasionMaxAcceptedThreshold",
            self.outcome.max_accepted_threshold,
        );
    }

    fn name(&self) -> &'static str {
        "InvasionPercolation"
    }
}

/// Typed rejections for invasion percolation.
#[derive(Debug)]
pub enum InvasionError {
    /// A graph with zero vertices.
    EmptyGraph,
    /// A zero-length recorded array.
    NoSteps,
    /// An empty seed set.
    EmptySeeds,
    /// A seed outside the graph's dense index range.
    SeedDomain { seed: usize, vertices: usize },
    /// The same seed listed twice.
    DuplicateSeed { seed: usize },
    /// A spanning endpoint outside the graph's dense index range.
    SpanningDomain { vertex: usize, vertices: usize },
    /// An exit (outlet) vertex outside the graph's dense index range.
    ExitDomain { vertex: usize, vertices: usize },
    /// The same exit vertex listed twice.
    DuplicateExit { vertex: usize },
    /// A trapping policy with an empty exit set.
    EmptyExitSet,
    /// A frozen threshold field with the wrong length.
    ThresholdLengthMismatch {
        /// The mode's required length (`V` site, `E` bond).
        expected: usize,
        actual: usize,
    },
    /// A NaN or infinite threshold (the heap ordering must be total).
    ThresholdNotFinite { index: usize, value: f64 },
    /// A buffer reservation failure.
    Capacity(ProcessError),
    /// The driven G1 substrate refused an operation (a process bug, not a
    /// caller error, once the domains validated).
    Substrate(String),
}

impl core::fmt::Display for InvasionError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::EmptyGraph => {
                write!(formatter, "invasion percolation needs at least one vertex")
            }
            Self::NoSteps => write!(formatter, "at least one recorded step is required"),
            Self::EmptySeeds => write!(formatter, "at least one seed vertex is required"),
            Self::SeedDomain { seed, vertices } => {
                write!(formatter, "seed {seed} outside the graph (0..{vertices})")
            }
            Self::DuplicateSeed { seed } => {
                write!(formatter, "seed {seed} listed more than once")
            }
            Self::SpanningDomain { vertex, vertices } => write!(
                formatter,
                "spanning endpoint {vertex} outside the graph (0..{vertices})"
            ),
            Self::ExitDomain { vertex, vertices } => write!(
                formatter,
                "exit vertex {vertex} outside the graph (0..{vertices})"
            ),
            Self::DuplicateExit { vertex } => {
                write!(formatter, "exit vertex {vertex} listed more than once")
            }
            Self::EmptyExitSet => {
                write!(formatter, "trapping invasion needs at least one exit vertex")
            }
            Self::ThresholdLengthMismatch { expected, actual } => write!(
                formatter,
                "threshold field has {actual} entries, expected {expected}"
            ),
            Self::ThresholdNotFinite { index, value } => write!(
                formatter,
                "threshold {value} at index {index} is not finite (non-finite values break the total heap order)"
            ),
            Self::Capacity(source) => {
                write!(formatter, "could not prepare invasion buffers: {source}")
            }
            Self::Substrate(reason) => {
                write!(formatter, "driven G1 substrate refused an operation: {reason}")
            }
        }
    }
}

impl std::error::Error for InvasionError {}
