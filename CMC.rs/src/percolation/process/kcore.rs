//! E3: k-core percolation — residual-degree peeling (roadmap §12).
//!
//! The k-core is the unique maximal induced subgraph in which every
//! surviving vertex `v` keeps at least `t_v` neighbors inside the subgraph
//! (`t_v = k` uniform, or one `u32` threshold per vertex). Two facts drive
//! the implementation:
//!
//! 1. **Order independence is mathematical, not implementational**: the
//!    family of subgraphs satisfying the min-residual-degree constraints is
//!    closed under union (`deg_{H1 | H2}(v) >= max(deg_H1(v), deg_H2(v)) >=
//!    t_v`), so a unique maximal element exists and every deletion order
//!    that only removes deficient vertices reaches it. The exhaustive and
//!    shuffled-relabeling gates in `tests/physics/kcore_percolation.rs`
//!    verify this empirically; the code does not enforce it.
//! 2. **Synchronous-round accounting is exact under FIFO peeling**: peeling
//!    with the E1 two-phase [`ProcessQueue`], popping each queue segment
//!    (the vertices enqueued during the previous round) as one round,
//!    reproduces the synchronous wave process — a vertex crossing below
//!    its threshold during round `r` deletions is deleted in round `r+1`,
//!    exactly as if all round-`r` deletions had happened simultaneously.
//!    `round_survivors()[r]` is therefore the active-vertex count after
//!    peeling round `r` (index `0` = the initial activity, before any
//!    deletion).
//!
//! Degrees follow the platform multigraph convention: one undirected
//! physical edge contributes one incidence at each endpoint, so a
//! self-loop contributes `2` to its vertex's residual degree and parallel
//! edges count separately.
//!
//! [`kcore_decompose`] is the reusable, RNG-free core analysis: generic
//! over [`UndirectedGraphView`], driven through an E1 [`KCoreWorkspace`]
//! (queue + [`NeighborCounters`] residual degrees + [`StampFlags`]
//! deleted/visited marks), allocation-free after prepare. [`KCoreMC`] is
//! the owned MonteCarlo adapter: one sweep draws the initial activity with
//! the F2 uniform site Bernoulli law (bitwise identical to
//! `StaticPercolationMC` site sampling — not a duplicate sampler) and
//! peels deterministically; one measure records the run as ONE C0 vector
//! sample (`KCoreSurvivorCount`, fixed recorded-round budget) plus
//! scalars. One realization is one vector, never a collection of
//! independent scalar samples.
//!
//! Physics anchor: k-core emergence on G(n, m) — Janson & Luczak, Random
//! Structures Algorithms 30, 50-62 (2007) (arXiv:math/0508453, archived in
//! `docs/ref/`), Theorem 1.1: with mean degree `lambda = 2m/n`, the k-core
//! (k >= 3) is empty whp below `lambda_k = min_{mu>0} mu/psi_{k-1}(mu)`
//! and giant above it with `v(Core_k)/n -> psi_k(mu_k(lambda))`,
//! `psi_j(mu) = P(Poisson(mu) >= j)` and `mu_k(lambda)` the largest
//! solution of `mu = lambda * psi_{k-1}(mu)`. Evidence in
//! `CMC.rs/PERCOLATION.md` (E3).

use std::collections::TryReserveError;

use carlo_rs::{Context, MonteCarlo};
use rand::Rng;
use rand_xoshiro::Xoshiro256PlusPlus;

use super::{NeighborCounters, ProcessError, ProcessQueue, StampFlags};
use crate::percolation::{
    Probability, ProbabilityError, ProbabilityField, SiteBernoulli, StaticConfiguration,
    VertexActivity,
};
use crate::{CsrLattice, GraphView, UndirectedGraphView, VertexId};

/// Core-membership thresholds: one `k` for every vertex, or one `u32` per
/// dense vertex index.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum CoreThresholds {
    /// Every vertex must keep at least `k` core neighbors.
    Uniform(u32),
    /// Vertex `i` must keep at least `thresholds[i]` core neighbors.
    PerVertex(Vec<u32>),
}

impl CoreThresholds {
    /// Uniform `k`; rejects values that do not fit the E1 `u32` counters.
    pub fn uniform(k: usize) -> Result<Self, KCoreError> {
        u32::try_from(k)
            .map(Self::Uniform)
            .map_err(|_| KCoreError::ThresholdTooLarge { requested: k })
    }
}

/// Result summary of one k-core decomposition.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct KCoreSummary {
    /// Vertices in the final core (active and never deleted).
    pub core_vertex_count: usize,
    /// Physical edges with both endpoints in the final core.
    pub core_edge_count: usize,
    /// Largest connected component inside the final core.
    pub largest_core_size: usize,
    /// Synchronous peeling rounds that deleted at least one vertex.
    pub peeling_rounds: usize,
    /// Vertices deleted by peeling (initially inactive vertices excluded).
    pub deleted_total: usize,
}

/// Reusable E-track workspace for k-core peeling.
///
/// Owns the E1 primitives plus per-run scratch (deletion waves, deletion
/// order, survivor curve, activity copy). `prepare` once; every
/// [`kcore_decompose`] call reuses the buffers without allocating.
#[derive(Default)]
pub struct KCoreWorkspace {
    queue: ProcessQueue,
    residual: NeighborCounters,
    deleted: StampFlags,
    visited: StampFlags,
    /// Per-vertex thresholds loaded for the current run.
    thresholds: Vec<u32>,
    /// Synchronous deletion round per vertex (valid for deleted vertices).
    wave: Vec<u32>,
    /// Dense vertex indices in deletion (pop) order.
    deletion_order: Vec<usize>,
    /// Active count after each peeling round; index `0` is the initial
    /// activity. Length `peeling_rounds + 1` after a run.
    round_survivors: Vec<f64>,
    /// Copy of the run's initial activity (post-run `in_core` queries).
    active: Vec<bool>,
}

impl KCoreWorkspace {
    /// Prepare every buffer for graphs up to `vertices` vertices.
    pub fn prepare(&mut self, vertices: usize) -> Result<(), KCoreError> {
        self.queue.reserve(vertices)?;
        self.residual.prepare(vertices)?;
        self.deleted.prepare(vertices)?;
        self.visited.prepare(vertices)?;
        reserve_vec(&mut self.thresholds, vertices)?;
        reserve_vec(&mut self.wave, vertices)?;
        reserve_vec(&mut self.deletion_order, vertices)?;
        // rounds <= vertices, plus the initial-count entry.
        reserve_vec(&mut self.round_survivors, vertices + 1)?;
        reserve_vec(&mut self.active, vertices)?;
        Ok(())
    }

    /// Deletion order of the last run: dense vertex indices in the order
    /// they were popped (deleted). Deterministic for a fixed graph,
    /// threshold set, and initial activity.
    pub fn deletion_order(&self) -> &[usize] {
        &self.deletion_order
    }

    /// Synchronous-round survivor curve of the last run: element `r` is the
    /// active-vertex count after peeling round `r` (`0` = initial
    /// activity). Length `peeling_rounds + 1`.
    pub fn round_survivors(&self) -> &[f64] {
        &self.round_survivors
    }

    /// Whether vertex `index` survived in the final core of the last run.
    pub fn in_core(&self, index: usize) -> bool {
        self.active.get(index).copied().unwrap_or(false) && !self.deleted.is_marked(index)
    }

    /// Synchronous round in which `index` was deleted (`None` for
    /// survivors and initially inactive vertices).
    pub fn deletion_round(&self, index: usize) -> Option<usize> {
        if index >= self.wave.len() || self.in_core(index) || !self.active[index] {
            return None;
        }
        Some(self.wave[index] as usize)
    }

    /// Residual degree of `index` at the end of the last run (for core
    /// vertices: the degree inside the core; for deleted vertices: the
    /// value at deletion).
    ///
    /// # Panics
    /// Panics when `index >= vertex_count` (review follow-up 2026-10-23:
    /// unlike `in_core`/`deletion_round` there is no meaningful default
    /// for a missing counter).
    pub fn residual_degree(&self, index: usize) -> u32 {
        self.residual.get(index)
    }
}

fn reserve_vec<T>(buffer: &mut Vec<T>, capacity: usize) -> Result<(), KCoreError> {
    if buffer.capacity() < capacity {
        buffer.clear();
        buffer.try_reserve_exact(capacity)?;
    }
    Ok(())
}

/// Peel the k-core of the subgraph induced by `initial_active` and return
/// its summary.
///
/// Deterministic: initial deficiency scan in ascending dense-index order,
/// neighbor iteration in CSR incidence order. The final core is a property
/// of the input, not of this order (see the module docs). Requires
/// `initial_active.len() == graph.vertex_count()`; heterogeneous
/// thresholds must match the vertex count.
pub fn kcore_decompose<G: UndirectedGraphView>(
    graph: &G,
    thresholds: &CoreThresholds,
    initial_active: &[bool],
    workspace: &mut KCoreWorkspace,
) -> Result<KCoreSummary, KCoreError> {
    let vertices = graph.vertex_count();
    if initial_active.len() != vertices {
        return Err(KCoreError::ActivityLength {
            expected: vertices,
            actual: initial_active.len(),
        });
    }
    if let CoreThresholds::PerVertex(per_vertex) = thresholds {
        if per_vertex.len() != vertices {
            return Err(KCoreError::ThresholdLength {
                expected: vertices,
                actual: per_vertex.len(),
            });
        }
    }
    // Wave and round counters are u32; one check covers both.
    u32::try_from(vertices).map_err(|_| KCoreError::CounterRange {
        value: vertices,
        what: "vertex count",
    })?;
    workspace.prepare(vertices)?;
    workspace.thresholds.clear();
    match thresholds {
        CoreThresholds::Uniform(k) => workspace.thresholds.resize(vertices, *k),
        CoreThresholds::PerVertex(per_vertex) => {
            workspace.thresholds.extend_from_slice(per_vertex);
        }
    }
    workspace.active.clear();
    workspace.active.extend_from_slice(initial_active);
    workspace.queue.reset();
    workspace.deletion_order.clear();
    workspace.round_survivors.clear();
    workspace.wave.clear();
    workspace.wave.resize(vertices, 0);
    workspace.deleted.new_generation();

    // Split the workspace into disjoint borrows for the peeling loop.
    let KCoreWorkspace {
        queue,
        residual,
        deleted,
        visited,
        thresholds,
        wave,
        deletion_order,
        round_survivors,
        active,
    } = workspace;
    let thresholds = &thresholds[..];
    let active = &active[..];

    // Residual degrees over the active subgraph (usize accumulation; one
    // checked conversion per vertex).
    let mut active_total = 0usize;
    for index in 0..vertices {
        if !active[index] {
            continue;
        }
        active_total += 1;
        let mut degree = 0usize;
        for incidence in graph.incidences(VertexId::from_valid_index(index)) {
            if active[incidence.neighbor.index()] {
                degree += 1;
            }
        }
        let degree = u32::try_from(degree).map_err(|_| KCoreError::CounterRange {
            value: degree,
            what: "degree",
        })?;
        residual.set(index, degree);
    }
    round_survivors.push(active_total as f64);

    // Initial deficiency scan (ascending dense index): round-1 deletions.
    for index in 0..vertices {
        if active[index] && residual.get(index) < thresholds[index] {
            deleted.mark(index);
            queue.push(index);
        }
    }

    // Peeling: one queue segment per synchronous round.
    let mut deleted_total = 0usize;
    let mut rounds = 0usize;
    while queue.pending() > 0 {
        rounds += 1;
        let round = rounds as u32;
        for _ in 0..queue.pending() {
            let victim = queue.pop().expect("pending() guards the pop");
            wave[victim] = round;
            deletion_order.push(victim);
            deleted_total += 1;
            for incidence in graph.incidences(VertexId::from_valid_index(victim)) {
                let neighbor = incidence.neighbor.index();
                if !active[neighbor] || deleted.is_marked(neighbor) {
                    continue;
                }
                residual.decrement(neighbor, 1);
                if residual.get(neighbor) < thresholds[neighbor] {
                    deleted.mark(neighbor);
                    queue.push(neighbor);
                }
            }
        }
        round_survivors.push((active_total - deleted_total) as f64);
    }

    // Core facts: vertex and edge counts plus the largest core component.
    let in_core = |index: usize| active[index] && !deleted.is_marked(index);
    let mut core_vertex_count = 0usize;
    let mut core_edge_count = 0usize;
    for edge_index in 0..graph.edge_count() {
        let edge = graph.edge_id(edge_index).expect("dense edge index");
        let [left, right] = graph.edge_endpoints(edge);
        if in_core(left.index()) && in_core(right.index()) {
            core_edge_count += 1;
        }
    }
    visited.new_generation();
    queue.reset();
    let mut largest_core_size = 0usize;
    for seed in 0..vertices {
        if !in_core(seed) || visited.is_marked(seed) {
            continue;
        }
        visited.mark(seed);
        queue.push(seed);
        let mut size = 0usize;
        while let Some(vertex) = queue.pop() {
            size += 1;
            for incidence in graph.incidences(VertexId::from_valid_index(vertex)) {
                let neighbor = incidence.neighbor.index();
                if in_core(neighbor) && !visited.is_marked(neighbor) {
                    visited.mark(neighbor);
                    queue.push(neighbor);
                }
            }
        }
        core_vertex_count += size;
        largest_core_size = largest_core_size.max(size);
    }
    debug_assert_eq!(deleted_total + core_vertex_count, active_total);

    Ok(KCoreSummary {
        core_vertex_count,
        core_edge_count,
        largest_core_size,
        peeling_rounds: rounds,
        deleted_total,
    })
}

/// Owned adapter: random initial activity plus deterministic k-core
/// peeling.
///
/// One `sweep()` draws the initial activity with the F2 uniform site
/// Bernoulli law (one `f64` per vertex in ascending dense order, active
/// iff draw `< keep`; exact `keep = 0` / `keep = 1` shortcuts, bitwise
/// identical to `StaticPercolationMC` site sampling on the same stream)
/// and peels to completion. One `measure()` records ONE C0 vector sample:
///
/// - `KCoreSurvivorCount[r]`: active count after peeling round `r`, for
///   `r = 0..=recorded_rounds` (fixed shape `recorded_rounds + 1`). The
///   peeling itself always runs to completion; if the process terminates
///   before the budget the trailing entries repeat the final core size,
///   and if it exceeds the budget only the recording is truncated (the
///   `KCoreRounds` scalar stays exact) — the D3 fixed-layer convention.
/// - Scalars `KCoreFinalSize`, `KCoreRounds`, `KCoreDeletedTotal`.
///
/// Researchers sweeping `k` (or the keep probability) run one adapter per
/// value: the realization ensemble differs per point, so no single vector
/// output could represent a sweep as one sample. `keep = 1` reduces the
/// adapter to the deterministic full-activity decomposition.
pub struct KCoreMC {
    lattice: CsrLattice,
    probability: Probability,
    law: SiteBernoulli<'static>,
    configuration: StaticConfiguration,
    workspace: KCoreWorkspace,
    thresholds: CoreThresholds,
    budget: usize,
    /// Fixed-length C0 buffer: survivors truncated/padded to the budget.
    survivors: Vec<f64>,
    /// Initial-activity copy handed to the analysis every sweep.
    initial: Vec<bool>,
    summary: KCoreSummary,
}

impl KCoreMC {
    /// Build the adapter over an owned lattice.
    ///
    /// `recorded_rounds >= 1` fixes the `KCoreSurvivorCount` shape.
    pub fn new(
        lattice: CsrLattice,
        keep: f64,
        thresholds: CoreThresholds,
        recorded_rounds: usize,
    ) -> Result<Self, KCoreError> {
        lattice
            .validate()
            .map_err(|error| KCoreError::Topology(error.to_string()))?;
        let vertices = lattice.vertex_count();
        if vertices == 0 {
            return Err(KCoreError::EmptyGraph);
        }
        if recorded_rounds == 0 {
            return Err(KCoreError::NoRecordedRounds);
        }
        if let CoreThresholds::PerVertex(per_vertex) = &thresholds {
            if per_vertex.len() != vertices {
                return Err(KCoreError::ThresholdLength {
                    expected: vertices,
                    actual: per_vertex.len(),
                });
            }
        }
        let probability = Probability::new(keep).map_err(KCoreError::Probability)?;
        let law = SiteBernoulli::new(ProbabilityField::Uniform(probability));
        let mut workspace = KCoreWorkspace::default();
        workspace.prepare(vertices)?;
        let edge_count = lattice.edge_count();
        Ok(Self {
            lattice,
            probability,
            law,
            configuration: StaticConfiguration::new(vertices, edge_count),
            workspace,
            thresholds,
            budget: recorded_rounds,
            survivors: vec![0.0; recorded_rounds + 1],
            initial: vec![false; vertices],
            summary: KCoreSummary::default(),
        })
    }

    /// The owned lattice under decomposition.
    pub const fn lattice(&self) -> &CsrLattice {
        &self.lattice
    }

    /// The initial-activity keep probability.
    pub const fn keep(&self) -> f64 {
        self.probability.get()
    }

    /// The threshold set.
    pub fn thresholds(&self) -> &CoreThresholds {
        &self.thresholds
    }

    /// Recorded-round budget (the `KCoreSurvivorCount` shape minus one).
    pub const fn recorded_rounds(&self) -> usize {
        self.budget
    }

    /// Summary of the last sweep.
    pub const fn summary(&self) -> KCoreSummary {
        self.summary
    }

    /// Fixed-length survivor curve of the last sweep (shape
    /// `recorded_rounds + 1`; see the type documentation for the
    /// truncation/padding convention).
    pub fn round_survivors(&self) -> &[f64] {
        &self.survivors
    }

    /// Deletion order of the last sweep.
    pub fn deletion_order(&self) -> &[usize] {
        self.workspace.deletion_order()
    }

    /// The initial activity drawn for the last sweep (the F2 site-law
    /// output; its bitwise draw parity with `StaticPercolationMC` is a
    /// tested gate).
    pub fn initial_activity(&self) -> &[bool] {
        &self.initial
    }

    /// Whether vertex `index` survived in the last sweep's core.
    pub fn in_core(&self, index: usize) -> bool {
        self.workspace.in_core(index)
    }

    /// Draw the initial activity and peel to completion.
    ///
    /// Exposed for tests and for replaying the documented draw contract;
    /// `measure` records what this produced.
    pub fn draw_and_peel<R: Rng + ?Sized>(&mut self, rng: &mut R) {
        self.law
            .sample(&self.lattice, &mut self.configuration, rng)
            .expect("constructor validated the configuration domains");
        for (index, slot) in self.initial.iter_mut().enumerate() {
            let vertex = VertexId::from_valid_index(index);
            *slot = self.configuration.vertex_active(vertex);
        }
        let summary = kcore_decompose(
            &self.lattice,
            &self.thresholds,
            &self.initial,
            &mut self.workspace,
        )
        .expect("constructor validated the analysis domains");
        let curve = self.workspace.round_survivors();
        let recorded = curve.len().min(self.survivors.len());
        self.survivors[..recorded].copy_from_slice(&curve[..recorded]);
        let final_size = summary.core_vertex_count as f64;
        for value in &mut self.survivors[recorded..] {
            *value = final_size;
        }
        self.summary = summary;
    }
}

impl MonteCarlo for KCoreMC {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        self.draw_and_peel(&mut ctx.rng);
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        ctx.measure_array("KCoreSurvivorCount", &self.survivors)
            .expect("recorded-round budget is fixed at construction");
        ctx.measure("KCoreFinalSize", self.summary.core_vertex_count as f64);
        ctx.measure("KCoreRounds", self.summary.peeling_rounds as f64);
        ctx.measure("KCoreDeletedTotal", self.summary.deleted_total as f64);
    }

    fn name(&self) -> &'static str {
        "KCore"
    }
}

/// Typed rejections for the k-core analysis and adapter.
#[derive(Debug)]
pub enum KCoreError {
    /// A buffer could not be reserved.
    Capacity(ProcessError),
    /// Heterogeneous thresholds do not match the vertex count.
    ThresholdLength { expected: usize, actual: usize },
    /// Initial activity does not match the vertex count.
    ActivityLength { expected: usize, actual: usize },
    /// Uniform threshold larger than the E1 `u32` counters.
    ThresholdTooLarge { requested: usize },
    /// A degree or round count exceeds the `u32` counter range.
    CounterRange { value: usize, what: &'static str },
    /// The adapter needs at least one vertex.
    EmptyGraph,
    /// The adapter needs at least one recorded round.
    NoRecordedRounds,
    /// The lattice failed validation.
    Topology(String),
    /// Keep probability outside `[0, 1]`.
    Probability(ProbabilityError),
}

impl From<ProcessError> for KCoreError {
    fn from(source: ProcessError) -> Self {
        Self::Capacity(source)
    }
}

impl From<TryReserveError> for KCoreError {
    fn from(source: TryReserveError) -> Self {
        Self::Capacity(ProcessError::from(source))
    }
}

impl core::fmt::Display for KCoreError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Capacity(source) => {
                write!(formatter, "could not reserve k-core buffers: {source}")
            }
            Self::ThresholdLength { expected, actual } => write!(
                formatter,
                "per-vertex core thresholds have length {actual}; expected {expected}"
            ),
            Self::ActivityLength { expected, actual } => write!(
                formatter,
                "initial activity has length {actual}; expected {expected}"
            ),
            Self::ThresholdTooLarge { requested } => write!(
                formatter,
                "uniform core threshold {requested} exceeds the u32 counter range"
            ),
            Self::CounterRange { value, what } => {
                write!(formatter, "{what} {value} exceeds the u32 counter range")
            }
            Self::EmptyGraph => write!(formatter, "k-core adapter needs at least one vertex"),
            Self::NoRecordedRounds => {
                write!(
                    formatter,
                    "k-core adapter needs at least one recorded round"
                )
            }
            Self::Topology(reason) => write!(formatter, "invalid lattice: {reason}"),
            Self::Probability(source) => write!(formatter, "invalid keep probability: {source}"),
        }
    }
}

impl std::error::Error for KCoreError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Capacity(source) => Some(source),
            Self::Probability(source) => Some(source),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Bond, BondType};

    fn path(vertices: usize) -> CsrLattice {
        CsrLattice::from_edges(
            vertices,
            (0..vertices - 1)
                .map(|site| Bond::new(site, site + 1, BondType::Generic, 1.0))
                .collect(),
        )
    }

    #[test]
    fn path_peeling_hand_checked_rounds() {
        // k = 2 on a 6-path: rounds delete {0,5}, {1,4}, {2,3}.
        let lattice = path(6);
        let mut workspace = KCoreWorkspace::default();
        let summary = kcore_decompose(
            &lattice,
            &CoreThresholds::uniform(2).unwrap(),
            &[true; 6],
            &mut workspace,
        )
        .unwrap();
        assert_eq!(summary.core_vertex_count, 0);
        assert_eq!(summary.peeling_rounds, 3);
        assert_eq!(summary.deleted_total, 6);
        assert_eq!(workspace.round_survivors(), &[6.0, 4.0, 2.0, 0.0]);
        assert_eq!(workspace.deletion_order(), &[0, 5, 1, 4, 2, 3]);
        assert_eq!(workspace.deletion_round(0), Some(1));
        assert_eq!(workspace.deletion_round(1), Some(2));
        assert_eq!(workspace.deletion_round(2), Some(3));
        assert!(workspace.deletion_round(9).is_none());
    }

    #[test]
    fn heterogeneous_thresholds_peel_in_waves() {
        // Thresholds [2,2,1,2] on a 4-path: leaves first, then 1, then 2.
        let lattice = path(4);
        let mut workspace = KCoreWorkspace::default();
        let summary = kcore_decompose(
            &lattice,
            &CoreThresholds::PerVertex(vec![2, 2, 1, 2]),
            &[true; 4],
            &mut workspace,
        )
        .unwrap();
        assert_eq!(summary.peeling_rounds, 3);
        assert_eq!(workspace.round_survivors(), &[4.0, 2.0, 1.0, 0.0]);
        assert_eq!(workspace.deletion_round(2), Some(3));
    }

    #[test]
    fn length_mismatches_are_typed_errors() {
        let lattice = path(4);
        let mut workspace = KCoreWorkspace::default();
        match kcore_decompose(
            &lattice,
            &CoreThresholds::PerVertex(vec![1; 3]),
            &[true; 4],
            &mut workspace,
        ) {
            Err(KCoreError::ThresholdLength {
                expected: 4,
                actual: 3,
            }) => {}
            other => panic!("expected ThresholdLength, got {other:?}"),
        }
        match kcore_decompose(
            &lattice,
            &CoreThresholds::uniform(1).unwrap(),
            &[true; 3],
            &mut workspace,
        ) {
            Err(KCoreError::ActivityLength {
                expected: 4,
                actual: 3,
            }) => {}
            other => panic!("expected ActivityLength, got {other:?}"),
        }
        assert!(matches!(
            CoreThresholds::uniform((u32::MAX as usize) + 1),
            Err(KCoreError::ThresholdTooLarge { .. })
        ));
    }
}
