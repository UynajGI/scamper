//! E2: bootstrap percolation (monotone activation-threshold closure).
//!
//! An inactive vertex activates once the number of its **active-neighbor
//! incidences** reaches a per-vertex threshold: homogeneous `r` or a
//! heterogeneous `u32` vector. Parallel edges count once per incidence (the
//! F3 physical-edge convention); a self-loop never helps its own vertex —
//! an inactive vertex is not its own active neighbor, and once active the
//! increment is skipped.
//!
//! Two explicit, non-interchangeable clocks (roadmap §12 E2: "specify
//! ... whether recorded time means synchronous rounds or queue closure
//! order"):
//!
//! - [`BootstrapDynamics::SynchronousRounds`]: the discrete parallel rounds
//!   of the literature (Holroyd's clock). `A_0` = seeds; every inactive
//!   vertex with `|N(v) ∩ A_t| >= r_v` activates in round `t + 1`; the
//!   trajectory is `[|A_0|, |A_1|, ...]` and `rounds` counts update rounds
//!   until no vertex activates.
//! - [`BootstrapDynamics::QueueClosure`]: an asynchronous FIFO monotone
//!   closure. Activation events apply immediately; `events` counts rule
//!   activations after the seeds and the trajectory is the active count
//!   after each event. The **final closure is independent of the update
//!   order** (monotone-process fixed point); only the event sequence
//!   depends on it.
//!
//! Threshold-0 vertices activate at the first update step in both clocks
//! (zero active neighbors already meet a zero requirement) and count as
//! rule activations, not seeds. Consequences: `r = 0` activates every
//! vertex; `r > max_degree` leaves exactly the seeds (plus any threshold-0
//! vertex, impossible for homogeneous `r > 0`).
//!
//! The core analysis ([`BootstrapWorkspace::synchronous_rounds`] /
//! [`BootstrapWorkspace::queue_closure`]) is RNG-free and generic over
//! `G: UndirectedGraphView` with seeds supplied as any
//! [`VertexActivity`] — researchers bring arbitrary graphs and
//! caller-sampled activity; a single O(V + E) closure replaces repeated
//! full scans. Buffers are the E1 primitives
//! ([`ProcessQueue`], [`NeighborCounters`], [`StampFlags`]) and are reused
//! across realizations with zero steady-state allocation; there is
//! deliberately no shared physical transition type.
//!
//! Physics anchor: 2D `r = 2` bootstrap on the L x L torus, Holroyd,
//! "Sharp metastability threshold for two-dimensional bootstrap
//! percolation", Probab. Theory Relat. Fields 125, 195-224 (2003)
//! (arXiv:math/0206132, archived in `docs/ref/`): `I(L, p) -> 1` iff
//! `p ln L > pi^2/18`, i.e. the finite-size seed threshold
//! `p_c(L) ~ pi^2 / (18 ln L)`. Evidence in `CMC.rs/PERCOLATION.md` (E2).

use carlo_rs::{Context, MonteCarlo};
use rand::{Rng, RngExt};
use rand_xoshiro::Xoshiro256PlusPlus;

use super::{NeighborCounters, ProcessError, ProcessQueue, StampFlags};
use crate::percolation::{Probability, ProbabilityError, VertexActivity};
use crate::{UndirectedGraphView, VertexId};

/// Activation thresholds.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum BootstrapThreshold {
    /// Every vertex activates once at least `value` of its neighbor
    /// incidences are active.
    Homogeneous {
        /// The common threshold `r`.
        value: u32,
    },
    /// One threshold per dense vertex index.
    Heterogeneous {
        /// Thresholds, indexed by dense vertex ID.
        values: Vec<u32>,
    },
}

/// How initial seeds are produced for one realization.
#[derive(Clone, Debug, PartialEq)]
pub enum BootstrapSeeds {
    /// Each vertex is an initial seed independently with probability `p`.
    ///
    /// RNG draw contract (bitwise reproducibility): one `f64` draw per
    /// vertex in ascending dense-index order, seed iff `draw < p`.
    Density { p: f64 },
    /// Explicit dense seed indices, identical in every realization. Pair
    /// with an externally sampled [`VertexActivity`] through the core
    /// workspace when more general seed laws are needed.
    Fixed { vertices: Vec<usize> },
}

/// Which dynamics the recorded time axis means.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BootstrapDynamics {
    /// Parallel discrete rounds: vertex `v` activates in round `t + 1` iff
    /// it is inactive after round `t` and has at least `r_v` active
    /// neighbors at the end of round `t`. Recorded steps are rounds
    /// (`|A_0|` = seeds at step 0).
    SynchronousRounds,
    /// Asynchronous FIFO monotone closure: a crossing vertex activates as
    /// soon as its queue turn arrives. Recorded steps are activation
    /// events (step 0 = seeds, step `i` = active count after the `i`-th
    /// rule activation). The final closure is order-independent.
    QueueClosure,
}

/// Summary of one closure realization.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BootstrapOutcome {
    /// Initially active vertices from the seed activity (threshold-0
    /// vertices are NOT counted here — they are rule activations).
    pub seed_count: usize,
    /// Active vertices in the final closure.
    pub final_active: usize,
    /// Rule activations beyond the seeds (`final_active - seed_count`),
    /// threshold-0 activations included.
    pub activated_total: usize,
    /// [`BootstrapDynamics::SynchronousRounds`]: update rounds until no
    /// vertex activated (0 when the seeds are already closed). Reads 0 in
    /// queue mode.
    pub rounds: usize,
    /// [`BootstrapDynamics::QueueClosure`]: rule-activation events after
    /// the seeds. Reads 0 in synchronous mode.
    pub events: usize,
}

/// Reusable closure workspace over the E1 primitives.
///
/// `prepare` once, then run either clock any number of times; every run re-zeroes the
/// counters (O(V), allocation-free) and reuses the retained capacity. The
/// caller owns the trajectory buffer so non-Carlo analyses never depend on
/// scheduler plumbing.
#[derive(Debug, Default)]
pub struct BootstrapWorkspace {
    /// Synchronous: current-round sources. Asynchronous: the closure queue.
    queue: ProcessQueue,
    /// Synchronous: next-round candidates (unused by queue closure).
    next: ProcessQueue,
    counters: NeighborCounters,
    flags: StampFlags,
    /// Flattened per-vertex thresholds, refilled per run.
    thresholds: Vec<u32>,
}

impl BootstrapWorkspace {
    /// An unprepared workspace.
    pub fn new() -> Self {
        Self::default()
    }

    /// Reserve every buffer for `vertices` vertices. Re-preparing with the
    /// same size only re-zeroes state (never allocates).
    pub fn prepare(&mut self, vertices: usize) -> Result<(), BootstrapError> {
        if vertices == 0 {
            return Err(BootstrapError::EmptyGraph);
        }
        self.queue
            .reserve(vertices)
            .map_err(BootstrapError::Capacity)?;
        self.next
            .reserve(vertices)
            .map_err(BootstrapError::Capacity)?;
        self.counters
            .prepare(vertices)
            .map_err(BootstrapError::Capacity)?;
        self.flags
            .prepare(vertices)
            .map_err(BootstrapError::Capacity)?;
        if self.thresholds.capacity() < vertices {
            self.thresholds.clear();
            self.thresholds
                .try_reserve_exact(vertices)
                .map_err(|source| BootstrapError::Capacity(ProcessError::from(source)))?;
        }
        self.thresholds.clear();
        self.thresholds.resize(vertices, 0);
        Ok(())
    }

    /// Whether `index` is active in the last run's closure.
    ///
    /// Valid until the next run on this workspace; indices are dense
    /// vertex positions.
    /// Whether `index` is active after the last run.
    ///
    /// Meaningful only after at least one closure run: before any run on a
    /// prepared workspace every vertex reads as active (slots hold stamp
    /// zero), and on an unprepared workspace the call panics (review
    /// follow-up 2026-10-23).
    ///
    /// # Panics
    /// Panics when the workspace was never prepared.
    pub fn is_active(&self, index: usize) -> bool {
        self.flags.is_marked(index)
    }

    /// Flatten the threshold rule into the cached per-vertex array.
    fn materialize(
        &mut self,
        thresholds: &BootstrapThreshold,
        vertices: usize,
    ) -> Result<(), BootstrapError> {
        match thresholds {
            BootstrapThreshold::Homogeneous { value } => self.thresholds.fill(*value),
            BootstrapThreshold::Heterogeneous { values } => {
                if values.len() != vertices {
                    return Err(BootstrapError::ThresholdLengthMismatch {
                        expected: vertices,
                        actual: values.len(),
                    });
                }
                self.thresholds.copy_from_slice(values);
            }
        }
        Ok(())
    }

    /// Mark the seeds, count them, and enqueue the threshold-0 vertices
    /// that qualify from the start. Returns the seed count.
    fn plant_seeds<V: VertexActivity>(
        &mut self,
        seeds: &V,
        vertices: usize,
    ) -> Result<usize, BootstrapError> {
        if seeds.vertex_count() != vertices {
            return Err(BootstrapError::ActivityDomain {
                expected: vertices,
                actual: seeds.vertex_count(),
            });
        }
        self.flags.new_generation();
        let mut seed_count = 0usize;
        for index in 0..vertices {
            let vertex = VertexId::from_valid_index(index);
            if seeds.vertex_active(vertex) {
                self.flags.mark(index);
                seed_count += 1;
            }
        }
        Ok(seed_count)
    }

    /// Synchronous-rounds closure (the literature's parallel clock).
    ///
    /// `trajectory` is cleared and refilled with `[|A_0|, ..., |A_R|]`
    /// (caller-owned buffer, capacity reused). Each vertex enters a round
    /// frontier at most once, so the trajectory has at most `V + 1`
    /// entries.
    pub fn synchronous_rounds<G, V>(
        &mut self,
        graph: &G,
        seeds: &V,
        thresholds: &BootstrapThreshold,
        trajectory: &mut Vec<usize>,
    ) -> Result<BootstrapOutcome, BootstrapError>
    where
        G: UndirectedGraphView,
        V: VertexActivity,
    {
        let vertices = graph.vertex_count();
        self.prepare(vertices)?;
        self.materialize(thresholds, vertices)?;
        let seed_count = self.plant_seeds(seeds, vertices)?;

        self.queue.reset();
        self.next.reset();
        for index in 0..vertices {
            if self.flags.is_marked(index) {
                self.queue.push(index);
            } else if self.thresholds[index] == 0 {
                // Zero threshold needs no neighbor information: the vertex
                // qualifies for the first update round from the start.
                self.next.push(index);
            }
        }

        let mut active = seed_count;
        let mut rounds = 0usize;
        trajectory.clear();
        trajectory.push(active);
        loop {
            // Propagate the current round's activations; first crossings
            // (counter == threshold, exact equality fires once) collect in
            // `next` as the next round's candidates.
            while let Some(source) = self.queue.pop() {
                for incidence in graph.incidences(VertexId::from_valid_index(source)) {
                    let neighbor = incidence.neighbor.index();
                    if !self.flags.is_marked(neighbor) {
                        self.counters.increment(neighbor, 1);
                        if self.counters.get(neighbor) == self.thresholds[neighbor] {
                            self.next.push(neighbor);
                        }
                    }
                }
            }
            if self.next.is_empty() {
                break;
            }
            // Promote: every candidate activates this round and becomes a
            // source of the next.
            rounds += 1;
            while let Some(candidate) = self.next.pop() {
                if !self.flags.is_marked(candidate) {
                    self.flags.mark(candidate);
                    active += 1;
                    self.queue.push(candidate);
                }
            }
            trajectory.push(active);
        }

        Ok(BootstrapOutcome {
            seed_count,
            final_active: active,
            activated_total: active - seed_count,
            rounds,
            events: 0,
        })
    }

    /// Asynchronous FIFO monotone closure.
    ///
    /// `trajectory` is cleared and refilled with `[seeds, active after
    /// event 1, ...]` (caller-owned buffer, capacity reused; at most
    /// `V + 1` entries). The initial queue is seeds in ascending order,
    /// then threshold-0 inactive vertices ascending; every later event is
    /// a first crossing. The final closure is independent of this order.
    pub fn queue_closure<G, V>(
        &mut self,
        graph: &G,
        seeds: &V,
        thresholds: &BootstrapThreshold,
        trajectory: &mut Vec<usize>,
    ) -> Result<BootstrapOutcome, BootstrapError>
    where
        G: UndirectedGraphView,
        V: VertexActivity,
    {
        let vertices = graph.vertex_count();
        self.prepare(vertices)?;
        self.materialize(thresholds, vertices)?;
        let seed_count = self.plant_seeds(seeds, vertices)?;

        self.queue.reset();
        for index in 0..vertices {
            if self.flags.is_marked(index) {
                self.queue.push(index);
            }
        }
        for index in 0..vertices {
            if !self.flags.is_marked(index) && self.thresholds[index] == 0 {
                self.queue.push(index);
            }
        }

        let mut active = seed_count;
        let mut events = 0usize;
        trajectory.clear();
        trajectory.push(active);
        while let Some(item) = self.queue.pop() {
            if !self.flags.is_marked(item) {
                self.flags.mark(item);
                active += 1;
                events += 1;
                trajectory.push(active);
            }
            for incidence in graph.incidences(VertexId::from_valid_index(item)) {
                let neighbor = incidence.neighbor.index();
                if !self.flags.is_marked(neighbor) {
                    self.counters.increment(neighbor, 1);
                    if self.counters.get(neighbor) == self.thresholds[neighbor] {
                        self.queue.push(neighbor);
                    }
                }
            }
        }

        Ok(BootstrapOutcome {
            seed_count,
            final_active: active,
            activated_total: active - seed_count,
            rounds: 0,
            events,
        })
    }
}

/// Caller-sampled seed set for the adapter's [`BootstrapSeeds::Density`]
/// law; implements [`VertexActivity`] so the core runs unchanged.
struct SampledSeeds {
    active: Vec<u8>,
}

impl VertexActivity for SampledSeeds {
    fn vertex_count(&self) -> usize {
        self.active.len()
    }

    fn vertex_active(&self, vertex: VertexId) -> bool {
        self.active[vertex.index()] != 0
    }
}

/// MonteCarlo adapter for bootstrap percolation on an arbitrary owned
/// undirected graph view (a `CsrLattice`, a borrowed CSR, or a
/// researcher's own `UndirectedGraphView` type).
///
/// One `sweep()` draws one complete realization (seed law plus closure to
/// the fixed point); the following `measure()` records it as ONE C0
/// element-preserving array sample plus closure scalars — never as
/// independent scalar samples per step:
///
/// - `BootstrapActiveCount[i]`, `i` in `0..steps`: active count at recorded
///   step `i` (synchronous round `i`; queue event `i`; step 0 is the seed
///   state in both clocks). Entries past the closure repeat the final
///   count; when closure takes more steps than `steps` the array truncates
///   while the scalars still report the full closure. The array length is
///   fixed by `steps`, so every realization keeps the same shape.
/// - `BootstrapFinalActive`: active vertices in the final closure.
/// - `BootstrapActivatedTotal`: rule activations beyond the seeds
///   (`final_active - seed_count`).
/// - `BootstrapRounds` (synchronous) or `BootstrapEvents` (queue): the
///   clock's closure time.
///
/// Configure `thermalization_sweeps = 0`; realizations are i.i.d.
pub struct BootstrapPercolationMC<G: UndirectedGraphView> {
    graph: G,
    thresholds: BootstrapThreshold,
    dynamics: BootstrapDynamics,
    seeds: BootstrapSeeds,
    steps: usize,
    workspace: BootstrapWorkspace,
    sampled: SampledSeeds,
    /// Full trajectory of the last realization (scratch, capacity V + 1).
    trajectory: Vec<usize>,
    /// Fixed-length recorded array of the last realization.
    history: Vec<f64>,
    outcome: BootstrapOutcome,
}

impl<G: UndirectedGraphView> BootstrapPercolationMC<G> {
    /// Build the adapter.
    ///
    /// `steps` is the recorded-array length (`>= 1`). All domains are
    /// validated here so `sweep` cannot fail.
    pub fn new(
        graph: G,
        thresholds: BootstrapThreshold,
        seeds: BootstrapSeeds,
        dynamics: BootstrapDynamics,
        steps: usize,
    ) -> Result<Self, BootstrapError> {
        let vertices = graph.vertex_count();
        if vertices == 0 {
            return Err(BootstrapError::EmptyGraph);
        }
        if steps == 0 {
            return Err(BootstrapError::NoSteps);
        }
        if let BootstrapThreshold::Heterogeneous { values } = &thresholds {
            if values.len() != vertices {
                return Err(BootstrapError::ThresholdLengthMismatch {
                    expected: vertices,
                    actual: values.len(),
                });
            }
        }
        match &seeds {
            BootstrapSeeds::Density { p } => {
                Probability::new(*p).map_err(BootstrapError::Probability)?;
            }
            BootstrapSeeds::Fixed { vertices: fixed } => {
                for &seed in fixed {
                    if seed >= vertices {
                        return Err(BootstrapError::SeedDomain { seed, vertices });
                    }
                }
            }
        }
        let mut workspace = BootstrapWorkspace::new();
        workspace.prepare(vertices)?;
        let sampled = SampledSeeds {
            active: vec![0; vertices],
        };
        let mut trajectory = Vec::new();
        trajectory
            .try_reserve_exact(vertices + 1)
            .map_err(|source| BootstrapError::Capacity(ProcessError::from(source)))?;
        let history = Vec::with_capacity(steps);
        Ok(Self {
            graph,
            thresholds,
            dynamics,
            seeds,
            steps,
            workspace,
            sampled,
            trajectory,
            history,
            outcome: BootstrapOutcome {
                seed_count: 0,
                final_active: 0,
                activated_total: 0,
                rounds: 0,
                events: 0,
            },
        })
    }

    /// The owned graph the process runs on.
    pub const fn graph(&self) -> &G {
        &self.graph
    }

    /// The activation-threshold rule.
    pub fn threshold(&self) -> &BootstrapThreshold {
        &self.thresholds
    }

    /// The seed law.
    pub fn seeds(&self) -> &BootstrapSeeds {
        &self.seeds
    }

    /// The recorded time convention.
    pub const fn dynamics(&self) -> BootstrapDynamics {
        self.dynamics
    }

    /// Recorded-array length (the C0 shape of `BootstrapActiveCount`).
    pub const fn steps(&self) -> usize {
        self.steps
    }

    /// Outcome of the last realized closure.
    pub const fn last_outcome(&self) -> BootstrapOutcome {
        self.outcome
    }

    /// The fixed-length recorded array of the last realization.
    pub fn last_history(&self) -> &[f64] {
        &self.history
    }

    /// The full per-round/per-event trajectory of the last realization.
    pub fn last_trajectory(&self) -> &[usize] {
        &self.trajectory
    }

    /// Draw the seeds and run one complete closure realization.
    ///
    /// Deterministic given the seed law and RNG stream; exposed for tests
    /// and RNG-free analysis paths that do not record measurements.
    pub fn realize<R: Rng + ?Sized>(&mut self, rng: &mut R) {
        match &self.seeds {
            BootstrapSeeds::Density { p } => {
                let p = Probability::new(*p)
                    .expect("constructor validated the seed density")
                    .get();
                for slot in self.sampled.active.iter_mut() {
                    *slot = u8::from(rng.random::<f64>() < p);
                }
            }
            BootstrapSeeds::Fixed { vertices } => {
                self.sampled.active.fill(0);
                for &seed in vertices {
                    self.sampled.active[seed] = 1;
                }
            }
        }
        let Self {
            graph,
            thresholds,
            dynamics,
            workspace,
            sampled,
            trajectory,
            ..
        } = self;
        let outcome = match dynamics {
            BootstrapDynamics::SynchronousRounds => {
                workspace.synchronous_rounds(graph, sampled, thresholds, trajectory)
            }
            BootstrapDynamics::QueueClosure => {
                workspace.queue_closure(graph, sampled, thresholds, trajectory)
            }
        }
        .expect("constructor validated every domain");
        self.outcome = outcome;
        // Compress the full trajectory into the fixed C0 shape: steps past
        // the closure repeat the final count.
        let trajectory = &self.trajectory;
        let final_active = trajectory[trajectory.len() - 1] as f64;
        self.history.clear();
        self.history.extend((0..self.steps).map(|step| {
            if step < trajectory.len() {
                trajectory[step] as f64
            } else {
                final_active
            }
        }));
    }
}

impl<G: UndirectedGraphView> MonteCarlo for BootstrapPercolationMC<G> {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        self.realize(&mut ctx.rng);
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        ctx.measure_array("BootstrapActiveCount", &self.history)
            .expect("array shape is fixed at construction");
        ctx.measure("BootstrapFinalActive", self.outcome.final_active as f64);
        ctx.measure(
            "BootstrapActivatedTotal",
            self.outcome.activated_total as f64,
        );
        match self.dynamics {
            BootstrapDynamics::SynchronousRounds => {
                ctx.measure("BootstrapRounds", self.outcome.rounds as f64);
            }
            BootstrapDynamics::QueueClosure => {
                ctx.measure("BootstrapEvents", self.outcome.events as f64);
            }
        }
    }

    fn name(&self) -> &'static str {
        "BootstrapPercolation"
    }
}

/// Typed rejections for bootstrap percolation.
#[derive(Debug)]
pub enum BootstrapError {
    /// A graph with zero vertices (the closure would be vacuous; the
    /// roadmap asks for typed semantics rather than a silent empty pass).
    EmptyGraph,
    /// A zero-length recorded array.
    NoSteps,
    /// Heterogeneous thresholds with the wrong length.
    ThresholdLengthMismatch { expected: usize, actual: usize },
    /// A fixed seed outside the graph's dense index range.
    SeedDomain { seed: usize, vertices: usize },
    /// Seed activity whose vertex count differs from the graph.
    ActivityDomain { expected: usize, actual: usize },
    /// A seed density outside `[0, 1]`.
    Probability(ProbabilityError),
    /// A buffer reservation failure.
    Capacity(ProcessError),
}

impl core::fmt::Display for BootstrapError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::EmptyGraph => write!(formatter, "bootstrap percolation needs at least one vertex"),
            Self::NoSteps => write!(formatter, "at least one recorded step is required"),
            Self::ThresholdLengthMismatch { expected, actual } => write!(
                formatter,
                "heterogeneous thresholds need one entry per vertex (expected {expected}, got {actual})"
            ),
            Self::SeedDomain { seed, vertices } => {
                write!(formatter, "seed {seed} outside the graph (0..{vertices})")
            }
            Self::ActivityDomain { expected, actual } => write!(
                formatter,
                "seed activity covers {actual} vertices, graph has {expected}"
            ),
            Self::Probability(source) => write!(formatter, "invalid seed density: {source}"),
            Self::Capacity(source) => write!(formatter, "could not prepare bootstrap buffers: {source}"),
        }
    }
}

impl std::error::Error for BootstrapError {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::build_chain;

    /// Dense byte-vector seed activity (test-side VertexActivity).
    struct MaskSeeds {
        active: Vec<u8>,
    }

    impl VertexActivity for MaskSeeds {
        fn vertex_count(&self) -> usize {
            self.active.len()
        }

        fn vertex_active(&self, vertex: VertexId) -> bool {
            self.active[vertex.index()] != 0
        }
    }

    fn mask(bits: &[u8]) -> MaskSeeds {
        MaskSeeds {
            active: bits.to_vec(),
        }
    }

    #[test]
    fn chain_run_off_never_activates() {
        let chain = build_chain(4, false);
        let mut workspace = BootstrapWorkspace::new();
        let seeds = mask(&[1, 0, 0, 0]);
        let mut trajectory = Vec::new();
        let outcome = workspace
            .synchronous_rounds(
                &chain,
                &seeds,
                &BootstrapThreshold::Homogeneous { value: 3 },
                &mut trajectory,
            )
            .unwrap();
        // r = 3 exceeds the chain's max degree 2: only the seed remains.
        assert_eq!(outcome.final_active, 1);
        assert_eq!(outcome.activated_total, 0);
        assert_eq!(outcome.rounds, 0);
        assert_eq!(trajectory, vec![1]);
    }

    #[test]
    fn zero_threshold_activates_everything_from_first_step() {
        let chain = build_chain(5, false);
        let mut workspace = BootstrapWorkspace::new();
        let seeds = mask(&[0, 1, 0, 0, 0]);
        let mut trajectory = Vec::new();
        for dynamics in [
            BootstrapDynamics::SynchronousRounds,
            BootstrapDynamics::QueueClosure,
        ] {
            let outcome = match dynamics {
                BootstrapDynamics::SynchronousRounds => workspace.synchronous_rounds(
                    &chain,
                    &seeds,
                    &BootstrapThreshold::Homogeneous { value: 0 },
                    &mut trajectory,
                ),
                BootstrapDynamics::QueueClosure => workspace.queue_closure(
                    &chain,
                    &seeds,
                    &BootstrapThreshold::Homogeneous { value: 0 },
                    &mut trajectory,
                ),
            }
            .unwrap();
            assert_eq!(outcome.final_active, 5, "r = 0 activates all vertices");
            assert_eq!(outcome.seed_count, 1);
            assert_eq!(outcome.activated_total, 4);
            assert_eq!(trajectory.first(), Some(&1));
            assert_eq!(trajectory.last(), Some(&5));
        }
    }

    #[test]
    fn heterogeneous_thresholds_hand_example() {
        // Path 0-1-2-3-4 with thresholds [2,1,1,2,3], seeds {0,4}:
        // round 1 activates 1; round 2 activates 2; round 3 activates 3.
        let chain = build_chain(5, false);
        let mut workspace = BootstrapWorkspace::new();
        let seeds = mask(&[1, 0, 0, 0, 1]);
        let thresholds = BootstrapThreshold::Heterogeneous {
            values: vec![2, 1, 1, 2, 3],
        };
        let mut trajectory = Vec::new();
        let outcome = workspace
            .synchronous_rounds(&chain, &seeds, &thresholds, &mut trajectory)
            .unwrap();
        assert_eq!(trajectory, vec![2, 3, 4, 5]);
        assert_eq!(outcome.rounds, 3);
        let queue_outcome = workspace
            .queue_closure(&chain, &seeds, &thresholds, &mut trajectory)
            .unwrap();
        assert_eq!(queue_outcome.final_active, 5);
        assert_eq!(queue_outcome.events, 3);
    }

    #[test]
    fn typed_rejections() {
        let chain = build_chain(3, false);
        let mut workspace = BootstrapWorkspace::new();
        let seeds = mask(&[1, 0, 0]);
        let mut trajectory = Vec::new();
        let mismatch = workspace.synchronous_rounds(
            &chain,
            &seeds,
            &BootstrapThreshold::Heterogeneous { values: vec![1, 1] },
            &mut trajectory,
        );
        assert!(matches!(
            mismatch,
            Err(BootstrapError::ThresholdLengthMismatch {
                expected: 3,
                actual: 2
            })
        ));
        let domain = workspace.synchronous_rounds(
            &chain,
            &mask(&[1, 0]),
            &BootstrapThreshold::Homogeneous { value: 1 },
            &mut trajectory,
        );
        assert!(matches!(
            domain,
            Err(BootstrapError::ActivityDomain {
                expected: 3,
                actual: 2
            })
        ));
        // `CsrLattice` cannot represent zero sites, so exercise the typed
        // empty-graph rejection through a minimal foreign view.
        struct NoVertices;
        impl crate::GraphView for NoVertices {
            fn vertex_count(&self) -> usize {
                0
            }
        }
        impl crate::UndirectedGraphView for NoVertices {
            fn edge_count(&self) -> usize {
                0
            }
            fn edge_endpoints(&self, _edge: crate::EdgeId) -> [crate::VertexId; 2] {
                panic!("no edges exist");
            }
            fn incidences(
                &self,
                _vertex: crate::VertexId,
            ) -> impl Iterator<Item = crate::Incidence> + '_ {
                std::iter::empty()
            }
        }
        let empty = BootstrapPercolationMC::new(
            NoVertices,
            BootstrapThreshold::Homogeneous { value: 2 },
            BootstrapSeeds::Density { p: 0.5 },
            BootstrapDynamics::SynchronousRounds,
            4,
        );
        assert!(matches!(empty, Err(BootstrapError::EmptyGraph)));
        let empty_run = workspace.queue_closure(
            &NoVertices,
            &mask(&[]),
            &BootstrapThreshold::Homogeneous { value: 2 },
            &mut trajectory,
        );
        assert!(matches!(empty_run, Err(BootstrapError::EmptyGraph)));
        let bad_seed = BootstrapPercolationMC::new(
            build_chain(3, false),
            BootstrapThreshold::Homogeneous { value: 2 },
            BootstrapSeeds::Fixed { vertices: vec![3] },
            BootstrapDynamics::QueueClosure,
            4,
        );
        assert!(matches!(
            bad_seed,
            Err(BootstrapError::SeedDomain {
                seed: 3,
                vertices: 3
            })
        ));
        let no_steps = BootstrapPercolationMC::new(
            build_chain(3, false),
            BootstrapThreshold::Homogeneous { value: 2 },
            BootstrapSeeds::Density { p: 0.5 },
            BootstrapDynamics::QueueClosure,
            0,
        );
        assert!(matches!(no_steps, Err(BootstrapError::NoSteps)));
    }
}
