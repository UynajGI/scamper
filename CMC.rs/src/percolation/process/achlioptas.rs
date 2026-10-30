//! G2: Achlioptas processes — rule-driven competitive edge addition (roadmap
//! §13).
//!
//! An Achlioptas process turns ordinary random edge addition into a
//! competition: every step proposes a fixed number of candidate edges, scores
//! each against the CURRENT component structure, and adds exactly one. The
//! canonical choice (Achlioptas 2000; "explosive percolation" literature)
//! scores by suppressing the giant — the candidate with the SMALLEST score
//! wins — with the product rule `|C(u)| * |C(v)|` and the sum rule
//! `|C(u)| + |C(v)|` as the two named rules. Riordan and Warnke,
//! "Achlioptas Process Phase Transitions Are Continuous", Science 335,
//! 1189-1191 (2012) (arXiv:1102.5306, archived in `docs/ref/`) proved that
//! such processes, including the product rule, have CONTINUOUS phase
//! transitions; the statistical anchor gate in
//! `tests/physics/achlioptas_process.rs` checks exactly that signature
//! (finite-cluster susceptibility growing with N, no O(1) jump in the
//! largest-component curve) on fixed sparse Erdos-Renyi graphs. The paper
//! itself describes both candidate-sampling conventions this module exposes:
//! edges "chosen ... uniformly at random from all possible edges [or] from
//! those edges not present" — with and without replacement.
//!
//! ## Physical parameters (all explicit, none defaultable away)
//!
//! - **Candidate count** `m >= 1`: candidates proposed per step.
//! - **Sampling** [`CandidateSampling`]: `WithoutReplacement` draws `m`
//!   distinct edges from the not-yet-added pool (partial Fisher-Yates);
//!   `WithReplacement` draws `m` iid uniform edges from the whole universe,
//!   so a candidate may repeat within a step or be already added.
//! - **Edge universe** [`EdgeUniverse`]: all physical edges, or an explicit
//!   subset of dense edge indices (sorted, deduplicated at construction).
//! - **Rule** [`AchlioptasRule`]: `Product` (`s_u * s_v`) or `Sum`
//!   (`s_u + s_v`), evaluated as `u128` over the two endpoint component
//!   sizes. A closure `score(size_a, size_b)` is deliberately insufficient
//!   (roadmap §13): sampling, universe, and tie policy are process
//!   parameters, not function arguments.
//! - **Tie policy** [`TiePolicy`]: a tie is two candidate edges with EQUAL
//!   score in one step. `LowestEdgeIndex` (default convention) selects the
//!   tied candidate with the lowest dense edge index; under
//!   `WithReplacement` a repeated edge is one candidate. `EarliestDraw`
//!   selects the tied candidate that appeared earliest in the step's draw
//!   order.
//!
//! ## Process definition and null steps
//!
//! Every vertex is active from step zero (this is an edge-addition process on
//! a fixed vertex set; the initial state is V isolates, largest component
//! size 1). Selection is always minimum-score. Under `WithoutReplacement`
//! every candidate is fresh, so every step adds an edge and the process ends
//! after exactly `|U|` steps (`U` = universe). Under `WithReplacement` an
//! already-added edge MAY win its step (its score `s*s` resp. `2*s` competes
//! with the others): such a **null step** consumes its draws, advances the
//! step counter, and adds nothing; the process terminates when every
//! permitted edge is added (a.s., coupon-collector time). The candidate pool
//! is NEVER filtered or redrawn — filtering would silently change the law
//! into a different process.
//!
//! The terminal state contains every permitted edge: the adapter runs to
//! complete occupation of the universe (roadmap gate).
//!
//! ## RNG draw contract (bitwise replay)
//!
//! Per realization, with `U = |universe|`, `L` = edges not yet added before
//! a step, and `k = m` (`WithReplacement`) resp. `k = min(m, L)`
//! (`WithoutReplacement`):
//!
//! 1. Initialization consumes NO randomness: vertices activate in ascending
//!    dense order `0..V`.
//! 2. Each step consumes exactly `k` draws, in candidate order
//!    `j = 0..k-1`, each one `rng.random_range(range)` call:
//!    - `WithReplacement`: `random_range(0..U)` every draw (`U` static for
//!      the whole run), mapped through the sorted universe table.
//!    - `WithoutReplacement`: the `j`-th draw is `random_range(0..L-j)`,
//!      selecting from the pool and swapping the selected edge to pool slot
//!      `L-1-j` (partial Fisher-Yates into the pool tail).
//! 3. After selection, the winner leaves the without-replacement pool by
//!    swapping ITS tail slot (`L-1-r` for draw position `r`) with the LAST
//!    pool slot `L-1` and shrinking the logical length to `L-1`; the
//!    `k-1` losers stay. (This compaction is part of the contract: an
//!    independent implementation that reproduces the draws and this
//!    compaction reproduces the sweep's accepted sequence bit-for-bit.)
//! 4. Nothing else consumes randomness.
//!
//! The pool is an unordered set from the draw's perspective (positions are
//! uniform), so pool ARRANGEMENT is internal: [`Self::step_with_candidates`]
//! replays a recorded candidate trace with identical accepted-edge
//! sequences, curves, and pool CONTENT, though the pool order may differ
//! from the drawing sweep that produced the trace.
//!
//! ## Output (C0)
//!
//! One realization is ONE vector sample over the occupation-count axis
//! (never a collection of independent per-step scalars):
//!
//! - `AchlioptasLargestSize[i]`: largest active component size at checkpoint
//!   `i` (with zero occupied edges this is 1 — the V active isolates).
//! - `AchlioptasOccupiedEdges[i]`: checkpoint coordinates (occupied-edge
//!   counts; constant across samples, so merged rows are self-describing).
//! - Scalars: `AchlioptasGiantReached` (0/1, always) and
//!   `AchlioptasStepsToGiant` (only when reached) — the percolation
//!   criterion is EXPLICIT: the first accepted edge whose summary satisfies
//!   `2 * largest_component_size > V` ("exceeds half of V"), recorded as the
//!   occupied-edge count at that acceptance. Plus `AchlioptasSteps` (proposal
//!   steps incl. null steps) and `AchlioptasNullSteps` (0 under
//!   `WithoutReplacement`).
//!
//! Provenance: callers record rule/sampling/tie/grid on the run's `Results`
//! metadata via the `as_provenance` strings.
//!
//! ## Efficiency
//!
//! The process drives the G1 [`IncrementalUnionWorkspace`] (prepare once per
//! topology; `reset` between realizations) and scores candidates through
//! `component_size` — the substrate's designated scoring primitive. All
//! buffers (pool, candidates, accepted order, curve, coordinates, added
//! flags) are reserved at construction; zero steady-state allocation
//! (asserted by the allocation probe).

use core::fmt;
use std::collections::TryReserveError;

use carlo_rs::{Context, MonteCarlo};
use rand::{Rng, RngExt};
use rand_xoshiro::Xoshiro256PlusPlus;

use super::{IncrementalError, IncrementalUnionWorkspace};
use crate::{CsrLattice, EdgeId, GraphView, UndirectedGraphView, VertexId};

/// The competitive scoring rule applied to a candidate edge `(u, v)`, where
/// `s_u`/`s_v` are the sizes of the endpoints' current components.
///
/// Selection is minimum-score in both variants (giant suppression).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AchlioptasRule {
    /// Score `s_u * s_v` (the classic "product rule").
    Product,
    /// Score `s_u + s_v` (the "sum rule").
    Sum,
}

impl AchlioptasRule {
    /// Stable provenance string for `Results` metadata attributes.
    pub const fn as_provenance(self) -> &'static str {
        match self {
            Self::Product => "product",
            Self::Sum => "sum",
        }
    }
}

/// Scoring rule for the two competing candidate edges.
///
/// The process scores every candidate from its endpoints' component sizes and
/// accepts the minimum score. Built-ins use exact `u128` arithmetic; a custom
/// rule must be deterministic (the process spends no RNG on scoring).
pub trait AchlioptasScore {
    /// Score a candidate joining components of sizes `size_a` and `size_b`.
    fn score(&self, size_a: usize, size_b: usize) -> u128;

    /// Stable provenance name.
    fn name(&self) -> &str;
}

impl AchlioptasScore for AchlioptasRule {
    #[inline]
    fn score(&self, size_a: usize, size_b: usize) -> u128 {
        match self {
            Self::Product => (size_a as u128) * (size_b as u128),
            Self::Sum => (size_a as u128) + (size_b as u128),
        }
    }

    fn name(&self) -> &str {
        self.as_provenance()
    }
}

/// How the `m` candidate edges of one step are drawn from the universe.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CandidateSampling {
    /// `m` distinct edges drawn from the not-yet-added pool (partial
    /// Fisher-Yates): every step adds a fresh edge, and the process ends
    /// after exactly `|U|` steps.
    WithoutReplacement,
    /// `m` iid uniform draws from the whole (static) universe: candidates may
    /// repeat within a step and may be already added; winning with an
    /// already-added edge is a null step (nothing added, draws consumed).
    WithReplacement,
}

impl CandidateSampling {
    /// Stable provenance string for `Results` metadata attributes.
    pub const fn as_provenance(self) -> &'static str {
        match self {
            Self::WithoutReplacement => "without-replacement",
            Self::WithReplacement => "with-replacement",
        }
    }
}

/// Resolution of a tie: two candidate edges with equal score in one step.
///
/// See the module documentation for the exact tie definition.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TiePolicy {
    /// Select the tied candidate with the lowest dense edge index
    /// (deterministic; identical winners for identical candidate sets
    /// regardless of draw order). The default convention.
    LowestEdgeIndex,
    /// Select the tied candidate that appeared earliest in the step's draw
    /// order.
    EarliestDraw,
}

impl TiePolicy {
    /// Stable provenance string for `Results` metadata attributes.
    pub const fn as_provenance(self) -> &'static str {
        match self {
            Self::LowestEdgeIndex => "lowest-edge-index",
            Self::EarliestDraw => "earliest-draw",
        }
    }
}

/// Tie-breaking rule when two candidates share the minimum score.
///
/// The process always compares a later candidate against the best-so-far
/// candidate, so a rule only decides whether the later one wins. Built-ins
/// reproduce the two documented reference orders.
pub trait AchlioptasTie {
    /// Whether the candidate at dense index `candidate` with `score` beats
    /// the best-so-far candidate `best` with `best_score`.
    fn beats(&self, candidate: usize, score: u128, best: usize, best_score: u128) -> bool;

    /// Stable provenance name.
    fn name(&self) -> &str;
}

impl AchlioptasTie for TiePolicy {
    #[inline]
    fn beats(&self, candidate: usize, score: u128, best: usize, best_score: u128) -> bool {
        match self {
            Self::LowestEdgeIndex => {
                score < best_score || (score == best_score && candidate < best)
            }
            // Strictly better only: on equal score the earlier draw stays
            // the winner.
            Self::EarliestDraw => score < best_score,
        }
    }

    fn name(&self) -> &str {
        self.as_provenance()
    }
}

/// The set of physical edges the process is permitted to add.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum EdgeUniverse {
    /// Every physical edge of the driving lattice (dense `0..E`).
    AllEdges,
    /// An explicit subset of dense edge indices. Sorted and deduplicated at
    /// construction; out-of-domain or duplicated indices are typed errors.
    Subset(Vec<usize>),
}

/// Where along the occupation-count axis a checkpoint is recorded (the
/// [`crate::percolation::AttackGrid`] convention, transplanted to edge
/// occupation).
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum AchlioptasGrid {
    /// A checkpoint after every accepted edge (plus the empty occupation):
    /// `|U| + 1` points.
    Full,
    /// Checkpoints at `0, step, 2*step, ...` occupied edges, always ending at
    /// `|U|` (the final interval may be shorter than `step`).
    Every { step: usize },
    /// Checkpoints at exactly these occupied-edge counts — no endpoints are
    /// forced. Values are deduplicated and sorted; every value must lie
    /// within `0..=|U|` and the grid must be non-empty.
    At { edges: Vec<usize> },
}

impl AchlioptasGrid {
    /// Checkpoint occupied-edge counts for a universe of `permitted` edges.
    fn checkpoints(&self, permitted: usize) -> Result<Vec<usize>, AchlioptasError> {
        match self {
            Self::Full => Ok((0..=permitted).collect()),
            Self::Every { step } => {
                if *step == 0 {
                    return Err(AchlioptasError::InvalidGrid("step must be >= 1"));
                }
                let mut points: Vec<usize> = (0..=permitted).step_by(*step).collect();
                if points.last() != Some(&permitted) {
                    points.push(permitted);
                }
                Ok(points)
            }
            Self::At { edges } => {
                if edges.is_empty() {
                    return Err(AchlioptasError::InvalidGrid(
                        "at least one checkpoint is required",
                    ));
                }
                let mut points = edges.clone();
                points.sort_unstable();
                points.dedup();
                if let Some(&edges) = points.last() {
                    if edges > permitted {
                        return Err(AchlioptasError::GridDomain { edges, permitted });
                    }
                }
                Ok(points)
            }
        }
    }
}

/// Result of driving one process step from a caller-supplied candidate list.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct AchlioptasStepRecord {
    /// Dense index of the selected (winning) candidate edge.
    pub selected_edge: usize,
    /// Whether the winner was fresh (an edge was added). `false` marks a null
    /// step: only possible under `WithReplacement`.
    pub added_edge: bool,
}

/// Owned MonteCarlo adapter for Achlioptas processes on an arbitrary
/// [`CsrLattice`].
///
/// One `sweep()` is one complete realization (to full occupation of the
/// universe); the following `measure()` records the curve as ONE C0 vector
/// sample plus the scalars. Configure `thermalization_sweeps = 0`;
/// realizations are independent given their seeds. Drive single steps from
/// an external candidate source (deterministic replay, traced diagnostics)
/// with [`Self::step_with_candidates`].
pub struct AchlioptasProcessMC<
    Rule: AchlioptasScore = AchlioptasRule,
    Tie: AchlioptasTie = TiePolicy,
> {
    lattice: CsrLattice,
    rule: Rule,
    sampling: CandidateSampling,
    candidates_per_step: usize,
    tie: Tie,
    /// Dense edge indices of the universe, ascending.
    universe: Vec<usize>,
    /// Sorted occupied-edge counts (the recording grid).
    checkpoints: Vec<usize>,
    workspace: IncrementalUnionWorkspace,
    /// Remaining-edge pool (dense indices, unordered); without-replacement
    /// draws and bookkeeping only.
    pool: Vec<usize>,
    /// Logical pool length (`pool[0..pool_len]` = not-yet-added edges).
    pool_len: usize,
    /// Candidate dense indices of the current sweep step.
    candidates: Vec<usize>,
    /// Added flag per physical edge (byte per edge; O(E) per realization).
    edge_added: Vec<u8>,
    /// Dense indices of accepted edges, in acceptance order (diagnostics,
    /// replay gate, terminal gate).
    accepted: Vec<usize>,
    /// Largest-component curve over the grid, refilled every realization.
    largest_curve: Vec<f64>,
    /// Checkpoint coordinates (occupied-edge counts), fixed at construction.
    coordinates: Vec<f64>,
    vertices: usize,
    /// Next grid index awaiting its checkpoint record.
    next_checkpoint: usize,
    /// Occupied-edge count of the current realization.
    added: usize,
    steps: usize,
    null_steps: usize,
    /// First occupied-edge count with `2 * largest > V` (the explicit
    /// percolation criterion), `None` until then.
    steps_to_giant: Option<usize>,
}

impl<Rule: AchlioptasScore, Tie: AchlioptasTie> AchlioptasProcessMC<Rule, Tie> {
    /// Build the adapter over an owned lattice.
    ///
    /// The trajectory budget of the underlying G1 workspace is `V + E`
    /// (V initialization activations plus one entry per accepted universe
    /// edge; `E >= |U|` always).
    pub fn new(
        lattice: CsrLattice,
        rule: Rule,
        sampling: CandidateSampling,
        candidates_per_step: usize,
        tie: Tie,
        universe: EdgeUniverse,
        grid: AchlioptasGrid,
    ) -> Result<Self, AchlioptasError> {
        lattice
            .validate()
            .map_err(|error| AchlioptasError::Topology(error.to_string()))?;
        let vertices = lattice.vertex_count();
        let edge_count = lattice.edge_count();
        if vertices == 0 {
            return Err(AchlioptasError::EmptyGraph);
        }
        if candidates_per_step == 0 {
            return Err(AchlioptasError::InvalidCandidateCount {
                requested: candidates_per_step,
            });
        }
        let universe = match universe {
            EdgeUniverse::AllEdges => (0..edge_count).collect(),
            EdgeUniverse::Subset(mut indices) => {
                if let Some(&edge) = indices.iter().find(|&&edge| edge >= edge_count) {
                    return Err(AchlioptasError::UniverseEdgeOutOfBounds {
                        requested: edge,
                        edges: edge_count,
                    });
                }
                indices.sort_unstable();
                if let Some(position) = indices.windows(2).position(|window| window[0] == window[1])
                {
                    return Err(AchlioptasError::DuplicateUniverseEdge {
                        edge: indices[position],
                    });
                }
                indices
            }
        };
        if universe.is_empty() {
            return Err(AchlioptasError::EmptyUniverse);
        }
        let permitted = universe.len();
        let checkpoints = grid.checkpoints(permitted)?;
        let mut workspace = IncrementalUnionWorkspace::new();
        workspace.prepare(vertices, edge_count)?;
        let mut pool = Vec::new();
        pool.try_reserve_exact(permitted)?;
        let mut candidates = Vec::new();
        candidates.try_reserve_exact(candidates_per_step)?;
        candidates.resize(candidates_per_step, 0);
        let mut edge_added = Vec::new();
        edge_added.try_reserve_exact(edge_count)?;
        edge_added.resize(edge_count, 0);
        let mut accepted = Vec::new();
        accepted.try_reserve_exact(permitted)?;
        let mut largest_curve = Vec::new();
        largest_curve.try_reserve_exact(checkpoints.len())?;
        largest_curve.resize(checkpoints.len(), 0.0);
        let coordinates = checkpoints.iter().map(|&edges| edges as f64).collect();
        let mut process = Self {
            lattice,
            rule,
            sampling,
            candidates_per_step,
            tie,
            universe,
            checkpoints,
            workspace,
            pool,
            pool_len: 0,
            candidates,
            edge_added,
            accepted,
            largest_curve,
            coordinates,
            vertices,
            next_checkpoint: 0,
            added: 0,
            steps: 0,
            null_steps: 0,
            steps_to_giant: None,
        };
        // Stand at step zero (V active isolates, empty occupation) so
        // `step_with_candidates` is legal on a fresh adapter.
        process.begin_realization();
        Ok(process)
    }

    /// The owned lattice the process runs on.
    pub const fn lattice(&self) -> &CsrLattice {
        &self.lattice
    }

    /// The competitive scoring rule.
    pub const fn rule(&self) -> &Rule {
        &self.rule
    }

    /// The candidate sampling mode.
    pub const fn sampling(&self) -> CandidateSampling {
        self.sampling
    }

    /// Candidates proposed per step.
    pub const fn candidates_per_step(&self) -> usize {
        self.candidates_per_step
    }

    /// The tie policy.
    pub const fn tie(&self) -> &Tie {
        &self.tie
    }

    /// Universe size (number of permitted edges).
    pub fn permitted_edge_count(&self) -> usize {
        self.universe.len()
    }

    /// Occupied-edge counts recorded per realization (the grid).
    pub fn checkpoints(&self) -> &[usize] {
        &self.checkpoints
    }

    /// Largest-component values of the last realization, aligned with
    /// [`Self::checkpoints`].
    pub fn largest_curve(&self) -> &[f64] {
        &self.largest_curve
    }

    /// Checkpoint coordinates (occupied-edge counts as `f64`; constant
    /// across realizations).
    pub fn occupied_coordinates(&self) -> &[f64] {
        &self.coordinates
    }

    /// Dense indices of accepted edges, in acceptance order, of the last
    /// realization.
    pub fn accepted_edges(&self) -> &[usize] {
        &self.accepted
    }

    /// Edges not yet added (the without-replacement pool). Empty under
    /// `WithReplacement` (that mode draws from the static universe; consult
    /// [`Self::accepted_edges`] and [`Self::added_edge_count`] instead).
    pub fn remaining_edges(&self) -> &[usize] {
        &self.pool[..self.pool_len]
    }

    /// Edges accepted so far in the current realization.
    pub const fn added_edge_count(&self) -> usize {
        self.added
    }

    /// Proposal steps consumed by the current realization (null steps
    /// included).
    pub const fn steps(&self) -> usize {
        self.steps
    }

    /// Steps of the current realization that added nothing (always 0 under
    /// `WithoutReplacement`).
    pub const fn null_steps(&self) -> usize {
        self.null_steps
    }

    /// First occupied-edge count whose summary satisfied `2 * largest > V`,
    /// or `None` when the giant was not reached (possible only when the
    /// universe's graph keeps the largest component at most `V / 2`).
    pub const fn steps_to_giant(&self) -> Option<usize> {
        self.steps_to_giant
    }

    /// The underlying G1 workspace (trajectory access for callers driving
    /// the adapter directly).
    pub const fn workspace(&self) -> &IncrementalUnionWorkspace {
        &self.workspace
    }

    /// Reset to the initial state (V active isolates, empty occupation) and
    /// run one complete realization: steps until every permitted edge is
    /// added, drawing candidates per the module's RNG contract.
    ///
    /// Deterministic given the RNG state. Exposed for tests and for callers
    /// running the process without a Carlo.rs context.
    pub fn run_realization<R: Rng + ?Sized>(&mut self, rng: &mut R) {
        self.begin_realization();
        while self.added < self.universe.len() {
            let count = self.draw_candidates(rng);
            let (winner, draw_position) = self.select_winner(count);
            self.commit_with_hint(winner, draw_position, count);
        }
        debug_assert_eq!(self.next_checkpoint, self.checkpoints.len());
        debug_assert_eq!(self.accepted.len(), self.universe.len());
    }

    /// Drive ONE process step from a caller-supplied candidate list,
    /// bypassing the RNG entirely (deterministic replay of a recorded trace,
    /// or external candidate sources).
    ///
    /// Contract: the list must be non-empty and no longer than
    /// [`Self::candidates_per_step`]; every candidate must be a universe
    /// member; under `WithoutReplacement` every candidate must additionally
    /// be not-yet-added (that mode never proposes occupied edges — its pool
    /// is exactly the remaining set). Returns the selected edge and whether
    /// it was added (null step under `WithReplacement`).
    ///
    /// The step consumes no randomness; pool maintenance removes the winner
    /// by search (O(pool)), so replayed realizations reproduce the drawing
    /// sweep's accepted sequence, curve, and pool content (see the module
    /// documentation for the arrangement caveat).
    pub fn step_with_candidates(
        &mut self,
        candidates: &[usize],
    ) -> Result<AchlioptasStepRecord, AchlioptasError> {
        if candidates.is_empty() {
            return Err(AchlioptasError::EmptyCandidateStep);
        }
        if candidates.len() > self.candidates_per_step {
            return Err(AchlioptasError::TooManyCandidates {
                offered: candidates.len(),
                maximum: self.candidates_per_step,
            });
        }
        for &edge in candidates {
            if self.universe.binary_search(&edge).is_err() {
                return Err(AchlioptasError::EdgeNotInUniverse { edge });
            }
            if matches!(self.sampling, CandidateSampling::WithoutReplacement)
                && self.edge_added[edge] != 0
            {
                return Err(AchlioptasError::EdgeAlreadyAdded { edge });
            }
        }
        let (winner, _) = self.select_winner_from(candidates);
        if matches!(self.sampling, CandidateSampling::WithoutReplacement) {
            self.remove_from_pool_by_search(winner)?;
        }
        let added_edge = self.finish_commit(winner);
        Ok(AchlioptasStepRecord {
            selected_edge: winner,
            added_edge,
        })
    }

    /// Reset every realization structure, re-activate the V vertices (dense
    /// ascending, step coordinate 0, no randomness), and record the empty-
    /// occupation checkpoint when the grid contains it.
    fn begin_realization(&mut self) {
        self.workspace.reset();
        self.edge_added.fill(0);
        self.accepted.clear();
        // The pool exists only for without-replacement draws; under
        // with-replacement the draws index the static universe table and
        // `remaining_edges()` stays empty (its documented contract).
        self.pool_len = match self.sampling {
            CandidateSampling::WithoutReplacement => {
                self.pool.clear();
                self.pool.extend_from_slice(&self.universe);
                self.universe.len()
            }
            CandidateSampling::WithReplacement => 0,
        };
        self.added = 0;
        self.steps = 0;
        self.null_steps = 0;
        self.steps_to_giant = None;
        for index in 0..self.vertices {
            self.workspace
                .activate_vertex(&self.lattice, VertexId::from_valid_index(index), 0)
                .expect("initialization activates each vertex exactly once after reset");
        }
        self.next_checkpoint = 0;
        if self.checkpoints.first() == Some(&0) {
            self.largest_curve[0] = self.workspace.summary().largest_component_size as f64;
            self.next_checkpoint = 1;
        }
    }

    /// Draw this step's candidates per the RNG contract into
    /// `self.candidates[0..count]`, returning `count`.
    ///
    /// `WithoutReplacement`: partial Fisher-Yates — draw `j` selects
    /// `random_range(0..L-j)` and swaps the selected edge to pool slot
    /// `L-1-j`, so after the draws the candidates sit in the pool tail in
    /// reverse draw order. `WithReplacement`: `m` iid
    /// `random_range(0..U)` draws through the static universe table.
    fn draw_candidates<R: Rng + ?Sized>(&mut self, rng: &mut R) -> usize {
        match self.sampling {
            CandidateSampling::WithoutReplacement => {
                let count = self.candidates_per_step.min(self.pool_len);
                for draw in 0..count {
                    let remaining = self.pool_len - draw;
                    let index = rng.random_range(0..remaining);
                    self.candidates[draw] = self.pool[index];
                    self.pool.swap(index, remaining - 1);
                }
                count
            }
            CandidateSampling::WithReplacement => {
                let permitted = self.universe.len();
                for draw in 0..self.candidates_per_step {
                    let position = rng.random_range(0..permitted);
                    self.candidates[draw] = self.universe[position];
                }
                self.candidates_per_step
            }
        }
    }

    /// Score and select over `self.candidates[0..count]` (the sweep path).
    fn select_winner(&self, count: usize) -> (usize, usize) {
        self.select_winner_from(&self.candidates[..count])
    }

    /// Score every candidate edge by the rule over its endpoints' component
    /// sizes (the G1 `component_size` scoring primitive) and select the
    /// minimum-score winner under the tie policy.
    ///
    /// Returns `(dense edge index, position in the candidate list)`.
    fn select_winner_from(&self, candidates: &[usize]) -> (usize, usize) {
        let mut best = 0_usize;
        let mut best_score = self.score_of(candidates[best]);
        for candidate in 1..candidates.len() {
            let score = self.score_of(candidates[candidate]);
            let better = self
                .tie
                .beats(candidates[candidate], score, candidates[best], best_score);
            if better {
                best = candidate;
                best_score = score;
            }
        }
        (candidates[best], best)
    }

    #[inline]
    fn score_of(&self, edge: usize) -> u128 {
        let [left, right] = self.lattice.edge_endpoints(EdgeId::from_valid_index(edge));
        let size_left = self
            .workspace
            .component_size(left)
            .expect("every vertex is active in an Achlioptas realization");
        let size_right = self
            .workspace
            .component_size(right)
            .expect("every vertex is active in an Achlioptas realization");
        self.rule.score(size_left, size_right)
    }

    /// Commit a sweep step whose winner sits in the pool tail (the partial
    /// Fisher-Yates arrangement: the candidate at draw position
    /// `draw_position` occupies pool slot `L - 1 - draw_position`). The
    /// winner leaves via a swap with the LAST pool slot `L - 1` — the one
    /// slot that exits the logical pool when it shrinks to `L - 1` — so
    /// every loser stays and the remaining set loses exactly the winner.
    fn commit_with_hint(&mut self, winner: usize, draw_position: usize, count: usize) {
        if matches!(self.sampling, CandidateSampling::WithoutReplacement) {
            debug_assert!(count >= 1);
            let winner_slot = self.pool_len - 1 - draw_position;
            debug_assert_eq!(self.pool[winner_slot], winner);
            let last = self.pool_len - 1;
            self.pool.swap(winner_slot, last);
            self.pool_len -= 1;
        }
        self.finish_commit(winner);
    }

    /// Remove `winner` from the without-replacement pool by search (the
    /// replay path: tail positions of a drawing sweep are unknown here).
    fn remove_from_pool_by_search(&mut self, winner: usize) -> Result<(), AchlioptasError> {
        let position = self.pool[..self.pool_len]
            .iter()
            .position(|&edge| edge == winner)
            .ok_or(AchlioptasError::PoolInconsistent { edge: winner })?;
        let last = self.pool_len - 1;
        self.pool.swap(position, last);
        self.pool_len -= 1;
        Ok(())
    }

    /// Activate the winner (or book a null step), update counters, the
    /// giant criterion, and every checkpoint the acceptance reaches.
    /// Returns whether an edge was added.
    fn finish_commit(&mut self, winner: usize) -> bool {
        self.steps += 1;
        if self.edge_added[winner] != 0 {
            self.null_steps += 1;
            return false;
        }
        let summary = self
            .workspace
            .activate_edge(
                &self.lattice,
                EdgeId::from_valid_index(winner),
                self.added + 1,
            )
            .expect("adapter activates only universe edges that are not yet added");
        self.edge_added[winner] = 1;
        self.added += 1;
        self.accepted.push(winner);
        if self.steps_to_giant.is_none() && 2 * summary.largest_component_size > self.vertices {
            self.steps_to_giant = Some(self.added);
        }
        while self.next_checkpoint < self.checkpoints.len()
            && self.checkpoints[self.next_checkpoint] <= self.added
        {
            self.largest_curve[self.next_checkpoint] = summary.largest_component_size as f64;
            self.next_checkpoint += 1;
        }
        true
    }
}

impl<Rule: AchlioptasScore, Tie: AchlioptasTie> MonteCarlo for AchlioptasProcessMC<Rule, Tie> {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        self.run_realization(&mut ctx.rng);
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        ctx.measure_array("AchlioptasLargestSize", &self.largest_curve)
            .expect("checkpoint count is fixed by the grid");
        ctx.measure_array("AchlioptasOccupiedEdges", &self.coordinates)
            .expect("checkpoint count is fixed by the grid");
        ctx.measure(
            "AchlioptasGiantReached",
            f64::from(u8::from(self.steps_to_giant.is_some())),
        );
        if let Some(steps) = self.steps_to_giant {
            ctx.measure("AchlioptasStepsToGiant", steps as f64);
        }
        ctx.measure("AchlioptasSteps", self.steps as f64);
        ctx.measure("AchlioptasNullSteps", self.null_steps as f64);
    }

    fn name(&self) -> &'static str {
        "Achlioptas"
    }
}

/// Failure of an Achlioptas-process operation.
#[derive(Debug)]
pub enum AchlioptasError {
    /// A buffer could not be reserved.
    Capacity(TryReserveError),
    /// The driving lattice failed validation.
    Topology(String),
    /// The lattice has no vertices.
    EmptyGraph,
    /// The edge universe is empty (no permitted edges).
    EmptyUniverse,
    /// The candidate count was zero.
    InvalidCandidateCount { requested: usize },
    /// A universe subset edge lies outside the lattice's edge domain.
    UniverseEdgeOutOfBounds { requested: usize, edges: usize },
    /// A universe subset listed the same edge twice.
    DuplicateUniverseEdge { edge: usize },
    /// The recording grid was malformed (empty, or `Every { step: 0 }`).
    InvalidGrid(&'static str),
    /// A grid checkpoint lies beyond the universe.
    GridDomain { edges: usize, permitted: usize },
    /// A caller-supplied candidate step was empty.
    EmptyCandidateStep,
    /// A caller-supplied candidate step exceeded the candidate count.
    TooManyCandidates { offered: usize, maximum: usize },
    /// A caller-supplied candidate is not a universe member.
    EdgeNotInUniverse { edge: usize },
    /// A without-replacement candidate step included an already-added edge.
    EdgeAlreadyAdded { edge: usize },
    /// The winner of a without-replacement step was absent from the pool (an
    /// internal invariant failure).
    PoolInconsistent { edge: usize },
}

impl From<TryReserveError> for AchlioptasError {
    fn from(source: TryReserveError) -> Self {
        Self::Capacity(source)
    }
}

impl From<IncrementalError> for AchlioptasError {
    fn from(source: IncrementalError) -> Self {
        match source {
            IncrementalError::Capacity(error) => Self::Capacity(error),
            // Workspace-level domain failures cannot occur from this
            // adapter (sizes validated at construction); surface them as
            // topology errors rather than silently re-labelling.
            other => Self::Topology(other.to_string()),
        }
    }
}

impl fmt::Display for AchlioptasError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Capacity(_) => {
                write!(formatter, "could not reserve Achlioptas-process buffers")
            }
            Self::Topology(message) => write!(formatter, "invalid lattice: {message}"),
            Self::EmptyGraph => {
                write!(formatter, "Achlioptas processes need at least one vertex")
            }
            Self::EmptyUniverse => {
                write!(formatter, "the edge universe is empty: no permitted edges")
            }
            Self::InvalidCandidateCount { requested } => {
                write!(formatter, "candidate count must be >= 1, got {requested}")
            }
            Self::UniverseEdgeOutOfBounds { requested, edges } => write!(
                formatter,
                "universe edge {requested} is outside the 0..{edges} edge domain"
            ),
            Self::DuplicateUniverseEdge { edge } => {
                write!(formatter, "universe lists edge {edge} more than once")
            }
            Self::InvalidGrid(message) => write!(formatter, "invalid recording grid: {message}"),
            Self::GridDomain { edges, permitted } => write!(
                formatter,
                "grid checkpoint {edges} exceeds the {permitted} permitted edges"
            ),
            Self::EmptyCandidateStep => {
                write!(formatter, "a candidate step needs at least one edge")
            }
            Self::TooManyCandidates { offered, maximum } => write!(
                formatter,
                "candidate step offered {offered} edges, more than the configured {maximum}"
            ),
            Self::EdgeNotInUniverse { edge } => {
                write!(formatter, "candidate edge {edge} is not in the universe")
            }
            Self::EdgeAlreadyAdded { edge } => write!(
                formatter,
                "edge {edge} is already added; without-replacement steps propose fresh edges only"
            ),
            Self::PoolInconsistent { edge } => write!(
                formatter,
                "internal inconsistency: winner {edge} missing from the remaining-edge pool"
            ),
        }
    }
}

impl std::error::Error for AchlioptasError {
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
    use crate::{build_chain, Bond, BondType};
    use rand::SeedableRng;

    /// 5 vertices, 5 edges incl. one self-loop, engineered for rule and tie
    /// exercises: (0,0) loop, (1,2), (3,4), (1,3), (2,4).
    fn toy_graph() -> CsrLattice {
        CsrLattice::try_from_edges(
            5,
            vec![
                Bond::new(0, 0, BondType::Generic, 1.0),
                Bond::new(1, 2, BondType::Generic, 1.0),
                Bond::new(3, 4, BondType::Generic, 1.0),
                Bond::new(1, 3, BondType::Generic, 1.0),
                Bond::new(2, 4, BondType::Generic, 1.0),
            ],
        )
        .expect("toy graph is valid")
    }

    fn process(m: usize, rule: AchlioptasRule, tie: TiePolicy) -> AchlioptasProcessMC {
        AchlioptasProcessMC::new(
            toy_graph(),
            rule,
            CandidateSampling::WithoutReplacement,
            m,
            tie,
            EdgeUniverse::AllEdges,
            AchlioptasGrid::Full,
        )
        .expect("toy process is valid")
    }

    #[test]
    fn one_step_scores_and_ties_are_hand_computable() {
        // Edges: 0=(0,0) loop, 1=(1,2), 2=(3,4), 3=(1,3), 4=(2,4).
        // First step from candidates {1, 2}: both join singleton pairs —
        // product 1*1 = 1, sum 1+1 = 2 for both. A TIE: LowestEdgeIndex
        // must pick edge 1 regardless of draw order; EarliestDraw picks
        // whichever candidate was listed first.
        let mut lowest = process(2, AchlioptasRule::Product, TiePolicy::LowestEdgeIndex);
        let record = lowest.step_with_candidates(&[2, 1]).unwrap();
        assert_eq!(record.selected_edge, 1, "tie -> lowest dense edge index");
        assert!(record.added_edge);

        let mut earliest = process(2, AchlioptasRule::Product, TiePolicy::EarliestDraw);
        let record = earliest.step_with_candidates(&[2, 1]).unwrap();
        assert_eq!(record.selected_edge, 2, "tie -> earliest draw position");

        let mut earliest_same = process(2, AchlioptasRule::Product, TiePolicy::EarliestDraw);
        let record = earliest_same.step_with_candidates(&[1, 2]).unwrap();
        assert_eq!(record.selected_edge, 1, "earliest draw of a swapped list");
    }

    #[test]
    fn product_and_sum_rules_disagree_by_hand() {
        // Hand-worked split: build components of sizes 5, 2, 3, 1, then
        // offer one candidate joining (5,1) and one joining (2,3).
        // Edge 7=(4,10): product 5*1 = 5, sum 5+1 = 6.
        // Edge 8=(6,7): product 2*3 = 6, sum 2+3 = 5.
        // The product rule accepts edge 7 (5 < 6); the sum rule accepts
        // edge 8 (5 < 6) — the same candidate set, different accepted edge.
        let lattice = CsrLattice::try_from_edges(
            11,
            vec![
                Bond::new(0, 1, BondType::Generic, 1.0), // 0  -> {0..4} size 5
                Bond::new(1, 2, BondType::Generic, 1.0), // 1
                Bond::new(2, 3, BondType::Generic, 1.0), // 2
                Bond::new(3, 4, BondType::Generic, 1.0), // 3
                Bond::new(5, 6, BondType::Generic, 1.0), // 4  -> {5,6} size 2
                Bond::new(7, 8, BondType::Generic, 1.0), // 5  -> {7,8,9} size 3
                Bond::new(8, 9, BondType::Generic, 1.0), // 6
                Bond::new(4, 10, BondType::Generic, 1.0), // 7: joins (5,1)
                Bond::new(6, 7, BondType::Generic, 1.0), // 8: joins (2,3)
            ],
        )
        .expect("rule-split graph is valid");
        let build = |rule| {
            AchlioptasProcessMC::new(
                lattice.clone(),
                rule,
                CandidateSampling::WithoutReplacement,
                2,
                TiePolicy::LowestEdgeIndex,
                EdgeUniverse::AllEdges,
                AchlioptasGrid::Full,
            )
            .unwrap()
        };
        let mut product = build(AchlioptasRule::Product);
        let mut sum = build(AchlioptasRule::Sum);
        for process in [&mut product, &mut sum] {
            for edge in 0..=6 {
                process.step_with_candidates(&[edge]).unwrap();
            }
            assert_eq!(process.workspace().summary().largest_component_size, 5);
        }
        let product_record = product.step_with_candidates(&[7, 8]).unwrap();
        let sum_record = sum.step_with_candidates(&[8, 7]).unwrap();
        assert_eq!(product_record.selected_edge, 7, "product: 5 < 6");
        assert_eq!(sum_record.selected_edge, 8, "sum: 5 < 6");
    }

    #[test]
    fn with_replacement_null_steps_are_recorded_not_filtered() {
        let mut process = AchlioptasProcessMC::new(
            build_chain(3, false),
            AchlioptasRule::Product,
            CandidateSampling::WithReplacement,
            2,
            TiePolicy::LowestEdgeIndex,
            EdgeUniverse::AllEdges,
            AchlioptasGrid::Full,
        )
        .unwrap();
        // Add edge 0, then propose it again alongside an equal-score fresh
        // edge 1: both candidates have score 1 here ONLY before edge 0 —
        // after edge 0=(0,1), candidate 0 scores (2,2): product 4, while
        // fresh edge 1=(1,2) scores (2,1): product 2. The fresh edge wins,
        // so force the null step by proposing the occupied edge alone.
        process.step_with_candidates(&[0]).unwrap();
        let record = process.step_with_candidates(&[0]).unwrap();
        assert_eq!(record.selected_edge, 0);
        assert!(!record.added_edge, "occupied winner is a null step");
        assert_eq!(process.null_steps(), 1);
        assert_eq!(process.steps(), 2);
        assert_eq!(process.added_edge_count(), 1);
        // A null step must leave the summary untouched.
        assert_eq!(process.workspace().summary().active_edge_count, 1);
    }

    #[test]
    fn constructor_rejects_malformed_configurations() {
        let graph = toy_graph();
        assert!(matches!(
            AchlioptasProcessMC::new(
                graph.clone(),
                AchlioptasRule::Product,
                CandidateSampling::WithoutReplacement,
                0,
                TiePolicy::LowestEdgeIndex,
                EdgeUniverse::AllEdges,
                AchlioptasGrid::Full,
            ),
            Err(AchlioptasError::InvalidCandidateCount { requested: 0 })
        ));
        assert!(matches!(
            AchlioptasProcessMC::new(
                graph.clone(),
                AchlioptasRule::Product,
                CandidateSampling::WithoutReplacement,
                1,
                TiePolicy::LowestEdgeIndex,
                EdgeUniverse::Subset(vec![5]),
                AchlioptasGrid::Full,
            ),
            Err(AchlioptasError::UniverseEdgeOutOfBounds {
                requested: 5,
                edges: 5
            })
        ));
        assert!(matches!(
            AchlioptasProcessMC::new(
                graph.clone(),
                AchlioptasRule::Product,
                CandidateSampling::WithoutReplacement,
                1,
                TiePolicy::LowestEdgeIndex,
                EdgeUniverse::Subset(vec![1, 1]),
                AchlioptasGrid::Full,
            ),
            Err(AchlioptasError::DuplicateUniverseEdge { edge: 1 })
        ));
        assert!(matches!(
            AchlioptasProcessMC::new(
                graph.clone(),
                AchlioptasRule::Product,
                CandidateSampling::WithoutReplacement,
                1,
                TiePolicy::LowestEdgeIndex,
                EdgeUniverse::AllEdges,
                AchlioptasGrid::At { edges: vec![6] },
            ),
            Err(AchlioptasError::GridDomain {
                edges: 6,
                permitted: 5
            })
        ));
        let mut fresh = AchlioptasProcessMC::new(
            build_chain(3, false),
            AchlioptasRule::Product,
            CandidateSampling::WithoutReplacement,
            2,
            TiePolicy::LowestEdgeIndex,
            EdgeUniverse::AllEdges,
            AchlioptasGrid::Full,
        )
        .unwrap();
        assert!(matches!(
            fresh.step_with_candidates(&[]),
            Err(AchlioptasError::EmptyCandidateStep)
        ));
        // Duplicate candidates are legal while the edge is fresh...
        assert_eq!(
            fresh.step_with_candidates(&[0, 0]).unwrap(),
            AchlioptasStepRecord {
                selected_edge: 0,
                added_edge: true
            }
        );
        // ...but a without-replacement step refuses already-added edges.
        assert!(matches!(
            fresh.step_with_candidates(&[0, 0]),
            Err(AchlioptasError::EdgeAlreadyAdded { edge: 0 })
        ));
        // And a candidate outside the universe is refused.
        assert!(matches!(
            fresh.step_with_candidates(&[1, 1]),
            Ok(AchlioptasStepRecord {
                selected_edge: 1,
                added_edge: true
            })
        ));
        let mut bounded = AchlioptasProcessMC::new(
            build_chain(3, false),
            AchlioptasRule::Product,
            CandidateSampling::WithReplacement,
            2,
            TiePolicy::LowestEdgeIndex,
            EdgeUniverse::AllEdges,
            AchlioptasGrid::Full,
        )
        .unwrap();
        assert!(matches!(
            bounded.step_with_candidates(&[0, 5]),
            Err(AchlioptasError::EdgeNotInUniverse { edge: 5 })
        ));
        assert!(matches!(
            bounded.step_with_candidates(&[0, 0, 1]),
            Err(AchlioptasError::TooManyCandidates {
                offered: 3,
                maximum: 2
            })
        ));
    }

    #[test]
    fn run_realization_reaches_the_complete_terminal_state() {
        let mut process = process(2, AchlioptasRule::Product, TiePolicy::LowestEdgeIndex);
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(0x00C2_4143);
        process.run_realization(&mut rng);
        let summary = process.workspace().summary();
        assert_eq!(summary.active_vertex_count, 5);
        assert_eq!(summary.active_edge_count, 5, "every permitted edge");
        // Vertex 0 only carries the self-loop: the terminal largest
        // component is {1,2,3,4}.
        assert_eq!(summary.largest_component_size, 4);
        let mut accepted = process.accepted_edges().to_vec();
        accepted.sort_unstable();
        assert_eq!(accepted, vec![0, 1, 2, 3, 4]);
        assert_eq!(process.steps(), 5, "without replacement: no null steps");
        assert_eq!(process.null_steps(), 0);
        assert_eq!(process.remaining_edges().len(), 0);
        // The largest curve starts at the V-isolate state and ends complete.
        assert_eq!(process.largest_curve().first(), Some(&1.0));
        assert_eq!(process.largest_curve().last(), Some(&4.0));
        // Giant criterion: the first acceptance with 2 * largest > 5
        // (largest >= 3). Recomputed exhaustively in the physics tests.
        let giant = process.steps_to_giant().expect("largest reaches 4 > 5/2");
        assert!((1..=5).contains(&giant));
    }
}
