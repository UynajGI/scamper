//! L1: long-range percolation power-law edge law on the G1 substrate
//! (roadmap v2 §4, stage L1).
//!
//! Anchor (archived
//! `docs/ref/md/Biskup_2003_long_range_chemical_distance.md`): Biskup,
//! *On the scaling of the chemical distance in long-range percolation*,
//! Ann. Probab. 32, 2938 (2004) — the edge law `p_xy = |x-y|^{-s+o(1)}`
//! over unordered pairs of distinct sites of Z^d, with `|·|` a norm
//! distance (Biskup takes ℓ2 "for definiteness; ... any other equivalent
//! norm will do", ~line 134). Stage L1 realizes the explicit power-law
//! member of that family on a finite box; the chemical-distance regimes
//! (L2) and diameter transitions (L3) build on this layer later.
//!
//! ## The law
//!
//! A box is a product of per-dimension side lengths (each >= 1); dense
//! vertex indices are the row-major lexicographic ranks of the coordinate
//! vectors. For explicit `beta >= 0` and `s >= 0` (0 is legal but
//! degenerate in each case: no edges, resp. a distance-independent
//! probability `min(1, beta)`), each unordered pair (x, y) of distinct box
//! vertices at distance r is occupied independently with
//!
//! ```text
//! p(r) = min(1, beta * r^(-s))    (r >= 1; the exact f64 expression is
//!                                  [`LongRangeBox::pair_probability`])
//! ```
//!
//! **Norm choice: ℓ1 (Manhattan).** Biskup's regime statements (s relative
//! to d and 2d) hold for any norm equivalent to ℓ2, and ℓ1 is one. ℓ1 is
//! integer-exact on Z^d — the distance is a `usize`, so `p(r)` depends on
//! floating rounding only through the single `powf` — and in d = 1 (the
//! stage-L1 target) ℓ1, ℓ2, and the coordinate difference coincide.
//!
//! **Boundary framing.** The pair universe is the finite-box restriction
//! of the infinite-measure law: every unordered pair {x, y} with BOTH
//! endpoints in the box and none outside — no edge leaves the box — each
//! independent with the same `p(r)`, i.e. the restriction of Biskup's
//! infinite product measure to the box's finite pair set. Boundary
//! effects are confined to the box surface: a vertex at depth k from the
//! boundary lacks exactly the partners the infinite lattice would offer
//! beyond it.
//!
//! ## Pair order and RNG draw contract
//!
//! **Consideration order**: lexicographic by dense (x, y) with x < y —
//! equivalently ascending pair index `x·n − x(x+1)/2 + (y − x − 1)`
//! ([`LongRangeBox::pair_index`]); dense order is coordinate-lexicographic
//! order, so this is the coordinate lexicographic order too. The draw
//! contract (bitwise reproducible on one platform): **one uniform `f64`
//! per pair, in ascending pair order** — saturated (p = 1) and
//! near-zero-probability pairs included, nothing skipped. The pair at
//! index k is occupied iff `uniforms[k] < p(r)`.
//!
//! **Saturation.** `p(r) = 1` exactly when `beta * r^{-s} >= 1`; in
//! particular every nearest-neighbor pair (r = 1) is surely occupied at
//! `beta >= 1` (`u < 1` holds for every `u` in `[0, 1)`), and a huge
//! `beta` saturates every pair up to radius `beta^{1/s}` — the complete
//! graph in the limit. For the accepted domains (finite `beta`, `s >= 0`,
//! `r >= 1`) the factor `r^{-s}` lies in `(0, 1]`, so `beta * r^{-s}`
//! never overflows past `beta` itself; `min(1, ·)` saturates everything
//! at or above 1. No (s, d) combination is rejected.
//!
//! ## G1 substrate drive (process coordinate = added-edge count)
//!
//! One realization is one incremental construction on the G1
//! [`IncrementalUnionWorkspace`] over a zero-materialization complete pair
//! view of the box: all `n` vertices are activated first (dense ascending,
//! step coordinate 0 — the G2 convention; the fixed vertex set is the
//! stage), then each occupied pair is activated in consideration order
//! with step coordinate **the 1-based added-edge count**. The realized
//! added-edge count is random (a sum of independent non-identical
//! Bernoulli draws); the trajectory steps are exactly `[0; n] ++ [1..=m]`
//! and each entry carries the substrate's online [`IncrementalSummary`].
//!
//! ## C0 recording (adapter)
//!
//! One realization is ONE element-preserving vector sample plus scalars
//! (the G3 padding convention for random-length growth): the
//! `LongRangeLargestSize` array has a fixed construction-time length
//! `steps`, entry `i` = largest component size after the (i+1)-th added
//! edge, and entries past the realized length repeat the final value
//! (1 when nothing was added — every vertex is active). The scalars carry
//! the true realized values: `LongRangeAddedEdges`,
//! `LongRangeComponentCount`, `LongRangeSecondMoment`. Configure
//! `thermalization_sweeps = 0`; realizations are i.i.d.

use std::collections::TryReserveError;

use carlo_rs::{Context, MonteCarlo};
use rand::{Rng, RngExt};
use rand_xoshiro::Xoshiro256PlusPlus;

use super::process::{IncrementalError, IncrementalUnionWorkspace, ProcessError};
use crate::{EdgeId, GraphView, Incidence, UndirectedGraphView, VertexId};

/// A finite box of Z^d: per-dimension side lengths and their derived
/// vertex/pair counts.
///
/// Dense vertex indices are the row-major lexicographic ranks of the
/// coordinate vectors; the pair universe is every unordered pair of
/// distinct vertices (the complete graph, never materialized — see the
/// module documentation for the law defined on it).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct LongRangeBox {
    sides: Vec<usize>,
    vertices: usize,
    pairs: usize,
}

impl LongRangeBox {
    /// Build the box; every side must be >= 1, and the vertex/pair counts
    /// must fit a `usize` (typed rejections otherwise).
    pub fn new(sides: &[usize]) -> Result<Self, LongRangeError> {
        if sides.is_empty() {
            return Err(LongRangeError::EmptyBox);
        }
        let mut product = 1_u128;
        for (axis, &side) in sides.iter().enumerate() {
            if side == 0 {
                return Err(LongRangeError::ZeroSide { axis });
            }
            product *= side as u128;
            if product > usize::MAX as u128 {
                return Err(LongRangeError::BoxTooLarge {
                    sides: sides.to_vec(),
                });
            }
        }
        let vertices = product as usize;
        let wide_pairs = vertices as u128 * (vertices as u128 - 1) / 2;
        if wide_pairs > usize::MAX as u128 {
            return Err(LongRangeError::PairUniverseTooLarge { vertices });
        }
        Ok(Self {
            sides: sides.to_vec(),
            vertices,
            pairs: wide_pairs as usize,
        })
    }

    /// Per-dimension side lengths.
    #[inline]
    pub fn sides(&self) -> &[usize] {
        &self.sides
    }

    /// The dimension d.
    #[inline]
    pub fn dimension(&self) -> usize {
        self.sides.len()
    }

    /// Number of vertices (product of the sides).
    #[inline]
    pub const fn vertices(&self) -> usize {
        self.vertices
    }

    /// Number of unordered vertex pairs (`vertices * (vertices - 1) / 2`).
    #[inline]
    pub const fn pairs(&self) -> usize {
        self.pairs
    }

    /// Dense vertex index of a coordinate vector, or `None` when the
    /// vector has the wrong dimension or an out-of-range coordinate.
    pub fn vertex_index(&self, coords: &[usize]) -> Option<usize> {
        if coords.len() != self.sides.len() {
            return None;
        }
        let mut index = 0_usize;
        for (&coordinate, &side) in coords.iter().zip(self.sides.iter()) {
            if coordinate >= side {
                return None;
            }
            index = index * side + coordinate;
        }
        Some(index)
    }

    /// Coordinate vector of a dense vertex index (row-major decomposition),
    /// or `None` outside `0..vertices`. Allocates — an analysis accessor,
    /// never on the realization path.
    pub fn coordinates(&self, index: usize) -> Option<Vec<usize>> {
        if index >= self.vertices {
            return None;
        }
        let mut coords = vec![0_usize; self.sides.len()];
        let mut remaining = index;
        for (slot, &side) in coords.iter_mut().zip(self.sides.iter().rev()) {
            *slot = remaining % side;
            remaining /= side;
        }
        Some(coords)
    }

    /// ℓ1 distance `sum_i |left_i - right_i|` between two dense vertex
    /// indices (see the module documentation for the norm choice; 0 for
    /// equal indices).
    pub fn distance(&self, left: usize, right: usize) -> usize {
        let (mut left, mut right) = (left, right);
        let mut total = 0_usize;
        for &side in self.sides.iter().rev() {
            total += (left % side).abs_diff(right % side);
            left /= side;
            right /= side;
        }
        total
    }

    /// Pair index of two distinct vertices in the documented consideration
    /// order (lexicographic (x, y), x < y), or `None` for equal vertices.
    pub fn pair_index(&self, left: usize, right: usize) -> Option<usize> {
        if left == right {
            return None;
        }
        let (low, high) = (left.min(right), left.max(right));
        let wide = low as u128 * self.vertices as u128 - low as u128 * (low as u128 + 1) / 2
            + (high - low - 1) as u128;
        Some(wide as usize)
    }

    /// Endpoints `(x, y)` with `x < y` of the pair at `pair`, or `None`
    /// outside `0..pairs` (the inverse of [`Self::pair_index`]).
    pub fn pair_endpoints(&self, pair: usize) -> Option<(usize, usize)> {
        if pair >= self.pairs {
            return None;
        }
        // base(a) = pairs whose first endpoint is below a; find the largest
        // a with base(a) <= pair by binary search (base is increasing on
        // 0..vertices-1), then the offset picks the second endpoint.
        let mut low = 0_usize;
        let mut high = self.vertices - 1;
        while low < high {
            let mid = (low + high).div_ceil(2);
            let base = mid as u128 * self.vertices as u128 - mid as u128 * (mid as u128 + 1) / 2;
            if base <= pair as u128 {
                low = mid;
            } else {
                high = mid - 1;
            }
        }
        let base = low as u128 * self.vertices as u128 - low as u128 * (low as u128 + 1) / 2;
        let second = low + 1 + (pair as u128 - base) as usize;
        Some((low, second))
    }

    /// The pair-occupation probability at ℓ1 distance `r >= 1` — the exact
    /// bitwise draw-contract expression (same-platform replay mirrors
    /// these bits; `min` saturates any product at or above 1).
    #[inline]
    pub fn pair_probability(distance: usize, beta: f64, s: f64) -> f64 {
        (beta * (distance as f64).powf(-s)).min(1.0)
    }
}

/// Zero-materialization complete pair view of a box: vertex k is the k-th
/// dense box vertex, edge k is the k-th pair in the documented
/// consideration order. Only the G1 substrate reads it (vertex activations
/// scan incidences for waiting edges; every edge endpoint lookup is O(1)).
struct PairUniverse<'a> {
    box_: &'a LongRangeBox,
}

impl GraphView for PairUniverse<'_> {
    fn vertex_count(&self) -> usize {
        self.box_.vertices()
    }
}

impl UndirectedGraphView for PairUniverse<'_> {
    fn edge_count(&self) -> usize {
        self.box_.pairs()
    }

    fn edge_endpoints(&self, edge: EdgeId) -> [VertexId; 2] {
        let (first, second) = self
            .box_
            .pair_endpoints(edge.index())
            .expect("edge id from this view");
        [
            VertexId::from_valid_index(first),
            VertexId::from_valid_index(second),
        ]
    }

    fn incidences(&self, vertex: VertexId) -> impl Iterator<Item = Incidence> + '_ {
        let me = vertex.index();
        (0..self.box_.vertices()).filter_map(move |other| {
            let pair = self.box_.pair_index(me, other)?;
            Some(Incidence {
                neighbor: VertexId::from_valid_index(other),
                edge: EdgeId::from_valid_index(pair),
            })
        })
    }
}

/// Scalar outcome of one long-range realization.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct LongRangeOutcome {
    /// Occupied pairs (the realized added-edge count; random per
    /// realization).
    pub added_edges: usize,
    /// Final component count (all box vertices active — isolates
    /// included).
    pub component_count: usize,
    /// Final largest-component size.
    pub largest_component_size: usize,
    /// Final raw second moment `sum(s_i^2)`.
    pub raw_second_moment: u128,
}

/// Reusable long-range workspace: the occupied-pair record over the driven
/// G1 substrate.
///
/// `prepare` once per box size (idempotent; every buffer reserved), `reset`
/// between realizations; zero steady-state allocation thereafter (the
/// trajectory budget `vertices + pairs` is exact — see the substrate).
#[derive(Debug, Default)]
pub struct LongRangeWorkspace {
    substrate: IncrementalUnionWorkspace,
    /// Occupied pair indices in consideration order (the substrate's
    /// activated-edge set, exposed as a list).
    occupied: Vec<usize>,
}

impl LongRangeWorkspace {
    /// An unprepared workspace.
    pub fn new() -> Self {
        Self::default()
    }

    /// Reserve every buffer for a box of exactly `vertices` vertices and
    /// `pairs` unordered pairs, resetting all state. Re-preparing with
    /// different sizes is supported and leaves no stale state.
    pub fn prepare(&mut self, vertices: usize, pairs: usize) -> Result<(), LongRangeError> {
        if self.occupied.capacity() < pairs {
            self.occupied.clear();
            self.occupied
                .try_reserve_exact(pairs)
                .map_err(LongRangeError::from)?;
        }
        self.occupied.clear();
        self.substrate
            .prepare(vertices, pairs)
            .map_err(|error| match error {
                IncrementalError::Capacity(source) => {
                    LongRangeError::Capacity(ProcessError::from(source))
                }
                other => LongRangeError::Substrate(other.to_string()),
            })?;
        Ok(())
    }

    /// Reset for the next realization, retaining all capacity.
    pub fn reset(&mut self) {
        self.substrate.reset();
        self.occupied.clear();
    }

    /// The driven G1 substrate (online component summaries along the
    /// added-edge count; trajectory steps are `[0; n] ++ [1..=m]`).
    #[inline]
    pub const fn substrate(&self) -> &IncrementalUnionWorkspace {
        &self.substrate
    }

    /// Occupied pair indices of the last run, in consideration order.
    #[inline]
    pub fn occupied_pairs(&self) -> &[usize] {
        &self.occupied
    }

    /// Size of the component containing dense box vertex `index`, or
    /// `None` outside `0..vertices` (every vertex is active after a run).
    /// The pair-universe view is module-private, so its `VertexId`s are
    /// minted here rather than by the caller.
    #[inline]
    pub fn component_size(&self, index: usize) -> Option<usize> {
        self.substrate
            .component_size(VertexId::from_valid_index(index))
    }
}

/// Run one complete long-range realization (the reusable, RNG-free core).
///
/// `uniforms` holds one `f64` per pair in the documented consideration
/// order (ascending pair index); every entry must be finite and in
/// `[0, 1]` — the pair at index k is occupied iff `uniforms[k] <
/// min(1, beta * r^{-s})` with r the ℓ1 distance of pair k's endpoints.
/// The workspace is re-prepared idempotently, so a correctly sized
/// workspace never allocates in steady state.
pub fn long_range_run(
    box_: &LongRangeBox,
    beta: f64,
    s: f64,
    uniforms: &[f64],
    workspace: &mut LongRangeWorkspace,
) -> Result<LongRangeOutcome, LongRangeError> {
    validate_parameters(beta, s)?;
    if uniforms.len() != box_.pairs() {
        return Err(LongRangeError::UniformLengthMismatch {
            expected: box_.pairs(),
            actual: uniforms.len(),
        });
    }
    for (index, &value) in uniforms.iter().enumerate() {
        if !value.is_finite() || value < 0.0 || value > 1.0 {
            return Err(LongRangeError::UniformInvalid { index, value });
        }
    }
    let vertices = box_.vertices();
    workspace.prepare(vertices, box_.pairs())?;
    let universe = PairUniverse { box_ };

    // Stage: all box vertices, dense ascending, step coordinate 0 (the G2
    // convention — the fixed vertex set is the stage).
    for index in 0..vertices {
        workspace
            .substrate
            .activate_vertex(&universe, VertexId::from_valid_index(index), 0)
            .map_err(substrate_error)?;
    }

    // One uniform per pair, ascending pair order; occupied pairs drive the
    // substrate with the 1-based added-edge count as the step coordinate.
    let mut added = 0_usize;
    let mut pair = 0_usize;
    for first in 0..vertices {
        for second in (first + 1)..vertices {
            let probability = LongRangeBox::pair_probability(box_.distance(first, second), beta, s);
            if uniforms[pair] < probability {
                added += 1;
                workspace
                    .substrate
                    .activate_edge(&universe, EdgeId::from_valid_index(pair), added)
                    .map_err(substrate_error)?;
                workspace.occupied.push(pair);
            }
            pair += 1;
        }
    }
    debug_assert_eq!(pair, box_.pairs());
    let summary = workspace.substrate.summary();
    Ok(LongRangeOutcome {
        added_edges: added,
        component_count: summary.component_count,
        largest_component_size: summary.largest_component_size,
        raw_second_moment: summary.raw_second_moment,
    })
}

/// Shared (beta, s) domain validation (core and adapter): both must be
/// finite and non-negative. Zero is legal but degenerate (no edges, resp.
/// a distance-independent probability) — documented, not rejected.
fn validate_parameters(beta: f64, s: f64) -> Result<(), LongRangeError> {
    if !beta.is_finite() || beta < 0.0 {
        return Err(LongRangeError::BetaInvalid { value: beta });
    }
    if !s.is_finite() || s < 0.0 {
        return Err(LongRangeError::ExponentInvalid { value: s });
    }
    Ok(())
}

fn substrate_error(error: IncrementalError) -> LongRangeError {
    LongRangeError::Substrate(error.to_string())
}

/// MonteCarlo adapter for the long-range power-law edge law on a box.
///
/// One `sweep()` draws the uniform field (the documented dense-pair-order
/// draw contract) and runs one complete realization on the G1 substrate;
/// `measure()` records it as ONE C0 element-preserving vector sample plus
/// scalars (see the module documentation for the padding convention).
/// Configure `thermalization_sweeps = 0`; realizations are i.i.d.
pub struct LongRangePercolationMC {
    box_: LongRangeBox,
    beta: f64,
    s: f64,
    steps: usize,
    workspace: LongRangeWorkspace,
    uniforms: Vec<f64>,
    history: Vec<f64>,
    outcome: LongRangeOutcome,
}

impl LongRangePercolationMC {
    /// Build the adapter. `steps` is the recorded-array length (`>= 1`);
    /// every domain is validated here so `sweep` cannot fail.
    pub fn new(
        box_: LongRangeBox,
        beta: f64,
        s: f64,
        steps: usize,
    ) -> Result<Self, LongRangeError> {
        validate_parameters(beta, s)?;
        if steps == 0 {
            return Err(LongRangeError::NoSteps);
        }
        let mut workspace = LongRangeWorkspace::new();
        workspace.prepare(box_.vertices(), box_.pairs())?;
        let mut uniforms = Vec::new();
        uniforms
            .try_reserve_exact(box_.pairs())
            .map_err(LongRangeError::from)?;
        uniforms.resize(box_.pairs(), 0.0);
        let mut history = Vec::new();
        history
            .try_reserve_exact(steps)
            .map_err(LongRangeError::from)?;
        Ok(Self {
            box_,
            beta,
            s,
            steps,
            workspace,
            uniforms,
            history,
            outcome: LongRangeOutcome::default(),
        })
    }

    /// The box the law runs on.
    pub const fn geometry(&self) -> &LongRangeBox {
        &self.box_
    }

    /// The intensity `beta`.
    pub const fn beta(&self) -> f64 {
        self.beta
    }

    /// The decay exponent `s`.
    pub const fn exponent(&self) -> f64 {
        self.s
    }

    /// Recorded-array length (the C0 shape of `LongRangeLargestSize`).
    pub const fn steps(&self) -> usize {
        self.steps
    }

    /// Outcome of the last realized box.
    pub const fn last_outcome(&self) -> LongRangeOutcome {
        self.outcome
    }

    /// The fixed-length recorded array of the last realization.
    pub fn last_history(&self) -> &[f64] {
        &self.history
    }

    /// The workspace of the last realization (occupied pairs and substrate
    /// trajectory).
    pub const fn workspace(&self) -> &LongRangeWorkspace {
        &self.workspace
    }

    /// Draw the uniform field and run one complete realization.
    ///
    /// Draw contract (bitwise reproducibility): one `f64` per pair in
    /// ascending pair order (the consideration order of the module
    /// documentation). Deterministic given the seed law and RNG stream;
    /// exposed for tests and RNG-free analysis paths.
    pub fn realize<R: Rng + ?Sized>(&mut self, rng: &mut R) {
        for slot in self.uniforms.iter_mut() {
            *slot = rng.random::<f64>();
        }
        self.outcome = {
            let Self {
                box_,
                beta,
                s,
                workspace,
                uniforms,
                ..
            } = self;
            long_range_run(box_, *beta, *s, uniforms, workspace)
                .expect("constructor validated every domain")
        };
        // Compress into the fixed C0 shape (G3 padding convention): entry
        // i is the largest size after the (i+1)-th added edge; entries past
        // the realized length repeat the final value.
        let summaries = self.workspace.substrate().trajectory().summaries();
        let base = self.box_.vertices();
        let final_largest = self.outcome.largest_component_size;
        self.history.clear();
        self.history.extend((0..self.steps).map(|step| {
            if step < self.outcome.added_edges {
                summaries[base + step].largest_component_size as f64
            } else {
                final_largest as f64
            }
        }));
    }
}

impl MonteCarlo for LongRangePercolationMC {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        self.realize(&mut ctx.rng);
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        ctx.measure_array("LongRangeLargestSize", &self.history)
            .expect("array shape is fixed at construction");
        ctx.measure("LongRangeAddedEdges", self.outcome.added_edges as f64);
        ctx.measure(
            "LongRangeComponentCount",
            self.outcome.component_count as f64,
        );
        ctx.measure(
            "LongRangeSecondMoment",
            self.outcome.raw_second_moment as f64,
        );
    }

    fn name(&self) -> &'static str {
        "LongRangePercolation"
    }
}

/// Typed rejections for the long-range edge law.
#[derive(Debug)]
pub enum LongRangeError {
    /// A box description with no dimensions.
    EmptyBox,
    /// A zero side length.
    ZeroSide { axis: usize },
    /// The vertex count does not fit a `usize`.
    BoxTooLarge { sides: Vec<usize> },
    /// The pair count does not fit a `usize`.
    PairUniverseTooLarge { vertices: usize },
    /// A negative or non-finite intensity.
    BetaInvalid { value: f64 },
    /// A negative or non-finite decay exponent.
    ExponentInvalid { value: f64 },
    /// A zero-length recorded array.
    NoSteps,
    /// A uniform field with the wrong length.
    UniformLengthMismatch { expected: usize, actual: usize },
    /// A uniform outside `[0, 1]` or not finite.
    UniformInvalid { index: usize, value: f64 },
    /// A buffer reservation failure.
    Capacity(ProcessError),
    /// The driven G1 substrate refused an operation (a process bug, not a
    /// caller error, once the domains validated).
    Substrate(String),
}

impl From<TryReserveError> for LongRangeError {
    fn from(source: TryReserveError) -> Self {
        Self::Capacity(ProcessError::from(source))
    }
}

impl core::fmt::Display for LongRangeError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::EmptyBox => {
                write!(formatter, "a long-range box needs at least one dimension")
            }
            Self::ZeroSide { axis } => {
                write!(formatter, "side length of axis {axis} is zero")
            }
            Self::BoxTooLarge { sides } => {
                write!(
                    formatter,
                    "box {sides:?} has more vertices than a usize can count"
                )
            }
            Self::PairUniverseTooLarge { vertices } => write!(
                formatter,
                "box with {vertices} vertices has more pairs than a usize can count"
            ),
            Self::BetaInvalid { value } => write!(
                formatter,
                "beta {value} is negative or not finite (zero is legal but degenerate)"
            ),
            Self::ExponentInvalid { value } => write!(
                formatter,
                "exponent s {value} is negative or not finite (zero is legal but degenerate)"
            ),
            Self::NoSteps => write!(formatter, "at least one recorded step is required"),
            Self::UniformLengthMismatch { expected, actual } => write!(
                formatter,
                "uniform field has {actual} entries, expected one per pair ({expected})"
            ),
            Self::UniformInvalid { index, value } => write!(
                formatter,
                "uniform {value} at pair {index} is outside [0, 1] or not finite"
            ),
            Self::Capacity(source) => {
                write!(formatter, "could not prepare long-range buffers: {source}")
            }
            Self::Substrate(reason) => {
                write!(
                    formatter,
                    "driven G1 substrate refused an operation: {reason}"
                )
            }
        }
    }
}

impl std::error::Error for LongRangeError {}
