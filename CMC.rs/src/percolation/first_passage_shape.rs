//! F3: first-passage shape theorem brackets, geodesic counting and a
//! finite-horizon geodesic-aging statistic on open Z^2 boxes (roadmap v2
//! §3 F3), over the F1 weighted shortest-path layer.
//!
//! # Anchor
//!
//! Anchor (archived
//! `docs/ref/md/Auffinger_2017_50yrs_first_passage_percolation.md`):
//! Auffinger-Damron-Hanson, *50 years of first-passage percolation*, AMS
//! ULS 68 (2017). The exact statements this layer brackets:
//!
//! - **Time constant** (Theorem 2.1, ~lines 392-405): under
//!   `E min[t_1, ..., t_{2d}] < infinity` there is a constant
//!   `mu(e_1) in [0, infinity)` with `T(0, n e_1) / n -> mu(e_1)` a.s.
//!   and in `L1`. The proof is subadditivity: a path from `0` to `n e_1`
//!   need not go through `m e_1`, so `T(0, n e_1) <= T(0, m e_1) +
//!   T(m e_1, n e_1)` (Figure 2, ~lines 442-464) and
//!   `lim E[T(0, n e_1)]/n = inf_n E[T(0, n e_1)]/n` (~lines 424-433) —
//!   convergence from ABOVE.
//! - **Bounds** (~lines 518-526): `0 <= mu(e_1) <= E tau_e` (the direct
//!   path plus the law of large numbers), with Theorem 2.4
//!   (Hammersley-Welsh): `mu(e_1) < E tau_e` whenever `F` is not a
//!   trivial distribution.
//! - **Shape theorem** (Theorem 2.16, Cox-Durrett, ~lines 745-755): for
//!   each weight law there is a deterministic, convex, compact `B_nu`
//!   with non-empty interior, symmetric about the axes, such that for
//!   every `epsilon > 0`, `(1 - epsilon) B_nu subset B(t)/t subset
//!   (1 + epsilon) B_nu` eventually almost surely. The finite boxes here
//!   measure directional time constants `mu(u)`; no closed-form shape is
//!   claimed.
//! - **Geodesics** (Definition 4.1, ~lines 3904-3914): a path `gamma`
//!   from `x` to `y` with `T(gamma) = T(x, y)` is a geodesic. With
//!   continuous passage times the metric is uniquely geodesic and the
//!   geodesic graph `T(x)` is a spanning tree (~lines 4104-4117). At
//!   two-level weights geodesics are non-unique; the distance-based
//!   predecessor DAG (all edges with `dist[u] + w = dist[v]`) enumerates
//!   them exactly.
//! - **Overlap** (Question 16, ~lines 2280-2283): "Show that in `d = 2`
//!   there exists `a > 0` such that the expected overlap of two geodesics
//!   from `0` to `nx` sampled independently is at least `n^a`." The
//!   survey gives no crisp finite-horizon aging law and does not use the
//!   word "aging"; the statistic below is EXPLORATORY and is gated for
//!   correctness only (see [`geodesic_aging`]).
//!
//! # Directional time constants
//!
//! [`ShapeGeometry`] builds one open `side x side` Z^2 box
//! ([`build_square`] row-major dense ids `x + y * side`) with a single
//! origin at the center and one target per `(Direction, distance)` pair:
//! `(1, 0)`, `(0, 1)`, `(1, 1)` and `(1, -1)` times `n`. All targets live
//! in the same box, so ONE F1 Dijkstra per realization yields the whole
//! directional ladder. The box reserves `margin` free cells beyond the
//! farthest target on every side, which bounds the boundary bias of the
//! finite-box measurement.
//!
//! [`ShapeSampler`] draws one weight per box edge under [`EdgeWeightLaw`]
//! (the F2 bitwise contract) and records the C0 passage-time profile —
//! one `f64` per target. [`ShapeEstimator`] accumulates per-target mean,
//! unbiased variance and standard error. The gates (validation layer)
//! require, for each direction, `mu_n = E[T(nu)]/n` to be positive,
//! bounded above by `E tau_e`, non-increasing within noise along
//! `n = 8, 16, 32`, and consistent with the norm inequalities
//! `mu(u) + mu(v) >= mu(u + v)` (paired per realization).
//!
//! # Geodesic counting
//!
//! [`geodesic_count_run`] runs one F1 pass and then a dynamic program over
//! the tight-edge DAG in increasing distance order:
//! `ways[source] = 1`, and for every tight edge `u -> v` (with
//! `dist[u] + w == dist[v]`), `ways[v] += ways[u]`. Counts are exact
//! `u128` with explicit saturation, parallel edges counted separately,
//! and zero weights rejected: a zero-weight edge makes the tight graph
//! cyclic and the geodesic count infinite (the F1 layer allows zero
//! weights; this counting layer does not).
//!
//! # Exploratory aging statistic
//!
//! [`geodesic_aging`] measures the fraction of the VERTEX sequence of the
//! geodesic to `n u` that survives as the prefix of the geodesic to `2 n u`
//! in the SAME environment and under F1's deterministic predecessor tie
//! policy. (The survey's overlap object is a set of edges; this vertex-
//! prefix statistic is a related but distinct exploratory proxy.) The
//! survey's finite-horizon concern is Question 30 (~lines 5044-5046):
//! "Consider the collection of all edges that belong to an infinite
//! number of `Gamma(v)` ... Does `f(r) -> infinity` as `r -> infinity`;
//! and, if so, how fast?" The statistic is a one-environment,
//! one-tie-broken-path analog of that persistence question, NOT the
//! survey's overlap object; no asymptotic claim is made. Correctness
//! (every reconstructed path is a true geodesic and simple) is gated, and
//! the ladder's mean persistence is reported with CIs as a descriptive
//! exploratory statistic only (fresh-seed ladders are not monotone; no
//! convergence claim is made).
//!
//! # Cost
//!
//! One realization is one F1 Dijkstra, `O((V + E) log E)` over the box
//! plus an `O(targets)` profile copy; the sampler reserves every buffer
//! once, so steady-state realizations allocate nothing. Geodesic counting
//! adds an `O(V log V)` distance sort and an `O(V + E)` DAG pass; the
//! workspace reserves the sort buffer and the way counts once.

use core::fmt;

use rand::{Rng, RngExt};

use super::first_passage::{
    first_passage_run, FirstPassageError, FirstPassageWorkspace, GeodesicTracking, NO_PREDECESSOR,
};
use super::first_passage_observables::{EdgeWeightLaw, PassageObservablesError};
use super::process::ProcessError;
use crate::{build_square, CsrLattice, UndirectedGraphView, VertexId};

/// Hard guard-rail on the open box side built by [`ShapeGeometry`]: the
/// all-target geometry is `side^2` vertices, so a side above `1024`
/// (`2^20` vertices) is a typed rejection rather than an allocation
/// attempt.
pub const SHAPE_SIDE_CAP: usize = 1024;

/// One of the four Z^2 directions the stage brackets.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Direction {
    /// `(1, 0)`.
    East,
    /// `(0, 1)`.
    North,
    /// `(1, 1)`.
    Northeast,
    /// `(1, -1)`.
    Southeast,
}

impl Direction {
    /// The unit lattice step `(dx, dy)`.
    pub const fn step(self) -> (i64, i64) {
        match self {
            Self::East => (1, 0),
            Self::North => (0, 1),
            Self::Northeast => (1, 1),
            Self::Southeast => (1, -1),
        }
    }

    /// A stable short label for reports.
    pub const fn label(self) -> &'static str {
        match self {
            Self::East => "east",
            Self::North => "north",
            Self::Northeast => "northeast",
            Self::Southeast => "southeast",
        }
    }
}

/// One target of a [`ShapeGeometry`]: a direction, a step count, and the
/// dense vertex id of `origin + distance * direction`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ShapeTarget {
    /// The lattice direction.
    pub direction: Direction,
    /// The number of unit steps from the origin.
    pub distance: usize,
    /// The dense vertex id of the target.
    pub vertex: usize,
}

/// Open Z^2 box with one centered origin and one target per requested
/// `(direction, distance)` pair.
///
/// The box has `side = 2 * (max_distance + margin) + 1` vertices per axis
/// and dense ids `x + y * side` (the [`build_square`] convention). Target
/// order is directions-major, distances-minor in the caller's order (the
/// profile coordinate order). The geometry owns its
/// [`CsrLattice`](crate::CsrLattice); a [`ShapeSampler`] takes it by
/// value (clone first if several samplers share one description).
#[derive(Clone, Debug)]
pub struct ShapeGeometry {
    lattice: CsrLattice,
    origin: usize,
    side: usize,
    center: usize,
    targets: Vec<ShapeTarget>,
}

impl ShapeGeometry {
    /// Build the geometry for every `(direction, distance)` combination.
    ///
    /// `directions` and `distances` must be non-empty, every distance must
    /// be at least 1, and the derived side must stay at or below
    /// [`SHAPE_SIDE_CAP`]. Rejections are typed and leave no partial state
    /// (the lattice is built only after all checks).
    pub fn new(
        directions: &[Direction],
        distances: &[usize],
        margin: usize,
    ) -> Result<Self, ShapeError> {
        if directions.is_empty() {
            return Err(ShapeError::NoDirections);
        }
        if distances.is_empty() {
            return Err(ShapeError::NoDistances);
        }
        for (index, &distance) in distances.iter().enumerate() {
            if distance == 0 {
                return Err(ShapeError::ZeroDistance { index });
            }
        }
        let max_distance = *distances.iter().max().expect("non-empty distances");
        let reach = max_distance
            .checked_add(margin)
            .ok_or(ShapeError::SideOverflow)?;
        let side = reach
            .checked_mul(2)
            .and_then(|value| value.checked_add(1))
            .ok_or(ShapeError::SideOverflow)?;
        if side > SHAPE_SIDE_CAP {
            return Err(ShapeError::SideAboveCap {
                side,
                cap: SHAPE_SIDE_CAP,
            });
        }
        let lattice = build_square(side, side, false);
        let center = max_distance + margin;
        let origin = center + center * side;
        let mut targets = Vec::new();
        targets
            .try_reserve_exact(directions.len() * distances.len())
            .map_err(|source| ShapeError::Capacity {
                buffer: "shape targets",
                source: ProcessError::from(source),
            })?;
        for &direction in directions {
            let (dx, dy) = direction.step();
            for &distance in distances {
                let x = center as i64 + dx * distance as i64;
                let y = center as i64 + dy * distance as i64;
                debug_assert!(x >= 0 && y >= 0);
                targets.push(ShapeTarget {
                    direction,
                    distance,
                    vertex: x as usize + (y as usize) * side,
                });
            }
        }
        Ok(Self {
            lattice,
            origin,
            side,
            center,
            targets,
        })
    }

    /// The box topology (all dense edges carry drawn weights).
    #[inline]
    pub fn lattice(&self) -> &CsrLattice {
        &self.lattice
    }

    /// Dense vertex id of the origin.
    #[inline]
    pub const fn origin(&self) -> usize {
        self.origin
    }

    /// The box side length.
    #[inline]
    pub const fn side(&self) -> usize {
        self.side
    }

    /// The origin's coordinate on both axes.
    #[inline]
    pub const fn center(&self) -> usize {
        self.center
    }

    /// The targets in profile order.
    #[inline]
    pub fn targets(&self) -> &[ShapeTarget] {
        &self.targets
    }

    /// Profile length (the target count).
    #[inline]
    pub fn boundary_size(&self) -> usize {
        self.targets.len()
    }

    /// Number of box edges (all drawn; there are no structural edges).
    #[inline]
    pub fn edge_count(&self) -> usize {
        self.lattice.n_edges()
    }
}

/// Reusable directional realization sampler: the weight slice, the F1
/// workspace, and the passage-time profile of the last realization.
///
/// `new` reserves every buffer once; [`Self::realize`] then draws one
/// uniform per box edge (dense [`build_square`] order, the F2 contract),
/// runs one F1 pass from the origin, and copies the target distances into
/// the profile. One realization is ONE C0 vector sample of fixed shape
/// (one entry per target). [`Self::realize_tracked`] additionally records
/// the F1 predecessor tree for geodesic reconstruction; `realize` leaves
/// predecessors off.
#[derive(Debug)]
pub struct ShapeSampler {
    geometry: ShapeGeometry,
    weights: Vec<f64>,
    profile: Vec<f64>,
    workspace: FirstPassageWorkspace,
}

impl ShapeSampler {
    /// Reserve all buffers for `geometry` (taken by value) and pre-prepare
    /// the F1 workspace.
    pub fn new(geometry: ShapeGeometry) -> Result<Self, ShapeError> {
        let edges = geometry.edge_count();
        let mut weights = Vec::new();
        weights
            .try_reserve_exact(edges)
            .map_err(|source| ShapeError::Capacity {
                buffer: "shape weights",
                source: ProcessError::from(source),
            })?;
        weights.resize(edges, 0.0);
        let mut profile = Vec::new();
        profile
            .try_reserve_exact(geometry.boundary_size())
            .map_err(|source| ShapeError::Capacity {
                buffer: "shape profile",
                source: ProcessError::from(source),
            })?;
        profile.resize(geometry.boundary_size(), f64::INFINITY);
        let mut workspace = FirstPassageWorkspace::new();
        workspace
            .prepare(geometry.lattice().n_sites, edges)
            .map_err(ShapeError::Passage)?;
        Ok(Self {
            geometry,
            weights,
            profile,
            workspace,
        })
    }

    /// One realization: draw one weight per box edge under `law`, run one
    /// F1 pass from the origin (no predecessor tracking), and return the
    /// target profile.
    pub fn realize<R: Rng + ?Sized>(
        &mut self,
        law: EdgeWeightLaw,
        rng: &mut R,
    ) -> Result<&[f64], ShapeError> {
        self.draw(law, rng, GeodesicTracking::Off)
    }

    /// [`Self::realize`] with predecessor tracking (for geodesic
    /// reconstruction through [`Self::workspace`]).
    pub fn realize_tracked<R: Rng + ?Sized>(
        &mut self,
        law: EdgeWeightLaw,
        rng: &mut R,
    ) -> Result<&[f64], ShapeError> {
        self.draw(law, rng, GeodesicTracking::Predecessors)
    }

    /// The shared draw + run path of both realize variants.
    fn draw<R: Rng + ?Sized>(
        &mut self,
        law: EdgeWeightLaw,
        rng: &mut R,
        tracking: GeodesicTracking,
    ) -> Result<&[f64], ShapeError> {
        law.validate().map_err(ShapeError::Law)?;
        for weight in self.weights.iter_mut() {
            *weight = law.draw(rng.random::<f64>());
        }
        self.run(tracking)
    }

    /// One realization on a caller-supplied weight slice (one `f64` per
    /// dense edge; the exact-enumeration path). A length mismatch is a
    /// typed rejection that leaves the last profile intact; negative and
    /// non-finite entries are the wrapped F1 rejections.
    pub fn realize_with_weights(
        &mut self,
        weights: &[f64],
        tracking: GeodesicTracking,
    ) -> Result<&[f64], ShapeError> {
        let expected = self.geometry.edge_count();
        if weights.len() != expected {
            return Err(ShapeError::WeightLength {
                expected,
                actual: weights.len(),
            });
        }
        self.weights.copy_from_slice(weights);
        self.run(tracking)
    }

    /// The single F1 pass and the profile copy.
    fn run(&mut self, tracking: GeodesicTracking) -> Result<&[f64], ShapeError> {
        first_passage_run(
            &self.geometry.lattice,
            &self.weights,
            self.geometry.origin,
            tracking,
            &mut self.workspace,
        )
        .map_err(ShapeError::Passage)?;
        let distances = self.workspace.distances();
        for (slot, target) in self.profile.iter_mut().zip(&self.geometry.targets) {
            *slot = distances[target.vertex];
        }
        Ok(&self.profile)
    }

    /// The geometry this sampler was built from.
    #[inline]
    pub fn geometry(&self) -> &ShapeGeometry {
        &self.geometry
    }

    /// The passage-time profile of the last realization (one `f64` per
    /// target, all [`f64::INFINITY`] before the first).
    #[inline]
    pub fn profile(&self) -> &[f64] {
        &self.profile
    }

    /// The weight slice of the last realization (one `f64` per dense
    /// edge) — the replay handle for the bitwise draw contract.
    #[inline]
    pub fn weights(&self) -> &[f64] {
        &self.weights
    }

    /// The F1 workspace of the last realization (distances and optional
    /// predecessor tree).
    #[inline]
    pub fn workspace(&self) -> &FirstPassageWorkspace {
        &self.workspace
    }
}

/// Accumulates passage-time profiles into per-target sample statistics.
///
/// All buffers are reserved at `new`, so [`Self::add`] allocates nothing;
/// the estimate's vectors are allocated at report time only.
#[derive(Clone, Debug)]
pub struct ShapeEstimator {
    targets: usize,
    samples: usize,
    sum: Vec<f64>,
    sum_sq: Vec<f64>,
}

impl ShapeEstimator {
    /// Reserve accumulators for profiles of `targets` entries (`>= 1`).
    pub fn new(targets: usize) -> Result<Self, ShapeError> {
        if targets == 0 {
            return Err(ShapeError::NoTargets);
        }
        let mut sum = Vec::new();
        sum.try_reserve_exact(targets)
            .map_err(|source| ShapeError::Capacity {
                buffer: "shape sums",
                source: ProcessError::from(source),
            })?;
        sum.resize(targets, 0.0);
        let mut sum_sq = Vec::new();
        sum_sq
            .try_reserve_exact(targets)
            .map_err(|source| ShapeError::Capacity {
                buffer: "shape sums of squares",
                source: ProcessError::from(source),
            })?;
        sum_sq.resize(targets, 0.0);
        Ok(Self {
            targets,
            samples: 0,
            sum,
            sum_sq,
        })
    }

    /// The reserved profile length.
    #[inline]
    pub fn targets(&self) -> usize {
        self.targets
    }

    /// Profiles accumulated so far.
    #[inline]
    pub const fn samples(&self) -> usize {
        self.samples
    }

    /// Add one realization's profile. Length mismatches are typed
    /// rejections; nothing is added on rejection.
    pub fn add(&mut self, profile: &[f64]) -> Result<(), ShapeError> {
        if profile.len() != self.targets {
            return Err(ShapeError::ProfileLength {
                expected: self.targets,
                actual: profile.len(),
            });
        }
        self.samples += 1;
        for ((sum, sum_sq), &value) in self.sum.iter_mut().zip(self.sum_sq.iter_mut()).zip(profile)
        {
            *sum += value;
            *sum_sq += value * value;
        }
        Ok(())
    }

    /// Sample statistics, or `None` before two realizations. Allocates the
    /// estimate's vectors (report-time, never per realization).
    pub fn estimate(&self) -> Option<ShapeEstimate> {
        if self.samples < 2 {
            return None;
        }
        let samples = self.samples as f64;
        let means: Vec<f64> = self.sum.iter().map(|&sum| sum / samples).collect();
        let variances: Vec<f64> = self
            .sum_sq
            .iter()
            .zip(self.sum.iter())
            .map(|(&sum_sq, &sum)| ((sum_sq - sum * sum / samples) / (samples - 1.0)).max(0.0))
            .collect();
        let standard_errors = variances
            .iter()
            .map(|&variance| (variance / samples).sqrt())
            .collect();
        Some(ShapeEstimate {
            samples: self.samples,
            means,
            variances,
            standard_errors,
        })
    }
}

/// Sample statistics over independent directional realizations.
#[derive(Clone, Debug, PartialEq)]
pub struct ShapeEstimate {
    /// Profiles accumulated.
    pub samples: usize,
    /// Per-target sample mean (target order).
    pub means: Vec<f64>,
    /// Per-target unbiased sample variance.
    pub variances: Vec<f64>,
    /// Per-target standard error of [`Self::means`].
    pub standard_errors: Vec<f64>,
}

/// Scalar outcome of one geodesic-count computation.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct GeodesicCountOutcome {
    /// The dense target vertex.
    pub target: usize,
    /// Exact geodesic count from the source to the target
    /// (`u128::MAX` when [`Self::saturated`]).
    pub geodesic_count: u128,
    /// Whether the exact count overflowed `u128` (saturating convention).
    pub saturated: bool,
    /// Vertices reached by the F1 pass (source included).
    pub reached: usize,
}

/// Reusable geodesic-count workspace: the F1 distances/predecessor
/// workspace, the per-vertex way counts, and the distance-sorted vertex
/// order.
///
/// `prepare` is implicit in [`geodesic_count_run`] (which re-prepares the
/// F1 workspace idempotently); buffers grow on first use and are reused,
/// so repeated calls on similarly sized graphs allocate nothing.
#[derive(Debug, Default)]
pub struct GeodesicCountWorkspace {
    passage: FirstPassageWorkspace,
    ways: Vec<u128>,
    order: Vec<usize>,
}

impl GeodesicCountWorkspace {
    /// An unprepared workspace.
    pub fn new() -> Self {
        Self::default()
    }

    /// Exact geodesic counts per vertex of the last run (`u128::MAX` at
    /// saturated vertices); empty before the first run.
    #[inline]
    pub fn ways(&self) -> &[u128] {
        &self.ways
    }

    /// Distances of the last run (the F1 workspace accessor).
    #[inline]
    pub fn distances(&self) -> &[f64] {
        self.passage.distances()
    }
}

/// Count geodesics from `source` to `target` on an undirected multigraph
/// with a caller-supplied non-negative weight field.
///
/// One F1 pass gives every distance; then, in increasing distance order,
/// every TIGHT edge `u -> v` (`dist[u] + w == dist[v]`, exact f64
/// equality) adds `ways[u]` to `ways[v]`. Positive weights are required:
/// a zero-weight edge creates a zero-cost cycle in the tight graph and
/// infinitely many geodesic walks, so `w == 0` is the typed
/// [`ShapeError::ZeroWeight`] rejection (F1 itself allows zero weights).
/// Counts are exact `u128` and saturate with an explicit flag; parallel
/// edges are distinct geodesics (edge sequences), self-loops never tight.
/// An unreachable target yields count 0 (not an error); `source == target`
/// yields 1 (the empty path).
///
/// Validation happens before the F1 pass is started, so typed rejections
/// leave the workspace's previous counts untouched.
pub fn geodesic_count_run<G: UndirectedGraphView>(
    graph: &G,
    weights: &[f64],
    source: usize,
    target: usize,
    workspace: &mut GeodesicCountWorkspace,
) -> Result<GeodesicCountOutcome, ShapeError> {
    let vertices = graph.vertex_count();
    let edges = graph.edge_count();
    if weights.len() != edges {
        return Err(ShapeError::WeightLength {
            expected: edges,
            actual: weights.len(),
        });
    }
    for (edge, &value) in weights.iter().enumerate() {
        if value == 0.0 {
            return Err(ShapeError::ZeroWeight { edge });
        }
    }
    if source >= vertices {
        return Err(ShapeError::SourceDomain { source, vertices });
    }
    if target >= vertices {
        return Err(ShapeError::TargetDomain { target, vertices });
    }
    first_passage_run(
        graph,
        weights,
        source,
        GeodesicTracking::Off,
        &mut workspace.passage,
    )
    .map_err(ShapeError::Passage)?;

    reserve_exact(&mut workspace.ways, vertices, 0_u128)?;
    reserve_exact(&mut workspace.order, vertices, 0_usize)?;
    let dist = workspace.passage.distances();
    for (index, slot) in workspace.order.iter_mut().enumerate() {
        *slot = index;
    }
    workspace.order.sort_unstable_by(|&left, &right| {
        dist[left]
            .total_cmp(&dist[right])
            .then_with(|| left.cmp(&right))
    });
    workspace.ways[source] = 1;
    let reached = dist.iter().filter(|distance| distance.is_finite()).count();
    for index in 0..workspace.order.len() {
        let u = workspace.order[index];
        let ways_u = workspace.ways[u];
        if ways_u == 0 || !dist[u].is_finite() {
            continue;
        }
        for incidence in graph.incidences(VertexId::from_valid_index(u)) {
            let v = incidence.neighbor.index();
            let edge = incidence.edge.index();
            if dist[v].is_finite() && dist[u] + weights[edge] == dist[v] && dist[v] > dist[u] {
                match workspace.ways[v].checked_add(ways_u) {
                    Some(value) => workspace.ways[v] = value,
                    None => {
                        workspace.ways[v] = u128::MAX;
                    }
                }
            }
        }
    }
    Ok(GeodesicCountOutcome {
        target,
        geodesic_count: workspace.ways[target],
        saturated: workspace.ways[target] == u128::MAX,
        reached,
    })
}

/// Scalar outcome of one geodesic-path reconstruction.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct GeodesicPathOutcome {
    /// The dense target vertex.
    pub target: usize,
    /// Vertices on the reconstructed path (the source and target included).
    pub vertices: usize,
    /// Edges on the path (`vertices - 1`).
    pub edges: usize,
}

/// Reconstruct the F1 predecessor-tree geodesic from `source` to `target`
/// into caller-owned buffers (cleared first; no allocation when their
/// capacity suffices).
///
/// The predecessor arrays must come from a [`GeodesicTracking::Predecessors`]
/// run on the same topology. The returned path is deterministic under the
/// F1 tie policy. With positive weights the path is simple (distances
/// strictly decrease along the predecessor chain); this function does not
/// re-verify simplicity or weights. Failures: equal-length arrays
/// expected, source/target in range, target reached (a
/// [`NO_PREDECESSOR`] chain end before the source), predecessor indices in
/// range, and a source not reachable from the chain within `vertices`
/// steps (cycle guard).
pub fn geodesic_path(
    predecessor_vertices: &[usize],
    predecessor_edges: &[usize],
    source: usize,
    target: usize,
    path: &mut Vec<usize>,
    path_edges: &mut Vec<usize>,
) -> Result<GeodesicPathOutcome, ShapeError> {
    let vertices = predecessor_vertices.len();
    if predecessor_edges.len() != vertices {
        return Err(ShapeError::PredecessorLength {
            expected: vertices,
            actual: predecessor_edges.len(),
        });
    }
    if source >= vertices {
        return Err(ShapeError::SourceDomain { source, vertices });
    }
    if target >= vertices {
        return Err(ShapeError::TargetDomain { target, vertices });
    }
    path.clear();
    path_edges.clear();
    path.push(target);
    let mut current = target;
    let mut steps = 0_usize;
    while current != source {
        if steps > vertices {
            return Err(ShapeError::PredecessorCycle { target });
        }
        let predecessor = predecessor_vertices[current];
        if predecessor == NO_PREDECESSOR {
            return Err(ShapeError::UnreachableTarget { target });
        }
        if predecessor >= vertices {
            return Err(ShapeError::PredecessorDomain {
                vertex: current,
                predecessor,
                vertices,
            });
        }
        let edge = predecessor_edges[current];
        if edge == NO_PREDECESSOR {
            return Err(ShapeError::PredecessorEdgeMissing { vertex: current });
        }
        path_edges.push(edge);
        current = predecessor;
        path.push(current);
        steps += 1;
    }
    path.reverse();
    path_edges.reverse();
    Ok(GeodesicPathOutcome {
        target,
        vertices: path.len(),
        edges: path_edges.len(),
    })
}

/// Finite-horizon geodesic-aging statistic: the shared leading PREFIX of
/// two paths (short = geodesic to `n u`, long = geodesic to `2 n u` in
/// the same environment).
///
/// Anchored (loosely, and labelled EXPLORATORY) in the survey's Question
/// 30 (~lines 5044-5046), which asks whether the edges belonging to
/// infinitely many geodesics `Gamma(v)` have a diverging count on circles
/// `|x| = r`; this statistic is the finite-`n`, single-tie-broken-path
/// analog. It is NOT the survey's independent-geodesic overlap of
/// Question 16, and no asymptotic law is claimed.
///
/// Returns `None` when either path is empty, the two paths start at
/// different vertices, or the short path has no edge. Otherwise
/// `persistence = shared_prefix_edges / short_edges in [0, 1]`, where
/// `shared_prefix_edges` counts the number of leading edges on which the
/// two paths agree (so `1` means the short geodesic is a prefix of the
/// long one).
pub fn geodesic_aging(short: &[usize], long: &[usize]) -> Option<GeodesicAging> {
    if short.len() < 2 || long.len() < 2 || short[0] != long[0] {
        return None;
    }
    let limit = short.len().min(long.len());
    let mut shared_vertices = 1_usize;
    while shared_vertices < limit && short[shared_vertices] == long[shared_vertices] {
        shared_vertices += 1;
    }
    let short_edges = short.len() - 1;
    let shared_prefix_edges = shared_vertices - 1;
    Some(GeodesicAging {
        short_edges,
        long_edges: long.len() - 1,
        shared_prefix_edges,
        persistence: shared_prefix_edges as f64 / short_edges as f64,
    })
}

/// Result of [`geodesic_aging`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GeodesicAging {
    /// Edges on the short (`n u`) path.
    pub short_edges: usize,
    /// Edges on the long (`2 n u`) path.
    pub long_edges: usize,
    /// Leading edges shared by the two paths.
    pub shared_prefix_edges: usize,
    /// `shared_prefix_edges / short_edges`, in `[0, 1]`.
    pub persistence: f64,
}

/// Reserve `size` logical entries (rewriting the fill value so a shrink
/// leaves no stale state, capacity retained).
fn reserve_exact<T: Clone>(storage: &mut Vec<T>, size: usize, fill: T) -> Result<(), ShapeError> {
    if storage.capacity() < size {
        storage.clear();
        storage
            .try_reserve_exact(size)
            .map_err(|source| ShapeError::Capacity {
                buffer: "geodesic-count buffers",
                source: ProcessError::from(source),
            })?;
    }
    storage.clear();
    storage.resize(size, fill);
    Ok(())
}

/// Typed rejections of the shape/geodesic layer.
#[derive(Debug)]
pub enum ShapeError {
    /// A geometry without directions.
    NoDirections,
    /// A geometry without distances.
    NoDistances,
    /// A zero step distance (the origin would be its own target).
    ZeroDistance { index: usize },
    /// The `2 * (max_distance + margin) + 1` side overflows `usize`.
    SideOverflow,
    /// The box side exceeds [`SHAPE_SIDE_CAP`].
    SideAboveCap { side: usize, cap: usize },
    /// An invalid weight law (the wrapped F2 rejection).
    Law(PassageObservablesError),
    /// A weight slice without exactly one entry per dense edge.
    WeightLength { expected: usize, actual: usize },
    /// A zero weight in the geodesic-count path (ill-defined count).
    ZeroWeight { edge: usize },
    /// A profile whose length differs from the estimator's target count.
    ProfileLength { expected: usize, actual: usize },
    /// A zero-target estimator.
    NoTargets,
    /// A source outside the graph's dense index range.
    SourceDomain { source: usize, vertices: usize },
    /// A target outside the graph's dense index range.
    TargetDomain { target: usize, vertices: usize },
    /// A target not reached by the predecessor run.
    UnreachableTarget { target: usize },
    /// Predecessor vertex/edge arrays of different lengths.
    PredecessorLength { expected: usize, actual: usize },
    /// A predecessor vertex outside the array.
    PredecessorDomain {
        vertex: usize,
        predecessor: usize,
        vertices: usize,
    },
    /// A reached vertex without a predecessor edge.
    PredecessorEdgeMissing { vertex: usize },
    /// A predecessor chain that did not reach the source within the
    /// vertex count (a cyclic or corrupt tree).
    PredecessorCycle { target: usize },
    /// A rejection raised by the underlying F1 pass.
    Passage(FirstPassageError),
    /// A buffer reservation failure.
    Capacity {
        buffer: &'static str,
        source: ProcessError,
    },
}

impl fmt::Display for ShapeError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoDirections => write!(formatter, "the shape geometry needs a direction"),
            Self::NoDistances => write!(formatter, "the shape geometry needs a distance"),
            Self::ZeroDistance { index } => write!(
                formatter,
                "distance 0 at position {index} would make the origin its own target"
            ),
            Self::SideOverflow => {
                write!(formatter, "the derived box side overflows usize")
            }
            Self::SideAboveCap { side, cap } => write!(
                formatter,
                "shape box side {side} exceeds the guard-rail {cap}"
            ),
            Self::Law(source) => write!(formatter, "weight law rejected: {source}"),
            Self::WeightLength { expected, actual } => write!(
                formatter,
                "weight field has {actual} entries, expected {expected} (one per dense edge)"
            ),
            Self::ZeroWeight { edge } => write!(
                formatter,
                "weight at edge {edge} is zero: the tight-edge graph then has zero-cost cycles, \
                 so the geodesic count is not finite (zero weights are legal in F1, not here)"
            ),
            Self::ProfileLength { expected, actual } => write!(
                formatter,
                "profile has {actual} entries, expected {expected} (the target count)"
            ),
            Self::NoTargets => write!(formatter, "the estimator needs at least one target"),
            Self::SourceDomain { source, vertices } => {
                write!(
                    formatter,
                    "source {source} outside the graph (0..{vertices})"
                )
            }
            Self::TargetDomain { target, vertices } => {
                write!(
                    formatter,
                    "target {target} outside the graph (0..{vertices})"
                )
            }
            Self::UnreachableTarget { target } => {
                write!(
                    formatter,
                    "target {target} was not reached by the predecessor run"
                )
            }
            Self::PredecessorLength { expected, actual } => write!(
                formatter,
                "predecessor edge array has {actual} entries, expected {expected}"
            ),
            Self::PredecessorDomain {
                vertex,
                predecessor,
                vertices,
            } => write!(
                formatter,
                "predecessor {predecessor} of vertex {vertex} is outside 0..{vertices}"
            ),
            Self::PredecessorEdgeMissing { vertex } => write!(
                formatter,
                "reached vertex {vertex} carries no predecessor edge"
            ),
            Self::PredecessorCycle { target } => write!(
                formatter,
                "the predecessor chain of target {target} did not reach the source"
            ),
            Self::Passage(source) => write!(formatter, "shortest-path pass rejected: {source}"),
            Self::Capacity { buffer, source } => {
                write!(formatter, "could not reserve the {buffer}: {source}")
            }
        }
    }
}

impl std::error::Error for ShapeError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Law(source) => Some(source),
            Self::Passage(source) => Some(source),
            Self::Capacity { source, .. } => Some(source),
            _ => None,
        }
    }
}
