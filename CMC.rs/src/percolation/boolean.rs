//! Boolean (continuum) percolation over overlapping grains (P track,
//! roadmap §14 P2).
//!
//! Initial scope (verbatim from the roadmap): orthorhombic cells;
//! disks/spheres; monodisperse and finite-polydisperse radii; binomial and
//! Poisson point processes. Arbitrary cell geometry, arbitrary grains, and
//! multilevel cell lists remain out of scope until measured workloads
//! justify them.
//!
//! # Composition on P1, not a second geometry
//!
//! P1 ([`crate::percolation::continuum`]) enumerates candidate pairs under
//! a single scalar cutoff. Boolean connectivity needs the per-pair contact
//! distance `r_i + r_j`, which varies with the radii. The composition is
//! the one prescribed by the P track design rule: enumerate candidates at
//! the maximum contact cutoff `2 * r_max` through P1, then keep exactly the
//! pairs whose minimum-image squared distance satisfies
//! `d^2 <= (r_i + r_j)^2` — boundary-inclusive, squared comparison, no
//! epsilon, the same discipline as P1's cutoff test.
//!
//! - **Monodisperse radii** never pay a filter: `r_i + r_j = 2 r` equals
//!   the enumeration cutoff for every pair, and `(r + r)^2` is bit-identical
//!   to P1's `(2r)^2` (doubling and quadrupling are exact), so the P1 pair
//!   list IS the overlap list and passes through unfiltered.
//! - **Per-particle radii** run one in-place `retain` over the P1 pair list.
//!   Each candidate pair's squared distance is recomputed with the particle
//!   module's own metric (`OrthorhombicCell::distance_squared`, identical
//!   fold order to the enumeration), so the filter adds one metric call per
//!   *candidate* pair and duplicates no geometry code. No true overlapping
//!   pair is missed: `r_i + r_j <= 2 r_max` is the enumeration cutoff.
//!
//! The P1 minimum-image ceiling therefore carries over: the contact cutoff
//! `2 r_max` must be at most half the shortest cell side (the reused
//! `CellList` rejects larger cutoffs with a typed error).
//!
//! # Ephemeral [`OverlapGraph`], not a persistent CSR
//!
//! The output follows the documented P1 transient pattern: an owned vertex
//! count plus a flat pair `Vec` that exists for one F3 analysis pass and is
//! then dropped. No adjacency index is built (`incidences` are derived by
//! scanning the pair list on demand); an equilibrium simulation measuring
//! connectivity every sweep must not pay for a persistent per-snapshot CSR.
//!
//! # Sampling versus measurement
//!
//! [`PointProcess`] (binomial with a fixed count, Poisson with a mean)
//! *samples* new configurations under a documented RNG draw contract. For
//! an existing equilibrium particle simulation the module does the opposite
//! — nothing is resampled: the caller supplies coordinates and radii
//! ([`overlap_pairs`], [`overlap_pairs_from_configuration`], or the frozen
//! [`BooleanPercolationMC::measuring`] adapter) and connectivity is a
//! *measurement* of that chain's accepted coordinates, not a second
//! geometry sampler.
//!
//! Realizations are i.i.d.: configure `thermalization_sweeps = 0` and read
//! every sweep as one independent sample. A Poisson realization may draw
//! zero particles — a legitimate vacuum realization, recorded as the
//! all-zero component summary rather than an error (a *request* to analyze
//! an empty coordinate buffer through [`overlap_pairs`] is still the typed
//! P1 rejection: that is a caller bug, not a sample).
//!
//! Known continuum thresholds (for example the 2D disk covered-fraction
//! threshold) are used only as statistical brackets with real error bars —
//! never labeled exact.
//!
//! # Steady-state allocation
//!
//! The adapter itself allocates nothing per realization after construction:
//! positions land in one retained buffer reserved for the count ceiling, the
//! F3 workspace and the query-outcome scratch are reused, and the observable
//! plan records scalars. The allocations that do recur per realization are
//! the documented P1 transient (coordinate copy, cell-list buckets, pair
//! list). The allocation probe pins exactly that steady-state subset: the
//! full realization path performs the same number of allocations as the
//! isolated transient construction replayed over identical coordinate
//! snapshots.

use carlo_rs::{Context, MonteCarlo};
use rand::RngExt;
use rand_xoshiro::Xoshiro256PlusPlus;

use super::continuum::{candidate_pairs_from_configuration, CandidatePairGraph, ContinuumError};
use super::{
    analyze, analyze_with_labels, ComponentWorkspace, EdgeActivity, ObservablePlan,
    ObservablePlanError, SampleView, VertexActivity,
};
use crate::particle::{OrthorhombicCell, ParticleConfiguration};
use crate::topology::{EdgeId, GraphView, Incidence, UndirectedGraphView, VertexId};

/// Chunk ceiling for the Poisson count draw; see [`PointProcess::draw_count`].
const POISSON_CHUNK: f64 = 512.0;

/// Largest accepted Poisson mean. The draw costs `O(mean)` uniform f64s per
/// realization, and a Boolean realization of more than a million grains is
/// beyond the measured P-track workloads; raise the cap when one justifies
/// it.
const MAX_POISSON_MEAN: f64 = 1.0e6;

/// Typed rejection of a Boolean-percolation request.
#[derive(Debug, Clone, PartialEq)]
pub enum BooleanError {
    /// The underlying P1 candidate-pair enumeration rejected the request
    /// (non-finite coordinates, cutoff beyond the minimum-image maximum,
    /// empty coordinates, ...); the P1 diagnostic passes through unchanged.
    Continuum(ContinuumError),
    /// The observable plan is invalid for the realizations this adapter
    /// analyzes (moment observables beyond the exact `f64` integer range).
    Observable(ObservablePlanError),
    /// A radius was non-finite or not strictly positive. `particle` is
    /// `None` for the monodisperse radius and `Some(index)` per particle.
    InvalidRadius {
        particle: Option<usize>,
        radius: f64,
    },
    /// A per-particle radius vector held no entries.
    EmptyRadii,
    /// The per-particle radius length does not match the particle count it
    /// must describe.
    RadiiLengthMismatch { expected: usize, actual: usize },
    /// A binomial point process was constructed with zero particles.
    InvalidParticleCount { count: usize },
    /// The Poisson mean was non-finite or not strictly positive.
    InvalidPoissonMean { mean: f64 },
    /// The Poisson mean exceeds the documented sampling cap.
    PoissonMeanTooLarge { mean: f64, maximum: f64 },
    /// Per-particle radii pin the particle count, so they cannot pair with
    /// a Poisson process whose count is random. Use a binomial process
    /// (count equal to the radius length) or monodisperse radii.
    PerParticleRadiiNeedFixedCount,
    /// Boundary queries carry a fixed vertex domain, but a Poisson
    /// realization has a random particle count. Use a binomial process or
    /// the frozen measuring adapter.
    BoundaryQueriesNeedFixedCount { queries: usize },
    /// A boundary query was built for a different vertex count than the
    /// realizations this adapter analyzes.
    QueryDomain {
        query: usize,
        expected: usize,
        actual: usize,
    },
    /// The plan requests winding observables, which need the continuum
    /// geometry contract of the W track; no such embedding exists here.
    WindingUnsupported,
    /// A retained buffer could not be reserved.
    Capacity,
}

impl core::fmt::Display for BooleanError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Continuum(source) => write!(formatter, "candidate-pair stage: {source}"),
            Self::Observable(source) => write!(formatter, "invalid observable plan: {source}"),
            Self::InvalidRadius { particle, radius } => match particle {
                Some(index) => write!(
                    formatter,
                    "radius of particle {index} must be finite and positive, got {radius}"
                ),
                None => write!(
                    formatter,
                    "monodisperse radius must be finite and positive, got {radius}"
                ),
            },
            Self::EmptyRadii => {
                formatter.write_str("per-particle radii must hold at least one entry")
            }
            Self::RadiiLengthMismatch { expected, actual } => write!(
                formatter,
                "per-particle radii hold {actual} entries but must describe {expected} particles"
            ),
            Self::InvalidParticleCount { count } => {
                write!(
                    formatter,
                    "binomial particle count must be positive, got {count}"
                )
            }
            Self::InvalidPoissonMean { mean } => {
                write!(
                    formatter,
                    "poisson mean must be finite and positive, got {mean}"
                )
            }
            Self::PoissonMeanTooLarge { mean, maximum } => {
                write!(
                    formatter,
                    "poisson mean {mean} exceeds the sampling cap {maximum}"
                )
            }
            Self::PerParticleRadiiNeedFixedCount => formatter.write_str(
                "per-particle radii pin the particle count; pair them with a binomial \
                 process or use monodisperse radii",
            ),
            Self::BoundaryQueriesNeedFixedCount { queries } => write!(
                formatter,
                "{queries} boundary queries need a fixed particle count; a Poisson \
                 realization draws its count at random",
            ),
            Self::QueryDomain {
                query,
                expected,
                actual,
            } => write!(
                formatter,
                "boundary query {query} has V={actual}; realizations have V={expected}"
            ),
            Self::WindingUnsupported => formatter.write_str(
                "winding observables require the continuum W-track contract; no periodic \
                 embedding exists for Boolean grains",
            ),
            Self::Capacity => formatter.write_str("could not reserve a retained buffer"),
        }
    }
}

impl std::error::Error for BooleanError {}

/// Grain radii of a Boolean model: one shared radius, or one radius per
/// particle (finite polydispersity — the caller assigns the species).
///
/// Every radius is validated finite and strictly positive by the
/// constructors; per-particle radii additionally pin the particle count they
/// describe. The enumeration entry points re-validate, so a `Radii` value
/// built directly as an enum literal cannot smuggle an invalid radius into
/// geometry.
#[derive(Clone, Debug, PartialEq)]
pub enum Radii {
    /// One disk/sphere radius shared by every particle.
    Monodisperse {
        /// The shared radius.
        radius: f64,
    },
    /// One radius per particle, particle-major.
    PerParticle {
        /// Radius of particle `i` at index `i`.
        radii: Vec<f64>,
    },
}

impl Radii {
    /// Monodisperse radii with a validated shared radius.
    pub fn monodisperse(radius: f64) -> Result<Self, BooleanError> {
        if !radius.is_finite() || radius <= 0.0 {
            return Err(BooleanError::InvalidRadius {
                particle: None,
                radius,
            });
        }
        Ok(Self::Monodisperse { radius })
    }

    /// Per-particle radii, each validated finite and strictly positive.
    pub fn per_particle(radii: impl IntoIterator<Item = f64>) -> Result<Self, BooleanError> {
        let radii = radii.into_iter().collect::<Vec<_>>();
        if radii.is_empty() {
            return Err(BooleanError::EmptyRadii);
        }
        for (index, &radius) in radii.iter().enumerate() {
            if !radius.is_finite() || radius <= 0.0 {
                return Err(BooleanError::InvalidRadius {
                    particle: Some(index),
                    radius,
                });
            }
        }
        Ok(Self::PerParticle { radii })
    }

    /// Radius of one particle. The caller guarantees `particle` is in range
    /// (the enumeration entry points enforce it against the particle count).
    #[inline]
    pub fn radius(&self, particle: usize) -> f64 {
        match self {
            Self::Monodisperse { radius } => *radius,
            Self::PerParticle { radii } => radii[particle],
        }
    }

    /// Largest radius; the pair-enumeration cutoff is twice this.
    #[inline]
    pub fn max_radius(&self) -> f64 {
        match self {
            Self::Monodisperse { radius } => *radius,
            Self::PerParticle { radii } => radii.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        }
    }

    /// The fixed particle count the radii describe, or `None` when
    /// monodisperse radii apply to any count.
    #[inline]
    pub const fn particle_count(&self) -> Option<usize> {
        match self {
            Self::Monodisperse { .. } => None,
            Self::PerParticle { radii } => Some(radii.len()),
        }
    }

    /// Pair-enumeration cutoff `2 * r_max` handed to P1: the smallest
    /// scalar cutoff that provably contains every overlapping pair.
    #[inline]
    pub fn contact_cutoff(&self) -> f64 {
        2.0 * self.max_radius()
    }

    /// Reject per-particle radii whose length does not describe
    /// `particles` particles.
    fn validate_count(&self, particles: usize) -> Result<(), BooleanError> {
        if let Some(actual) = self.particle_count() {
            if actual != particles {
                return Err(BooleanError::RadiiLengthMismatch {
                    expected: particles,
                    actual,
                });
            }
        }
        Ok(())
    }

    /// Re-validate every radius (guards direct enum-literal construction).
    fn ensure_valid(&self) -> Result<(), BooleanError> {
        match self {
            Self::Monodisperse { radius } => {
                if !radius.is_finite() || *radius <= 0.0 {
                    return Err(BooleanError::InvalidRadius {
                        particle: None,
                        radius: *radius,
                    });
                }
            }
            Self::PerParticle { radii } => {
                if radii.is_empty() {
                    return Err(BooleanError::EmptyRadii);
                }
                for (index, &radius) in radii.iter().enumerate() {
                    if !radius.is_finite() || radius <= 0.0 {
                        return Err(BooleanError::InvalidRadius {
                            particle: Some(index),
                            radius,
                        });
                    }
                }
            }
        }
        Ok(())
    }
}

/// Point process that *samples* Boolean-model configurations.
///
/// Positions are uniform in the primary cell, drawn particle-major with one
/// `random_range(0.0 .. L_axis)` per coordinate (axis ascending); no
/// wrapping is needed because the draws already land in the cell. The count
/// draw is part of each variant's contract below. Realizations are i.i.d. —
/// nothing depends on the previous one.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PointProcess {
    /// Exactly `count` uniform particles (the fixed-`n` binomial process).
    /// The count draw consumes no randomness; only the positions vary.
    Binomial {
        /// Fixed particle count of every realization.
        count: usize,
    },
    /// A Poisson-distributed particle count with mean `mean`.
    Poisson {
        /// Mean particle count of one realization.
        mean: f64,
    },
}

impl PointProcess {
    /// Binomial process with a positive fixed count. Zero particles is a
    /// typed rejection: a deterministically empty process is not a
    /// percolation realization (a Poisson vacuum draw, by contrast, is a
    /// legitimate sample — see the module docs).
    pub fn binomial(count: usize) -> Result<Self, BooleanError> {
        if count == 0 {
            return Err(BooleanError::InvalidParticleCount { count });
        }
        Ok(Self::Binomial { count })
    }

    /// Poisson process with a finite positive mean within the sampling cap.
    pub fn poisson(mean: f64) -> Result<Self, BooleanError> {
        if !mean.is_finite() || mean <= 0.0 {
            return Err(BooleanError::InvalidPoissonMean { mean });
        }
        if mean > MAX_POISSON_MEAN {
            return Err(BooleanError::PoissonMeanTooLarge {
                mean,
                maximum: MAX_POISSON_MEAN,
            });
        }
        Ok(Self::Poisson { mean })
    }

    /// Fixed particle count, or `None` when the count is random (Poisson).
    #[inline]
    pub const fn fixed_count(&self) -> Option<usize> {
        match self {
            Self::Binomial { count } => Some(*count),
            Self::Poisson { .. } => None,
        }
    }

    /// Mean particle count of one realization.
    #[inline]
    pub const fn expected_count(&self) -> f64 {
        match self {
            Self::Binomial { count } => *count as f64,
            Self::Poisson { mean } => *mean,
        }
    }

    /// Draw the particle count of one realization.
    ///
    /// **Binomial:** returns `count`; no RNG consumption.
    ///
    /// **Poisson:** Knuth's multiplication method, chunked. The mean is
    /// split into `ceil(mean / 512)` equal parts `mean / chunks` (one
    /// rounding, so the exact law is `Poisson(chunks * fl(mean / chunks))`,
    /// a sub-ulp relative bias). Each part is sampled by multiplying
    /// standard uniform f64s — one `rng.random::<f64>()` draw per step,
    /// sequential — until the running product first drops to or below
    /// `exp(-part)`, returning the number of draws before that step; the
    /// chunk counts are summed. Poisson additivity makes the sum exact, and
    /// the chunk ceiling keeps `exp(-part) >= e^-512`, far above subnormal
    /// underflow, for every accepted mean.
    pub fn draw_count<R: RngExt + ?Sized>(&self, rng: &mut R) -> usize {
        match *self {
            Self::Binomial { count } => count,
            Self::Poisson { mean } => {
                let chunks = (mean / POISSON_CHUNK).ceil().max(1.0) as usize;
                let part = mean / chunks as f64;
                (0..chunks).map(|_| poisson_knuth(part, rng)).sum()
            }
        }
    }
}

/// Knuth's multiplication method for one Poisson part `<= POISSON_CHUNK`.
fn poisson_knuth<R: RngExt + ?Sized>(part: f64, rng: &mut R) -> usize {
    let threshold = (-part).exp();
    let mut product = 1.0;
    let mut draws = 0;
    loop {
        draws += 1;
        product *= rng.random::<f64>();
        if product <= threshold {
            return draws - 1;
        }
    }
}

/// Ephemeral owned snapshot of one Boolean overlap analysis (the P1
/// transient pattern; see the module docs).
///
/// A vertex per particle and one physical edge per overlapping pair
/// `(i, j)`, `i < j`, in P1's canonical cell-major order (the polydisperse
/// filter removes pairs without reordering). It exists for the duration of
/// one F3 [`crate::percolation::analyze`] pass and is then dropped; no
/// adjacency index is materialized, `incidences` scan the pair list on
/// demand.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OverlapGraph {
    vertex_count: usize,
    pairs: Vec<[usize; 2]>,
}

impl OverlapGraph {
    /// Number of overlapping pairs (physical edges).
    #[inline]
    pub const fn pair_count(&self) -> usize {
        self.pairs.len()
    }

    /// Overlapping pairs in P1's canonical order.
    #[inline]
    pub fn pairs(&self) -> &[[usize; 2]] {
        &self.pairs
    }

    /// Consume the snapshot into its flat pair list.
    #[inline]
    pub fn into_pairs(self) -> Vec<[usize; 2]> {
        self.pairs
    }
}

impl GraphView for OverlapGraph {
    fn vertex_count(&self) -> usize {
        self.vertex_count
    }
}

impl UndirectedGraphView for OverlapGraph {
    fn edge_count(&self) -> usize {
        self.pairs.len()
    }

    fn edge_endpoints(&self, edge: EdgeId) -> [VertexId; 2] {
        let [left, right] = self.pairs[edge.index()];
        [
            VertexId::from_valid_index(left),
            VertexId::from_valid_index(right),
        ]
    }

    /// Derived by scanning the pair list — no adjacency index (same
    /// contract and cost as P1's `PairGraph::incidences`).
    fn incidences(&self, vertex: VertexId) -> impl Iterator<Item = Incidence> + '_ {
        let target = vertex.index();
        self.pairs
            .iter()
            .enumerate()
            .filter_map(move |(index, &[left, right])| {
                if left == target {
                    Some(Incidence {
                        neighbor: VertexId::from_valid_index(right),
                        edge: EdgeId::from_valid_index(index),
                    })
                } else if right == target {
                    Some(Incidence {
                        neighbor: VertexId::from_valid_index(left),
                        edge: EdgeId::from_valid_index(index),
                    })
                } else {
                    None
                }
            })
    }
}

/// Enumerate the overlapping pairs of a flat-coordinate Boolean
/// configuration: particles `i < j` with minimum-image
/// `d(i, j)^2 <= (r_i + r_j)^2`, boundary-inclusive.
///
/// Composition: P1 candidate enumeration at `2 r_max` (see the module docs),
/// with the per-pair contact filter applied in place for per-particle radii
/// and skipped for monodisperse radii (where the P1 list is already exact).
/// The output order is P1's canonical cell-major sequence. Zero particles,
/// non-finite coordinates, and a contact cutoff beyond the minimum-image
/// maximum are the P1 typed rejections, wrapped in
/// [`BooleanError::Continuum`].
pub fn overlap_pairs<const D: usize>(
    coordinates: &[f64],
    cell: &OrthorhombicCell<D>,
    radii: &Radii,
) -> Result<OverlapGraph, BooleanError> {
    // Radii are validated first so radius bugs report as radius errors, not
    // as the downstream cutoff/geometry rejections they trigger.
    radii.ensure_valid()?;
    let adapter = CandidatePairGraph::new(coordinates, cell, radii.contact_cutoff())
        .map_err(BooleanError::Continuum)?;
    radii.validate_count(adapter.particle_count())?;
    let pairs = filter_or_pass(
        adapter
            .pairs()
            .map_err(BooleanError::Continuum)?
            .into_pairs(),
        |first, second| {
            let mut left = [0.0; D];
            let mut right = [0.0; D];
            left.copy_from_slice(&coordinates[first * D..(first + 1) * D]);
            right.copy_from_slice(&coordinates[second * D..(second + 1) * D]);
            (left, right)
        },
        cell,
        radii,
    );
    Ok(OverlapGraph {
        vertex_count: adapter.particle_count(),
        pairs,
    })
}

/// Enumerate overlapping pairs directly from the particle module's accepted
/// configuration (zero-copy equilibrium-chain path).
///
/// Equivalent to [`overlap_pairs`] on the flattened positions; see the
/// module docs for why an equilibrium simulation measures instead of
/// resampling.
pub fn overlap_pairs_from_configuration<const D: usize>(
    configuration: &ParticleConfiguration<D>,
    radii: &Radii,
) -> Result<OverlapGraph, BooleanError> {
    radii.ensure_valid()?;
    radii.validate_count(configuration.len())?;
    let positions = configuration.positions();
    let cell = configuration.cell();
    let pairs = filter_or_pass(
        candidate_pairs_from_configuration(configuration, radii.contact_cutoff())
            .map_err(BooleanError::Continuum)?
            .into_pairs(),
        |first, second| (positions[first], positions[second]),
        cell,
        radii,
    );
    Ok(OverlapGraph {
        vertex_count: positions.len(),
        pairs,
    })
}

/// Apply the per-pair contact filter for per-particle radii, or pass the P1
/// list through untouched for monodisperse radii.
///
/// `positions_of` decodes one pair's endpoints for the two entry points
/// (flat slices versus `&[[f64; D]]`).
fn filter_or_pass<const D: usize>(
    mut pairs: Vec<[usize; 2]>,
    positions_of: impl Fn(usize, usize) -> ([f64; D], [f64; D]),
    cell: &OrthorhombicCell<D>,
    radii: &Radii,
) -> Vec<[usize; 2]> {
    if matches!(radii, Radii::Monodisperse { .. }) {
        return pairs;
    }
    pairs.retain(|&[first, second]| {
        let (left, right) = positions_of(first, second);
        let contact = radii.radius(first) + radii.radius(second);
        cell.distance_squared(&left, &right) <= contact * contact
    });
    pairs
}

/// All-active site mask over one snapshot (every grain exists).
struct AllGrains(usize);

impl VertexActivity for AllGrains {
    fn vertex_count(&self) -> usize {
        self.0
    }

    fn vertex_active(&self, _vertex: VertexId) -> bool {
        true
    }
}

/// All-active bond mask: every overlapping pair connects.
struct AllOverlaps(usize);

impl EdgeActivity for AllOverlaps {
    fn edge_count(&self) -> usize {
        self.0
    }

    fn edge_active(&self, _edge: EdgeId) -> bool {
        true
    }
}

/// Carlo.rs adapter for Boolean percolation (P2).
///
/// One `sweep()` produces one realization and the following `measure()`
/// analyzes it:
///
/// - **Sampling adapters** ([`Self::sampling`]) draw the point process —
///   count, then uniform positions into the retained buffer — every sweep.
/// - **Measuring adapters** ([`Self::measuring`]) own an immutable copy of
///   caller-supplied coordinates and never resample: `sweep()` is a
///   documented no-op and every `measure()` analyzes the identical frozen
///   configuration (the equilibrium-chain measurement contract). A live
///   chain should instead call [`overlap_pairs_from_configuration`] per
///   accepted snapshot, keeping its own coordinates.
///
/// `measure()` runs the F3 component analysis with every grain active and
/// every overlap connecting, and records the F4 [`ObservablePlan`] scalars
/// (`ObservablePlan::all` covers GiantFraction, susceptibility, ...).
/// There is deliberately no degree-array observable: degree support is
/// unbounded for a continuum contact graph, so a fixed-shape C0 array would
/// need an arbitrary cap; the moment scalars carry the same connectivity
/// information. Winding observables are rejected (no continuum embedding
/// exists; that is the W-track contract). Configure
/// `thermalization_sweeps = 0` — samples are i.i.d., and for the measuring
/// adapter thermalization is meaningless by construction.
pub struct BooleanPercolationMC<const D: usize> {
    cell: OrthorhombicCell<D>,
    radii: Radii,
    /// `None` for the frozen measuring adapter.
    process: Option<PointProcess>,
    /// Retained particle-major coordinate buffer: resampled by sampling
    /// sweeps inside the reserved capacity, immutable when frozen.
    coordinates: Vec<f64>,
    plan: ObservablePlan,
    workspace: ComponentWorkspace,
    /// Boundary outcomes of a vacuum realization (always false); reused so
    /// the empty path allocates nothing.
    vacuum_outcomes: Vec<bool>,
}

impl<const D: usize> BooleanPercolationMC<D> {
    /// Sampling adapter: every sweep draws a fresh point-process
    /// realization (binomial or Poisson) of `radii` grains in `cell`.
    pub fn sampling(
        cell: OrthorhombicCell<D>,
        radii: Radii,
        process: PointProcess,
        plan: ObservablePlan,
    ) -> Result<Self, BooleanError> {
        if process.fixed_count().is_none() {
            // Poisson: per-particle radii and boundary queries both need a
            // fixed particle count.
            if radii.particle_count().is_some() {
                return Err(BooleanError::PerParticleRadiiNeedFixedCount);
            }
            if !plan.boundary_queries().is_empty() {
                return Err(BooleanError::BoundaryQueriesNeedFixedCount {
                    queries: plan.boundary_queries().len(),
                });
            }
        }
        let ceiling = count_ceiling(&process);
        Self::from_parts(cell, radii, Some(process), plan, Vec::new(), ceiling)
    }

    /// Measuring adapter: analyze one frozen caller-supplied configuration
    /// (flat particle-major coordinates) identically on every measure; no
    /// resampling ever happens.
    pub fn measuring(
        cell: OrthorhombicCell<D>,
        coordinates: Vec<f64>,
        radii: Radii,
        plan: ObservablePlan,
    ) -> Result<Self, BooleanError> {
        if !coordinates.len().is_multiple_of(D) {
            return Err(BooleanError::Continuum(
                ContinuumError::CoordinateCountMismatch {
                    length: coordinates.len(),
                    dimension: D,
                },
            ));
        }
        if coordinates.is_empty() {
            return Err(BooleanError::Continuum(ContinuumError::EmptyCoordinates {
                dimension: D,
            }));
        }
        for (flat, &value) in coordinates.iter().enumerate() {
            if !value.is_finite() {
                return Err(BooleanError::Continuum(
                    ContinuumError::NonFiniteCoordinate {
                        particle: flat / D,
                        axis: flat % D,
                    },
                ));
            }
        }
        let count = coordinates.len() / D;
        Self::from_parts(cell, radii, None, plan, coordinates, count)
    }

    fn from_parts(
        cell: OrthorhombicCell<D>,
        radii: Radii,
        process: Option<PointProcess>,
        plan: ObservablePlan,
        coordinates: Vec<f64>,
        ceiling: usize,
    ) -> Result<Self, BooleanError> {
        if plan.requests_winding() {
            return Err(BooleanError::WindingUnsupported);
        }
        radii.ensure_valid()?;
        let fixed_count = match process {
            Some(process) => process.fixed_count(),
            None => Some(coordinates.len() / D),
        };
        if let Some(expected) = fixed_count {
            radii.validate_count(expected)?;
            for (index, query) in plan.boundary_queries().iter().enumerate() {
                if query.vertex_count() != expected {
                    return Err(BooleanError::QueryDomain {
                        query: index,
                        expected,
                        actual: query.vertex_count(),
                    });
                }
            }
        }
        // Mirror the reused CellList minimum-image gate so the infallible
        // per-realization path can never hit it at enumeration time.
        let maximum = 0.5 * cell.minimum_length();
        let cutoff = radii.contact_cutoff();
        if cutoff > maximum * (1.0 + 16.0 * f64::EPSILON) {
            return Err(BooleanError::Continuum(ContinuumError::CutoffTooLarge {
                cutoff,
                maximum,
            }));
        }
        plan.validate(ceiling).map_err(BooleanError::Observable)?;
        let mut coordinates = coordinates;
        coordinates
            .try_reserve_exact(ceiling.saturating_mul(D).saturating_sub(coordinates.len()))
            .map_err(|_| BooleanError::Capacity)?;
        let mut workspace = ComponentWorkspace::new();
        workspace
            .prepare(
                ceiling,
                plan.boundary_queries().len(),
                /* with_labels = */ false,
            )
            .map_err(|_| BooleanError::Capacity)?;
        let vacuum_outcomes = vec![false; plan.boundary_queries().len()];
        Ok(Self {
            cell,
            radii,
            process,
            coordinates,
            plan,
            workspace,
            vacuum_outcomes,
        })
    }

    /// The periodic cell every realization lives in.
    pub const fn cell(&self) -> &OrthorhombicCell<D> {
        &self.cell
    }

    /// The grain radii.
    pub const fn radii(&self) -> &Radii {
        &self.radii
    }

    /// The sampled point process, or `None` for the frozen measuring
    /// adapter.
    pub const fn process(&self) -> Option<PointProcess> {
        self.process
    }

    /// The current realization's flat particle-major coordinates (the
    /// frozen buffer for measuring adapters; the last drawn sample for
    /// sampling adapters).
    pub fn coordinates(&self) -> &[f64] {
        &self.coordinates
    }

    /// The scalar F4 observable plan recorded per realization.
    pub const fn observable_plan(&self) -> &ObservablePlan {
        &self.plan
    }

    /// Draw one realization's positions into the retained buffer.
    ///
    /// Draw contract: particle-major, axis-ascending, one
    /// `random_range(0.0 .. L_axis)` uniform per coordinate.
    fn draw_positions<R: RngExt + ?Sized>(&mut self, count: usize, rng: &mut R) {
        self.coordinates.clear();
        let lengths = self.cell.lengths();
        for _ in 0..count {
            for &length in lengths {
                self.coordinates.push(rng.random_range(0.0..length));
            }
        }
    }
}

/// Validation ceiling for the plan and the retained coordinate buffer: the
/// fixed binomial count, or a Poisson excursion guard (`mean + 1024`).
/// A Poisson tail draw beyond the guard grows the retained buffer once and
/// is re-validated against the exact-moment limit by `measure` itself.
fn count_ceiling(process: &PointProcess) -> usize {
    match *process {
        PointProcess::Binomial { count } => count,
        PointProcess::Poisson { mean } => mean as usize + 1024,
    }
}

impl<const D: usize> MonteCarlo for BooleanPercolationMC<D> {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        if let Some(process) = self.process {
            let count = process.draw_count(&mut ctx.rng);
            self.draw_positions(count, &mut ctx.rng);
        }
        // The frozen measuring adapter resamples nothing by contract.
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        let count = self.coordinates.len() / D;
        if count == 0 {
            // Poisson vacuum realization: record the all-zero summary —
            // GiantFraction and the susceptibility scalar are undefined and
            // absent for this sample (the defined indicator records 0).
            self.plan
                .record(&SampleView::vacuum(&self.vacuum_outcomes), ctx);
            return;
        }
        let graph = overlap_pairs(&self.coordinates, &self.cell, &self.radii)
            .expect("constructor validated radii domains, finiteness, and the cutoff ceiling");
        self.plan
            .validate(graph.vertex_count())
            .expect("constructor validated the plan at the count ceiling");
        let result = if self.plan.requests_labels() {
            analyze_with_labels(
                &graph,
                &AllGrains(graph.vertex_count()),
                &AllOverlaps(graph.pair_count()),
                self.plan.boundary_queries(),
                &mut self.workspace,
            )
        } else {
            analyze(
                &graph,
                &AllGrains(graph.vertex_count()),
                &AllOverlaps(graph.pair_count()),
                self.plan.boundary_queries(),
                &mut self.workspace,
            )
        }
        .expect("constructor prepared matching analysis domains");
        let sample = SampleView::new(
            graph.vertex_count(),
            result.summary(),
            result.query_outcomes(),
            None,
            result.labels(),
        );
        self.plan.record(&sample, ctx);
    }

    fn name(&self) -> &'static str {
        "BooleanPercolation"
    }
}

#[cfg(test)]
mod tests {
    use super::super::StaticObservable;
    use super::*;
    use rand::SeedableRng;

    fn flat<const D: usize>(positions: &[[f64; D]]) -> Vec<f64> {
        positions.iter().flat_map(|p| p.iter().copied()).collect()
    }

    fn context() -> Context<Xoshiro256PlusPlus> {
        Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(11), 0, 1)
    }

    #[test]
    fn radii_validation_is_typed() {
        for radius in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                Radii::monodisperse(radius),
                Err(BooleanError::InvalidRadius { particle: None, .. })
            ));
        }
        assert!(matches!(
            Radii::per_particle(Vec::new()),
            Err(BooleanError::EmptyRadii)
        ));
        assert!(matches!(
            Radii::per_particle([0.5, 0.0]),
            Err(BooleanError::InvalidRadius {
                particle: Some(1),
                ..
            })
        ));
        // Direct enum-literal construction is re-validated by the
        // enumeration entry points, not silently trusted.
        let invalid = Radii::Monodisperse { radius: f64::NAN };
        let cell = OrthorhombicCell::new([4.0, 4.0]).unwrap();
        assert!(matches!(
            overlap_pairs(&[1.0, 1.0], &cell, &invalid),
            Err(BooleanError::InvalidRadius { .. })
        ));

        let monodisperse = Radii::monodisperse(0.5).unwrap();
        assert_eq!(monodisperse.radius(7), 0.5);
        assert_eq!(monodisperse.contact_cutoff(), 1.0);
        assert_eq!(monodisperse.particle_count(), None);
        let per_particle = Radii::per_particle([0.3, 1.2, 0.9]).unwrap();
        assert_eq!(per_particle.max_radius(), 1.2);
        assert_eq!(per_particle.contact_cutoff(), 2.4);
        assert_eq!(per_particle.particle_count(), Some(3));
        assert_eq!(per_particle.radius(2), 0.9);
    }

    #[test]
    fn point_process_validation_is_typed() {
        assert!(matches!(
            PointProcess::binomial(0),
            Err(BooleanError::InvalidParticleCount { count: 0 })
        ));
        for mean in [0.0, -3.0, f64::NAN] {
            assert!(matches!(
                PointProcess::poisson(mean),
                Err(BooleanError::InvalidPoissonMean { .. })
            ));
        }
        assert!(matches!(
            PointProcess::poisson(MAX_POISSON_MEAN * 2.0),
            Err(BooleanError::PoissonMeanTooLarge { .. })
        ));
        assert_eq!(PointProcess::binomial(10).unwrap().fixed_count(), Some(10));
        assert_eq!(PointProcess::poisson(10.0).unwrap().fixed_count(), None);
        // A mean at the cap still constructs and draws (chunked path).
        let capped = PointProcess::poisson(MAX_POISSON_MEAN).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(1);
        let drawn = capped.draw_count(&mut rng);
        assert!(
            (drawn as f64 - MAX_POISSON_MEAN).abs() < 5.0e4,
            "draw {drawn} vs mean {MAX_POISSON_MEAN}"
        );
    }

    #[test]
    fn monodisperse_pass_through_matches_filtered_all_equal_radii() {
        // Both composition paths on one geometry: the unfiltered
        // monodisperse pass-through and the per-particle filter with every
        // radius equal must agree exactly.
        let cell = OrthorhombicCell::new([8.0, 8.0]).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x33);
        let positions: Vec<[f64; 2]> = (0..80)
            .map(|_| [rng.random_range(0.0..8.0), rng.random_range(0.0..8.0)])
            .collect();
        let flat = flat(&positions);
        let monodisperse = overlap_pairs(&flat, &cell, &Radii::monodisperse(0.9).unwrap()).unwrap();
        let filtered =
            overlap_pairs(&flat, &cell, &Radii::per_particle([0.9; 80]).unwrap()).unwrap();
        assert!(!monodisperse.pairs().is_empty());
        assert_eq!(monodisperse.pairs(), filtered.pairs());
    }

    #[test]
    fn adapter_rejects_incompatible_domains() {
        let cell = OrthorhombicCell::new([10.0, 10.0]).unwrap();
        let plan = ObservablePlan::all(Vec::new());
        let rejection = |error: Result<BooleanPercolationMC<2>, BooleanError>| {
            error.err().expect("domain mismatch must be rejected")
        };
        // Poisson + per-particle radii.
        assert!(matches!(
            rejection(BooleanPercolationMC::<2>::sampling(
                cell,
                Radii::per_particle([0.3; 10]).unwrap(),
                PointProcess::poisson(10.0).unwrap(),
                plan.clone(),
            )),
            BooleanError::PerParticleRadiiNeedFixedCount
        ));
        // Radii length mismatch against the binomial count.
        assert_eq!(
            rejection(BooleanPercolationMC::<2>::sampling(
                cell,
                Radii::per_particle([0.3; 9]).unwrap(),
                PointProcess::binomial(10).unwrap(),
                plan.clone(),
            )),
            BooleanError::RadiiLengthMismatch {
                expected: 10,
                actual: 9
            }
        );
        // Winding observables.
        let winding = ObservablePlan::new([StaticObservable::Wrapping], Vec::new());
        assert!(matches!(
            rejection(BooleanPercolationMC::<2>::sampling(
                cell,
                Radii::monodisperse(0.3).unwrap(),
                PointProcess::binomial(10).unwrap(),
                winding,
            )),
            BooleanError::WindingUnsupported
        ));
        // Non-finite frozen coordinates.
        assert!(matches!(
            rejection(BooleanPercolationMC::<2>::measuring(
                cell,
                vec![0.0, f64::NAN],
                Radii::monodisperse(0.3).unwrap(),
                plan.clone(),
            )),
            BooleanError::Continuum(ContinuumError::NonFiniteCoordinate { .. })
        ));
        // Radius too large for the minimum-image ceiling.
        assert!(matches!(
            rejection(BooleanPercolationMC::<2>::sampling(
                cell,
                Radii::monodisperse(3.0).unwrap(),
                PointProcess::binomial(10).unwrap(),
                plan,
            )),
            BooleanError::Continuum(ContinuumError::CutoffTooLarge { .. })
        ));
    }

    #[test]
    fn vacuum_realization_records_zero_summary() {
        // Poisson mean small enough that empty realizations occur: they
        // must record zero rows, never fail.
        let cell = OrthorhombicCell::new([20.0, 20.0]).unwrap();
        let mut adapter = BooleanPercolationMC::<2>::sampling(
            cell,
            Radii::monodisperse(0.2).unwrap(),
            PointProcess::poisson(0.7).unwrap(),
            ObservablePlan::new(
                [
                    StaticObservable::ActiveVertexCount,
                    StaticObservable::ComponentCount,
                    StaticObservable::FiniteClusterSusceptibility,
                ],
                Vec::new(),
            ),
        )
        .unwrap();
        let mut ctx = context();
        let mut empties = 0;
        for _ in 0..400 {
            adapter.sweep(&mut ctx);
            if adapter.coordinates().is_empty() {
                empties += 1;
            }
            adapter.measure(&mut ctx);
        }
        assert!(empties > 0, "fixture must draw vacuum realizations");
        let values = ctx.finalize_measurements();
        // Every sample — vacuum included — contributes its row.
        assert_eq!(values["ActiveVertexCount"].n_bins, 400);
        assert_eq!(values["ComponentCount"].n_bins, 400);
        // Expected count 0.7: mean inside a generous 5-sigma interval.
        let mean = values["ActiveVertexCount"].mean;
        assert!((mean - 0.7).abs() < 5.0 * (0.7_f64 / 400.0).sqrt());
    }
}
