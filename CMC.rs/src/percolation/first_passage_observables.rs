//! F2: first-passage percolation observables on Z^2 over the F1 weighted
//! shortest-path layer (roadmap v2 §3 F2).
//!
//! # Geometries and the passage-time observable
//!
//! [`PassageGeometry`] builds the two F2 source/target layouts on the open
//! Z^2 box (dense vertex id `x + y * width`, the [`build_square`]
//! row-major convention):
//!
//! - **Corridor** ([`PassageGeometry::corridor`]): a `length` r >= 1 by
//!   `width` w >= 1 corridor; the **source set** is the left column x = 0
//!   and the **target set** the right column x = r (covered distance r).
//!   This is the fluctuation geometry of the validation gates.
//! - **Segment** ([`PassageGeometry::segment`]): a single **origin**
//!   vertex and an explicit list of target points along the +x row, s
//!   steps apart, inside a box with symmetric longitudinal `margin` and
//!   `half_height` transversal room. One profile entry per target — the
//!   point-to-point passage time `T(0, x)` of the literature, measured
//!   several distances at once in one box.
//!
//! A realization draws one weight per box edge and produces the
//! **passage-time profile**: one `f64` per target vertex,
//! `T_y = min over paths from any source-set vertex to target y` of the
//! sum of edge weights. One realization is ONE C0 vector sample of fixed
//! shape (w for corridors, one entry per listed step for segments) — never
//! the scalar minimum alone; the scalar passage time `T_r = min_y T_y` is
//! the derived accessor [`PassageSampler::passage_time`].
//!
//! The minimum over the source set costs ONE F1 pass: the topology carries
//! one extra **super-source** vertex (the last dense vertex) joined to
//! every source-set vertex by structural zero-weight edges (the last dense
//! edges, never assigned drawn weights), so
//! `dist(super-source, y) = min_{u in source set} T(u, y)` exactly. Zero
//! weights are legal Dijkstra inputs by the F1 contract.
//!
//! Path degeneracy never reaches this observable: `T_y` is a minimum over
//! path WEIGHTS, a unique value however many optimal paths realize it.
//! The scalar minimum may be attained at several targets — the value is
//! still unique, and no entry is ever singled out (the profile keeps them
//! all). Geodesic reconstruction (`GeodesicTracking`) is deliberately not
//! requested here; that is F3's layer.
//!
//! # Weight laws and the bitwise draw contract
//!
//! [`EdgeWeightLaw`] carries the two standard FPP weight distributions,
//! both applied to a uniform draw `u in [0, 1)`:
//! [`EdgeWeightLaw::Exponential`] `w = -ln(1-u)/rate` (rate 1 is the
//! classical FPP model) and [`EdgeWeightLaw::Uniform`] `w = u`. The draw
//! contract, replayable bitwise from any seed: **one `rng.random::<f64>()`
//! per box edge, dense edge order** ([`build_square`]'s per-site X-then-Y
//! order), then — consuming no draws — the super-source edges receive
//! structural `0.0`. [`PassageSampler::weights`] exposes the slice of the
//! last run so the contract can be verified against an independent
//! replay, and [`PassageSampler::realize_with_weights`] runs production on
//! a caller-supplied slice (the exact-enumeration path: the weights are
//! data, exactly as in F1).
//!
//! # Fluctuation and sample machinery
//!
//! The fluctuation observable is the sample variance of the scalar passage
//! time across realizations; the survey anchor (Auffinger-Damron-Hanson,
//! *50 years of first-passage percolation*, AMS ULS 68 (2017), archived in
//! `docs/ref/`) states `Var T(0,x) ~ ||x||^{2*chi}` with the rigorous bound
//! `0 <= chi <= 1/2` (Kesten) and the physicists' prediction `chi = 1/3`
//! in d = 2. This layer only MEASURES: [`PassageEstimator`] accumulates
//! mean / unbiased variance / standard error of the scalar `T_r` plus
//! per-target mean and variance of the profile (vector-mean accumulation);
//! any exponent or time-constant claim is a statistical bracket produced
//! by the validation layer, never an output here.

use core::fmt;

use rand::{Rng, RngExt};

use super::first_passage::{
    first_passage_run, FirstPassageError, FirstPassageWorkspace, GeodesicTracking,
};
use super::process::ProcessError;
use crate::{build_square, Bond, BondType, CsrLattice};

/// Source/target layout on the open Z^2 box plus the F2 super-source
/// topology: the corridor (column to column) and the segment (origin to
/// explicit points).
///
/// See the module documentation for the geometry, boundary
/// conventions, and dense edge order. The geometry owns its
/// [`CsrLattice`](crate::CsrLattice); a [`PassageSampler`] takes it by
/// value (clone first if several samplers share one description).
#[derive(Clone, Debug)]
pub struct PassageGeometry {
    lattice: CsrLattice,
    /// Largest source-to-target distance (corridor length; last segment step).
    length: usize,
    /// Dense ids of the target vertices — the profile coordinates.
    targets: Vec<usize>,
    /// Box edges carrying drawn weights (dense edges `0..drawn_edges`).
    drawn_edges: usize,
}

impl PassageGeometry {
    /// Column-to-column corridor of `length` r >= 1 passage steps and
    /// `width` w >= 1 rows: source set the left column, target set the
    /// right column (w profile entries, row order).
    ///
    /// Panics on `length == 0` or `width == 0` — the same
    /// programmatic-validity convention as the crate's lattice builders.
    pub fn corridor(length: usize, width: usize) -> Self {
        assert!(length >= 1, "corridor length must be at least 1");
        assert!(width >= 1, "corridor width must be at least 1");
        let row_stride = length + 1;
        let attachments: Vec<usize> = (0..width).map(|row| row * row_stride).collect();
        let targets: Vec<usize> = (0..width).map(|row| length + row * row_stride).collect();
        Self::with_super_source(length, length + 1, width, &attachments, targets)
    }

    /// Origin-to-points segment: an open box of `2*half_height + 1` rows
    /// and `steps.last() + 2*margin + 1` columns; the super-source joins
    /// the single origin vertex (`margin` columns from the left on the
    /// vertically centered row), and the target set holds one vertex per
    /// entry of `steps` (ascending, unique, each >= 1), that many steps to
    /// the right of the origin on the same row — one profile entry per
    /// step.
    ///
    /// With `steps = [r, 2r]` the box is symmetric under the reflection
    /// about its vertical center line, which swaps the origin and the `2r`
    /// target and fixes the `r` target; the i.i.d. weight law is invariant
    /// under that reflection, so `E T(0, r e) = E T(r e, 2r e)` exactly and
    /// path concatenation gives the exact-in-expectation subadditivity
    /// `E T_{2r} <= 2 E T_r` — the validation gates' anchor inequality.
    /// (The corridor observable instead carries entry/exit surface
    /// corrections and converges to its cylinder constant from below; the
    /// two families bracket complementary physics and both are exposed.)
    ///
    /// Panics on empty `steps`, non-ascending or repeated steps, or a
    /// zero step.
    pub fn segment(steps: &[usize], margin: usize, half_height: usize) -> Self {
        assert!(!steps.is_empty(), "segment needs at least one step");
        assert!(
            steps.windows(2).all(|pair| pair[0] < pair[1]),
            "segment steps must be strictly ascending"
        );
        assert!(steps[0] >= 1, "segment steps start at 1");
        let length = steps[steps.len() - 1];
        let box_width = length + 2 * margin + 1;
        let origin = margin + half_height * box_width;
        let targets: Vec<usize> = steps.iter().map(|&step| origin + step).collect();
        Self::with_super_source(length, box_width, 2 * half_height + 1, &[origin], targets)
    }

    /// The shared construction: the open `box_width x box_height` Z^2 box
    /// from [`build_square`], then one super-source vertex joined to every
    /// attachment vertex by a structural zero-weight edge.
    fn with_super_source(
        length: usize,
        box_width: usize,
        box_height: usize,
        attachments: &[usize],
        targets: Vec<usize>,
    ) -> Self {
        let box_lattice = build_square(box_width, box_height, false);
        let drawn_edges = box_lattice.n_edges();
        let super_source = box_lattice.n_sites;
        let mut bonds = box_lattice.edges;
        for &vertex in attachments {
            // Structural edges: their weight lives in the sampler's slice
            // (always 0.0); the Bond field mirrors that for consistency.
            bonds.push(Bond::new(super_source, vertex, BondType::Generic, 0.0));
        }
        let lattice = CsrLattice::from_edges(super_source + 1, bonds);
        debug_assert_eq!(lattice.n_edges(), drawn_edges + attachments.len());
        Self {
            lattice,
            length,
            targets,
            drawn_edges,
        }
    }

    /// The geometry's topology (box plus super-source).
    #[inline]
    pub fn lattice(&self) -> &CsrLattice {
        &self.lattice
    }

    /// Largest source-to-target distance (corridor length; last segment step).
    #[inline]
    pub const fn length(&self) -> usize {
        self.length
    }

    /// Dense vertex id of the super-source (one past the last box vertex).
    #[inline]
    pub fn super_source(&self) -> usize {
        self.lattice.n_sites - 1
    }

    /// Dense ids of the target vertices (the profile coordinates); corridor
    /// row order, segment step order.
    #[inline]
    pub fn targets(&self) -> &[usize] {
        &self.targets
    }

    /// Profile length: the target count (the corridor width, or the number
    /// of segment steps).
    #[inline]
    pub fn boundary_size(&self) -> usize {
        self.targets.len()
    }

    /// Box edges that carry drawn weights (dense edges `0..drawn_edges`).
    #[inline]
    pub const fn drawn_edge_count(&self) -> usize {
        self.drawn_edges
    }

    /// All dense edges: drawn box edges plus the structural zero-weight
    /// super-source edges.
    #[inline]
    pub fn total_edge_count(&self) -> usize {
        self.lattice.n_edges()
    }
}

/// Open extension point for a first-passage edge-weight law.
///
/// A law is a pure transform of one uniform draw: it validates its own
/// parameters and maps `u in [0, 1)` to a non-negative weight. Keeping the
/// transform separate from the realization sampler preserves the bitwise
/// draw contract: the sampler owns the RNG, the law owns only the mapping.
pub trait WeightLaw {
    /// Law-specific rejection type.
    type Error;

    /// Reject invalid parameters before any realization is drawn.
    fn validate(&self) -> Result<(), Self::Error>;

    /// The weight for one uniform draw.
    fn draw(&self, uniform: f64) -> f64;
}

impl WeightLaw for EdgeWeightLaw {
    type Error = PassageObservablesError;

    fn validate(&self) -> Result<(), PassageObservablesError> {
        EdgeWeightLaw::validate(self)
    }

    fn draw(&self, uniform: f64) -> f64 {
        EdgeWeightLaw::draw(self, uniform)
    }
}

/// Edge-weight law for first-passage realizations: a pure transform of one
/// uniform `u in [0, 1)` draw.
///
/// Both members are standard FPP choices; see the module documentation
/// for the bitwise draw contract.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum EdgeWeightLaw {
    /// Exponential with rate `rate` (mean `1/rate`):
    /// `w = -ln(1-u)/rate`. The classical FPP model (Richardson 1973).
    /// `rate` must be finite and strictly positive.
    Exponential { rate: f64 },
    /// Uniform on [0, 1): `w = u`, the raw draw.
    Uniform,
}

impl EdgeWeightLaw {
    /// Exponential law (mean `1/rate`). Also constructible as
    /// [`EdgeWeightLaw::Exponential`]` { rate }`; invalid rates are
    /// rejected at realization time with the offending value.
    pub const fn exponential(rate: f64) -> Self {
        Self::Exponential { rate }
    }

    /// Validate the law's parameters.
    pub fn validate(&self) -> Result<(), PassageObservablesError> {
        match *self {
            Self::Exponential { rate } if rate.is_finite() && rate > 0.0 => Ok(()),
            Self::Exponential { rate } => Err(PassageObservablesError::InvalidRate { rate }),
            Self::Uniform => Ok(()),
        }
    }

    /// The weight for one uniform draw `u in [0, 1)` (the documented
    /// transform; total for validated laws). Written as
    /// `0.0 - ln(1-u)/rate` so the draw `u = 0` yields `+0.0`, never
    /// `-0.0` (a negative zero weight would pass F1 validation but
    /// bitwise-compare unequally against a naive `-x` replay).
    #[inline]
    pub fn draw(&self, uniform: f64) -> f64 {
        match *self {
            Self::Exponential { rate } => 0.0 - (1.0 - uniform).ln() / rate,
            Self::Uniform => uniform,
        }
    }
}

/// Reusable realization sampler: the weight slice, the F1 workspace, and
/// the profile buffer of the last realization.
///
/// `new` reserves every buffer once; [`PassageSampler::realize`] then
/// performs zero steady-state allocation per realization (draws land in
/// the reused weight slice, the F1 workspace is pre-prepared). The
/// returned profile of each call is invalidated by the next call on this
/// sampler — copy it ([`PassageSampler::profile`]) or feed an estimator
/// before the next realization.
#[derive(Debug)]
pub struct PassageSampler {
    geometry: PassageGeometry,
    /// One `f64` per dense edge: drawn box weights then structural zeros.
    weights: Vec<f64>,
    /// One `f64` per target vertex (the C0 vector sample).
    profile: Vec<f64>,
    workspace: FirstPassageWorkspace,
}

impl PassageSampler {
    /// Reserve all buffers for `geometry` (taken by value) and pre-prepare
    /// the F1 workspace.
    pub fn new(geometry: PassageGeometry) -> Result<Self, PassageObservablesError> {
        let edges = geometry.total_edge_count();
        let mut weights = Vec::new();
        weights
            .try_reserve_exact(edges)
            .map_err(|source| PassageObservablesError::Capacity(ProcessError::from(source)))?;
        weights.resize(edges, 0.0);
        let mut profile = Vec::new();
        profile
            .try_reserve_exact(geometry.boundary_size())
            .map_err(|source| PassageObservablesError::Capacity(ProcessError::from(source)))?;
        profile.resize(geometry.boundary_size(), f64::INFINITY);
        let mut workspace = FirstPassageWorkspace::new();
        workspace
            .prepare(geometry.lattice().n_sites, edges)
            .map_err(PassageObservablesError::Passage)?;
        Ok(Self {
            geometry,
            weights,
            profile,
            workspace,
        })
    }

    /// One realization: draw the box weights under `law` from `rng` (the
    /// documented bitwise contract), write the structural zeros, run one
    /// F1 pass from the super-source, and return the passage-time profile
    /// along the target set.
    pub fn realize<R: Rng + ?Sized>(
        &mut self,
        law: EdgeWeightLaw,
        rng: &mut R,
    ) -> Result<&[f64], PassageObservablesError> {
        law.validate()?;
        let drawn = self.geometry.drawn_edge_count();
        for weight in &mut self.weights[..drawn] {
            *weight = law.draw(rng.random::<f64>());
        }
        for weight in &mut self.weights[drawn..] {
            *weight = 0.0;
        }
        self.run()
    }

    /// One realization on a caller-supplied weight slice: one `f64` per
    /// dense edge (drawn box edges in dense order, then the super-source
    /// entries — zero under the documented geometries). The
    /// exact-enumeration validation path; weights are data exactly as in
    /// F1, so non-finite or negative entries are typed rejections from the
    /// underlying pass and a rejected call leaves the last profile intact.
    pub fn realize_with_weights(
        &mut self,
        weights: &[f64],
    ) -> Result<&[f64], PassageObservablesError> {
        let expected = self.geometry.total_edge_count();
        if weights.len() != expected {
            return Err(PassageObservablesError::WeightLength {
                expected,
                actual: weights.len(),
            });
        }
        self.weights.copy_from_slice(weights);
        // On a downstream rejection the sampler's internal weight buffer
        // retains these values (the caller's slice is never mutated).
        self.run()
    }

    /// The single F1 pass and the profile copy; both realize paths share it.
    fn run(&mut self) -> Result<&[f64], PassageObservablesError> {
        let super_source = self.geometry.super_source();
        first_passage_run(
            &self.geometry.lattice,
            &self.weights,
            super_source,
            GeodesicTracking::Off,
            &mut self.workspace,
        )
        .map_err(PassageObservablesError::Passage)?;
        let distances = self.workspace.distances();
        for (slot, &target) in self.profile.iter_mut().zip(self.geometry.targets()) {
            *slot = distances[target];
        }
        Ok(&self.profile)
    }

    /// The geometry this sampler was built from.
    #[inline]
    pub fn geometry(&self) -> &PassageGeometry {
        &self.geometry
    }

    /// The passage-time profile of the last realization (one `f64` per
    /// target vertex). All [`f64::INFINITY`] before the first realization.
    #[inline]
    pub fn profile(&self) -> &[f64] {
        &self.profile
    }

    /// The weight slice of the last realization (one `f64` per dense edge:
    /// drawn box weights, then the structural zeros) — the replay handle
    /// for the bitwise draw contract.
    #[inline]
    pub fn weights(&self) -> &[f64] {
        &self.weights
    }

    /// The scalar passage time `T_r = min_y T_y` of the last realization —
    /// a derived accessor over the profile, never the primary sample.
    /// [`f64::INFINITY`] before the first realization.
    pub fn passage_time(&self) -> f64 {
        self.profile.iter().copied().fold(f64::INFINITY, f64::min)
    }
}

/// Accumulates independent passage-time profiles into sample statistics.
///
/// The scalar statistics (mean, unbiased variance, standard error) describe
/// the scalar passage time `T_r = min_y T_y` computed per profile by
/// [`Self::add`]; the vector accumulators keep the per-target mean and
/// variance of the profile. All buffers are reserved at `new`, so `add`
/// allocates nothing per realization.
#[derive(Clone, Debug)]
pub struct PassageEstimator {
    samples: usize,
    sum: f64,
    sum_sq: f64,
    profile_sum: Vec<f64>,
    profile_sum_sq: Vec<f64>,
}

impl PassageEstimator {
    /// Reserve the accumulators for profiles of `boundary_size` entries.
    pub fn new(boundary_size: usize) -> Result<Self, PassageObservablesError> {
        let mut profile_sum = Vec::new();
        profile_sum
            .try_reserve_exact(boundary_size)
            .map_err(|source| PassageObservablesError::Capacity(ProcessError::from(source)))?;
        profile_sum.resize(boundary_size, 0.0);
        let mut profile_sum_sq = Vec::new();
        profile_sum_sq
            .try_reserve_exact(boundary_size)
            .map_err(|source| PassageObservablesError::Capacity(ProcessError::from(source)))?;
        profile_sum_sq.resize(boundary_size, 0.0);
        Ok(Self {
            samples: 0,
            sum: 0.0,
            sum_sq: 0.0,
            profile_sum,
            profile_sum_sq,
        })
    }

    /// Realizations accumulated so far.
    #[inline]
    pub const fn samples(&self) -> usize {
        self.samples
    }

    /// The reserved boundary size.
    #[inline]
    pub fn boundary_size(&self) -> usize {
        self.profile_sum.len()
    }

    /// Add one realization's profile (the scalar passage time is its
    /// minimum, computed here). Length mismatches are typed rejections;
    /// nothing is added on rejection.
    pub fn add(&mut self, profile: &[f64]) -> Result<(), PassageObservablesError> {
        if profile.len() != self.profile_sum.len() {
            return Err(PassageObservablesError::ProfileLength {
                expected: self.profile_sum.len(),
                actual: profile.len(),
            });
        }
        let passage = profile.iter().copied().fold(f64::INFINITY, f64::min);
        self.samples += 1;
        self.sum += passage;
        self.sum_sq += passage * passage;
        for ((sum, sum_sq), &value) in self
            .profile_sum
            .iter_mut()
            .zip(self.profile_sum_sq.iter_mut())
            .zip(profile)
        {
            *sum += value;
            *sum_sq += value * value;
        }
        Ok(())
    }

    /// Sample statistics, or `None` until at least two realizations.
    ///
    /// Allocates the two profile vectors of the returned estimate (a
    /// report-time cost, never per realization).
    pub fn estimate(&self) -> Option<PassageEstimate> {
        if self.samples < 2 {
            return None;
        }
        let samples = self.samples as f64;
        let mean = self.sum / samples;
        let variance = ((self.sum_sq - self.sum * self.sum / samples) / (samples - 1.0)).max(0.0);
        let profile_mean = self.profile_sum.iter().map(|&sum| sum / samples).collect();
        let profile_variance = self
            .profile_sum_sq
            .iter()
            .zip(self.profile_sum.iter())
            .map(|(&sum_sq, &sum)| ((sum_sq - sum * sum / samples) / (samples - 1.0)).max(0.0))
            .collect();
        Some(PassageEstimate {
            samples: self.samples,
            mean,
            variance,
            standard_error: (variance / samples).sqrt(),
            profile_mean,
            profile_variance,
        })
    }
}

/// Sample statistics over independent first-passage realizations.
#[derive(Clone, Debug, PartialEq)]
pub struct PassageEstimate {
    /// Realizations accumulated.
    pub samples: usize,
    /// Sample mean of the scalar passage time `T_r`.
    pub mean: f64,
    /// Unbiased sample variance of `T_r` — the fluctuation observable
    /// (bracketed against `Var T_r ~ r^{2*chi}`, never claimed exactly).
    pub variance: f64,
    /// Standard error of [`Self::mean`] across realizations.
    pub standard_error: f64,
    /// Per-target sample mean of the profile (target order).
    pub profile_mean: Vec<f64>,
    /// Per-target unbiased sample variance of the profile.
    pub profile_variance: Vec<f64>,
}

/// Typed rejections of the first-passage observables layer.
#[derive(Debug)]
pub enum PassageObservablesError {
    /// An exponential rate that is not finite and strictly positive.
    InvalidRate { rate: f64 },
    /// A weight slice without exactly one entry per dense edge.
    WeightLength { expected: usize, actual: usize },
    /// A profile whose length differs from the estimator's boundary size.
    ProfileLength { expected: usize, actual: usize },
    /// A buffer reservation failure.
    Capacity(ProcessError),
    /// A rejection raised by the underlying F1 pass.
    Passage(FirstPassageError),
}

impl fmt::Display for PassageObservablesError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidRate { rate } => write!(
                formatter,
                "exponential rate {rate} must be finite and strictly positive"
            ),
            Self::WeightLength { expected, actual } => write!(
                formatter,
                "weight field has {actual} entries, expected {expected} (one per dense edge)"
            ),
            Self::ProfileLength { expected, actual } => write!(
                formatter,
                "profile has {actual} entries, expected {expected} (the boundary size)"
            ),
            Self::Capacity(source) => {
                write!(
                    formatter,
                    "could not prepare first-passage observable buffers: {source}"
                )
            }
            Self::Passage(source) => write!(formatter, "passage-time pass rejected: {source}"),
        }
    }
}

impl std::error::Error for PassageObservablesError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Capacity(source) => Some(source),
            Self::Passage(source) => Some(source),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn weight_law_trait_is_bitwise_identical_to_the_inherent_transform() {
        for law in [EdgeWeightLaw::Uniform, EdgeWeightLaw::exponential(2.5)] {
            <EdgeWeightLaw as WeightLaw>::validate(&law).expect("valid law");
            for step in 0..=1000 {
                let uniform = f64::from(step) / 1000.0;
                assert_eq!(
                    <EdgeWeightLaw as WeightLaw>::draw(&law, uniform).to_bits(),
                    law.draw(uniform).to_bits()
                );
            }
        }
    }

    #[test]
    fn weight_law_rejects_invalid_rates_before_drawing() {
        let law = EdgeWeightLaw::exponential(0.0);
        assert!(matches!(
            <EdgeWeightLaw as WeightLaw>::validate(&law),
            Err(PassageObservablesError::InvalidRate { rate: 0.0 })
        ));
    }
}
