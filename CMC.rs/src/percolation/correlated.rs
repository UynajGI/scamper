//! C1: power-law correlated Gaussian activity fields (roadmap v2 §5, stage C1).
//!
//! The stage realizes the Schrenk/Pruessner Fourier-filtering generator for
//! site occupation on a periodic hypercubic lattice: a Gaussian random field
//! whose two-point correlation decays as a power law is thresholded at a
//! fixed quantile, giving a configuration whose per-site marginal is exactly
//! the requested occupation probability `p` while its sites are correlated.
//!
//! ## Anchors and conventions (quoted exactly)
//!
//! - Schrenk et al., PRE 88, 052102 (2013), archived
//!   `docs/ref/md/Schrenk_2013_correlated_disorder_percolation.md`: the
//!   Fourier-filtering method imposes the power spectrum
//!   `S(f) ~ |f|^{-beta_c} = (f_1^2 + f_2^2)^{-beta_c/2}` (their Eq. 1,
//!   md lines ~142-146) with `beta_c = 2(H + 1)`, which by the
//!   Wiener-Khinchin theorem gives the height correlation
//!   `c(r) = <h(x) h(x+r)>_x ~ |r|^{2H}` (their Eq. 2, md lines ~160-163),
//!   and their footnote 33 (md line ~1422) maps their `H` to the
//!   correlation exponent `a` used here by `a = -2H`, i.e.
//!   `c(r) ~ |r|^{-a}`. Schrenk's correlated range is `-1 <= H <= 0`,
//!   equivalently `0 <= a <= 2`.
//! - Pruessner et al., PRE 97, 062141 (2018), archived
//!   `docs/ref/md/Pruessner_2018_correlated_thresholds.md`: the correlation
//!   parameter is defined by `C(r) ~ |r|^{-a}` (their Eq. 1, md line 62);
//!   their discrete FFM (Sec. II and Appendix A, md lines ~219-344 and
//!   ~1293-1390) takes a symmetric target `C_r`, transforms it to the
//!   discrete spectral density `S_k` (Eq. 8), sets negative `S_k` to zero
//!   ("zero-cutoff", Appendix B, md line ~1407), draws
//!   `Phi_k^{re/im} = sqrt(S_k) U` with `U` Gaussian of mean 0 and variance
//!   `N` (Eq. 11), and inverse-transforms. Their target is
//!   `C(r) = (1 + r^2)^{-a/2}` with `sigma_tau^2 = C_0 = 1` (Eq. 12,
//!   md lines ~380-390). This module uses exactly that target family.
//! - Saberi, EPJ B 89, 253 (2016), archived
//!   `docs/ref/md/Saberi_2015_percolation_review.md`, md lines ~239-241:
//!   with long-range correlations "percolation clusters become more
//!   compact", the qualitative trend C2's fractal-dimension observable
//!   tracks.
//!
//! ## The field contract (bitwise reproducible draws)
//!
//! The lattice is a product of power-of-two sides `L_0, ..., L_{d-1}`
//! (radix-2 axes), with dense vertex rank
//! `x = sum_i coord_i * stride_i`, `stride_0 = 1`,
//! `stride_i = prod_{j<i} L_j` — the row-major order of
//! [`crate::build_hypercubic`] (axis 0 is the fastest index). The
//! displacement squared between two sites is the periodic minimum-image
//! Euclidean norm.
//!
//! 1. **Design spectrum.** `C_r = (1 + r^2)^{-a/2}` on every periodic
//!    displacement `r` (including `r = 0`, where `C_0 = 1`), transformed
//!    with the forward DFT `S_k = sum_r C_r exp(-2 pi i k.r / L)`.
//! 2. **Zero-cutoff and normalization.** `S_k <- max(S_k, 0)`; with
//!    `v = (1/N) sum_k S_k`, the per-mode complex power is
//!    `P_k = N S_k / v`. Then `sum_k P_k = N^2` exactly (up to rounding),
//!    so every site of the realized field has variance exactly `1`.
//! 3. **Mode draws (Hermitian symmetry).** Modes are visited in ascending
//!    dense rank `k`; `-k` is the componentwise negation modulo `L`.
//!    A self-conjugate mode (`k = -k`, every component `0` or `L_i/2`)
//!    draws one standard normal `g` and sets
//!    `X_k = sqrt(P_k) g` (real). A conjugate pair draws one Box-Muller
//!    pair `(g1, g2)` at its canonical (lower-rank) member and sets
//!    `X_k = sqrt(P_k/2) (g1 + i g2)`, `X_{-k} = conj(X_k)`. The lower-rank
//!    member is always visited first because dense rank is lexicographic.
//! 4. **Inverse transform.** `f_x = (1/N) sum_k X_k exp(+2 pi i k.x / L)`;
//!    the imaginary residue (roundoff, order `1e-16`) is discarded.
//!    [`CorrelatedGaussianField::sample`] consumes exactly one uniform per
//!    Box-Muller variate (two per conjugate pair, one per self-conjugate
//!    mode), using `u1 = 1 - U[0,1)` so `ln u1` is finite.
//!
//! The draws are reproducible bitwise for a fixed toolchain and platform
//! (the same qualification as the L1 draw contract): the transform order,
//! the mode order, and the Box-Muller formula are fixed here, but last-bit
//! `sin`/`cos` differences across libm implementations are out of scope.
//!
//! ## Occupation law and the `a -> infinity` reduction
//!
//! [`CorrelatedSiteLaw`] occupies site `x` iff `f_x > Phi^{-1}(1 - p)`.
//! Because every site marginal is exactly standard normal, the occupation
//! probability is exactly `p` in the ensemble; the inverse-CDF rational
//! approximation contributes at most `1.2e-9` absolute probability error
//! (Acklam's approximation, documented on
//! [`inverse_normal_cdf`]). Every physical edge is active, i.e. a site
//! configuration.
//!
//! **`a = +infinity` is the stage's iid reference by contract.** The FFM
//! limit `a -> infinity` gives `C_r = delta_{r0}` (white noise), which the
//! platform defines as iid Bernoulli and implements by delegating to
//! [`SiteBernoulli`] with the same RNG stream and no preliminary draws —
//! bitwise identical activity to `StaticLaw::site(0)` under the same seed,
//! mirroring the R2 `q = 1` delegation pattern. This is the
//! uncorrelated reference row of Pruessner's table (`a = infinity`,
//! 0.592746), so that row is exact rather than merely statistically white.
//! By contrast the opposite limit `a = 0` gives `C_r = 1` for every `r`: a
//! CONSTANT field, not white noise. The thresholded law is then
//! all-occupied or all-empty and carries no percolation transition; the
//! constant field remains available from [`CorrelatedGaussianField`] for
//! completeness, and no iid claim is attached to it. For
//! `0 < a < infinity` the field path runs.
//!
//! ## Adapter and allocation contract
//!
//! [`CorrelatedPercolationMC`] is the owned scheduler adapter: one `sweep()`
//! draws one field and thresholds it into its
//! [`StaticConfiguration`]; the following `measure()` runs the unchanged F3
//! [`analyze`] pass and records the existing [`ObservablePlan`] scalars.
//! Samples are i.i.d. — configure `thermalization_sweeps = 0`. Every buffer
//! (field, spectrum, masks, union-find workspace) is allocated in the
//! constructor; steady-state sweeps and measurements allocate nothing.
//! Winding observables are not part of this adapter (as in the R2 adapter);
//! use [`analyze_winding`](super::analyze_winding) on
//! [`CorrelatedPercolationMC::configuration`] with an embedding for those.
//!
//! No FFT crate is used: the transform is a fixed-size radix-2
//! Cooley-Tukey over power-of-two axes, written here so the draw contract
//! (mode order, Hermitian pairing, normalization) stays explicit and the
//! workspace can be reused without allocation.

#![allow(clippy::neg_cmp_op_on_partial_ord)]

use core::f64::consts::PI;
use core::fmt;

use carlo_rs::{Context, MonteCarlo};
use rand::{Rng, RngExt};
use rand_xoshiro::Xoshiro256PlusPlus;

use super::recording::record_component_analysis;
use super::{
    ComponentWorkspace, ObservablePlan, ObservablePlanError, OccupationLaw, Probability,
    ProbabilityError, SamplingError, SiteBernoulli, StaticConfiguration,
};
use crate::{CsrLattice, GraphView, UndirectedGraphView};

/// Failure to build a correlated field, a correlated law, or its adapter.
#[derive(Clone, Debug, PartialEq)]
pub enum CorrelatedError {
    /// The side vector is empty (a lattice needs at least one axis).
    EmptySides,
    /// A side is smaller than two; periodicity with one site per axis has
    /// no translation structure to correlate.
    SideTooSmall { axis: usize, side: usize },
    /// A side is not a power of two (the radix-2 transform contract).
    SideNotPowerOfTwo { axis: usize, side: usize },
    /// The vertex count does not fit a `usize`.
    LatticeTooLarge { sides: Vec<usize> },
    /// The exponent is negative or `NaN` (`+inf` is the white-noise limit
    /// and only the occupation law accepts it).
    InvalidExponent { exponent: f64 },
    /// A required buffer could not be reserved.
    Capacity,
    /// The graph vertex count differs from the field geometry's vertex count.
    GeometryMismatch { expected: usize, actual: usize },
    /// A configuration domain differs from the graph's domain.
    Sampling(SamplingError),
    /// The occupation probability is outside `[0, 1]` or not finite.
    Probability(ProbabilityError),
    /// The observable plan failed its own validation.
    Observable(ObservablePlanError),
    /// Winding observables need an embedding and are outside this adapter.
    WindingUnsupported,
    /// The adapter topology failed validation.
    Topology(String),
    /// A boundary query was built for a different vertex domain.
    QueryDomain {
        query: usize,
        expected: usize,
        actual: usize,
    },
    /// The requested displacement does not match the field dimension.
    DisplacementDimension { expected: usize, actual: usize },
}

impl fmt::Display for CorrelatedError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptySides => write!(formatter, "correlated field needs at least one axis"),
            Self::SideTooSmall { axis, side } => write!(
                formatter,
                "correlated field axis {axis} has side {side}; every side must be at least 2"
            ),
            Self::SideNotPowerOfTwo { axis, side } => write!(
                formatter,
                "correlated field axis {axis} has side {side}; every side must be a power of two"
            ),
            Self::LatticeTooLarge { sides } => write!(
                formatter,
                "correlated field dimensions {sides:?} exceed the addressable vertex count"
            ),
            Self::InvalidExponent { exponent } => write!(
                formatter,
                "correlation exponent must be finite and non-negative (or +inf), got {exponent}"
            ),
            Self::Capacity => write!(formatter, "correlated field workspace allocation failed"),
            Self::GeometryMismatch { expected, actual } => write!(
                formatter,
                "graph has {actual} vertices; correlated field geometry has {expected}"
            ),
            Self::Sampling(error) => error.fmt(formatter),
            Self::Probability(error) => error.fmt(formatter),
            Self::Observable(error) => error.fmt(formatter),
            Self::WindingUnsupported => write!(
                formatter,
                "the correlated adapter does not record winding observables; analyze its \
                 configuration with an embedding instead"
            ),
            Self::Topology(message) => write!(formatter, "correlated adapter topology: {message}"),
            Self::QueryDomain {
                query,
                expected,
                actual,
            } => write!(
                formatter,
                "boundary query {query} covers {actual} vertices; adapter has {expected}"
            ),
            Self::DisplacementDimension { expected, actual } => write!(
                formatter,
                "displacement has {actual} components; the field has {expected}"
            ),
        }
    }
}

impl std::error::Error for CorrelatedError {}

/// Geometry and exponent of a periodic power-law correlated Gaussian field.
///
/// Sides must be powers of two (the radix-2 contract) and at least two. The
/// exponent `a` is the target correlation exponent,
/// `C(r) = (1 + r^2)^{-a/2} ~ |r|^{-a}`; `a = 0` is accepted by the raw
/// field (constant field, documented) but the occupation law redirects
/// `a = +inf` to its iid delegate; `a = 0` is the constant-field limit.
#[derive(Clone, Debug, PartialEq)]
pub struct CorrelatedFieldSpec {
    sides: Vec<usize>,
    exponent: f64,
    vertices: usize,
    strides: Vec<usize>,
}

impl CorrelatedFieldSpec {
    /// Validate a geometry and exponent; no buffer is allocated.
    pub fn new(sides: &[usize], exponent: f64) -> Result<Self, CorrelatedError> {
        if sides.is_empty() {
            return Err(CorrelatedError::EmptySides);
        }
        if !(exponent >= 0.0) {
            // Rejects NaN and negative values; +inf is the white-noise limit.
            return Err(CorrelatedError::InvalidExponent { exponent });
        }
        let mut strides = Vec::with_capacity(sides.len());
        let mut product = 1_u128;
        for (axis, &side) in sides.iter().enumerate() {
            if side < 2 {
                return Err(CorrelatedError::SideTooSmall { axis, side });
            }
            if !side.is_power_of_two() {
                return Err(CorrelatedError::SideNotPowerOfTwo { axis, side });
            }
            strides.push(product as usize);
            product *= side as u128;
            if product > usize::MAX as u128 {
                return Err(CorrelatedError::LatticeTooLarge {
                    sides: sides.to_vec(),
                });
            }
        }
        Ok(Self {
            sides: sides.to_vec(),
            exponent,
            vertices: product as usize,
            strides,
        })
    }

    /// Per-axis side lengths.
    #[inline]
    pub fn sides(&self) -> &[usize] {
        &self.sides
    }

    /// The correlation exponent `a`.
    #[inline]
    pub const fn exponent(&self) -> f64 {
        self.exponent
    }

    /// The dimension `d`.
    #[inline]
    pub fn dimension(&self) -> usize {
        self.sides.len()
    }

    /// Total vertex count `N = prod_i L_i`.
    #[inline]
    pub const fn vertices(&self) -> usize {
        self.vertices
    }

    /// Dense-rank strides (axis 0 fastest), mirroring [`crate::build_hypercubic`].
    #[inline]
    pub fn strides(&self) -> &[usize] {
        &self.strides
    }

    /// Decode a dense rank into integer coordinates.
    pub fn coordinates(&self, index: usize) -> Option<Vec<usize>> {
        if index >= self.vertices {
            return None;
        }
        let mut remaining = index;
        let mut coordinates = vec![0; self.sides.len()];
        for axis in (0..self.sides.len()).rev() {
            coordinates[axis] = remaining / self.strides[axis];
            remaining %= self.strides[axis];
        }
        Some(coordinates)
    }
}

fn try_filled<T: Copy>(length: usize, value: T) -> Result<Vec<T>, CorrelatedError> {
    let mut vector = Vec::new();
    vector
        .try_reserve_exact(length)
        .map_err(|_| CorrelatedError::Capacity)?;
    vector.resize(length, value);
    Ok(vector)
}

/// Per-axis radix-2 tables: bit reversal and the twiddle table
/// `exp(-2 pi i k / L)` for `k in 0..L/2`.
struct AxisTransform {
    length: usize,
    reversed: Vec<usize>,
    cos: Vec<f64>,
    sin: Vec<f64>,
}

impl AxisTransform {
    fn new(length: usize) -> Result<Self, CorrelatedError> {
        let bits = length.trailing_zeros();
        let mut reversed = try_filled(length, 0_usize)?;
        for (index, slot) in reversed.iter_mut().enumerate() {
            *slot = index.reverse_bits() >> (usize::BITS - bits);
        }
        let half = length / 2;
        let mut cos = try_filled(half, 0.0)?;
        let mut sin = try_filled(half, 0.0)?;
        for (k, (cosine, sine)) in cos.iter_mut().zip(sin.iter_mut()).enumerate() {
            let angle = 2.0 * PI * (k as f64) / (length as f64);
            *cosine = angle.cos();
            *sine = angle.sin();
        }
        Ok(Self {
            length,
            reversed,
            cos,
            sin,
        })
    }
}

/// In-place radix-2 Cooley-Tukey transform of one axis.
///
/// The forward sign is `exp(-i...)`, the inverse `exp(+i...)` followed by a
/// `1/L` scaling. The caller gathers/scatters; this works on contiguous
/// slices of exactly `axis.length`.
fn transform_1d(re: &mut [f64], im: &mut [f64], axis: &AxisTransform, inverse: bool) {
    let n = axis.length;
    debug_assert_eq!(re.len(), n);
    debug_assert_eq!(im.len(), n);
    for index in 0..n {
        let reversed = axis.reversed[index];
        if index < reversed {
            re.swap(index, reversed);
            im.swap(index, reversed);
        }
    }
    let sign = if inverse { 1.0 } else { -1.0 };
    let mut length = 2;
    while length <= n {
        let half = length / 2;
        let step = n / length;
        let mut start = 0;
        while start < n {
            for offset in 0..half {
                let twiddle = offset * step;
                let (cos, sin) = (axis.cos[twiddle], sign * axis.sin[twiddle]);
                let even = start + offset;
                let odd = even + half;
                let odd_re = re[odd] * cos - im[odd] * sin;
                let odd_im = re[odd] * sin + im[odd] * cos;
                let even_re = re[even];
                let even_im = im[even];
                re[even] = even_re + odd_re;
                im[even] = even_im + odd_im;
                re[odd] = even_re - odd_re;
                im[odd] = even_im - odd_im;
            }
            start += length;
        }
        length <<= 1;
    }
    if inverse {
        let scale = 1.0 / n as f64;
        for value in re.iter_mut() {
            *value *= scale;
        }
        for value in im.iter_mut() {
            *value *= scale;
        }
    }
}

/// Draw one standard normal by Box-Muller with `u1 = 1 - U[0,1)`.
#[inline]
fn standard_normal<R: Rng + ?Sized>(rng: &mut R) -> f64 {
    let u1 = 1.0 - rng.random::<f64>();
    let u2 = rng.random::<f64>();
    (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
}

/// Draw a Box-Muller pair: two independent standard normals from two uniforms.
#[inline]
fn standard_normal_pair<R: Rng + ?Sized>(rng: &mut R) -> (f64, f64) {
    let u1 = 1.0 - rng.random::<f64>();
    let u2 = rng.random::<f64>();
    let radius = (-2.0 * u1.ln()).sqrt();
    let angle = 2.0 * PI * u2;
    (radius * angle.cos(), radius * angle.sin())
}

/// Reusable periodic power-law correlated Gaussian field.
///
/// Construction builds the design spectrum once (one forward transform);
/// [`Self::sample`] overwrites the reused spectrum and field buffers, so
/// steady-state sampling allocates nothing. See the module documentation
/// for the exact draw contract.
pub struct CorrelatedGaussianField {
    spec: CorrelatedFieldSpec,
    /// `P_k = N S_k / v`, the target `E|X_k|^2` of each mode.
    power: Vec<f64>,
    /// Dense rank of the conjugate partner `-k`, precomputed.
    partner: Vec<usize>,
    spectrum_re: Vec<f64>,
    spectrum_im: Vec<f64>,
    values: Vec<f64>,
    axes: Vec<AxisTransform>,
    scratch_re: Vec<f64>,
    scratch_im: Vec<f64>,
}

impl CorrelatedGaussianField {
    /// Build the design spectrum for `spec` (one forward transform).
    pub fn new(spec: CorrelatedFieldSpec) -> Result<Self, CorrelatedError> {
        let vertices = spec.vertices();
        let dimension = spec.dimension();
        let mut axes = Vec::with_capacity(dimension);
        for &side in spec.sides() {
            axes.push(AxisTransform::new(side)?);
        }
        let scratch = spec.sides().iter().copied().max().unwrap_or(0);
        let power = try_filled(vertices, 0.0)?;
        let partner = try_filled(vertices, 0_usize)?;
        let mut spectrum_re = try_filled(vertices, 0.0)?;
        let spectrum_im = try_filled(vertices, 0.0)?;
        let values = try_filled(vertices, 0.0)?;
        let scratch_re = try_filled(scratch, 0.0)?;
        let scratch_im = try_filled(scratch, 0.0)?;

        // 1. Design correlation C_r = (1 + r^2)^(-a/2) on the torus.
        let mut coordinates = vec![0_usize; dimension];
        for (index, slot) in spectrum_re.iter_mut().enumerate() {
            decode(&spec, index, &mut coordinates);
            let mut squared = 0.0;
            for (axis, &coordinate) in coordinates.iter().enumerate() {
                let side = spec.sides()[axis];
                let minimum_image = coordinate.min(side - coordinate);
                squared += (minimum_image as f64) * (minimum_image as f64);
            }
            *slot = (1.0 + squared).powf(-0.5 * spec.exponent());
        }

        let mut field = Self {
            spec,
            power,
            partner,
            spectrum_re,
            spectrum_im,
            values,
            axes,
            scratch_re,
            scratch_im,
        };
        field.build_partner_map();
        field.forward_transform();
        field.build_power();
        field.spectrum_re.fill(0.0);
        field.spectrum_im.fill(0.0);
        Ok(field)
    }

    /// The geometry and exponent this field realizes.
    #[inline]
    pub const fn spec(&self) -> &CorrelatedFieldSpec {
        &self.spec
    }

    /// The most recently sampled field (length `N`).
    #[inline]
    pub fn values(&self) -> &[f64] {
        &self.values
    }

    /// The design correlation `C(delta)/(v)` of the cutoff spectrum,
    /// `(1/N^2) sum_k P_k exp(2 pi i k.delta / L)`.
    ///
    /// This is the ensemble covariance the sampler targets after the
    /// zero-cutoff and normalization steps; compare sampled correlations
    /// against it. Returns `None` for a displacement with the wrong
    /// dimension.
    pub fn design_correlation(&self, displacement: &[isize]) -> Option<f64> {
        if displacement.len() != self.spec.dimension() {
            return None;
        }
        let mut coordinates = vec![0_usize; self.spec.dimension()];
        let mut total = 0.0;
        for (index, &power) in self.power.iter().enumerate() {
            decode(&self.spec, index, &mut coordinates);
            let mut phase = 0.0;
            for (axis, &coordinate) in coordinates.iter().enumerate() {
                let side = self.spec.sides()[axis] as isize;
                let product = (coordinate as isize) * displacement[axis];
                phase += 2.0 * PI * (product.rem_euclid(side) as f64) / (side as f64);
            }
            total += power * phase.cos();
        }
        let normalization = (self.spec.vertices() as f64) * (self.spec.vertices() as f64);
        Some(total / normalization)
    }

    /// Draw one independent field realization.
    pub fn sample<R: Rng + ?Sized>(&mut self, rng: &mut R) {
        for index in 0..self.spec.vertices() {
            let partner = self.partner[index];
            if partner == index {
                let draw = standard_normal(rng);
                self.spectrum_re[index] = self.power[index].sqrt() * draw;
                self.spectrum_im[index] = 0.0;
            } else if index < partner {
                let (real, imaginary) = standard_normal_pair(rng);
                let amplitude = (0.5 * self.power[index]).sqrt();
                self.spectrum_re[index] = amplitude * real;
                self.spectrum_im[index] = amplitude * imaginary;
            } else {
                self.spectrum_re[index] = self.spectrum_re[partner];
                self.spectrum_im[index] = -self.spectrum_im[partner];
            }
        }
        self.inverse_transform();
        self.values.copy_from_slice(&self.spectrum_re);
        debug_assert!(
            self.spectrum_im.iter().all(|value| value.abs() < 1e-9),
            "Hermitian inverse transform left a non-negligible imaginary residue"
        );
    }

    /// Precompute the conjugate partner of every dense mode rank.
    fn build_partner_map(&mut self) {
        let dimension = self.spec.dimension();
        let mut coordinates = vec![0_usize; dimension];
        for index in 0..self.spec.vertices() {
            decode(&self.spec, index, &mut coordinates);
            let mut partner = 0_usize;
            for (axis, coordinate) in coordinates.iter_mut().enumerate() {
                let side = self.spec.sides()[axis];
                let negated = (side - *coordinate) % side;
                partner += negated * self.spec.strides()[axis];
            }
            self.partner[index] = partner;
        }
    }

    /// Zero-cutoff, normalize, and cache the per-mode power.
    fn build_power(&mut self) {
        let vertices = self.spec.vertices();
        let mut sum = 0.0;
        for value in &mut self.spectrum_re {
            if *value < 0.0 {
                *value = 0.0;
            }
            sum += *value;
        }
        debug_assert!(sum > 0.0, "C_0 = 1 keeps the DC mode positive");
        let variance = sum / vertices as f64;
        for (power, &spectral) in self.power.iter_mut().zip(self.spectrum_re.iter()) {
            *power = vertices as f64 * spectral / variance;
        }
    }

    fn forward_transform(&mut self) {
        self.transform(false);
    }

    fn inverse_transform(&mut self) {
        self.transform(true);
    }

    /// Apply the 1D radix-2 transform along every axis.
    fn transform(&mut self, inverse: bool) {
        let vertices = self.spec.vertices();
        for (axis, plan) in self.axes.iter().enumerate() {
            let side = plan.length;
            let stride = self.spec.strides()[axis];
            let block = side * stride;
            for outer in 0..vertices / block {
                for offset in 0..stride {
                    let base = outer * block + offset;
                    for index in 0..side {
                        let source = base + index * stride;
                        self.scratch_re[index] = self.spectrum_re[source];
                        self.scratch_im[index] = self.spectrum_im[source];
                    }
                    transform_1d(&mut self.scratch_re, &mut self.scratch_im, plan, inverse);
                    for index in 0..side {
                        let target = base + index * stride;
                        self.spectrum_re[target] = self.scratch_re[index];
                        self.spectrum_im[target] = self.scratch_im[index];
                    }
                }
            }
        }
    }
}

/// Decode a dense rank into `coordinates` (axis 0 fastest).
fn decode(spec: &CorrelatedFieldSpec, index: usize, coordinates: &mut [usize]) {
    let mut remaining = index;
    for axis in (0..spec.dimension()).rev() {
        coordinates[axis] = remaining / spec.strides()[axis];
        remaining %= spec.strides()[axis];
    }
}

/// Acklam's rational approximation to the standard normal quantile.
///
/// `p` is clamped to `[0, 1]`; `0` and `1` map to `-inf`/`+inf`. The
/// absolute probability error is at most `1.15e-9` (Acklam 2003), so a
/// threshold at this quantile moves the marginal by at most `1.2e-9` —
/// far below any sampling resolution. No Halley refinement is applied
/// (that would need an `erf` implementation).
pub fn inverse_normal_cdf(p: f64) -> f64 {
    #[allow(clippy::excessive_precision)]
    const A: [f64; 6] = [
        -3.969_683_028_665_38e1,
        2.209_460_984_245_205e2,
        -2.759_285_104_469_687e2,
        1.383_577_518_672_690e2,
        -3.066_479_806_614_716e1,
        2.506_628_277_459_239,
    ];
    const B: [f64; 5] = [
        -5.447_609_879_822_406e1,
        1.615_858_368_580_409e2,
        -1.556_989_798_598_866e2,
        6.680_131_188_771_972e1,
        -1.328_068_155_288_572e1,
    ];
    const C: [f64; 6] = [
        -7.784_894_002_430_293e-3,
        -3.223_964_580_411_365e-1,
        -2.400_758_277_161_838,
        -2.549_732_539_343_734,
        4.374_664_141_464_968,
        2.938_163_982_698_783,
    ];
    const D: [f64; 4] = [
        7.784_695_709_041_462e-3,
        3.224_671_290_700_398e-1,
        2.445_134_137_142_996,
        3.754_408_661_907_416,
    ];
    const LOW: f64 = 0.024_25;
    const HIGH: f64 = 1.0 - LOW;

    if p.is_nan() {
        return f64::NAN;
    }
    if p <= 0.0 {
        return f64::NEG_INFINITY;
    }
    if p >= 1.0 {
        return f64::INFINITY;
    }
    if p < LOW {
        let q = (-2.0 * p.ln()).sqrt();
        return (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0);
    }
    if p <= HIGH {
        let q = p - 0.5;
        let r = q * q;
        return (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q
            / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0);
    }
    let q = (-2.0 * (1.0 - p).ln()).sqrt();
    -(((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
        / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
}

/// Correlated site occupation law with an exact per-site marginal.
///
/// The field geometry must match the graph's dense vertex order: use a
/// periodic hypercubic builder ([`crate::build_square`],
/// [`crate::build_hypercubic`]) with the same sides. The law validates the
/// vertex count only; a graph with `N` vertices in a different order is a
/// caller contract violation (documented, not detectable cheaply).
///
/// `a = +inf` delegates bitwise to [`SiteBernoulli`] (see the
/// module documentation); every other `0 < a < inf` draws one field per
/// [`Self::sample`].
pub struct CorrelatedSiteLaw {
    spec: CorrelatedFieldSpec,
    probability: Probability,
    threshold: f64,
    field: Option<CorrelatedGaussianField>,
}

impl CorrelatedSiteLaw {
    /// Build the law; allocates the reusable field when `0 < a < inf`.
    pub fn new(sides: &[usize], exponent: f64, probability: f64) -> Result<Self, CorrelatedError> {
        let spec = CorrelatedFieldSpec::new(sides, exponent)?;
        let probability = Probability::new(probability).map_err(CorrelatedError::Probability)?;
        let threshold = inverse_normal_cdf(1.0 - probability.get());
        let field = if exponent.is_infinite() {
            None
        } else {
            Some(CorrelatedGaussianField::new(spec.clone())?)
        };
        Ok(Self {
            spec,
            probability,
            threshold,
            field,
        })
    }

    /// The field geometry and exponent.
    #[inline]
    pub const fn spec(&self) -> &CorrelatedFieldSpec {
        &self.spec
    }

    /// The occupation probability `p`.
    #[inline]
    pub const fn probability(&self) -> Probability {
        self.probability
    }

    /// The threshold `Phi^{-1}(1 - p)` applied to the field.
    #[inline]
    pub const fn threshold(&self) -> f64 {
        self.threshold
    }

    /// Whether this law takes the iid delegate (`a = +inf`).
    #[inline]
    pub const fn is_iid_reduction(&self) -> bool {
        self.field.is_none()
    }

    /// The most recently sampled field, or `None` under the iid delegate.
    #[inline]
    pub fn field(&self) -> Option<&CorrelatedGaussianField> {
        self.field.as_ref()
    }

    /// Draw one occupation configuration.
    ///
    /// Domain validation happens before any RNG draw, so a rejected call
    /// leaves both the configuration and the RNG stream untouched.
    pub fn sample<G: UndirectedGraphView, R: Rng + ?Sized>(
        &mut self,
        graph: &G,
        configuration: &mut StaticConfiguration,
        rng: &mut R,
    ) -> Result<(), CorrelatedError> {
        if graph.vertex_count() != self.spec.vertices() {
            return Err(CorrelatedError::GeometryMismatch {
                expected: self.spec.vertices(),
                actual: graph.vertex_count(),
            });
        }
        let Some(field) = &mut self.field else {
            // a = +inf: the iid delegate, same stream, no draws here.
            return SiteBernoulli::new(self.probability.into())
                .sample(graph, configuration, rng)
                .map_err(CorrelatedError::Sampling);
        };
        if configuration.vertex_count() != self.spec.vertices()
            || configuration.edge_count() != graph.edge_count()
        {
            return Err(CorrelatedError::Sampling(
                SamplingError::ConfigurationSize {
                    expected_vertices: self.spec.vertices(),
                    actual_vertices: configuration.vertex_count(),
                    expected_edges: graph.edge_count(),
                    actual_edges: configuration.edge_count(),
                },
            ));
        }
        field.sample(rng);
        let threshold = self.threshold;
        for (active, &value) in configuration
            .vertex_mask_mut()
            .dense_mut()
            .iter_mut()
            .zip(field.values())
        {
            *active = u8::from(value > threshold);
        }
        configuration.edge_mask_mut().fill(true);
        Ok(())
    }
}

impl<G: UndirectedGraphView> OccupationLaw<G> for CorrelatedSiteLaw {
    type Error = CorrelatedError;

    fn sample<R: Rng + ?Sized>(
        &mut self,
        graph: &G,
        configuration: &mut StaticConfiguration,
        rng: &mut R,
    ) -> Result<(), CorrelatedError> {
        CorrelatedSiteLaw::sample(self, graph, configuration, rng)
    }
}

/// Owned scheduler adapter for correlated site percolation.
///
/// One `sweep()` draws one correlated field and thresholds it into the
/// owned [`StaticConfiguration`]; `measure()` runs the unchanged F3
/// [`analyze`] pass and records the [`ObservablePlan`] scalars (boundary
/// queries included). Samples are i.i.d.: set `thermalization_sweeps = 0`.
/// Winding observables are rejected (use the core analyzer with an
/// embedding). Construction materializes every buffer, so steady-state
/// sweeps and measurements allocate nothing.
pub struct CorrelatedPercolationMC {
    lattice: CsrLattice,
    law: CorrelatedSiteLaw,
    configuration: StaticConfiguration,
    workspace: ComponentWorkspace,
    plan: ObservablePlan,
}

impl CorrelatedPercolationMC {
    /// Build the adapter; validates topology, geometry, plan, and queries.
    pub fn new(
        lattice: CsrLattice,
        law: CorrelatedSiteLaw,
        plan: ObservablePlan,
    ) -> Result<Self, CorrelatedError> {
        if plan.requests_winding() {
            return Err(CorrelatedError::WindingUnsupported);
        }
        lattice
            .validate()
            .map_err(|error| CorrelatedError::Topology(error.to_string()))?;
        if lattice.vertex_count() != law.spec().vertices() {
            return Err(CorrelatedError::GeometryMismatch {
                expected: law.spec().vertices(),
                actual: lattice.vertex_count(),
            });
        }
        plan.validate(lattice.vertex_count())
            .map_err(CorrelatedError::Observable)?;
        for (index, query) in plan.boundary_queries().iter().enumerate() {
            if query.vertex_count() != lattice.vertex_count() {
                return Err(CorrelatedError::QueryDomain {
                    query: index,
                    expected: lattice.vertex_count(),
                    actual: query.vertex_count(),
                });
            }
        }
        let mut configuration =
            StaticConfiguration::new(lattice.vertex_count(), lattice.edge_count());
        // Materialize the dense mask now so steady-state sweeps only write.
        configuration.vertex_mask_mut().dense_mut();
        let mut workspace = ComponentWorkspace::new();
        workspace
            .prepare(lattice.vertex_count(), plan.boundary_queries().len(), false)
            .map_err(|_| CorrelatedError::Capacity)?;
        Ok(Self {
            lattice,
            law,
            configuration,
            workspace,
            plan,
        })
    }

    /// The adapter topology.
    #[inline]
    pub const fn lattice(&self) -> &CsrLattice {
        &self.lattice
    }

    /// The correlated occupation law.
    #[inline]
    pub const fn law(&self) -> &CorrelatedSiteLaw {
        &self.law
    }

    /// The current configuration (the last drawn realization).
    #[inline]
    pub const fn configuration(&self) -> &StaticConfiguration {
        &self.configuration
    }

    /// The explicit observable plan.
    #[inline]
    pub const fn observable_plan(&self) -> &ObservablePlan {
        &self.plan
    }
}

impl MonteCarlo for CorrelatedPercolationMC {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        self.law
            .sample(&self.lattice, &mut self.configuration, &mut ctx.rng)
            .expect("constructor validated the geometry and configuration domains");
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        record_component_analysis(
            &self.lattice,
            &self.configuration,
            &self.configuration,
            &self.plan,
            &mut self.workspace,
            None,
            ctx,
        );
    }

    fn name(&self) -> &'static str {
        "CorrelatedPercolation"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{build_chain, build_square, VertexActivity};
    use rand::SeedableRng;

    fn field(sides: &[usize], exponent: f64) -> CorrelatedGaussianField {
        CorrelatedGaussianField::new(CorrelatedFieldSpec::new(sides, exponent).unwrap()).unwrap()
    }

    fn flags(graph: &impl UndirectedGraphView, configuration: &StaticConfiguration) -> Vec<bool> {
        graph
            .vertex_ids()
            .map(|vertex| configuration.vertex_active(vertex))
            .collect()
    }

    #[test]
    fn radix_two_transform_round_trips() {
        let plan = AxisTransform::new(16).unwrap();
        let original: Vec<f64> = (0..16).map(|index| (index * 7 % 13) as f64 - 6.0).collect();
        let mut re = original.clone();
        let mut im = vec![0.0; 16];
        transform_1d(&mut re, &mut im, &plan, false);
        transform_1d(&mut re, &mut im, &plan, true);
        for (index, &value) in original.iter().enumerate() {
            assert!((re[index] - value).abs() < 1e-12, "index {index}");
            assert!(im[index].abs() < 1e-12, "index {index}");
        }

        // A known forward transform: the constant vector maps to (L, 0, ...).
        let mut re = vec![1.0; 8];
        let mut im = vec![0.0; 8];
        transform_1d(&mut re, &mut im, &AxisTransform::new(8).unwrap(), false);
        assert!((re[0] - 8.0).abs() < 1e-12);
        assert!(re[1..].iter().all(|value| value.abs() < 1e-12));
    }

    #[test]
    fn unit_variance_and_design_correlation_hold_in_the_ensemble() {
        let mut field = field(&[32, 32], 1.5);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(11);
        let replicas = 400;
        let side = 32_usize;
        let mut mean = 0.0;
        let mut square = 0.0;
        let mut covariance = [0.0; 9];
        for _ in 0..replicas {
            field.sample(&mut rng);
            let values = field.values();
            for &value in values {
                mean += value;
                square += value * value;
            }
            for (lag, slot) in covariance.iter_mut().enumerate() {
                let lag = lag as isize;
                for y in 0..side {
                    for x in 0..side {
                        let here = values[x + side * y];
                        let shifted = (x as isize + lag).rem_euclid(side as isize) as usize;
                        *slot += here * values[shifted + side * y];
                    }
                }
            }
        }
        let sites = (replicas * side * side) as f64;
        let average = mean / sites;
        let variance = square / sites - average * average;
        assert!(average.abs() < 0.03, "mean {average}");
        assert!((variance - 1.0).abs() < 0.08, "variance {variance}");
        for (lag, slot) in covariance.iter_mut().enumerate() {
            *slot /= (replicas * side * side) as f64;
            let design = field
                .design_correlation(&[lag as isize, 0])
                .expect("two-dimensional displacement");
            assert!(
                (*slot - design).abs() < 0.05,
                "lag {lag}: sampled {} vs design {design}",
                *slot
            );
        }
    }

    #[test]
    fn marginal_is_p_and_the_fraction_fluctuates_beyond_binomial() {
        let graph = build_square(16, 16, true);
        let probability = 0.55;
        let mut law = CorrelatedSiteLaw::new(&[16, 16], 1.5, probability).unwrap();
        let iid = SiteBernoulli::new(Probability::new(probability).unwrap().into());
        let mut configuration = StaticConfiguration::new(256, graph.edge_count());
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(21);
        let replicas = 500;
        let mut correlated = 0.0;
        let mut correlated_square = 0.0;
        let mut independent = 0.0;
        let mut independent_square = 0.0;
        for _ in 0..replicas {
            law.sample(&graph, &mut configuration, &mut rng).unwrap();
            let occupied = flags(&graph, &configuration)
                .into_iter()
                .filter(|&active| active)
                .count() as f64
                / 256.0;
            correlated += occupied;
            correlated_square += occupied * occupied;

            iid.sample(&graph, &mut configuration, &mut rng).unwrap();
            let occupied = flags(&graph, &configuration)
                .into_iter()
                .filter(|&active| active)
                .count() as f64
                / 256.0;
            independent += occupied;
            independent_square += occupied * occupied;
        }
        let mean = correlated / replicas as f64;
        let variance = correlated_square / replicas as f64 - mean * mean;
        let standard_error = (variance / replicas as f64).sqrt();
        assert!(
            (mean - probability).abs() <= 3.0 * standard_error.max(0.002),
            "occupation fraction {mean} vs p {probability} (se {standard_error})"
        );
        let iid_mean = independent / replicas as f64;
        let iid_variance = independent_square / replicas as f64 - iid_mean * iid_mean;
        assert!(
            variance > 1.5 * iid_variance,
            "correlated fraction variance {variance} is not above the iid {iid_variance}"
        );
    }

    #[test]
    fn a_zero_is_a_constant_field() {
        // a = 0 gives C_r = 1: one global Gaussian draw, not white noise.
        let spec = CorrelatedFieldSpec::new(&[32], 0.0).unwrap();
        let mut field = CorrelatedGaussianField::new(spec).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x4331_5a45_524f);
        field.sample(&mut rng);
        let first = field.values()[0].to_bits();
        assert!(
            field.values().iter().all(|value| value.to_bits() == first),
            "a = 0 must realize a spatially constant field"
        );
        // Thresholding a constant field is all-or-nothing.
        let graph = build_chain(32, false);
        let mut law = CorrelatedSiteLaw::new(&[32], 0.0, 0.37).unwrap();
        assert!(!law.is_iid_reduction());
        let mut configuration = StaticConfiguration::new(32, 31);
        law.sample(
            &graph,
            &mut configuration,
            &mut Xoshiro256PlusPlus::seed_from_u64(7),
        )
        .unwrap();
        let flags = flags(&graph, &configuration);
        assert!(
            flags.iter().all(|&flag| flag == flags[0]),
            "a = 0 occupation must be all-or-nothing"
        );
    }

    #[test]
    fn infinite_exponent_is_the_iid_delegate_too() {
        let graph = build_chain(32, false);
        let probability = 0.61;
        let mut correlated = StaticConfiguration::new(32, 31);
        let mut iid = StaticConfiguration::new(32, 31);
        let mut correlated_rng = Xoshiro256PlusPlus::seed_from_u64(99);
        let mut iid_rng = Xoshiro256PlusPlus::seed_from_u64(99);
        let mut law = CorrelatedSiteLaw::new(&[32], f64::INFINITY, probability).unwrap();
        assert!(law.is_iid_reduction());
        law.sample(&graph, &mut correlated, &mut correlated_rng)
            .unwrap();
        SiteBernoulli::new(Probability::new(probability).unwrap().into())
            .sample(&graph, &mut iid, &mut iid_rng)
            .unwrap();
        assert_eq!(flags(&graph, &correlated), flags(&graph, &iid));
    }

    #[test]
    fn field_draws_are_bitwise_deterministic() {
        let spec = CorrelatedFieldSpec::new(&[16, 16], 0.75).unwrap();
        let mut left = CorrelatedGaussianField::new(spec.clone()).unwrap();
        let mut right = CorrelatedGaussianField::new(spec).unwrap();
        let mut left_rng = Xoshiro256PlusPlus::seed_from_u64(2024);
        let mut right_rng = Xoshiro256PlusPlus::seed_from_u64(2024);
        for _ in 0..4 {
            left.sample(&mut left_rng);
            right.sample(&mut right_rng);
            assert!(left
                .values()
                .iter()
                .zip(right.values())
                .all(|(a, b)| a.to_bits() == b.to_bits()));
        }
        let mut other =
            CorrelatedGaussianField::new(CorrelatedFieldSpec::new(&[16, 16], 0.75).unwrap())
                .unwrap();
        other.sample(&mut Xoshiro256PlusPlus::seed_from_u64(2025));
        assert!(other
            .values()
            .iter()
            .zip(left.values())
            .any(|(a, b)| a.to_bits() != b.to_bits()));
    }

    #[test]
    fn rejected_samples_are_atomic() {
        let graph = build_square(8, 8, true);
        let mut law = CorrelatedSiteLaw::new(&[8, 8], 1.0, 0.5).unwrap();
        let mut configuration = StaticConfiguration::new(8, graph.edge_count());
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(3);
        let mut untouched = rng.clone();
        assert!(matches!(
            law.sample(&graph, &mut configuration, &mut rng),
            Err(CorrelatedError::Sampling(_))
        ));
        for _ in 0..8 {
            assert_eq!(rng.random::<u64>(), untouched.random::<u64>());
        }

        let wrong = build_chain(8, false);
        let mut configuration = StaticConfiguration::new(8, wrong.edge_count());
        assert!(matches!(
            law.sample(&wrong, &mut configuration, &mut rng),
            Err(CorrelatedError::GeometryMismatch { .. })
        ));
    }

    #[test]
    fn typed_rejections_cover_the_domain() {
        for (sides, exponent) in [
            (vec![], 1.0),
            (vec![1], 1.0),
            (vec![3], 1.0),
            (vec![4, 6], 1.0),
            (vec![4], -0.5),
            (vec![4], f64::NAN),
        ] {
            assert!(
                CorrelatedFieldSpec::new(&sides, exponent).is_err(),
                "accepted sides {sides:?} exponent {exponent}"
            );
        }
        assert!(matches!(
            CorrelatedFieldSpec::new(&[4, 4], -1.0),
            Err(CorrelatedError::InvalidExponent { .. })
        ));
        assert!(matches!(
            CorrelatedSiteLaw::new(&[8], 1.0, 1.5),
            Err(CorrelatedError::Probability(_))
        ));
        assert!(CorrelatedSiteLaw::new(&[8], f64::NEG_INFINITY, 0.5).is_err());
        assert_eq!(inverse_normal_cdf(0.0), f64::NEG_INFINITY);
        assert_eq!(inverse_normal_cdf(1.0), f64::INFINITY);
        assert!(inverse_normal_cdf(f64::NAN).is_nan());
        assert!(inverse_normal_cdf(0.5).abs() < 1e-12);
        // Acklam's approximation: |Phi(x) - p| <= 1.2e-9 is not re-checked
        // here without erf; the published bound is documented instead.
    }

    #[test]
    fn adapter_matches_direct_analysis_and_rejects_winding() {
        use super::super::{ObservablePlan, StaticObservable};
        let lattice = build_square(16, 16, true);
        let law = CorrelatedSiteLaw::new(&[16, 16], 2.0, 0.55).unwrap();
        let plan = ObservablePlan::new(
            [
                StaticObservable::LargestSize,
                StaticObservable::FiniteClusterSusceptibility,
            ],
            Vec::new(),
        );
        let mut adapter = CorrelatedPercolationMC::new(lattice.clone(), law, plan).unwrap();
        let mut ctx = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(77), 0, 1);
        adapter.sweep(&mut ctx);
        adapter.measure(&mut ctx);
        let values = ctx.finalize_measurements();
        let mut direct_workspace = ComponentWorkspace::new();
        let direct = super::super::analyze(
            &lattice,
            adapter.configuration(),
            adapter.configuration(),
            &[],
            &mut direct_workspace,
        )
        .unwrap();
        assert_eq!(
            values["LargestSize"].mean,
            direct.summary().largest_component_size as f64
        );

        let winding = ObservablePlan::new([StaticObservable::Wrapping], Vec::new());
        let law = CorrelatedSiteLaw::new(&[16, 16], 2.0, 0.55).unwrap();
        assert!(matches!(
            CorrelatedPercolationMC::new(lattice, law, winding),
            Err(CorrelatedError::WindingUnsupported)
        ));
    }
}
