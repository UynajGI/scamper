//! C2: shifted percolation thresholds under power-law correlated disorder
//! (roadmap v2 §5, stage C2).
//!
//! The stage measures `p_c(a)` on periodic `L x L` square lattices with the
//! C1 Fourier-filtering field and extrapolates the finite-size
//! pseudo-critical points to `L -> infinity`.
//!
//! ## Anchors
//!
//! - Pruessner et al., PRE 97, 062141 (2018), archived
//!   `docs/ref/md/Pruessner_2018_correlated_thresholds.md`, Table I
//!   (md lines ~596-624, 2D square, Euclidean metric) lists
//!   `a = inf: 0.592746` (the uncorrelated reference), `a = 3: 0.561406(4)`,
//!   `2.5: 0.556214(4)`, `2: 0.550143(5)`, `1.75: 0.546717(7)`,
//!   `1.5: 0.54299(1)`, `1.25: 0.53895(2)`, `1: 0.53452(4)`,
//!   `0.75: 0.5296(1)`, `0.5: 0.5239(3)`. Their fits use
//!   `|p_c - p_c^L| ~ L^{-1/\nu_a}` over `L >= 64` (md lines ~580-590) and
//!   they generated `10^5` replicas per size `L = 64 ... 8192`.
//! - The correlation target and generator are C1's
//!   `C(r) = (1 + r^2)^{-a/2}` discrete FFM (Pruessner Eq. 12 and
//!   Appendix A/B).
//! - Saberi, EPJ B 89, 253 (2016), archived
//!   `docs/ref/md/Saberi_2015_percolation_review.md`, md lines ~239-241:
//!   correlated clusters become more compact. The secondary observable here
//!   is the finite-size fractal dimension `d_f(a)` extracted from the
//!   largest-cluster size at the pseudo-critical point.
//!
//! **No four-digit claim is made.** The table values are used only as
//! tolerance-band centers: this stage runs a deliberately small ladder
//! (`L` up to 64/128 and a few hundred replicas per point) on a periodic
//! square lattice, so finite-size bias and the zero-cutoff distortion of
//! the FFM at small `L` are expected to exceed the quoted table
//! uncertainties. [`ThresholdLadder::pc_extrapolation`] reports the fit
//! residuals and the per-size bias so the comparison is explicit.
//!
//! ## Estimator (defined before measurement)
//!
//! One realization draws one C1 field, ranks all `N = L^2` sites by
//! descending field value, and activates them in that order on an
//! [`IncrementalUnionWorkspace`] (all physical edges pre-activated, i.e.
//! site percolation). Because `{field > z_p}` is exactly the first `n_p`
//! ranked sites, one ranking yields the whole `p` grid in one union-find
//! pass — the Newman-Ziff incremental trick applied to the *correlated*
//! ranking. For each grid point the reduction records two susceptibilities:
//!
//! ```text
//! chi_op(p)  = N Var(S_max / N) = Var(S_max) / N     (primary)
//! chi_fin(p) = (M2 - S_max^2) / (V_active - S_max)  (secondary, F3)
//! ```
//!
//! `chi_op` is the order-parameter susceptibility; it has the standard
//! finite-size form `chi_op = L^{gamma/nu} X((p - p_c) L^{1/nu})`, so its
//! peak defines the pseudo-critical point `p_c(L)`. `chi_fin` is the F3
//! finite-cluster susceptibility (mean cluster size with the largest
//! cluster excluded); it is recorded for context but is noisier and biased
//! in the C3 calibration, so it is not used for `nu`. The secondary
//! "large-cluster" fraction curve is monotone and is used only for its
//! crossing.
//!
//! The raw peak is refined by a parabola through the three grid points
//! around the maximum. That three-point estimator is adequate for C2's
//! wide windows but has a seed-to-seed scatter of several `1e-3` at
//! `L >= 512` (the peak top is broad and the variance estimator carries a
//! `sqrt(2/R)` common-mode fluctuation). C3 therefore provides two
//! add-ons, both additive to C2's draw stream:
//!
//! - [`ThresholdLadder::pool`] averages independent ladders and reports
//!   the between-ladder scatter of the peak, which the within-ladder
//!   bootstrap cannot see.
//! - [`ThresholdLadder::nu_from_pc_shift`] refines each peak with a local
//!   log-quadratic fit over [`PEAK_FIT_RADIUS`] and fits
//!   `p_c(inf) - p_c(L) = A L^{-1/nu}` with `p_c(inf)` fixed.
//!
//! The older [`ThresholdLadder::nu_from_log_derivative`] is retained for
//! continuity, but it is **not** a reliable `nu` estimator: the maximum of
//! the central difference of `ln chi_op` is selected on a noisy flank with
//! a one-grid-spacing difference window, so it is not scale covariant. On
//! the exact iid control it returns `nu` between 1.2 and 3.0 depending on
//! the window and smoothing (see the C3 validation record). Do not use it
//! to certify the Weinrib-Halperin branches.
//!
//! Statistical intervals from [`ThresholdLadder::run`] are 95% percentile
//! bootstrap intervals over realizations (200 resamples, fixed seeds
//! derived from the spec seed); they quantify sampling noise only. Pooled
//! ladders replace them with the between-ladder scatter. The extrapolation
//! additionally reports the reduced chi-square of the `L^{-1/nu}` fit,
//! which exposes the finite-size/systematic scatter that the bootstrap
//! cannot see. Callers must vary `minimum_side` to see the fit-window
//! sensitivity; this module never hides it behind a single number.
//!
//! ## Draw and reduction contract
//!
//! The field draws are C1's (bitwise reproducible per platform). `a = 0` is
//! rejected: the C1 `a = 0` field is constant, so there is no percolation
//! transition. `a = +inf` is the FFM white-noise limit and runs through the
//! same field path (exactly iid Bernoulli in law).

#![allow(clippy::neg_cmp_op_on_partial_ord)]

use core::fmt;

use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

use super::process::IncrementalUnionWorkspace;
use super::{inverse_normal_cdf, CorrelatedError, CorrelatedFieldSpec, CorrelatedGaussianField};
use crate::{build_square, EdgeId, GraphView, UndirectedGraphView, VertexId};

/// Number of bootstrap resamples behind every reported confidence interval.
pub const BOOTSTRAP_RESAMPLES: usize = 200;

/// Default half-window (in occupation probability) of the local
/// log-quadratic refinement of a susceptibility peak.
///
/// The finite-size peak is fitted with `ln chi = a0 + a1 x + a2 x^2`,
/// `x = p - p_max`, over the grid points within this distance of the
/// largest measured value. The radius is a physical probability scale, not
/// a number of grid points, so the fit follows the peak as it narrows with
/// increasing `L`; use a grid spacing well below it (C3 uses `2e-4`).
pub const PEAK_FIT_RADIUS: f64 = 0.004;

/// Two-sided 95% percentile bootstrap interval plus its point estimate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IntervalEstimate {
    pub value: f64,
    pub lower: f64,
    pub upper: f64,
    pub standard_error: f64,
}

/// Failure to configure or run a threshold ladder.
#[derive(Clone, Debug, PartialEq)]
pub enum CorrelatedThresholdError {
    /// The field geometry or exponent is invalid.
    Field(CorrelatedError),
    /// The exponent must be strictly positive (or `+inf`).
    Exponent { exponent: f64 },
    /// Fewer than two sides, or sides not strictly ascending.
    Sides { sides: Vec<usize> },
    /// The probability window is invalid.
    ProbabilityWindow {
        minimum: f64,
        maximum: f64,
        points: usize,
    },
    /// At least two samples are required (bootstrap intervals need more).
    Samples { samples: usize },
    /// A rung is missing for the requested fit window.
    FitWindow {
        minimum_side: usize,
        available: Vec<usize>,
    },
    /// The susceptibility peak fell on the window boundary.
    PeakOutsideWindow { side: usize },
    /// A crossing request named a side not present in the ladder.
    MissingSide { side: usize },
    /// A crossing request needs two distinct sides.
    SameSide { side: usize },
    /// The susceptibility curve had no finite value.
    UndefinedCurve { side: usize },
    /// The correlation-length exponent must be finite and positive.
    FitExponent { nu: f64 },
    /// The infinite-size threshold must be a finite probability.
    InfiniteThreshold { pc_infinity: f64 },
    /// The infinite-size threshold must lie strictly above every measured
    /// pseudo-critical point.
    ThresholdOrder {
        pc_infinity: f64,
        pseudo_critical: f64,
    },
    /// The peak-refinement half-window must be finite and positive.
    PeakRadius { radius: f64 },
    /// Ladders cannot be pooled unless exponent, grid and sizes agree.
    PoolMismatch,
    /// A local log-quadratic peak fit needs enough interior points.
    PeakFit { side: usize },
    /// A buffer could not be reserved.
    Capacity,
    /// The incremental union workspace rejected an activation.
    Incremental(String),
}

impl fmt::Display for CorrelatedThresholdError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Field(error) => error.fmt(formatter),
            Self::Exponent { exponent } => write!(
                formatter,
                "threshold ladders need a positive correlation exponent (or +inf), got {exponent}"
            ),
            Self::Sides { sides } => write!(
                formatter,
                "threshold ladder sides must be at least two strictly ascending values, got {sides:?}"
            ),
            Self::ProbabilityWindow {
                minimum,
                maximum,
                points,
            } => write!(
                formatter,
                "probability window [{minimum}, {maximum}] with {points} points is invalid \
                 (need 0 < min < max < 1 and at least five points)"
            ),
            Self::Samples { samples } => write!(
                formatter,
                "threshold ladder needs at least two samples per point, got {samples}"
            ),
            Self::FitWindow {
                minimum_side,
                available,
            } => write!(
                formatter,
                "no fit window: minimum side {minimum_side} leaves {available:?}"
            ),
            Self::PeakOutsideWindow { side } => write!(
                formatter,
                "susceptibility peak for L={side} lies on the probability-window boundary"
            ),
            Self::UndefinedCurve { side } => {
                write!(formatter, "susceptibility curve for L={side} is undefined")
            }
            Self::FitExponent { nu } => {
                write!(formatter, "fit exponent must be finite and positive, got {nu}")
            }
            Self::InfiniteThreshold { pc_infinity } => write!(
                formatter,
                "infinite-size threshold must be a finite probability, got {pc_infinity}"
            ),
            Self::ThresholdOrder {
                pc_infinity,
                pseudo_critical,
            } => write!(
                formatter,
                "infinite-size threshold {pc_infinity} must exceed the measured \
                 pseudo-critical point {pseudo_critical}"
            ),
            Self::PeakRadius { radius } => write!(
                formatter,
                "peak-refinement radius must be finite and positive, got {radius}"
            ),
            Self::PoolMismatch => write!(
                formatter,
                "pooled ladders must share exponent, probability grid and sizes"
            ),
            Self::PeakFit { side } => write!(
                formatter,
                "L={side} has too few positive susceptibility points for a local peak fit"
            ),
            Self::MissingSide { side } => write!(
                formatter,
                "side {side} is not present in this threshold ladder"
            ),
            Self::SameSide { side } => write!(
                formatter,
                "a crossing needs two distinct sides, got {side} twice"
            ),
            Self::Capacity => write!(formatter, "threshold ladder allocation failed"),
            Self::Incremental(message) => write!(formatter, "incremental union: {message}"),
        }
    }
}

impl std::error::Error for CorrelatedThresholdError {}

/// Ladder configuration.
///
/// `exponent` is the C1 correlation exponent `a`; `sides` are the linear
/// sizes `L`; the probability window is shared by every rung so bootstrap
/// curves and finite-size comparisons use one grid. Use
/// [`Self::new`] for the default window `[0.49, 0.60]` with 89 points
/// (`Delta p = 0.00125`, covering every Table I threshold from
/// `a = 0.5` to `a = inf`), or [`Self::with_probability_window`] to narrow
/// it.
#[derive(Clone, Debug, PartialEq)]
pub struct ThresholdLadderSpec {
    exponent: f64,
    sides: Vec<usize>,
    probability_min: f64,
    probability_max: f64,
    probability_points: usize,
    samples: usize,
    seed: u64,
}

impl ThresholdLadderSpec {
    /// Default probability window `[0.49, 0.60]`, 89 points.
    pub fn new(
        exponent: f64,
        sides: &[usize],
        samples: usize,
        seed: u64,
    ) -> Result<Self, CorrelatedThresholdError> {
        Self::with_probability_window(exponent, sides, 0.49, 0.60, 89, samples, seed)
    }

    /// Explicit probability window.
    pub fn with_probability_window(
        exponent: f64,
        sides: &[usize],
        probability_min: f64,
        probability_max: f64,
        probability_points: usize,
        samples: usize,
        seed: u64,
    ) -> Result<Self, CorrelatedThresholdError> {
        if !(exponent > 0.0) {
            // Rejects zero, negative, and NaN; +inf is the white-noise limit.
            return Err(CorrelatedThresholdError::Exponent { exponent });
        }
        CorrelatedFieldSpec::new(sides, exponent).map_err(CorrelatedThresholdError::Field)?;
        if sides.len() < 2 || sides.windows(2).any(|pair| pair[0] >= pair[1]) {
            return Err(CorrelatedThresholdError::Sides {
                sides: sides.to_vec(),
            });
        }
        if !(probability_min > 0.0
            && probability_min < probability_max
            && probability_max < 1.0
            && probability_min.is_finite()
            && probability_max.is_finite()
            && probability_points >= 5)
        {
            return Err(CorrelatedThresholdError::ProbabilityWindow {
                minimum: probability_min,
                maximum: probability_max,
                points: probability_points,
            });
        }
        if samples < 2 {
            return Err(CorrelatedThresholdError::Samples { samples });
        }
        Ok(Self {
            exponent,
            sides: sides.to_vec(),
            probability_min,
            probability_max,
            probability_points,
            samples,
            seed,
        })
    }

    /// The correlation exponent `a`.
    #[inline]
    pub const fn exponent(&self) -> f64 {
        self.exponent
    }

    /// The requested linear sizes.
    #[inline]
    pub fn sides(&self) -> &[usize] {
        &self.sides
    }

    /// Replicas per probability point.
    #[inline]
    pub const fn samples(&self) -> usize {
        self.samples
    }

    /// Base seed (bootstrap streams are derived from it deterministically).
    #[inline]
    pub const fn seed(&self) -> u64 {
        self.seed
    }

    /// The shared probability grid.
    pub fn probabilities(&self) -> Vec<f64> {
        probability_grid(
            self.probability_min,
            self.probability_max,
            self.probability_points,
        )
    }
}

fn probability_grid(minimum: f64, maximum: f64, points: usize) -> Vec<f64> {
    if points < 2 {
        return vec![minimum];
    }
    let step = (maximum - minimum) / (points - 1) as f64;
    (0..points)
        .map(|index| minimum + step * index as f64)
        .collect()
}

/// Per-size measurements of one ladder rung.
#[derive(Clone, Debug)]
pub struct RungObservables {
    /// Linear size `L`.
    pub side: usize,
    /// Susceptibility peak location (the primary `p_c(L)` estimator).
    pub pseudo_critical: IntervalEstimate,
    /// Peak value of the ensemble-averaged finite-cluster susceptibility.
    pub susceptibility_peak: f64,
    /// Ensemble-averaged largest-cluster fraction at the peak location
    /// (linearly interpolated between grid points).
    pub largest_fraction_at_peak: f64,
    /// Ensemble-averaged order-parameter susceptibility `N Var(S_max/N)`
    /// on the shared grid (the primary estimator).
    pub susceptibility: Vec<f64>,
    /// Ensemble-averaged F3 finite-cluster susceptibility
    /// `(M2 - S_max^2)/(V_active - S_max)` on the shared grid (secondary;
    /// the definition is reused from the F3 observable).
    pub finite_cluster_susceptibility: Vec<f64>,
    /// Ensemble-averaged `S_max/N` on the shared grid.
    pub largest_fraction: Vec<f64>,
    /// Fraction of samples with a defined finite-cluster susceptibility.
    pub defined_fraction: Vec<f64>,
}

/// Linear fit `y = intercept + slope x`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinearFit {
    pub intercept: f64,
    pub slope: f64,
    pub intercept_standard_error: f64,
    pub slope_standard_error: f64,
    pub reduced_chi_square: f64,
    pub degrees_of_freedom: usize,
}

/// Extrapolation `p_c(L) = p_c(inf) + A L^{-1/nu}`.
#[derive(Clone, Debug, PartialEq)]
pub struct PcExtrapolation {
    /// The intercept `p_c(inf)` with its 95% interval.
    pub pc_infinity: IntervalEstimate,
    /// The fitted amplitude `A`.
    pub amplitude: f64,
    /// The exponent used in `L^{-1/nu}`.
    pub nu: f64,
    /// Sizes entering the fit.
    pub used_sides: Vec<usize>,
    /// The raw finite-size bias `p_c(L) - p_c(inf)` per rung.
    pub finite_size_bias: Vec<(usize, f64)>,
    /// Reduced chi-square of the weighted fit.
    pub reduced_chi_square: f64,
    pub degrees_of_freedom: usize,
}

/// One rung's log-derivative slope with its bootstrap interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RungSlope {
    pub side: usize,
    pub slope: IntervalEstimate,
}

/// Correlation-length exponent from the finite-size scaling of the
/// susceptibility log-derivative.
#[derive(Clone, Debug, PartialEq)]
pub struct NuEstimate {
    /// `nu` with its 95% interval.
    pub nu: IntervalEstimate,
    /// The fitted `1/nu` (the slope of `ln(max d ln chi/dp)` vs `ln L`).
    pub inverse_nu: f64,
    pub inverse_nu_standard_error: f64,
    /// Per-rung log-derivative values (the fitted data).
    pub rung_slopes: Vec<RungSlope>,
    pub used_sides: Vec<usize>,
    pub reduced_chi_square: f64,
    pub degrees_of_freedom: usize,
}

/// One rung of a [`PcShiftNuEstimate`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PcShiftRung {
    /// Linear size `L`.
    pub side: usize,
    /// Refined pseudo-critical point `p_c(L)`.
    pub pseudo_critical: f64,
    /// Refined susceptibility peak height at `p_c(L)`.
    pub susceptibility_peak: f64,
    /// `p_c(infinity) - p_c(L)`, the fitted finite-size bias.
    pub bias: f64,
}

/// Correlation-length exponent from the finite-size shift of the
/// pseudo-critical point with `p_c(infinity)` held fixed.
#[derive(Clone, Debug, PartialEq)]
pub struct PcShiftNuEstimate {
    /// `nu` with its 95% interval.
    pub nu: IntervalEstimate,
    /// The fitted `1/nu`.
    pub inverse_nu: f64,
    pub inverse_nu_standard_error: f64,
    /// The fixed infinite-size threshold.
    pub pc_infinity: f64,
    /// The fitted amplitude `A` in `p_c(inf) - p_c(L) = A L^{-1/nu}`.
    pub amplitude: f64,
    /// Per-size refined peaks (the fitted data).
    pub rungs: Vec<PcShiftRung>,
    pub used_sides: Vec<usize>,
    pub reduced_chi_square: f64,
    pub degrees_of_freedom: usize,
}

/// Fractal dimension of the largest cluster at the pseudo-critical point.
#[derive(Clone, Debug, PartialEq)]
pub struct FractalDimensionEstimate {
    pub dimension: IntervalEstimate,
    pub used_sides: Vec<usize>,
    pub reduced_chi_square: f64,
    pub degrees_of_freedom: usize,
}

#[derive(Clone, Debug)]
struct RungSamples {
    slope_value: f64,
    slope_bootstrap: Vec<f64>,
}

/// One measured threshold ladder (`p_c(a,L)` curves plus their analyses).
#[derive(Clone, Debug)]
pub struct ThresholdLadder {
    exponent: f64,
    probabilities: Vec<f64>,
    rungs: Vec<RungObservables>,
    samples: Vec<RungSamples>,
    replicas: usize,
}

impl ThresholdLadder {
    /// Run the ladder: for every side, `replicas` independent fields on the
    /// shared probability grid, one incremental union-find pass each.
    pub fn run(spec: &ThresholdLadderSpec) -> Result<Self, CorrelatedThresholdError> {
        let probabilities = spec.probabilities();
        let thresholds: Vec<f64> = probabilities
            .iter()
            .map(|&p| inverse_normal_cdf(1.0 - p))
            .collect();
        let mut rungs = Vec::with_capacity(spec.sides.len());
        let mut samples = Vec::with_capacity(spec.sides.len());
        for (rung_index, &side) in spec.sides.iter().enumerate() {
            let (rung, rung_samples) =
                run_rung(spec, side, &probabilities, &thresholds, rung_index)?;
            rungs.push(rung);
            samples.push(rung_samples);
        }
        Ok(Self {
            exponent: spec.exponent,
            probabilities,
            rungs,
            samples,
            replicas: spec.samples,
        })
    }

    /// Pool independent ladders that share one spec into a single ladder.
    ///
    /// The per-`p` curves of the input ladders are averaged (the finite-
    /// cluster susceptibility is weighted by its defined-sample fraction),
    /// which is exactly the ensemble average of the pooled replicas. The
    /// rung's [`RungObservables::pseudo_critical`] is refined with
    /// [`PEAK_FIT_RADIUS`] on both the pooled curve and each input ladder;
    /// its standard error is the scatter of those per-ladder refined peaks.
    /// That between-ladder error includes the common-mode fluctuation that
    /// the within-ladder bootstrap cannot see. The pooled
    /// [`RungSamples::slope_bootstrap`] holds the per-ladder log-derivative
    /// slopes, so the existing estimators use the between-ladder scatter.
    ///
    /// Pooling is the supported way to reach the replica counts the C3
    /// `nu` fit needs without changing [`Self::run`]'s draw stream (the C2
    /// gate depends on it).
    pub fn pool(ladders: &[ThresholdLadder]) -> Result<Self, CorrelatedThresholdError> {
        if ladders.len() < 2 {
            return Err(CorrelatedThresholdError::Samples {
                samples: ladders.len(),
            });
        }
        let first = &ladders[0];
        for ladder in ladders {
            let same_sides = ladder.rungs.len() == first.rungs.len()
                && ladder
                    .rungs
                    .iter()
                    .zip(&first.rungs)
                    .all(|(left, right)| left.side == right.side);
            if ladder.exponent != first.exponent
                || ladder.probabilities != first.probabilities
                || !same_sides
                || ladder.replicas != first.replicas
            {
                return Err(CorrelatedThresholdError::PoolMismatch);
            }
        }
        let points = first.probabilities.len();
        let count = ladders.len() as f64;
        let mut rungs = Vec::with_capacity(first.rungs.len());
        let mut samples = Vec::with_capacity(first.rungs.len());
        for (rung_index, reference) in first.rungs.iter().enumerate() {
            let side = reference.side;
            let mut susceptibility = vec![0.0_f64; points];
            let mut finite_susceptibility = vec![0.0_f64; points];
            let mut finite_weight = vec![0.0_f64; points];
            let mut largest_fraction = vec![0.0_f64; points];
            let mut defined_fraction = vec![0.0_f64; points];
            let mut peaks = Vec::with_capacity(ladders.len());
            let mut slopes = Vec::with_capacity(ladders.len());
            for ladder in ladders {
                let rung = &ladder.rungs[rung_index];
                for point in 0..points {
                    susceptibility[point] += rung.susceptibility[point];
                    largest_fraction[point] += rung.largest_fraction[point];
                    defined_fraction[point] += rung.defined_fraction[point];
                    let weight = rung.defined_fraction[point];
                    let finite = rung.finite_cluster_susceptibility[point];
                    if weight > 0.0 && finite.is_finite() {
                        finite_susceptibility[point] += weight * finite;
                        finite_weight[point] += weight;
                    }
                }
                peaks.push(
                    log_quadratic_peak(
                        &first.probabilities,
                        &rung.susceptibility,
                        side,
                        PEAK_FIT_RADIUS,
                    )
                    .map_or(rung.pseudo_critical.value, |(location, _)| location),
                );
                let slope = ladder.samples[rung_index].slope_value;
                if slope.is_finite() && slope > 0.0 {
                    slopes.push(slope);
                }
            }
            for point in 0..points {
                susceptibility[point] /= count;
                largest_fraction[point] /= count;
                defined_fraction[point] /= count;
                finite_susceptibility[point] = if finite_weight[point] > 0.0 {
                    finite_susceptibility[point] / finite_weight[point]
                } else {
                    f64::NAN
                };
            }
            let (location, height) =
                log_quadratic_peak(&first.probabilities, &susceptibility, side, PEAK_FIT_RADIUS)?;
            let sorted_peaks = sorted(&peaks);
            let pseudo_critical = IntervalEstimate {
                value: location,
                lower: percentile(&sorted_peaks, 0.025).min(location),
                upper: percentile(&sorted_peaks, 0.975).max(location),
                standard_error: standard_error(&peaks),
            };
            let largest_fraction_at_peak =
                interpolate(&first.probabilities, &largest_fraction, location);
            let slope_value = maximum_log_derivative(&first.probabilities, &susceptibility);
            if !(slope_value > 0.0) || slopes.len() < 2 {
                return Err(CorrelatedThresholdError::UndefinedCurve { side });
            }
            rungs.push(RungObservables {
                side,
                pseudo_critical,
                susceptibility_peak: height,
                largest_fraction_at_peak,
                susceptibility,
                finite_cluster_susceptibility: finite_susceptibility,
                largest_fraction,
                defined_fraction,
            });
            samples.push(RungSamples {
                slope_value,
                slope_bootstrap: slopes,
            });
        }
        Ok(Self {
            exponent: first.exponent,
            probabilities: first.probabilities.clone(),
            rungs,
            samples,
            replicas: first.replicas.saturating_mul(ladders.len()),
        })
    }

    /// The correlation exponent `a`.
    #[inline]
    pub const fn exponent(&self) -> f64 {
        self.exponent
    }

    /// The shared probability grid.
    #[inline]
    pub fn probabilities(&self) -> &[f64] {
        &self.probabilities
    }

    /// Per-size measurements.
    #[inline]
    pub fn rungs(&self) -> &[RungObservables] {
        &self.rungs
    }

    /// Linear-interpolated crossing of the largest-cluster-fraction curves
    /// of two sizes on the shared probability grid.
    ///
    /// This is the standard finite-size crossing estimator for `p_c`; unlike
    /// the susceptibility-peak location it is calibrated against a white-noise
    /// ladder (which reproduces the exact iid threshold). `None` means the
    /// curves do not cross inside the measured window.
    pub fn largest_fraction_crossing(
        &self,
        left_side: usize,
        right_side: usize,
    ) -> Result<Option<f64>, CorrelatedThresholdError> {
        if left_side == right_side {
            return Err(CorrelatedThresholdError::SameSide { side: left_side });
        }
        let left = self
            .rungs
            .iter()
            .find(|rung| rung.side == left_side)
            .ok_or(CorrelatedThresholdError::MissingSide { side: left_side })?;
        let right = self
            .rungs
            .iter()
            .find(|rung| rung.side == right_side)
            .ok_or(CorrelatedThresholdError::MissingSide { side: right_side })?;
        for index in 1..self.probabilities.len() {
            let (a, b) = (index - 1, index);
            let left_a = left.largest_fraction[a] - right.largest_fraction[a];
            let left_b = left.largest_fraction[b] - right.largest_fraction[b];
            if left_a == 0.0 {
                return Ok(Some(self.probabilities[a]));
            }
            if left_a.signum() != left_b.signum() {
                let weight = left_a / (left_a - left_b);
                return Ok(Some(
                    self.probabilities[a]
                        + weight * (self.probabilities[b] - self.probabilities[a]),
                ));
            }
        }
        Ok(None)
    }

    /// Replicas per probability point.
    #[inline]
    pub const fn replicas(&self) -> usize {
        self.replicas
    }

    /// Extrapolate `p_c(L) = p_c(inf) + A L^{-1/nu}` over `side >= minimum_side`.
    ///
    /// Weights are the inverse bootstrap variances of the per-size peak
    /// locations; the reported interval is the intercept standard error
    /// inflated by the reduced chi-square when the fit is poor
    /// (`sqrt(max(1, chi2/dof))`, the usual conservative practice). The
    /// per-size bias is returned so the extrapolation can be audited.
    pub fn pc_extrapolation(
        &self,
        minimum_side: usize,
        nu: f64,
    ) -> Result<PcExtrapolation, CorrelatedThresholdError> {
        if !(nu.is_finite() && nu > 0.0) {
            return Err(CorrelatedThresholdError::FitExponent { nu });
        }
        let used: Vec<usize> = self
            .rungs
            .iter()
            .enumerate()
            .filter(|(_, rung)| rung.side >= minimum_side)
            .map(|(index, _)| index)
            .collect();
        if used.len() < 3 {
            return Err(CorrelatedThresholdError::FitWindow {
                minimum_side,
                available: self.rungs.iter().map(|rung| rung.side).collect(),
            });
        }
        let x: Vec<f64> = used
            .iter()
            .map(|&index| (self.rungs[index].side as f64).powf(-1.0 / nu))
            .collect();
        let y: Vec<f64> = used
            .iter()
            .map(|&index| self.rungs[index].pseudo_critical.value)
            .collect();
        let sigma: Vec<f64> = used
            .iter()
            .map(|&index| self.rungs[index].pseudo_critical.standard_error.max(1e-6))
            .collect();
        let fit = weighted_line(&x, &y, &sigma).expect("three points define a line");
        let scale = fit.reduced_chi_square.max(1.0).sqrt();
        let pc_infinity = IntervalEstimate {
            value: fit.intercept,
            lower: fit.intercept - 1.96 * fit.intercept_standard_error * scale,
            upper: fit.intercept + 1.96 * fit.intercept_standard_error * scale,
            standard_error: fit.intercept_standard_error * scale,
        };
        let finite_size_bias = used
            .iter()
            .map(|&index| {
                (
                    self.rungs[index].side,
                    self.rungs[index].pseudo_critical.value - fit.intercept,
                )
            })
            .collect();
        Ok(PcExtrapolation {
            pc_infinity,
            amplitude: fit.slope,
            nu,
            used_sides: used.iter().map(|&index| self.rungs[index].side).collect(),
            finite_size_bias,
            reduced_chi_square: fit.reduced_chi_square,
            degrees_of_freedom: fit.degrees_of_freedom,
        })
    }

    /// Estimate `nu` from `max_p d ln chi / dp ~ L^{1/nu}`.
    ///
    /// The derivative is a central difference of the ensemble-averaged
    /// `ln chi(p)` on the shared grid; its maximum is refined by a parabola
    /// through the three grid points around the largest value. Bootstrap
    /// resampling propagates sampling noise into each rung's slope, and a
    /// weighted line through `ln(slope)` vs `ln(L)` gives `1/nu`.
    pub fn nu_from_log_derivative(
        &self,
        minimum_side: usize,
    ) -> Result<NuEstimate, CorrelatedThresholdError> {
        let used: Vec<usize> = self
            .rungs
            .iter()
            .enumerate()
            .filter(|(_, rung)| rung.side >= minimum_side)
            .map(|(index, _)| index)
            .collect();
        if used.len() < 3 {
            return Err(CorrelatedThresholdError::FitWindow {
                minimum_side,
                available: self.rungs.iter().map(|rung| rung.side).collect(),
            });
        }
        let mut x = Vec::with_capacity(used.len());
        let mut y = Vec::with_capacity(used.len());
        let mut sigma = Vec::with_capacity(used.len());
        let mut rung_slopes = Vec::with_capacity(used.len());
        for &index in &used {
            let direct = self.samples[index].slope_value;
            let values: Vec<f64> = self.samples[index]
                .slope_bootstrap
                .iter()
                .copied()
                .filter(|value| value.is_finite() && *value > 0.0)
                .collect();
            if !(direct > 0.0) || values.len() < 2 {
                return Err(CorrelatedThresholdError::UndefinedCurve {
                    side: self.rungs[index].side,
                });
            }
            let error = standard_error(&values).max(1e-9);
            let side = self.rungs[index].side;
            let sorted_values = sorted(&values);
            rung_slopes.push(RungSlope {
                side,
                slope: IntervalEstimate {
                    value: direct,
                    lower: percentile(&sorted_values, 0.025),
                    upper: percentile(&sorted_values, 0.975),
                    standard_error: error,
                },
            });
            x.push((side as f64).ln());
            y.push(direct.ln());
            sigma.push(error / direct);
        }
        let fit = weighted_line(&x, &y, &sigma).expect("three points define a line");
        let inverse_nu = fit.slope;
        if !(inverse_nu > 0.0) {
            return Err(CorrelatedThresholdError::FitExponent { nu: inverse_nu });
        }
        let nu = 1.0 / inverse_nu;
        let nu_error = fit.slope_standard_error / (inverse_nu * inverse_nu);
        let scale = fit.reduced_chi_square.max(1.0).sqrt();
        Ok(NuEstimate {
            nu: IntervalEstimate {
                value: nu,
                lower: nu - 1.96 * nu_error * scale,
                upper: nu + 1.96 * nu_error * scale,
                standard_error: nu_error * scale,
            },
            inverse_nu,
            inverse_nu_standard_error: fit.slope_standard_error,
            rung_slopes,
            used_sides: used.iter().map(|&index| self.rungs[index].side).collect(),
            reduced_chi_square: fit.reduced_chi_square,
            degrees_of_freedom: fit.degrees_of_freedom,
        })
    }

    /// Estimate `nu` from `p_c(infinity) - p_c(L) = A L^{-1/nu}` with the
    /// infinite-size threshold held fixed.
    ///
    /// Each rung's pseudo-critical point is first refined by a local
    /// log-quadratic fit of the order-parameter susceptibility over
    /// `peak_radius` (use [`PEAK_FIT_RADIUS`] unless auditing the
    /// sensitivity). The fit is linear in `ln(p_c(inf) - p_c(L))` versus
    /// `ln L`; weights are the rung peak standard errors (for a pooled
    /// ladder these are the between-ladder scatter of the refined peaks).
    /// This is the estimator calibrated in the C3 validation: on an
    /// `a = +infinity` control ladder it recovers `nu = 4/3` within its
    /// interval, unlike the log-derivative estimator.
    ///
    /// `p_c(infinity)` must be supplied by the caller because a three- or
    /// four-point free fit of `(p_c(inf), A, nu)` is not stable at C3's
    /// sizes (the free fit biases `nu` high by ~10-20% even on the exact
    /// iid control). For `a = +infinity` use the exact iid site threshold
    /// `0.592746`; for the correlated rows use the table value under test.
    pub fn nu_from_pc_shift(
        &self,
        minimum_side: usize,
        pc_infinity: f64,
        peak_radius: f64,
    ) -> Result<PcShiftNuEstimate, CorrelatedThresholdError> {
        if !(pc_infinity.is_finite() && pc_infinity > 0.0 && pc_infinity < 1.0) {
            return Err(CorrelatedThresholdError::InfiniteThreshold { pc_infinity });
        }
        if !(peak_radius.is_finite() && peak_radius > 0.0) {
            return Err(CorrelatedThresholdError::PeakRadius {
                radius: peak_radius,
            });
        }
        let used: Vec<usize> = self
            .rungs
            .iter()
            .enumerate()
            .filter(|(_, rung)| rung.side >= minimum_side)
            .map(|(index, _)| index)
            .collect();
        if used.len() < 3 {
            return Err(CorrelatedThresholdError::FitWindow {
                minimum_side,
                available: self.rungs.iter().map(|rung| rung.side).collect(),
            });
        }
        let mut x = Vec::with_capacity(used.len());
        let mut y = Vec::with_capacity(used.len());
        let mut sigma = Vec::with_capacity(used.len());
        let mut rungs = Vec::with_capacity(used.len());
        for &index in &used {
            let rung = &self.rungs[index];
            let (location, height) = log_quadratic_peak(
                &self.probabilities,
                &rung.susceptibility,
                rung.side,
                peak_radius,
            )?;
            let bias = pc_infinity - location;
            if !(bias > 0.0) {
                return Err(CorrelatedThresholdError::ThresholdOrder {
                    pc_infinity,
                    pseudo_critical: location,
                });
            }
            let error = rung.pseudo_critical.standard_error.max(1e-12);
            x.push((rung.side as f64).ln());
            y.push(bias.ln());
            sigma.push((error / bias).max(1e-12));
            rungs.push(PcShiftRung {
                side: rung.side,
                pseudo_critical: location,
                susceptibility_peak: height,
                bias,
            });
        }
        let fit = weighted_line(&x, &y, &sigma).expect("three points define a line");
        let inverse_nu = -fit.slope;
        if !(inverse_nu > 0.0) {
            return Err(CorrelatedThresholdError::FitExponent { nu: inverse_nu });
        }
        let nu = 1.0 / inverse_nu;
        let nu_error = fit.slope_standard_error / (inverse_nu * inverse_nu);
        let scale = fit.reduced_chi_square.max(1.0).sqrt();
        Ok(PcShiftNuEstimate {
            nu: IntervalEstimate {
                value: nu,
                lower: nu - 1.96 * nu_error * scale,
                upper: nu + 1.96 * nu_error * scale,
                standard_error: nu_error * scale,
            },
            inverse_nu,
            inverse_nu_standard_error: fit.slope_standard_error,
            pc_infinity,
            amplitude: fit.intercept.exp(),
            rungs,
            used_sides: used.iter().map(|&index| self.rungs[index].side).collect(),
            reduced_chi_square: fit.reduced_chi_square,
            degrees_of_freedom: fit.degrees_of_freedom,
        })
    }

    /// Fractal dimension `S_max(p_c(L)) ~ L^{d_f}` at the pseudo-critical
    /// points (secondary compactness observable).
    pub fn fractal_dimension(
        &self,
        minimum_side: usize,
    ) -> Result<FractalDimensionEstimate, CorrelatedThresholdError> {
        let used: Vec<usize> = self
            .rungs
            .iter()
            .enumerate()
            .filter(|(_, rung)| rung.side >= minimum_side)
            .map(|(index, _)| index)
            .collect();
        if used.len() < 3 {
            return Err(CorrelatedThresholdError::FitWindow {
                minimum_side,
                available: self.rungs.iter().map(|rung| rung.side).collect(),
            });
        }
        let x: Vec<f64> = used
            .iter()
            .map(|&index| (self.rungs[index].side as f64).ln())
            .collect();
        let y: Vec<f64> = used
            .iter()
            .map(|&index| {
                let rung = &self.rungs[index];
                let size = rung.largest_fraction_at_peak * (rung.side * rung.side) as f64;
                size.max(1.0).ln()
            })
            .collect();
        let sigma = vec![1.0; used.len()];
        let fit = weighted_line(&x, &y, &sigma).expect("three points define a line");
        let scale = fit.reduced_chi_square.max(1.0).sqrt();
        Ok(FractalDimensionEstimate {
            dimension: IntervalEstimate {
                value: fit.slope,
                lower: fit.slope - 1.96 * fit.slope_standard_error * scale,
                upper: fit.slope + 1.96 * fit.slope_standard_error * scale,
                standard_error: fit.slope_standard_error * scale,
            },
            used_sides: used.iter().map(|&index| self.rungs[index].side).collect(),
            reduced_chi_square: fit.reduced_chi_square,
            degrees_of_freedom: fit.degrees_of_freedom,
        })
    }
}

/// One rung: sample `replicas` rankings and reduce them to curves plus
/// bootstrap distributions.
fn run_rung(
    spec: &ThresholdLadderSpec,
    side: usize,
    probabilities: &[f64],
    thresholds: &[f64],
    rung_index: usize,
) -> Result<(RungObservables, RungSamples), CorrelatedThresholdError> {
    let lattice = build_square(side, side, true);
    let vertices = lattice.vertex_count();
    let edges = lattice.edge_count();
    let field_spec = CorrelatedFieldSpec::new(&[side, side], spec.exponent)
        .map_err(CorrelatedThresholdError::Field)?;
    let mut field =
        CorrelatedGaussianField::new(field_spec).map_err(CorrelatedThresholdError::Field)?;
    let mut workspace = IncrementalUnionWorkspace::new();
    workspace
        .prepare(vertices, edges)
        .map_err(|_| CorrelatedThresholdError::Capacity)?;

    let points = probabilities.len();
    let mut order: Vec<usize> = Vec::new();
    order
        .try_reserve_exact(vertices)
        .map_err(|_| CorrelatedThresholdError::Capacity)?;
    order.extend(0..vertices);
    let mut counts = vec![0_usize; points];
    let mut chi = vec![f64::NAN; spec.samples * points];
    let mut largest = vec![0.0_f64; spec.samples * points];
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(
        spec.seed ^ 0x9e37_79b9_7f4a_7c15_u64.wrapping_mul(rung_index as u64 + 1),
    );

    for replica in 0..spec.samples {
        field.sample(&mut rng);
        let values = field.values();
        order.sort_unstable_by(|&left, &right| {
            values[right]
                .total_cmp(&values[left])
                .then(left.cmp(&right))
        });

        // Occupation counts n_p = #{field > z_p}; thresholds descend as p
        // ascends, so one advancing pointer computes every count.
        let mut count = 0;
        for (point, &threshold) in thresholds.iter().enumerate() {
            while count < vertices && values[order[count]] > threshold {
                count += 1;
            }
            counts[point] = count;
        }

        workspace.reset();
        for edge in 0..edges {
            workspace
                .activate_edge(&lattice, EdgeId::from_valid_index(edge), 0)
                .map_err(|error| CorrelatedThresholdError::Incremental(error.to_string()))?;
        }
        let mut active = 0;
        for point in 0..points {
            while active < counts[point] {
                workspace
                    .activate_vertex(&lattice, VertexId::from_valid_index(order[active]), 0)
                    .map_err(|error| CorrelatedThresholdError::Incremental(error.to_string()))?;
                active += 1;
            }
            let summary = workspace.summary();
            let largest_size = summary.largest_component_size;
            largest[replica * points + point] = if vertices == 0 {
                0.0
            } else {
                largest_size as f64 / vertices as f64
            };
            let finite = summary.active_vertex_count.saturating_sub(largest_size);
            if finite != 0 {
                let largest_squared = (largest_size as u128) * (largest_size as u128);
                let finite_moment = summary.raw_second_moment.saturating_sub(largest_squared);
                chi[replica * points + point] = finite_moment as f64 / finite as f64;
            }
        }
    }

    let mut susceptibility = vec![f64::NAN; points];
    let mut finite_susceptibility = vec![f64::NAN; points];
    let mut largest_fraction = vec![0.0; points];
    let mut defined_fraction = vec![0.0; points];
    for point in 0..points {
        let mut defined = 0;
        let mut finite_sum = 0.0;
        let mut largest_sum = 0.0;
        let mut largest_square = 0.0;
        for replica in 0..spec.samples {
            let value = chi[replica * points + point];
            if value.is_finite() {
                defined += 1;
                finite_sum += value;
            }
            let fraction = largest[replica * points + point];
            largest_sum += fraction;
            largest_square += fraction * fraction;
        }
        defined_fraction[point] = defined as f64 / spec.samples as f64;
        if defined != 0 {
            finite_susceptibility[point] = finite_sum / defined as f64;
        }
        let mean_fraction = largest_sum / spec.samples as f64;
        let variance =
            (largest_square / spec.samples as f64 - mean_fraction * mean_fraction).max(0.0);
        largest_fraction[point] = mean_fraction;
        // Order-parameter susceptibility N Var(S_max/N) = Var(S_max)/N.
        susceptibility[point] = vertices as f64 * variance;
    }

    let pseudo_critical = peak_estimate(probabilities, &susceptibility, side)?;
    let largest_fraction_at_peak =
        interpolate(probabilities, &largest_fraction, pseudo_critical.value);
    let susceptibility_peak = interpolate(probabilities, &susceptibility, pseudo_critical.value);

    // Bootstrap: resample replicas, recompute the variance susceptibility
    // curve, its peak, and its maximum log-derivative.
    let mut bootstrap_pc = Vec::with_capacity(BOOTSTRAP_RESAMPLES);
    let mut bootstrap_slope = Vec::with_capacity(BOOTSTRAP_RESAMPLES);
    let mut bootstrap_rng = Xoshiro256PlusPlus::seed_from_u64(
        spec.seed ^ 0xd1b5_4a32_d192_ed03_u64.wrapping_mul(rung_index as u64 + 1),
    );
    let mut indices = vec![0_usize; spec.samples];
    let mut curvature = vec![f64::NAN; points];
    let mut fractions = vec![0.0_f64; points];
    for _ in 0..BOOTSTRAP_RESAMPLES {
        for index in &mut indices {
            *index = bootstrap_rng.random_range(0..spec.samples);
        }
        fractions.fill(0.0);
        for &replica in &indices {
            let base = replica * points;
            for point in 0..points {
                fractions[point] += largest[base + point];
            }
        }
        let mut sums = vec![0.0_f64; points];
        for &replica in &indices {
            let base = replica * points;
            for point in 0..points {
                let value = largest[base + point];
                sums[point] += value * value;
            }
        }
        for point in 0..points {
            let mean = fractions[point] / spec.samples as f64;
            let variance = (sums[point] / spec.samples as f64 - mean * mean).max(0.0);
            curvature[point] = vertices as f64 * variance;
        }
        match peak_estimate(probabilities, &curvature, side) {
            Ok(estimate) => bootstrap_pc.push(estimate.value),
            Err(_) => bootstrap_pc.push(f64::NAN),
        }
        bootstrap_slope.push(maximum_log_derivative(probabilities, &curvature));
    }

    let finite_pc: Vec<f64> = bootstrap_pc
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .collect();
    if finite_pc.len() < 2 {
        return Err(CorrelatedThresholdError::UndefinedCurve { side });
    }
    let sorted_pc = sorted(&finite_pc);
    let pseudo_critical = IntervalEstimate {
        value: pseudo_critical.value,
        lower: percentile(&sorted_pc, 0.025),
        upper: percentile(&sorted_pc, 0.975),
        standard_error: standard_deviation(&finite_pc),
    };
    let slope_values: Vec<f64> = bootstrap_slope
        .iter()
        .copied()
        .filter(|value| value.is_finite() && *value > 0.0)
        .collect();
    if slope_values.len() < 2 {
        return Err(CorrelatedThresholdError::UndefinedCurve { side });
    }
    let slope_value = maximum_log_derivative(probabilities, &susceptibility);
    if !(slope_value > 0.0) {
        return Err(CorrelatedThresholdError::UndefinedCurve { side });
    }

    Ok((
        RungObservables {
            side,
            pseudo_critical,
            susceptibility_peak,
            largest_fraction_at_peak,
            susceptibility,
            finite_cluster_susceptibility: finite_susceptibility,
            largest_fraction,
            defined_fraction,
        },
        RungSamples {
            slope_value,
            slope_bootstrap: slope_values,
        },
    ))
}

/// Parabolic peak of `curve` on an equally spaced grid; the peak must be
/// interior.
fn peak_estimate(
    probabilities: &[f64],
    curve: &[f64],
    side: usize,
) -> Result<IntervalEstimate, CorrelatedThresholdError> {
    let mut best_index = usize::MAX;
    let mut best_value = f64::NEG_INFINITY;
    for (index, &value) in curve.iter().enumerate() {
        if value.is_finite() && value > best_value {
            best_index = index;
            best_value = value;
        }
    }
    if best_index == usize::MAX {
        return Err(CorrelatedThresholdError::UndefinedCurve { side });
    }
    let (index, value) = (best_index, best_value);
    if index == 0 || index + 1 >= curve.len() {
        return Err(CorrelatedThresholdError::PeakOutsideWindow { side });
    }
    let (left, center, right) = (curve[index - 1], value, curve[index + 1]);
    if !(left.is_finite() && right.is_finite()) {
        return Err(CorrelatedThresholdError::UndefinedCurve { side });
    }
    let denominator = left - 2.0 * center + right;
    let step = probabilities[index + 1] - probabilities[index];
    let offset = if denominator < 0.0 {
        0.5 * (left - right) / denominator
    } else {
        0.0
    };
    let location = probabilities[index] + step * offset;
    Ok(IntervalEstimate {
        value: location,
        lower: location,
        upper: location,
        standard_error: 0.0,
    })
}

/// Local log-quadratic peak of `curve`: returns `(location, height)`.
///
/// The fit is `ln chi = a0 + a1 x + a2 x^2` with `x = p - p_max` over the
/// grid points within `radius` of the largest measured value. Refining the
/// peak over a physical probability window (rather than three adjacent
/// grid points) suppresses the point-to-point noise that makes the raw
/// peak location jump by several grid steps between independent replica
/// pools.
fn log_quadratic_peak(
    probabilities: &[f64],
    curve: &[f64],
    side: usize,
    radius: f64,
) -> Result<(f64, f64), CorrelatedThresholdError> {
    let mut best_index = usize::MAX;
    let mut best_value = f64::NEG_INFINITY;
    for (index, &value) in curve.iter().enumerate() {
        if value.is_finite() && value > best_value {
            best_index = index;
            best_value = value;
        }
    }
    if best_index == usize::MAX {
        return Err(CorrelatedThresholdError::UndefinedCurve { side });
    }
    let center = probabilities[best_index];
    let (mut s0, mut s1, mut s2, mut s3, mut s4) = (0.0, 0.0, 0.0, 0.0, 0.0);
    let (mut t0, mut t1, mut t2) = (0.0, 0.0, 0.0);
    for (&probability, &value) in probabilities.iter().zip(curve) {
        let offset = probability - center;
        if offset.abs() <= radius && value > 0.0 && value.is_finite() {
            let y = value.ln();
            let x2 = offset * offset;
            s0 += 1.0;
            s1 += offset;
            s2 += x2;
            s3 += x2 * offset;
            s4 += x2 * x2;
            t0 += y;
            t1 += offset * y;
            t2 += x2 * y;
        }
    }
    if s0 < 5.0 {
        return Err(CorrelatedThresholdError::PeakFit { side });
    }
    // 3x3 normal equations, solved by Gaussian elimination with pivoting.
    let mut matrix = [[s0, s1, s2], [s1, s2, s3], [s2, s3, s4]];
    let mut rhs = [t0, t1, t2];
    for column in 0..3 {
        let pivot = (column..3)
            .max_by(|&left, &right| {
                matrix[left][column]
                    .abs()
                    .total_cmp(&matrix[right][column].abs())
            })
            .expect("non-empty pivot range");
        matrix.swap(column, pivot);
        rhs.swap(column, pivot);
        if !(matrix[column][column].abs() > 0.0) {
            return Err(CorrelatedThresholdError::PeakFit { side });
        }
        for row in column + 1..3 {
            let factor = matrix[row][column] / matrix[column][column];
            let pivot_row = matrix[column];
            for (entry, pivot) in matrix[row][column..].iter_mut().zip(&pivot_row[column..]) {
                *entry -= factor * pivot;
            }
            rhs[row] -= factor * rhs[column];
        }
    }
    let mut coefficients = [0.0_f64; 3];
    for row in (0..3).rev() {
        let known: f64 = (row + 1..3)
            .map(|column| matrix[row][column] * coefficients[column])
            .sum();
        coefficients[row] = (rhs[row] - known) / matrix[row][row];
    }
    let [a0, a1, a2] = coefficients;
    if !(a2 < 0.0) {
        return Err(CorrelatedThresholdError::PeakFit { side });
    }
    let location = center - a1 / (2.0 * a2);
    let first = probabilities
        .first()
        .copied()
        .ok_or(CorrelatedThresholdError::UndefinedCurve { side })?;
    let last = probabilities
        .last()
        .copied()
        .ok_or(CorrelatedThresholdError::UndefinedCurve { side })?;
    if !(location >= first && location <= last) {
        return Err(CorrelatedThresholdError::PeakOutsideWindow { side });
    }
    let height = (a0 - a1 * a1 / (4.0 * a2)).exp();
    Ok((location, height))
}

/// Maximum central-difference derivative of `ln chi(p)`, parabola-refined;
/// `NaN` when the curve has no interior maximum.
fn maximum_log_derivative(probabilities: &[f64], curve: &[f64]) -> f64 {
    let points = curve.len();
    if points < 3 {
        return f64::NAN;
    }
    let mut derivative = vec![f64::NAN; points];
    for index in 1..points - 1 {
        let (left, center, right) = (curve[index - 1], curve[index], curve[index + 1]);
        if left > 0.0 && center > 0.0 && right > 0.0 {
            derivative[index] =
                (right.ln() - left.ln()) / (probabilities[index + 1] - probabilities[index - 1]);
        }
    }
    let mut best_index = usize::MAX;
    let mut best_value = f64::NEG_INFINITY;
    for (index, &value) in derivative.iter().enumerate() {
        if value.is_finite() && value > best_value {
            best_index = index;
            best_value = value;
        }
    }
    if best_index == usize::MAX {
        return f64::NAN;
    }
    let (index, value) = (best_index, best_value);
    if index == 0 || index + 1 >= points {
        return value;
    }
    let (left, right) = (derivative[index - 1], derivative[index + 1]);
    if !(left.is_finite() && right.is_finite()) {
        return value;
    }
    let denominator = left - 2.0 * value + right;
    if denominator >= 0.0 {
        return value;
    }
    let offset = 0.5 * (left - right) / denominator;
    // Vertex of the parabola through the three derivative samples.
    value - 0.25 * offset * (left - right)
}

fn interpolate(x: &[f64], y: &[f64], target: f64) -> f64 {
    if x.is_empty() || y.is_empty() {
        return f64::NAN;
    }
    if target <= x[0] {
        return y[0];
    }
    for index in 1..x.len() {
        if target <= x[index] {
            let width = x[index] - x[index - 1];
            if width <= 0.0 {
                return y[index];
            }
            let weight = (target - x[index - 1]) / width;
            return y[index - 1] * (1.0 - weight) + y[index] * weight;
        }
    }
    y[y.len() - 1]
}

fn mean(values: &[f64]) -> f64 {
    values.iter().sum::<f64>() / values.len() as f64
}

fn standard_error(values: &[f64]) -> f64 {
    if values.len() < 2 {
        return 0.0;
    }
    let average = mean(values);
    let variance = values
        .iter()
        .map(|value| (value - average) * (value - average))
        .sum::<f64>()
        / (values.len() - 1) as f64;
    (variance / values.len() as f64).sqrt()
}

/// Sample standard deviation of bootstrap replicates. This IS the bootstrap
/// standard-error estimate for a single-ladder statistic; dividing it by
/// `sqrt(B)` would understate the CI. (`standard_error` above is used only
/// for the mean of independent ladder peaks in `pool`, where the /sqrt(n)
/// is the correct between-ladder error.)
fn standard_deviation(values: &[f64]) -> f64 {
    if values.len() < 2 {
        return 0.0;
    }
    let average = mean(values);
    let variance = values
        .iter()
        .map(|value| (value - average) * (value - average))
        .sum::<f64>()
        / (values.len() - 1) as f64;
    variance.sqrt()
}

fn sorted(values: &[f64]) -> Vec<f64> {
    let mut sorted = values.to_vec();
    sorted.sort_unstable_by(f64::total_cmp);
    sorted
}

fn percentile(sorted: &[f64], quantile: f64) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    let position = quantile * (sorted.len() - 1) as f64;
    let lower = position.floor() as usize;
    let upper = position.ceil() as usize;
    let weight = position - lower as f64;
    sorted[lower] * (1.0 - weight) + sorted[upper] * weight
}

/// Weighted least-squares line through at least three points.
fn weighted_line(x: &[f64], y: &[f64], sigma: &[f64]) -> Option<LinearFit> {
    if x.len() < 3 || x.len() != y.len() || x.len() != sigma.len() {
        return None;
    }
    let mut sum_w = 0.0;
    let mut sum_x = 0.0;
    let mut sum_y = 0.0;
    let mut sum_xx = 0.0;
    let mut sum_xy = 0.0;
    for ((&x, &y), &sigma) in x.iter().zip(y).zip(sigma) {
        if !(x.is_finite() && y.is_finite() && sigma.is_finite() && sigma > 0.0) {
            return None;
        }
        let weight = 1.0 / (sigma * sigma);
        sum_w += weight;
        sum_x += weight * x;
        sum_y += weight * y;
        sum_xx += weight * x * x;
        sum_xy += weight * x * y;
    }
    let determinant = sum_w * sum_xx - sum_x * sum_x;
    if !(determinant > 0.0) {
        return None;
    }
    let slope = (sum_w * sum_xy - sum_x * sum_y) / determinant;
    let intercept = (sum_y - slope * sum_x) / sum_w;
    let slope_variance = sum_w / determinant;
    let intercept_variance = sum_xx / determinant;
    let mut chi_square = 0.0;
    for ((&x, &y), &sigma) in x.iter().zip(y).zip(sigma) {
        let residual = y - (intercept + slope * x);
        chi_square += (residual / sigma) * (residual / sigma);
    }
    let degrees_of_freedom = x.len() - 2;
    Some(LinearFit {
        intercept,
        slope,
        intercept_standard_error: intercept_variance.sqrt(),
        slope_standard_error: slope_variance.sqrt(),
        reduced_chi_square: chi_square / degrees_of_freedom as f64,
        degrees_of_freedom,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn probability_grid_and_spec_validation() {
        let spec = ThresholdLadderSpec::new(1.5, &[16, 32], 32, 7).unwrap();
        let grid = spec.probabilities();
        assert_eq!(grid.len(), 89);
        assert!((grid[0] - 0.49).abs() < 1e-12);
        assert!((grid[grid.len() - 1] - 0.60).abs() < 1e-12);

        assert!(matches!(
            ThresholdLadderSpec::new(0.0, &[16], 32, 7),
            Err(CorrelatedThresholdError::Exponent { .. })
        ));
        assert!(matches!(
            ThresholdLadderSpec::new(f64::NAN, &[16], 32, 7),
            Err(CorrelatedThresholdError::Exponent { .. })
        ));
        assert!(matches!(
            ThresholdLadderSpec::new(1.0, &[32, 16], 32, 7),
            Err(CorrelatedThresholdError::Sides { .. })
        ));
        assert!(matches!(
            ThresholdLadderSpec::with_probability_window(1.0, &[16, 32], 0.6, 0.5, 10, 32, 7),
            Err(CorrelatedThresholdError::ProbabilityWindow { .. })
        ));
        assert!(matches!(
            ThresholdLadderSpec::new(1.0, &[16, 32], 1, 7),
            Err(CorrelatedThresholdError::Samples { .. })
        ));
    }

    #[test]
    fn weighted_line_recovers_an_exact_affine_relation() {
        let x = [0.1, 0.2, 0.5, 0.9];
        let y: Vec<f64> = x.iter().map(|x| 0.4 + 2.5 * x).collect();
        let sigma = [0.01, 0.02, 0.01, 0.03];
        let fit = weighted_line(&x, &y, &sigma).unwrap();
        assert!((fit.intercept - 0.4).abs() < 1e-10);
        assert!((fit.slope - 2.5).abs() < 1e-10);
        assert!(fit.reduced_chi_square < 1e-12);
    }

    #[test]
    fn percentile_matches_hand_values() {
        let values = sorted(&[1.0, 2.0, 3.0, 4.0, 5.0]);
        assert_eq!(percentile(&values, 0.0), 1.0);
        assert_eq!(percentile(&values, 1.0), 5.0);
        assert_eq!(percentile(&values, 0.5), 3.0);
    }

    /// A ladder whose susceptibility curves are exact Gaussians centred on
    /// `p_c(L) = pc_inf - amplitude L^{-1/nu}`, so the pc-shift estimator
    /// has an exact answer.
    fn synthetic_ladder(
        exponent: f64,
        pc_inf: f64,
        nu: f64,
        amplitude: f64,
        shift: f64,
    ) -> ThresholdLadder {
        let sides = [64_usize, 128, 256, 512];
        let probabilities: Vec<f64> = (0..=400)
            .map(|index| 0.47 + 0.0002 * index as f64)
            .collect();
        let mut rungs = Vec::new();
        let mut samples = Vec::new();
        for side in sides {
            let location = pc_inf - amplitude * (side as f64).powf(-1.0 / nu) + shift;
            let susceptibility: Vec<f64> = probabilities
                .iter()
                .map(|&p| 100.0 * (-0.5 * ((p - location) / 0.01).powi(2)).exp())
                .collect();
            let peak = susceptibility
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max);
            rungs.push(RungObservables {
                side,
                pseudo_critical: IntervalEstimate {
                    value: location,
                    lower: location,
                    upper: location,
                    standard_error: 1e-4,
                },
                susceptibility_peak: peak,
                largest_fraction_at_peak: 0.0,
                susceptibility,
                finite_cluster_susceptibility: vec![0.0; probabilities.len()],
                largest_fraction: vec![0.0; probabilities.len()],
                defined_fraction: vec![1.0; probabilities.len()],
            });
            samples.push(RungSamples {
                slope_value: 1.0,
                slope_bootstrap: vec![1.0, 1.1],
            });
        }
        ThresholdLadder {
            exponent,
            probabilities,
            rungs,
            samples,
            replicas: 1,
        }
    }

    #[test]
    fn log_quadratic_peak_recovers_a_gaussian_peak() {
        let probabilities: Vec<f64> = (0..=200)
            .map(|index| 0.48 + 0.0002 * index as f64)
            .collect();
        let (location, width, height) = (0.5001, 0.006, 123.0);
        let curve: Vec<f64> = probabilities
            .iter()
            .map(|&p| height * (-0.5 * ((p - location) / width).powi(2)).exp())
            .collect();
        let (found, found_height) = log_quadratic_peak(&probabilities, &curve, 512, 0.004).unwrap();
        assert!((found - location).abs() < 1e-9, "{found} vs {location}");
        assert!((found_height - height).abs() < 1e-6);
    }

    #[test]
    fn pc_shift_recovers_the_synthetic_exponent() {
        let ladder = synthetic_ladder(1.6, 0.543, 1.6, 0.8, 0.0);
        let estimate = ladder
            .nu_from_pc_shift(64, 0.543, 0.004)
            .expect("pc-shift fit");
        assert!(
            (estimate.nu.value - 1.6).abs() < 1e-6,
            "nu={}",
            estimate.nu.value
        );
        assert!(estimate.reduced_chi_square < 1e-6);
        assert_eq!(estimate.used_sides, vec![64, 128, 256, 512]);
        assert!((estimate.amplitude - 0.8).abs() < 1e-9);
    }

    #[test]
    fn pc_shift_rejects_a_threshold_below_a_peak() {
        let ladder = synthetic_ladder(1.6, 0.543, 1.6, 0.8, 0.0);
        assert!(matches!(
            ladder.nu_from_pc_shift(64, 0.50, 0.004),
            Err(CorrelatedThresholdError::ThresholdOrder { .. })
        ));
        assert!(matches!(
            ladder.nu_from_pc_shift(64, 0.543, 0.0),
            Err(CorrelatedThresholdError::PeakRadius { .. })
        ));
        assert!(matches!(
            ladder.nu_from_pc_shift(64, f64::NAN, 0.004),
            Err(CorrelatedThresholdError::InfiniteThreshold { .. })
        ));
    }

    #[test]
    fn pool_averages_the_curves_and_uses_the_between_ladder_scatter() {
        let left = synthetic_ladder(1.6, 0.543, 1.6, 0.8, 0.0);
        let right = synthetic_ladder(1.6, 0.543, 1.6, 0.8, 0.002);
        let pooled = ThresholdLadder::pool(&[left, right]).expect("pool");
        assert_eq!(pooled.replicas(), 2);
        let rung = &pooled.rungs()[3];
        assert_eq!(rung.side, 512);
        // The pooled Gaussian centre is the midpoint of the two shifts.
        let midpoint = 0.543 - 0.8 * 512.0_f64.powf(-1.0 / 1.6) + 0.001;
        assert!((rung.pseudo_critical.value - midpoint).abs() < 1e-6);
        assert!(rung.pseudo_critical.standard_error > 0.0);
        assert!(rung.pseudo_critical.lower < rung.pseudo_critical.value);
        assert!(rung.pseudo_critical.upper > rung.pseudo_critical.value);
    }

    #[test]
    fn pool_rejects_mismatched_ladders() {
        let left = synthetic_ladder(1.6, 0.543, 1.6, 0.8, 0.0);
        let right = synthetic_ladder(1.25, 0.543, 1.6, 0.8, 0.0);
        assert!(matches!(
            ThresholdLadder::pool(&[left, right]),
            Err(CorrelatedThresholdError::PoolMismatch)
        ));
        let single = synthetic_ladder(1.6, 0.543, 1.6, 0.8, 0.0);
        assert!(matches!(
            ThresholdLadder::pool(&[single]),
            Err(CorrelatedThresholdError::Samples { .. })
        ));
    }
}
