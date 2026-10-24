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
//! ranking. For each grid point the summary gives
//!
//! ```text
//! chi(p) = (M2 - S_max^2) / (V_active - S_max)
//! ```
//!
//! the F3 per-sample finite-cluster susceptibility (mean cluster size with
//! one largest cluster excluded), and `S_max / N`. `chi(p)` is averaged over
//! realizations (undefined samples are skipped and their fraction is
//! reported). The pseudo-critical point `p_c(L)` is the **peak of the
//! ensemble-averaged `chi(p)`**, refined by a parabola through the three
//! grid points around the maximum. The secondary "large-cluster" curve is
//! reported for context; the susceptibility peak is the primary estimator
//! because the largest-cluster fraction is monotone rather than peaked.
//!
//! Statistical intervals are 95% percentile bootstrap intervals over
//! realizations (200 resamples, fixed seeds derived from the spec seed).
//! They quantify sampling noise only. The extrapolation additionally
//! reports the reduced chi-square of the `L^{-1/nu}` fit, which exposes the
//! finite-size/systematic scatter that the bootstrap cannot see. Callers
//! must vary `minimum_side` to see the fit-window sensitivity; this module
//! never hides it behind a single number.
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
    /// The susceptibility curve had no finite value.
    UndefinedCurve { side: usize },
    /// The correlation-length exponent must be finite and positive.
    FitExponent { nu: f64 },
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
        standard_error: standard_error(&finite_pc),
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
}
