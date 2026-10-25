//! Generic finite-size-scaling toolkit (v3 A5).
//!
//! The tools here consume *reported* per-realization quantities, never a
//! particular adapter. A5a delivers:
//!
//! - [`BinderEstimator`]: the Binder cumulant of the order parameter
//!   `m = S_max / V` from independent realizations,
//!   `U = 1 - <m^4> / (3 <m^2>^2)`, with a first-order (delta-method)
//!   standard error. Binder's definitional convention; no literature value
//!   is claimed by this layer.
//! - [`crossing`]: the raw linear crossing of two finite-size curves sampled
//!   on the same probability grid. This is a finite-size estimate, not a
//!   thermodynamic-limit extrapolation.
//!
//! Everything is deterministic, allocation-free per sample, and returns
//! typed rejections for undefined statistics. Sample-variance conventions
//! match [`crate::Moments`] (unbiased, clamped at zero).

use super::estimator::Moments;

/// Finite-size observation of one scalar at one probability.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FssPoint {
    /// Linear lattice size or another strictly ordered size label.
    pub size: usize,
    /// Occupation probability (or the tuned control parameter).
    pub probability: f64,
    /// Measured value at this point.
    pub value: f64,
    /// Standard error of `value` across independent realizations.
    pub standard_error: f64,
}

impl FssPoint {
    /// Validate one point's domain: non-zero size, probability in `[0, 1]`,
    /// and a finite value with a finite, non-negative error.
    pub fn new(
        size: usize,
        probability: f64,
        value: f64,
        standard_error: f64,
    ) -> Result<Self, FssError> {
        if size == 0 {
            return Err(FssError::ZeroSize);
        }
        if !(0.0..=1.0).contains(&probability) {
            return Err(FssError::Probability { probability });
        }
        if !value.is_finite() {
            return Err(FssError::Value { value });
        }
        if !standard_error.is_finite() || standard_error < 0.0 {
            return Err(FssError::StandardError { standard_error });
        }
        Ok(Self {
            size,
            probability,
            value,
            standard_error,
        })
    }
}

/// A finite-size curve: one size sampled on a strictly ascending
/// probability grid.
#[derive(Clone, Debug, PartialEq)]
pub struct FssCurve {
    points: Vec<FssPoint>,
}

impl FssCurve {
    /// Validate a curve: non-empty, one size, strictly ascending
    /// probabilities, and every point valid by [`FssPoint::new`].
    pub fn new(points: impl IntoIterator<Item = FssPoint>) -> Result<Self, FssError> {
        let points: Vec<FssPoint> = points.into_iter().collect();
        let Some(first) = points.first() else {
            return Err(FssError::EmptyCurve);
        };
        for (index, point) in points.iter().enumerate() {
            FssPoint::new(
                point.size,
                point.probability,
                point.value,
                point.standard_error,
            )?;
            if point.size != first.size {
                return Err(FssError::CurveSize {
                    index,
                    expected: first.size,
                    actual: point.size,
                });
            }
            if index > 0 && point.probability <= points[index - 1].probability {
                return Err(FssError::ProbabilityOrder { index });
            }
        }
        Ok(Self { points })
    }

    pub fn size(&self) -> usize {
        self.points[0].size
    }

    pub fn points(&self) -> &[FssPoint] {
        &self.points
    }
}

/// Linear-crossing estimate with a first-order interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CrossingEstimate {
    /// Interpolated crossing probability.
    pub probability: f64,
    /// First-order standard error from the four endpoint errors.
    pub standard_error: f64,
    /// `probability - 1.96 * standard_error`.
    pub lower: f64,
    /// `probability + 1.96 * standard_error`.
    pub upper: f64,
    /// Size of the lower curve.
    pub lower_size: usize,
    /// Size of the upper curve.
    pub upper_size: usize,
}

/// First raw crossing of two curves sharing a probability grid.
///
/// The curves must have different sizes and identical probabilities; the
/// difference `upper - lower` is linearly interpolated between the first
/// adjacent grid pair where it changes sign. Endpoint measurements are
/// treated as independent and the two grid points are treated as
/// independent (only their reported errors are available), so the interval
/// is a first-order bracket, not a rigorous confidence statement. Curves
/// that never cross are a typed rejection, not an empty estimate.
pub fn crossing(lower: &FssCurve, upper: &FssCurve) -> Result<CrossingEstimate, FssError> {
    if upper.size() <= lower.size() {
        return Err(FssError::SizeOrder {
            lower: lower.size(),
            upper: upper.size(),
        });
    }
    for (index, (lower_point, upper_point)) in lower.points().iter().zip(upper.points()).enumerate()
    {
        if lower_point.probability.to_bits() != upper_point.probability.to_bits() {
            return Err(FssError::GridMismatch {
                index,
                lower: lower_point.probability,
                upper: upper_point.probability,
            });
        }
    }
    if lower.points().len() != upper.points().len() {
        return Err(FssError::GridLength {
            lower: lower.points().len(),
            upper: upper.points().len(),
        });
    }
    for index in 1..lower.points().len() {
        let first = upper.points()[index - 1].value - lower.points()[index - 1].value;
        let second = upper.points()[index].value - lower.points()[index].value;
        if first.signum() == second.signum() {
            continue;
        }
        let span = upper.points()[index].probability - upper.points()[index - 1].probability;
        let difference = first - second;
        let probability = upper.points()[index - 1].probability + span * first / difference;
        // p* = p1 + span * d1/(d1 - d2); derivatives with respect to the
        // two differences (signs cancel in the ratio).
        let derivative_first = -span * second / (difference * difference);
        let derivative_second = span * first / (difference * difference);
        let first_error = upper.points()[index - 1].standard_error.powi(2)
            + lower.points()[index - 1].standard_error.powi(2);
        let second_error = upper.points()[index].standard_error.powi(2)
            + lower.points()[index].standard_error.powi(2);
        let standard_error = (derivative_first.powi(2) * first_error
            + derivative_second.powi(2) * second_error)
            .sqrt();
        return Ok(CrossingEstimate {
            probability,
            standard_error,
            lower: probability - 1.96 * standard_error,
            upper: probability + 1.96 * standard_error,
            lower_size: lower.size(),
            upper_size: upper.size(),
        });
    }
    Err(FssError::NoCrossing)
}

/// Binder cumulant and its delta-method standard error.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinderCumulant {
    /// `1 - <m^4> / (3 <m^2>^2)`.
    pub mean: f64,
    /// First-order standard error (bootstrap is a later A5 slice).
    pub standard_error: f64,
}

/// Running Binder cumulant of `m = S_max / V` over independent samples.
///
/// Accumulates the first moments of `m^2`, `m^4`, and their product; the
/// product supplies the covariance term of the delta-method error. The
/// error is first-order and becomes undefined when `<m^2> = 0` (all
/// realizations empty) or with fewer than two samples.
///
/// ```
/// use cmc_rs::BinderEstimator;
///
/// let mut binder = BinderEstimator::new();
/// binder.add(4, 4).expect("full occupancy");
/// binder.add(0, 4).expect("empty realization");
/// let cumulant = binder.estimate().expect("two samples");
/// assert!((cumulant.mean - 1.0 / 3.0).abs() < 1e-12);
/// ```
#[derive(Clone, Copy, Debug, Default)]
pub struct BinderEstimator {
    samples: usize,
    m2: Moments,
    m4: Moments,
    m2_m4: Moments,
}

impl BinderEstimator {
    /// An empty accumulator.
    pub const fn new() -> Self {
        Self {
            samples: 0,
            m2: Moments::new(),
            m4: Moments::new(),
            m2_m4: Moments::new(),
        }
    }

    /// Realizations accumulated.
    pub const fn samples(&self) -> usize {
        self.samples
    }

    /// Add one realization's largest-cluster size and topology size.
    ///
    /// Rejects `total == 0` and `largest > total`: both are adapter
    /// contract violations, not data to clamp.
    pub fn add(&mut self, largest: usize, total: usize) -> Result<(), FssError> {
        if total == 0 {
            return Err(FssError::EmptyTopology);
        }
        if largest > total {
            return Err(FssError::LargestExceedsTotal { largest, total });
        }
        let m = largest as f64 / total as f64;
        let m2 = m * m;
        let m4 = m2 * m2;
        self.samples += 1;
        self.m2.add(m2);
        self.m4.add(m4);
        self.m2_m4.add(m2 * m4);
        Ok(())
    }

    /// The cumulative Binder cumulant, or `None` before two samples or
    /// when every realization is empty (`<m^2> = 0`).
    pub fn estimate(&self) -> Option<BinderCumulant> {
        let mean2 = self.m2.mean()?;
        let mean4 = self.m4.mean()?;
        if mean2 == 0.0 {
            return None;
        }
        let mean = 1.0 - mean4 / (3.0 * mean2 * mean2);
        let variance2 = self.m2.sample_variance()?;
        let variance4 = self.m4.sample_variance()?;
        // cov(m^2, m^4) from E[m^6] - E[m^2]E[m^4], with the unbiased
        // sample correction n / (n - 1).
        let samples = self.samples as f64;
        let covariance = (self.m2_m4.mean()? - mean2 * mean4) * samples / (samples - 1.0);
        // U = 1 - A/(3B^2) with A = <m^4>, B = <m^2>.
        let variance = (variance4 / samples) / (9.0 * mean2.powi(4))
            + (4.0 * mean4 * mean4 * (variance2 / samples)) / (9.0 * mean2.powi(6))
            - (4.0 * mean4 * (covariance / samples)) / (9.0 * mean2.powi(5));
        Some(BinderCumulant {
            mean,
            standard_error: variance.max(0.0).sqrt(),
        })
    }
}

/// Typed finite-size-tool rejection.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FssError {
    ZeroSize,
    Probability {
        probability: f64,
    },
    Value {
        value: f64,
    },
    StandardError {
        standard_error: f64,
    },
    EmptyTopology,
    LargestExceedsTotal {
        largest: usize,
        total: usize,
    },
    EmptyCurve,
    CurveSize {
        index: usize,
        expected: usize,
        actual: usize,
    },
    ProbabilityOrder {
        index: usize,
    },
    SizeOrder {
        lower: usize,
        upper: usize,
    },
    GridLength {
        lower: usize,
        upper: usize,
    },
    GridMismatch {
        index: usize,
        lower: f64,
        upper: f64,
    },
    NoCrossing,
}

impl core::fmt::Display for FssError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::ZeroSize => write!(formatter, "finite-size labels must be non-zero"),
            Self::Probability { probability } => {
                write!(formatter, "probability {probability} outside [0, 1]")
            }
            Self::Value { value } => write!(formatter, "value {value} is not finite"),
            Self::StandardError { standard_error } => {
                write!(
                    formatter,
                    "standard error {standard_error} is not finite and non-negative"
                )
            }
            Self::EmptyTopology => write!(formatter, "Binder sample needs a non-empty topology"),
            Self::LargestExceedsTotal { largest, total } => write!(
                formatter,
                "largest component {largest} exceeds topology size {total}"
            ),
            Self::EmptyCurve => write!(formatter, "finite-size curve has no points"),
            Self::CurveSize {
                index,
                expected,
                actual,
            } => write!(
                formatter,
                "point {index} has size {actual}, curve size is {expected}"
            ),
            Self::ProbabilityOrder { index } => {
                write!(formatter, "point {index} does not increase the probability")
            }
            Self::SizeOrder { lower, upper } => write!(
                formatter,
                "crossing needs increasing sizes, got {lower} and {upper}"
            ),
            Self::GridLength { lower, upper } => write!(
                formatter,
                "curves have different grid lengths ({lower} vs {upper})"
            ),
            Self::GridMismatch {
                index,
                lower,
                upper,
            } => write!(formatter, "grid point {index} differs ({lower} vs {upper})"),
            Self::NoCrossing => write!(formatter, "the two curves never cross"),
        }
    }
}

impl std::error::Error for FssError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn binder_of_binary_occupancy_is_one_third() {
        let mut binder = BinderEstimator::new();
        for largest in [0, 4, 0, 4] {
            binder.add(largest, 4).expect("sample");
        }
        let estimate = binder.estimate().expect("cumulant");
        assert!((estimate.mean - 1.0 / 3.0).abs() < 1e-15);
        assert!(estimate.standard_error.is_finite());
    }

    #[test]
    fn binder_of_a_constant_order_parameter_is_two_thirds() {
        let mut binder = BinderEstimator::new();
        for _ in 0..10 {
            binder.add(2, 4).expect("sample");
        }
        let estimate = binder.estimate().expect("cumulant");
        assert!((estimate.mean - 2.0 / 3.0).abs() < 1e-15);
        assert_eq!(estimate.standard_error, 0.0);
    }

    #[test]
    fn binder_is_undefined_before_two_samples_or_without_occupancy() {
        let mut binder = BinderEstimator::new();
        assert_eq!(binder.estimate(), None);
        binder.add(0, 4).expect("sample");
        assert_eq!(binder.estimate(), None);
        binder.add(0, 4).expect("sample");
        assert_eq!(binder.estimate(), None);
        assert!(matches!(
            binder.add(5, 4),
            Err(FssError::LargestExceedsTotal {
                largest: 5,
                total: 4
            })
        ));
        assert!(matches!(binder.add(0, 0), Err(FssError::EmptyTopology)));
    }

    fn linear_curves() -> (FssCurve, FssCurve) {
        // lower(p) = p; upper(p) = 2p - 0.5; crossing at p = 0.5.
        let lower = FssCurve::new([
            FssPoint::new(8, 0.4, 0.4, 0.01).expect("point"),
            FssPoint::new(8, 0.6, 0.6, 0.01).expect("point"),
        ])
        .expect("curve");
        let upper = FssCurve::new([
            FssPoint::new(16, 0.4, 0.3, 0.01).expect("point"),
            FssPoint::new(16, 0.6, 0.7, 0.01).expect("point"),
        ])
        .expect("curve");
        (lower, upper)
    }

    #[test]
    fn linear_crossing_interpolates_between_the_bracketing_points() {
        let (lower, upper) = linear_curves();
        let estimate = crossing(&lower, &upper).expect("crossing");
        assert!((estimate.probability - 0.5).abs() < 1e-15);
        assert_eq!(estimate.lower_size, 8);
        assert_eq!(estimate.upper_size, 16);
        assert!(estimate.standard_error > 0.0);
        assert!(estimate.lower < estimate.probability && estimate.probability < estimate.upper);
    }

    #[test]
    fn crossing_rejections_are_typed() {
        let (lower, _) = linear_curves();
        let same_size = FssCurve::new([
            FssPoint::new(8, 0.4, 0.3, 0.01).expect("point"),
            FssPoint::new(8, 0.6, 0.7, 0.01).expect("point"),
        ])
        .expect("curve");
        assert!(matches!(
            crossing(&lower, &same_size),
            Err(FssError::SizeOrder { .. })
        ));

        let never = FssCurve::new([
            FssPoint::new(16, 0.4, 0.5, 0.01).expect("point"),
            FssPoint::new(16, 0.6, 0.8, 0.01).expect("point"),
        ])
        .expect("curve");
        assert!(matches!(
            crossing(&lower, &never),
            Err(FssError::NoCrossing)
        ));

        let shifted = FssCurve::new([
            FssPoint::new(16, 0.41, 0.3, 0.01).expect("point"),
            FssPoint::new(16, 0.6, 0.7, 0.01).expect("point"),
        ])
        .expect("curve");
        assert!(matches!(
            crossing(&lower, &shifted),
            Err(FssError::GridMismatch { index: 0, .. })
        ));

        assert!(matches!(
            FssCurve::new(Vec::new()),
            Err(FssError::EmptyCurve)
        ));
        assert!(matches!(
            FssCurve::new([
                FssPoint::new(8, 0.6, 0.3, 0.01).expect("point"),
                FssPoint::new(8, 0.4, 0.7, 0.01).expect("point"),
            ]),
            Err(FssError::ProbabilityOrder { index: 1 })
        ));
    }
}
