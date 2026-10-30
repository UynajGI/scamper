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

use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

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

/// Finite-size scaling window `(p_c, nu, exponent)`.
///
/// `coordinate` is `(p - p_c) L^{1/nu}`; `scaled_value` is
/// `value * L^{-exponent}`. Exponent fitting is deliberately caller-driven:
/// score candidate windows with [`collapse_residual`] instead of hiding a
/// nonlinear optimizer here.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScalingWindow {
    pc: f64,
    nu: f64,
    exponent: f64,
}

impl ScalingWindow {
    pub fn new(pc: f64, nu: f64, exponent: f64) -> Result<Self, FssError> {
        if !pc.is_finite() || !(0.0..=1.0).contains(&pc) {
            return Err(FssError::Window {
                reason: "p_c must be finite and inside [0, 1]",
            });
        }
        if !nu.is_finite() || nu <= 0.0 {
            return Err(FssError::Window {
                reason: "nu must be finite and strictly positive",
            });
        }
        if !exponent.is_finite() {
            return Err(FssError::Window {
                reason: "the value exponent must be finite",
            });
        }
        Ok(Self { pc, nu, exponent })
    }

    pub fn coordinate(&self, size: usize, probability: f64) -> f64 {
        (probability - self.pc) * (size as f64).powf(1.0 / self.nu)
    }

    pub fn scaled_value(&self, size: usize, value: f64) -> f64 {
        value * (size as f64).powf(-self.exponent)
    }
}

/// One measurement to place on a collapse plot.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CollapsePoint {
    pub size: usize,
    pub probability: f64,
    pub value: f64,
    pub standard_error: f64,
}

/// Pairwise collapse consistency over cross-size points.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CollapseQuality {
    /// Cross-size pairs whose scaling coordinates are within tolerance.
    pub pairs: usize,
    /// Root-mean-square of `(y_i - y_j) / sqrt(se_i^2 + se_j^2)`.
    pub rms_z: f64,
}

/// Score a candidate scaling window by cross-size collapse consistency.
///
/// Every pair of points from *different* sizes whose `x` coordinates lie
/// within `coordinate_tolerance` contributes
/// `z = (y_i - y_j) / sqrt(se_i^2 + se_j^2)`; the returned `rms_z` is the
/// root mean square over all such pairs (lower is better). Zero combined
/// error with a nonzero scaled-value difference is a typed rejection: the
/// pair claims infinite disagreement, which is not a score. Pairs are not
/// independent (each point may contribute several times), so this is a
/// comparative quality metric for selecting a window, not a confidence
/// test.
pub fn collapse_residual(
    points: &[CollapsePoint],
    window: &ScalingWindow,
    coordinate_tolerance: f64,
) -> Result<CollapseQuality, FssError> {
    if points.is_empty() {
        return Err(FssError::EmptyCollapse);
    }
    if !coordinate_tolerance.is_finite() || coordinate_tolerance <= 0.0 {
        return Err(FssError::Tolerance {
            tolerance: coordinate_tolerance,
        });
    }
    let mut scaled = Vec::with_capacity(points.len());
    for point in points {
        FssPoint::new(
            point.size,
            point.probability,
            point.value,
            point.standard_error,
        )?;
        scaled.push((
            point.size,
            window.coordinate(point.size, point.probability),
            window.scaled_value(point.size, point.value),
            point.standard_error,
        ));
    }
    let mut pairs = 0usize;
    let mut sum_squares = 0.0;
    for first in 0..scaled.len() {
        for second in (first + 1)..scaled.len() {
            let (size_first, x_first, y_first, se_first) = scaled[first];
            let (size_second, x_second, y_second, se_second) = scaled[second];
            if size_first == size_second || (x_first - x_second).abs() > coordinate_tolerance {
                continue;
            }
            let difference = y_first - y_second;
            let error = (se_first * se_first + se_second * se_second).sqrt();
            if error == 0.0 {
                if difference == 0.0 {
                    continue;
                }
                return Err(FssError::DegeneratePair { first, second });
            }
            let z = difference / error;
            sum_squares += z * z;
            pairs += 1;
        }
    }
    if pairs == 0 {
        return Err(FssError::NoCollapsePairs);
    }
    Ok(CollapseQuality {
        pairs,
        rms_z: (sum_squares / pairs as f64).sqrt(),
    })
}

/// Percentile bootstrap interval over independent samples.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BootstrapInterval {
    /// Mean of the resampled statistics.
    pub mean: f64,
    /// Sample standard deviation of the resampled statistics.
    pub standard_error: f64,
    /// 2.5th percentile (nearest rank).
    pub lower: f64,
    /// 97.5th percentile (nearest rank).
    pub upper: f64,
    /// Number of bootstrap resamples.
    pub resamples: usize,
}

/// Percentile bootstrap of an arbitrary statistic over independent samples.
///
/// Draw contract: a fresh `Xoshiro256PlusPlus::seed_from_u64(seed)` draws
/// `samples.len()` indices with `random_range(0..len)` per resample, in
/// order, into one reused buffer; the statistic is then called on that
/// buffer. Identical inputs therefore reproduce identical intervals
/// bitwise. Allocations (the index-free buffer and the resample vector) are
/// reserved before the first draw. Non-finite inputs or statistic values
/// are typed rejections, never silently skipped.
pub fn bootstrap<F>(
    samples: &[f64],
    resamples: usize,
    seed: u64,
    mut statistic: F,
) -> Result<BootstrapInterval, FssError>
where
    F: FnMut(&[f64]) -> f64,
{
    if samples.is_empty() {
        return Err(FssError::EmptySamples);
    }
    if resamples < 2 {
        return Err(FssError::Resamples { resamples });
    }
    if let Some(value) = samples.iter().copied().find(|value| !value.is_finite()) {
        return Err(FssError::Value { value });
    }
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut buffer = vec![0.0; samples.len()];
    let mut values = Vec::new();
    values
        .try_reserve_exact(resamples)
        .map_err(|_| FssError::Capacity)?;
    for _ in 0..resamples {
        for slot in &mut buffer {
            *slot = samples[rng.random_range(0..samples.len())];
        }
        let value = statistic(&buffer);
        if !value.is_finite() {
            return Err(FssError::Statistic { value });
        }
        values.push(value);
    }
    let mut moments = Moments::new();
    for &value in &values {
        moments.add(value);
    }
    let mean = moments.mean().expect("resamples >= 2");
    let standard_error = moments.standard_error().expect("resamples >= 2");
    values.sort_by(f64::total_cmp);
    let last = values.len() - 1;
    Ok(BootstrapInterval {
        mean,
        standard_error,
        lower: values[(0.025 * last as f64).floor() as usize],
        upper: values[(0.975 * last as f64).ceil() as usize],
        resamples,
    })
}

/// One measured two-point correlation value at a distance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CorrelationPoint {
    /// Distance in lattice spacings.
    pub distance: f64,
    /// Non-negative correlation value `C(r)`.
    pub value: f64,
}

/// Second-moment correlation length from a truncated correlation function.
///
/// Definitional convention: with `d` spatial dimensions and non-negative
/// `C(r)`,
/// `xi^2 = sum_r r^2 C(r) / (2 d sum_r C(r))`,
/// the standard second-moment estimator (Amit & Martin-Mayor). Distances
/// are in lattice spacings; the caller is responsible for truncating the
/// sum where the signal ends (negative values are rejected rather than
/// silently clipping). `C(0)` contributes zero to the numerator, so
/// including or omitting it changes only the denominator.
pub fn second_moment_correlation_length(
    points: &[CorrelationPoint],
    dimensions: usize,
) -> Result<f64, FssError> {
    if dimensions == 0 {
        return Err(FssError::Dimensions { dimensions });
    }
    if points.is_empty() {
        return Err(FssError::EmptyCorrelation);
    }
    let mut numerator = 0.0;
    let mut denominator = 0.0;
    for (index, point) in points.iter().enumerate() {
        if !point.distance.is_finite() || point.distance < 0.0 {
            return Err(FssError::Distance {
                index,
                distance: point.distance,
            });
        }
        if !point.value.is_finite() || point.value < 0.0 {
            return Err(FssError::CorrelationValue {
                index,
                value: point.value,
            });
        }
        numerator += point.distance * point.distance * point.value;
        denominator += point.value;
    }
    if denominator <= 0.0 {
        return Err(FssError::UndefinedCorrelation);
    }
    Ok((numerator / (2.0 * dimensions as f64 * denominator)).sqrt())
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
    EmptyCollapse,
    Window {
        reason: &'static str,
    },
    Tolerance {
        tolerance: f64,
    },
    NoCollapsePairs,
    DegeneratePair {
        first: usize,
        second: usize,
    },
    EmptySamples,
    EmptyCorrelation,
    Dimensions {
        dimensions: usize,
    },
    Distance {
        index: usize,
        distance: f64,
    },
    CorrelationValue {
        index: usize,
        value: f64,
    },
    UndefinedCorrelation,
    Resamples {
        resamples: usize,
    },
    Statistic {
        value: f64,
    },
    Capacity,
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
            Self::EmptyCollapse => write!(formatter, "collapse needs at least one point"),
            Self::Window { reason } => write!(formatter, "invalid scaling window: {reason}"),
            Self::Tolerance { tolerance } => {
                write!(
                    formatter,
                    "coordinate tolerance {tolerance} must be finite and positive"
                )
            }
            Self::NoCollapsePairs => {
                write!(
                    formatter,
                    "no cross-size pair falls inside the coordinate tolerance"
                )
            }
            Self::DegeneratePair { first, second } => write!(
                formatter,
                "points {first} and {second} differ with zero combined uncertainty"
            ),
            Self::EmptySamples => write!(formatter, "bootstrap needs at least one sample"),
            Self::EmptyCorrelation => {
                write!(formatter, "correlation length needs at least one point")
            }
            Self::Dimensions { dimensions } => {
                write!(formatter, "dimensions must be at least 1, got {dimensions}")
            }
            Self::Distance { index, distance } => {
                write!(
                    formatter,
                    "point {index}: distance {distance} is not finite and non-negative"
                )
            }
            Self::CorrelationValue { index, value } => {
                write!(
                    formatter,
                    "point {index}: correlation {value} is not finite and non-negative"
                )
            }
            Self::UndefinedCorrelation => {
                write!(
                    formatter,
                    "correlation length is undefined with zero total weight"
                )
            }
            Self::Resamples { resamples } => {
                write!(
                    formatter,
                    "bootstrap needs at least two resamples, got {resamples}"
                )
            }
            Self::Statistic { value } => {
                write!(
                    formatter,
                    "bootstrap statistic returned a non-finite value {value}"
                )
            }
            Self::Capacity => write!(formatter, "could not reserve bootstrap storage"),
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

    fn mean(values: &[f64]) -> f64 {
        values.iter().sum::<f64>() / values.len() as f64
    }

    #[test]
    fn bootstrap_of_constant_samples_is_degenerate() {
        let interval = bootstrap(&[2.0; 8], 64, 0x4135_0001, mean).expect("bootstrap");
        assert_eq!(interval.mean, 2.0);
        assert_eq!(interval.standard_error, 0.0);
        assert_eq!(interval.lower, 2.0);
        assert_eq!(interval.upper, 2.0);
        assert_eq!(interval.resamples, 64);
    }

    #[test]
    fn bootstrap_mean_replays_bitwise_and_centers_on_the_sample_mean() {
        let samples: Vec<f64> = (0..10).map(f64::from).collect();
        let first = bootstrap(&samples, 1_000, 0x4135_0002, mean).expect("bootstrap");
        let second = bootstrap(&samples, 1_000, 0x4135_0002, mean).expect("bootstrap");
        assert_eq!(first.mean.to_bits(), second.mean.to_bits());
        assert_eq!(
            first.standard_error.to_bits(),
            second.standard_error.to_bits()
        );
        assert_eq!(first.lower.to_bits(), second.lower.to_bits());
        assert_eq!(first.upper.to_bits(), second.upper.to_bits());
        assert!((first.mean - 4.5).abs() < 0.15);
        assert!(first.lower < first.mean && first.mean < first.upper);

        let maximum = bootstrap(&samples, 256, 0x4135_0003, |values| {
            values.iter().copied().fold(f64::NEG_INFINITY, f64::max)
        })
        .expect("maximum bootstrap");
        assert!(maximum.lower >= 0.0 && maximum.upper <= 9.0);
    }

    #[test]
    fn bootstrap_rejections_are_typed() {
        assert!(matches!(
            bootstrap(&[], 10, 0, mean),
            Err(FssError::EmptySamples)
        ));
        assert!(matches!(
            bootstrap(&[1.0], 1, 0, mean),
            Err(FssError::Resamples { resamples: 1 })
        ));
        assert!(matches!(
            bootstrap(&[1.0, f64::NAN], 10, 0, mean),
            Err(FssError::Value { .. })
        ));
        assert!(matches!(
            bootstrap(&[1.0, 2.0], 10, 0, |_| f64::NAN),
            Err(FssError::Statistic { .. })
        ));
    }

    #[test]
    fn correlation_length_matches_the_flat_correlation_arithmetic() {
        let points: Vec<CorrelationPoint> = (0..=10)
            .map(|r| CorrelationPoint {
                distance: f64::from(r),
                value: 1.0,
            })
            .collect();
        let xi = second_moment_correlation_length(&points, 1).expect("correlation length");
        // sum r^2 = 385 over 11 points, divided by 2*d*11.
        assert!((xi - (385.0_f64 / 22.0).sqrt()).abs() < 1e-15);
    }

    #[test]
    fn correlation_length_recovers_an_exponential_decay() {
        let decay = 2.0_f64;
        let points: Vec<CorrelationPoint> = (0..=6_000)
            .map(|step| {
                let r = step as f64 * 0.01;
                CorrelationPoint {
                    distance: r,
                    value: (-r / decay).exp(),
                }
            })
            .collect();
        let xi = second_moment_correlation_length(&points, 1).expect("correlation length");
        assert!((xi - decay).abs() < 0.01, "xi = {xi} vs {decay}");
    }

    #[test]
    fn correlation_length_rejections_are_typed() {
        let zero = [CorrelationPoint {
            distance: 0.0,
            value: 0.0,
        }];
        assert!(matches!(
            second_moment_correlation_length(&zero, 1),
            Err(FssError::UndefinedCorrelation)
        ));
        assert!(matches!(
            second_moment_correlation_length(&zero, 0),
            Err(FssError::Dimensions { dimensions: 0 })
        ));
        assert!(matches!(
            second_moment_correlation_length(&[], 1),
            Err(FssError::EmptyCorrelation)
        ));
        assert!(matches!(
            second_moment_correlation_length(
                &[CorrelationPoint {
                    distance: -1.0,
                    value: 1.0
                }],
                1
            ),
            Err(FssError::Distance { index: 0, .. })
        ));
        assert!(matches!(
            second_moment_correlation_length(
                &[CorrelationPoint {
                    distance: 1.0,
                    value: -0.1
                }],
                1
            ),
            Err(FssError::CorrelationValue { index: 0, .. })
        ));
    }

    #[test]
    fn exact_scaling_collapses_and_a_perturbation_worsens_the_score() {
        let window = ScalingWindow::new(0.5, 2.0, 0.5).expect("window");
        let mut points = Vec::new();
        for size in [16usize, 32, 64] {
            for x in [-1.0_f64, 0.0, 1.0] {
                let probability = 0.5 + x / (size as f64).powf(0.5);
                let value = (1.0 / (1.0 + x * x)) * (size as f64).powf(0.5);
                points.push(CollapsePoint {
                    size,
                    probability,
                    value,
                    standard_error: 0.01,
                });
            }
        }
        let exact = collapse_residual(&points, &window, 1e-9).expect("pairs");
        assert_eq!(exact.pairs, 9);
        assert!(exact.rms_z < 1e-12, "rms_z = {}", exact.rms_z);

        points[0].value += 0.1;
        let perturbed = collapse_residual(&points, &window, 1e-9).expect("pairs");
        assert!(perturbed.rms_z > exact.rms_z);
    }

    #[test]
    fn collapse_rejections_are_typed() {
        let window = ScalingWindow::new(0.5, 2.0, 0.5).expect("window");
        let single = [CollapsePoint {
            size: 8,
            probability: 0.5,
            value: 1.0,
            standard_error: 0.1,
        }];
        assert!(matches!(
            collapse_residual(&single, &window, 1e-9),
            Err(FssError::NoCollapsePairs)
        ));
        assert!(matches!(
            collapse_residual(&[], &window, 1e-9),
            Err(FssError::EmptyCollapse)
        ));
        assert!(matches!(
            collapse_residual(&single, &window, 0.0),
            Err(FssError::Tolerance { .. })
        ));
        assert!(matches!(
            ScalingWindow::new(0.5, 0.0, 0.5),
            Err(FssError::Window { .. })
        ));

        let degenerate = [
            CollapsePoint {
                size: 8,
                probability: 0.5,
                value: 1.0,
                standard_error: 0.0,
            },
            CollapsePoint {
                size: 16,
                probability: 0.5,
                value: 2.0,
                standard_error: 0.0,
            },
        ];
        assert!(matches!(
            collapse_residual(&degenerate, &window, 1.0),
            Err(FssError::DegeneratePair {
                first: 0,
                second: 1
            })
        ));
    }

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
