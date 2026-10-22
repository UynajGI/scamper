use core::fmt;
use std::collections::TryReserveError;
#[cfg(feature = "allocation-probe")]
use std::mem::size_of;

use crate::Probability;

use super::{ScanCurve, ScanPoint};

/// Scalar observable requested from a scan curve.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ScanObservable {
    ComponentCount,
    ActiveEdgeCount,
    LargestComponentSize,
    /// `largest_component_size / total_vertices`; undefined at `V = 0`.
    GiantFraction,
    /// Raw `sum(s_i^2)`; reweighting converts the exact `u128` to `f64`.
    RawSecondMoment,
}

impl ScanObservable {
    pub const fn name(self) -> &'static str {
        match self {
            Self::ComponentCount => "ComponentCount",
            Self::ActiveEdgeCount => "ActiveEdgeCount",
            Self::LargestComponentSize => "LargestComponentSize",
            Self::GiantFraction => "GiantFraction",
            Self::RawSecondMoment => "RawSecondMoment",
        }
    }

    fn value(self, point: &ScanPoint, total_vertices: usize) -> f64 {
        match self {
            Self::ComponentCount => point.component_count as f64,
            Self::ActiveEdgeCount => point.active_edge_count as f64,
            Self::LargestComponentSize => point.largest_component_size as f64,
            Self::GiantFraction => point.largest_component_size as f64 / total_vertices as f64,
            Self::RawSecondMoment => point.raw_second_moment as f64,
        }
    }
}

/// Reusable buffer for binomial reweighting.
///
/// Call [`Self::prepare`] once for a curve length; repeated reweighting at
/// different probabilities then reuses the mode-relative weight buffer
/// without allocating.
#[derive(Clone, Debug, Default)]
pub struct ReweightWorkspace {
    weights: Vec<f64>,
}

impl ReweightWorkspace {
    pub const fn new() -> Self {
        Self {
            weights: Vec::new(),
        }
    }

    /// Reserve and initialize normalized weights for a curve with `trials`
    /// activation steps (`trials + 1` points).
    pub fn prepare(&mut self, trials: usize) -> Result<(), TryReserveError> {
        let len = trials.saturating_add(1);
        reserve_len(&mut self.weights, len)?;
        self.weights.resize(len, 0.0);
        Ok(())
    }

    #[inline]
    pub fn weight_capacity(&self) -> usize {
        self.weights.capacity()
    }

    /// Audit-only retained-capacity accounting for the allocation probe.
    #[cfg(feature = "allocation-probe")]
    #[doc(hidden)]
    pub fn capacity_audit(&self) -> ReweightCapacityAudit {
        ReweightCapacityAudit {
            weights: vector_capacity_bytes(&self.weights),
        }
    }

    fn fill_weights(&mut self, trials: usize, probability: f64) {
        // Start at the binomial mode and recurse outward. The mode weight is
        // the largest, so all relative weights are at most one and the common
        // normalization cancels in the weighted mean. This avoids one `exp`
        // per occupation count.
        let complement = 1.0 - probability;
        let mode = ((trials.saturating_add(1) as f64) * probability).floor() as usize;
        let mode = mode.min(trials);
        self.weights[mode] = 1.0;

        let upward = probability / complement;
        let mut weight = 1.0;
        for occupation in mode..trials {
            weight *= (trials - occupation) as f64 / (occupation + 1) as f64 * upward;
            self.weights[occupation + 1] = weight;
        }

        let downward = complement / probability;
        let mut weight = 1.0;
        for occupation in (1..=mode).rev() {
            weight *= occupation as f64 / (trials - occupation + 1) as f64 * downward;
            self.weights[occupation - 1] = weight;
        }
    }
}

/// Reweight one scan curve to a uniform Bernoulli probability.
///
/// A uniformly random activation order makes the size-`k` prefix uniformly
/// distributed over all `C(n,k)` subsets, so the binomial weights are exact
/// for site scans (`n = V`) and bond scans (`n = E`).
///
/// The returned scalar is the weighted mean of the requested observable over
/// that one curve. Its scan-to-scan uncertainty is estimated by feeding
/// independent scans into [`ScanEstimator`]; the binomial weights themselves
/// are not the Monte Carlo error.
pub fn reweight(
    curve: &ScanCurve<'_>,
    probability: Probability,
    observable: ScanObservable,
    workspace: &mut ReweightWorkspace,
) -> Result<f64, ReweightError> {
    let trials = curve.trials();
    let points = curve.points();
    let expected = trials.saturating_add(1);
    if points.len() != expected {
        return Err(ReweightError::CurveLength {
            expected,
            actual: points.len(),
        });
    }
    let total_vertices = curve.total_vertices();
    if observable == ScanObservable::GiantFraction && total_vertices == 0 {
        return Err(ReweightError::UndefinedObservable(observable));
    }
    if observable == ScanObservable::RawSecondMoment {
        let vertices = total_vertices as u128;
        if vertices.saturating_mul(vertices) > (1_u128 << 53) {
            return Err(ReweightError::InexactMomentConversion { total_vertices });
        }
    }

    let probability = probability.get();
    if probability == 0.0 {
        return Ok(observable.value(&points[0], total_vertices));
    }
    if probability == 1.0 {
        return Ok(observable.value(&points[trials], total_vertices));
    }

    workspace.prepare(trials).map_err(ReweightError::Capacity)?;
    workspace.fill_weights(trials, probability);
    let mut weight_sum = 0.0;
    let mut weighted_sum = 0.0;
    for (occupation, &weight) in workspace.weights[..=trials].iter().enumerate() {
        weight_sum += weight;
        weighted_sum += weight * observable.value(&points[occupation], total_vertices);
    }
    Ok(weighted_sum / weight_sum)
}

/// Accumulates independent per-scan reweighted values into a sample estimate.
///
/// `variance` and `standard_error` describe the scan-to-scan Monte Carlo
/// uncertainty; they are not the binomial occupation-count variance.
#[derive(Clone, Copy, Debug, Default)]
pub struct ScanEstimator {
    count: usize,
    sum: f64,
    sum_sq: f64,
}

impl ScanEstimator {
    pub const fn new() -> Self {
        Self {
            count: 0,
            sum: 0.0,
            sum_sq: 0.0,
        }
    }

    pub const fn samples(&self) -> usize {
        self.count
    }

    pub fn add(&mut self, value: f64) {
        self.count += 1;
        self.sum += value;
        self.sum_sq += value * value;
    }

    /// Sample mean, unbiased sample variance and standard error, or `None`
    /// until at least two independent scans have been added.
    pub fn estimate(&self) -> Option<ReweightedEstimate> {
        if self.count < 2 {
            return None;
        }
        let samples = self.count as f64;
        let mean = self.sum / samples;
        let variance = ((self.sum_sq - self.sum * self.sum / samples) / (samples - 1.0)).max(0.0);
        Some(ReweightedEstimate {
            mean,
            variance,
            standard_error: (variance / samples).sqrt(),
            samples: self.count,
        })
    }
}

/// Sample estimate over independent Newman-Ziff scans.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReweightedEstimate {
    pub mean: f64,
    pub variance: f64,
    pub standard_error: f64,
    pub samples: usize,
}

/// Failure of scan-curve reweighting.
#[derive(Debug)]
pub enum ReweightError {
    Capacity(TryReserveError),
    CurveLength { expected: usize, actual: usize },
    UndefinedObservable(ScanObservable),
    InexactMomentConversion { total_vertices: usize },
}

impl fmt::Display for ReweightError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Capacity(_) => write!(formatter, "could not reserve reweighting workspace"),
            Self::CurveLength { expected, actual } => write!(
                formatter,
                "scan curve has {actual} points; expected {expected}"
            ),
            Self::UndefinedObservable(observable) => write!(
                formatter,
                "observable {} is undefined for this scan curve",
                observable.name()
            ),
            Self::InexactMomentConversion { total_vertices } => write!(
                formatter,
                "raw second moment for {total_vertices} vertices exceeds exact f64 integer range"
            ),
        }
    }
}

impl std::error::Error for ReweightError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Capacity(source) => Some(source),
            _ => None,
        }
    }
}

#[cfg(feature = "allocation-probe")]
#[doc(hidden)]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ReweightCapacityAudit {
    pub weights: usize,
}

#[cfg(feature = "allocation-probe")]
impl ReweightCapacityAudit {
    pub const fn total_bytes(self) -> usize {
        self.weights
    }
}

fn reserve_len<T>(storage: &mut Vec<T>, len: usize) -> Result<(), TryReserveError> {
    if storage.capacity() < len {
        storage.try_reserve_exact(len - storage.len())?;
    }
    Ok(())
}

#[cfg(feature = "allocation-probe")]
fn vector_capacity_bytes<T>(storage: &Vec<T>) -> usize {
    storage.capacity().saturating_mul(size_of::<T>())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::percolation::scan::{ScanMode, ScanPoint};

    fn point(component_count: usize) -> ScanPoint {
        ScanPoint {
            occupation: 0,
            active_edge_count: 0,
            component_count,
            largest_component_size: component_count,
            raw_second_moment: (component_count as u128) * (component_count as u128),
        }
    }

    fn scan_curve<'a>(points: &'a [ScanPoint], total_vertices: usize) -> ScanCurve<'a> {
        ScanCurve::from_parts(
            ScanMode::Bond,
            points.len() - 1,
            total_vertices,
            &[],
            points,
        )
    }

    #[test]
    fn binomial_reweighting_matches_hand_values() {
        let points = [point(0), point(1), point(4)];
        let curve = scan_curve(&points, 2);
        let mut workspace = ReweightWorkspace::new();

        let half = reweight(
            &curve,
            Probability::new(0.5).unwrap(),
            ScanObservable::ComponentCount,
            &mut workspace,
        )
        .unwrap();
        assert!((half - 1.5).abs() < 1e-12);

        let zero = reweight(
            &curve,
            Probability::new(0.0).unwrap(),
            ScanObservable::ComponentCount,
            &mut workspace,
        )
        .unwrap();
        assert_eq!(zero, 0.0);

        let one = reweight(
            &curve,
            Probability::new(1.0).unwrap(),
            ScanObservable::ComponentCount,
            &mut workspace,
        )
        .unwrap();
        assert_eq!(one, 4.0);
    }

    #[test]
    fn giant_fraction_uses_total_vertices_and_rejects_empty_graph() {
        let points = [point(0), point(1), point(2)];
        let curve = scan_curve(&points, 2);
        let mut workspace = ReweightWorkspace::new();
        let fraction = reweight(
            &curve,
            Probability::new(0.5).unwrap(),
            ScanObservable::GiantFraction,
            &mut workspace,
        )
        .unwrap();
        assert!((fraction - 0.5).abs() < 1e-12);

        let empty = scan_curve(&points[..1], 0);
        let error = reweight(
            &empty,
            Probability::new(0.5).unwrap(),
            ScanObservable::GiantFraction,
            &mut workspace,
        )
        .unwrap_err();
        assert!(matches!(
            error,
            ReweightError::UndefinedObservable(ScanObservable::GiantFraction)
        ));
    }

    #[test]
    fn raw_second_moment_rejects_inexact_f64_conversion() {
        let points = [point(0)];
        let curve = scan_curve(&points, 100_000_000);
        let mut workspace = ReweightWorkspace::new();
        let error = reweight(
            &curve,
            Probability::new(0.5).unwrap(),
            ScanObservable::RawSecondMoment,
            &mut workspace,
        )
        .unwrap_err();
        assert!(matches!(
            error,
            ReweightError::InexactMomentConversion {
                total_vertices: 100_000_000
            }
        ));
    }

    #[test]
    fn reweight_rejects_curve_length_mismatch() {
        let points = [point(0)];
        let curve = ScanCurve::from_parts(ScanMode::Bond, 2, 2, &[], &points);
        let mut workspace = ReweightWorkspace::new();
        let error = reweight(
            &curve,
            Probability::new(0.5).unwrap(),
            ScanObservable::ComponentCount,
            &mut workspace,
        )
        .unwrap_err();
        assert!(matches!(
            error,
            ReweightError::CurveLength {
                expected: 3,
                actual: 1
            }
        ));
    }

    #[test]
    fn estimator_reports_sample_variance_and_standard_error() {
        let mut estimator = ScanEstimator::new();
        assert!(estimator.estimate().is_none());
        estimator.add(1.0);
        assert!(estimator.estimate().is_none());
        estimator.add(3.0);
        let estimate = estimator.estimate().unwrap();
        assert_eq!(estimate.mean, 2.0);
        assert_eq!(estimate.variance, 2.0);
        assert_eq!(estimate.standard_error, 1.0);
        assert_eq!(estimate.samples, 2);
    }

    #[test]
    fn workspace_is_reusable_across_probabilities() {
        let points = [point(0), point(1), point(4)];
        let curve = scan_curve(&points, 2);
        let mut workspace = ReweightWorkspace::new();
        for probability in [0.5, 0.25, 0.75, 0.5] {
            let mean = reweight(
                &curve,
                Probability::new(probability).unwrap(),
                ScanObservable::ComponentCount,
                &mut workspace,
            )
            .unwrap();
            assert!(mean.is_finite());
        }
        assert!(workspace.weight_capacity() >= 3);
    }
}
