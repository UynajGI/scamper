use core::fmt;
use std::collections::TryReserveError;
#[cfg(feature = "allocation-probe")]
use std::mem::size_of;

use crate::{Probability, ProbabilityField};

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

/// Reusable buffers for heterogeneous (Poisson-binomial) reweighting.
///
/// Call [`Self::prepare`] once per trial count; repeated reweighting at
/// different fields of that length then reuses the log-weight buffer and the
/// `ln(k)` table (the `log C(N,n)` recurrence consumes every integer up to
/// `N`, so the table is filled once instead of per call) without allocating.
/// Uniform fields delegate to the embedded [`ReweightWorkspace`].
#[derive(Clone, Debug, Default)]
pub struct HetReweightWorkspace {
    log_weights: Vec<f64>,
    log_integers: Vec<f64>,
    uniform: ReweightWorkspace,
    prepared_integers: usize,
}

impl HetReweightWorkspace {
    pub const fn new() -> Self {
        Self {
            log_weights: Vec::new(),
            log_integers: Vec::new(),
            uniform: ReweightWorkspace::new(),
            prepared_integers: 0,
        }
    }

    /// Reserve and size the buffers for a curve with `trials` activation
    /// steps, refilling the `ln(k)` table only when the trial count changes.
    pub fn prepare(&mut self, trials: usize) -> Result<(), TryReserveError> {
        let len = trials.saturating_add(1);
        reserve_len(&mut self.log_weights, len)?;
        self.log_weights.resize(len, 0.0);
        if self.prepared_integers != trials {
            reserve_len(&mut self.log_integers, trials)?;
            self.log_integers.resize(trials, 0.0);
            for (index, slot) in self.log_integers.iter_mut().enumerate() {
                *slot = ((index + 1) as f64).ln();
            }
            self.prepared_integers = trials;
        }
        Ok(())
    }

    #[inline]
    pub fn weight_capacity(&self) -> usize {
        self.log_weights.capacity()
    }

    /// Audit-only retained-capacity accounting for the allocation probe.
    #[cfg(feature = "allocation-probe")]
    #[doc(hidden)]
    pub fn capacity_audit(&self) -> HetReweightCapacityAudit {
        HetReweightCapacityAudit {
            log_weights: vector_capacity_bytes(&self.log_weights),
            log_integers: vector_capacity_bytes(&self.log_integers),
            uniform: self.uniform.capacity_audit().total_bytes(),
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
    let trials = validate_curve(curve, observable)?;
    let points = curve.points();
    let total_vertices = curve.total_vertices();

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

/// Curve and observable checks shared by the uniform and heterogeneous paths.
/// Returns the trial count on success.
fn validate_curve(
    curve: &ScanCurve<'_>,
    observable: ScanObservable,
) -> Result<usize, ReweightError> {
    let trials = curve.trials();
    let expected = trials.saturating_add(1);
    if curve.points().len() != expected {
        return Err(ReweightError::CurveLength {
            expected,
            actual: curve.points().len(),
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
    Ok(trials)
}

/// Per-scan heterogeneous reweighting output: the raw estimator plus the
/// exact-unit-mass diagnostic.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HetReweightResult {
    /// `sum_n C(N,n) w(A_n) f_n` over this one scan curve — an unbiased
    /// per-trajectory estimate of the heterogeneous-Bernoulli expectation of
    /// the observable. Feed independent scans into [`ScanEstimator`].
    pub value: f64,
    /// `sum_n C(N,n) w(A_n)` over this one scan curve. Its expectation over
    /// scans is exactly 1 (the total Poisson-binomial mass), so the scan
    /// average of `unit_mass` must agree with 1 within its standard error —
    /// a built-in correctness diagnostic for both the weights and the scan.
    pub unit_mass: f64,
}

/// Reweight one scan curve to a heterogeneous (Poisson-binomial) Bernoulli
/// field with one probability per universe element — vertices for site scans,
/// physical edges for bond scans.
///
/// # Derivation
///
/// A uniformly random activation order makes the size-`n` prefix `A_n` a
/// uniform `n`-subset of the `N`-element universe. For the heterogeneous
/// product measure with set weight `w(A) = prod_{i in A} p_i *
/// prod_{j not in A} (1 - p_j)`, splitting by subset size gives
///
/// ```text
/// E_het[f] = sum_n sum_{|A|=n} w(A) f(A)
///          = sum_n C(N,n) * E_unif,n[ w * f ],
/// ```
///
/// and `w(A_n) f_n` (the curve value at occupation `n`) is an unbiased
/// estimator of `E_unif,n[w f]` for every `n` at once, so
/// `sum_n C(N,n) w(A_n) f_n` is unbiased by linearity — the dependence
/// between prefixes of one permutation never enters. Taking `f = 1` exactly
/// turns the estimator into `sum_n C(N,n) w(A_n)`, whose expectation is the
/// total Poisson-binomial mass, exactly 1; that quantity is returned as
/// [`HetReweightResult::unit_mass`].
///
/// # Raw estimator, not a weighted average
///
/// Unlike the uniform [`reweight`], the returned `value` is the raw
/// (unnormalized) sum. For a heterogeneous field the normalization
/// `sum_n C(N,n) w(A_n)` fluctuates from scan to scan, so dividing by it
/// would turn the estimator into a biased ratio estimator; the uniform case
/// can normalize only because its weights depend on `n` alone, making the
/// denominator trajectory-independent. A uniform [`ProbabilityField`]
/// therefore delegates to [`reweight`] outright, with `unit_mass = 1` by
/// that normalization.
///
/// # Numerics
///
/// Everything runs in log space. The prefix log-likelihood starts at
/// `l_0 = sum_j ln(1 - p_j)`; entering element `e` adds
/// `ln p_e - ln(1 - p_e)` along `curve.activation_order()`.
/// `log C(N,n)` grows by the `O(N)` recurrence
/// `logC(N,n) = logC(N,n-1) + ln(N-n+1) - ln(n)` from `logC(N,0) = 0`
/// (no `lgamma`, no new dependencies). Terms are combined with a
/// running-max log-sum-exp rescale: `u_n = exp(l_n - l_max)` lies in
/// `(0, 1]`, so the internal accumulators can neither overflow nor — all
/// terms being positive — suffer catastrophic cancellation; relative error
/// stays at the accumulated-rounding scale of the log recurrences. Individual
/// `u_n` that underflow to zero are negligible against the running maximum by
/// construction. The final `exp(l_max)` restores the linear scale. Outside
/// the variance envelope below, the reachable failure mode is underflow: a
/// near-degenerate field can push the best single-scan log-mass below
/// `-ln(f64 range)`, making `exp(l_max)` — and hence `value` and
/// `unit_mass` — exactly `0.0` (observed at a half-0.01/half-0.99 field on
/// `N = 10^4`). The mirror side, saturation to infinity past
/// `ln(f64::MAX) ~ 709.8`, needs a positive excursion that large and has
/// not been observed even at `N = 2 * 10^5`.
///
/// # Variance caveat
///
/// The per-scan `value`s are heavy-tailed once the field spread times the
/// universe size is large: the log-weight of a size-`N/2` prefix fluctuates
/// across scans with standard deviation about
/// `sqrt((N/4) * Var(logit p))`, so scan values behave roughly lognormally
/// with that parameter and their mean converges far slower than a
/// light-tailed average. As a rule of thumb the reweighted scan path is
/// competitive while `sqrt((N/4) * Var(logit p))` stays of order one or
/// below; for spread fields on large universes (for example a
/// half-0.2/half-0.9 field at `N = 4096`) prefer the N1 direct heterogeneous
/// sampler ([`crate::SiteBernoulli`] and friends), which pays one draw per
/// element but keeps light-tailed estimates at any spread.
///
/// # Boundary probabilities
///
/// Entries equal to `0.0` or `1.0` are rejected with the offending index
/// ([`ReweightError::BoundaryProbability`]). A boundary entry makes the
/// target measure degenerate in that coordinate, and the honest coverage is
/// the N1 direct-sampling path, which handles deterministic entries exactly;
/// no universe-splitting is attempted here.
pub fn reweight_heterogeneous(
    curve: &ScanCurve<'_>,
    field: &ProbabilityField<'_>,
    observable: ScanObservable,
    workspace: &mut HetReweightWorkspace,
) -> Result<HetReweightResult, ReweightError> {
    let trials = validate_curve(curve, observable)?;

    if let Some(probability) = field.uniform() {
        let value = reweight(curve, probability, observable, &mut workspace.uniform)?;
        return Ok(HetReweightResult {
            value,
            unit_mass: 1.0,
        });
    }

    let actual = field.len().expect("heterogeneous fields carry a length");
    if actual != trials {
        return Err(ReweightError::ProbabilityFieldLength {
            expected: trials,
            actual,
        });
    }
    for index in 0..trials {
        let probability = field.at(index).get();
        if probability == 0.0 || probability == 1.0 {
            return Err(ReweightError::BoundaryProbability {
                index,
                value: probability,
            });
        }
    }

    workspace.prepare(trials).map_err(ReweightError::Capacity)?;
    let HetReweightWorkspace {
        log_weights,
        log_integers,
        ..
    } = workspace;
    let order = curve.activation_order();
    let points = curve.points();
    let total_vertices = curve.total_vertices();

    let mut log_likelihood = 0.0;
    for index in 0..trials {
        log_likelihood += (1.0 - field.at(index).get()).ln();
    }
    log_weights[0] = log_likelihood;

    let mut log_binomial = 0.0;
    for occupation in 1..=trials {
        let element = order[occupation - 1];
        let probability = field.at(element).get();
        // One ln per entry: the odds p/(1-p) collapse the two log terms of
        // the derivation into a single transcendental per element.
        log_likelihood += (probability / (1.0 - probability)).ln();
        log_binomial += log_integers[trials - occupation] - log_integers[occupation - 1];
        log_weights[occupation] = log_likelihood + log_binomial;
    }

    let log_max = log_weights
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let mut mass = 0.0;
    let mut weighted = 0.0;
    for (occupation, &log_weight) in log_weights.iter().enumerate() {
        let scaled = (log_weight - log_max).exp();
        mass += scaled;
        weighted += scaled * observable.value(&points[occupation], total_vertices);
    }
    let scale = log_max.exp();
    Ok(HetReweightResult {
        value: scale * weighted,
        unit_mass: scale * mass,
    })
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
    ProbabilityFieldLength { expected: usize, actual: usize },
    BoundaryProbability { index: usize, value: f64 },
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
            Self::ProbabilityFieldLength { expected, actual } => write!(
                formatter,
                "probability field has {actual} entries; expected {expected}"
            ),
            Self::BoundaryProbability { index, value } => write!(
                formatter,
                "probability field entry {index} is {value}; heterogeneous reweighting requires \
                 probabilities strictly between 0 and 1 (boundary fields belong to direct \
                 sampling)"
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

#[cfg(feature = "allocation-probe")]
#[doc(hidden)]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct HetReweightCapacityAudit {
    pub log_weights: usize,
    pub log_integers: usize,
    pub uniform: usize,
}

#[cfg(feature = "allocation-probe")]
impl HetReweightCapacityAudit {
    pub const fn total_bytes(self) -> usize {
        self.log_weights
            .saturating_add(self.log_integers)
            .saturating_add(self.uniform)
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

    fn het_curve<'a>(points: &'a [ScanPoint], order: &'a [usize]) -> ScanCurve<'a> {
        ScanCurve::from_parts(
            ScanMode::Site,
            points.len() - 1,
            points.len(),
            order,
            points,
        )
    }

    fn heterogeneous_field(values: &[f64]) -> ProbabilityField<'_> {
        ProbabilityField::Owned(
            values
                .iter()
                .map(|&value| Probability::new(value).expect("test probability"))
                .collect(),
        )
    }

    #[test]
    fn heterogeneous_hand_weights_match_both_activation_orders() {
        // Field [0.25, 0.75], observable values [0, 1, 4].
        // Order [1, 0]: unit mass 0.1875 + 1.125 + 0.1875 = 1.5 and
        // value 0.1875*0 + 1.125*1 + 0.1875*4 = 1.875.
        // Order [0, 1]: unit mass 0.5, value 0.875. Both orders average to
        // the exact heterogeneous mean 1.375 = 0.0625 + 0.5625 + 0.1875*4.
        let points = [point(0), point(1), point(4)];
        let field = heterogeneous_field(&[0.25, 0.75]);
        let mut workspace = HetReweightWorkspace::new();

        let first = het_curve(&points, &[1, 0]);
        let result = reweight_heterogeneous(
            &first,
            &field,
            ScanObservable::ComponentCount,
            &mut workspace,
        )
        .unwrap();
        assert!((result.unit_mass - 1.5).abs() < 1e-12);
        assert!((result.value - 1.875).abs() < 1e-12);

        let second = het_curve(&points, &[0, 1]);
        let result = reweight_heterogeneous(
            &second,
            &field,
            ScanObservable::ComponentCount,
            &mut workspace,
        )
        .unwrap();
        assert!((result.unit_mass - 0.5).abs() < 1e-12);
        assert!((result.value - 0.875).abs() < 1e-12);
        assert!((((result.value + 1.875) / 2.0) - 1.375).abs() < 1e-12);
    }

    #[test]
    fn uniform_field_delegates_to_binomial_reweight() {
        let points = [point(0), point(1), point(4), point(5)];
        let curve = scan_curve(&points, 3);
        let mut workspace = HetReweightWorkspace::new();
        let mut reference = ReweightWorkspace::new();
        for probability in [0.1, 0.5, 0.5927, 0.9, 0.0, 1.0] {
            let probability = Probability::new(probability).unwrap();
            let het = reweight_heterogeneous(
                &curve,
                &probability.into(),
                ScanObservable::ComponentCount,
                &mut workspace,
            )
            .unwrap();
            let binomial = reweight(
                &curve,
                probability,
                ScanObservable::ComponentCount,
                &mut reference,
            )
            .unwrap();
            assert_eq!(het.value, binomial, "p = {}", probability.get());
            assert_eq!(het.unit_mass, 1.0);
        }
    }

    #[test]
    fn constant_heterogeneous_field_matches_binomial_reweight() {
        let points: Vec<_> = (0..=10).map(|n| point(n * n)).collect();
        let order: Vec<_> = (0..10).collect();
        let curve = het_curve(&points, &order);
        let mut workspace = HetReweightWorkspace::new();
        let mut reference = ReweightWorkspace::new();
        for probability in [0.1, 0.25, 0.5, 0.5927, 0.75, 0.9] {
            let values = vec![probability; 10];
            let field = heterogeneous_field(&values);
            let het = reweight_heterogeneous(
                &curve,
                &field,
                ScanObservable::ComponentCount,
                &mut workspace,
            )
            .unwrap();
            let binomial = reweight(
                &curve,
                Probability::new(probability).unwrap(),
                ScanObservable::ComponentCount,
                &mut reference,
            )
            .unwrap();
            // Both estimate the same weighted mean; only the log-space
            // rounding of the Poisson-binomial path differs.
            let relative = ((het.value - binomial) / binomial).abs();
            assert!(relative < 1e-12, "p = {probability}: {relative}");
            assert!((het.unit_mass - 1.0).abs() < 1e-12, "p = {probability}");
        }
    }

    #[test]
    fn heterogeneous_rejects_length_mismatch_and_boundary_entries() {
        let points = [point(0), point(1), point(4), point(5)];
        let curve = het_curve(&points, &[0, 1, 2]);
        let mut workspace = HetReweightWorkspace::new();

        let short = heterogeneous_field(&[0.5, 0.5]);
        assert!(matches!(
            reweight_heterogeneous(
                &curve,
                &short,
                ScanObservable::ComponentCount,
                &mut workspace
            )
            .unwrap_err(),
            ReweightError::ProbabilityFieldLength {
                expected: 3,
                actual: 2
            }
        ));

        let zero = heterogeneous_field(&[0.5, 0.0, 0.5]);
        assert!(matches!(
            reweight_heterogeneous(
                &curve,
                &zero,
                ScanObservable::ComponentCount,
                &mut workspace
            )
            .unwrap_err(),
            ReweightError::BoundaryProbability {
                index: 1,
                value: 0.0
            }
        ));

        let one = heterogeneous_field(&[0.5, 0.5, 1.0]);
        assert!(matches!(
            reweight_heterogeneous(&curve, &one, ScanObservable::ComponentCount, &mut workspace)
                .unwrap_err(),
            ReweightError::BoundaryProbability {
                index: 2,
                value: 1.0
            }
        ));
    }

    #[test]
    fn heterogeneous_workspace_reuses_across_fields_and_growing_lengths() {
        let mut workspace = HetReweightWorkspace::new();
        let small = [point(0), point(1), point(4)];
        let curve = het_curve(&small, &[1, 0]);
        for values in [[0.25, 0.75], [0.6, 0.4]] {
            let field = heterogeneous_field(&values);
            let result = reweight_heterogeneous(
                &curve,
                &field,
                ScanObservable::ComponentCount,
                &mut workspace,
            )
            .unwrap();
            assert!(result.value.is_finite() && result.unit_mass.is_finite());
        }

        let large: Vec<_> = (0..=6).map(point).collect();
        let large_order: Vec<_> = (0..6).collect();
        let large_curve = het_curve(&large, &large_order);
        let field = heterogeneous_field(&[0.3, 0.4, 0.5, 0.6, 0.7, 0.8]);
        let grown = reweight_heterogeneous(
            &large_curve,
            &field,
            ScanObservable::ComponentCount,
            &mut workspace,
        )
        .unwrap();
        let mut fresh = HetReweightWorkspace::new();
        let reference = reweight_heterogeneous(
            &large_curve,
            &field,
            ScanObservable::ComponentCount,
            &mut fresh,
        )
        .unwrap();
        assert_eq!(grown.value, reference.value);
        assert_eq!(grown.unit_mass, reference.unit_mass);
        assert!(workspace.weight_capacity() >= 7);
    }

    #[test]
    fn heterogeneous_empty_curve_matches_first_point() {
        let points = [point(3)];
        let curve = het_curve(&points, &[]);
        let field = heterogeneous_field(&[]);
        let mut workspace = HetReweightWorkspace::new();
        let result = reweight_heterogeneous(
            &curve,
            &field,
            ScanObservable::ComponentCount,
            &mut workspace,
        )
        .unwrap();
        assert_eq!(result.value, 3.0);
        assert_eq!(result.unit_mass, 1.0);
    }
}
