use core::fmt;
use std::sync::Arc;

use carlo_rs::Context;
use rand_xoshiro::Xoshiro256PlusPlus;

use super::{BoundaryQuery, ComponentLabel, ComponentSummary, WindingSummary};

/// Scalar observables recorded for each independent static realization.
///
/// All counts are per-sample integers recorded as dimensionless f64 scalars;
/// the ensemble mean of each name is the reported estimate.
///
/// - `ActiveVertexCount` / `ActiveEdgeCount`: number of occupied sites /
///   occupied physical edges in the configuration.
/// - `ComponentCount`: number of tracked components among occupied sites.
/// - `LargestSize`: `S_max`, the largest component size (`0` when nothing is
///   occupied).
/// - `GiantFraction`: `S_max / V` with the denominator fixed to the total
///   topology vertex count `V`, occupied or not. At `V = 0` the observable is
///   undefined and silently absent from that sample's Results row.
/// - `RawSecondMoment`: `M2 = sum_i s_i^2` over all tracked components,
///   including the largest one, with no normalization.
/// - `FiniteClusterSusceptibility`: per sample, one largest component (equal
///   sizes: the canonical lowest-ID one) is excluded, then
///   `chi = (M2 - S_max^2) / (V_active - S_max)`. The exclusion happens before
///   averaging, so this cannot be reconstructed from separately aggregated
///   means. When `V_active = S_max` the denominator is zero and the sample is
///   undefined: the chi scalar is not written for that sample (its Results row
///   therefore has `n_bins` below the sweep count), while the
///   `FiniteClusterSusceptibilityDefined` indicator records `0`/`1` for every
///   sample.
/// - `Wrapping`: `1` when at least one tracked component has winding rank
///   greater than zero, `0` otherwise (W2). Requires an embedded (periodic)
///   topology; the adapter rejects winding observables without one.
/// - `WindingRank`: winding rank of the canonical largest component, an
///   integer in `0..=D` (W2). The mean is the ensemble-averaged rank, not a
///   probability; per-axis generator vectors remain analyzer-level output
///   until Carlo.rs grows vector estimates (C0).
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub enum StaticObservable {
    ActiveVertexCount,
    ActiveEdgeCount,
    ComponentCount,
    LargestSize,
    GiantFraction,
    RawSecondMoment,
    FiniteClusterSusceptibility,
    Wrapping,
    WindingRank,
}

impl StaticObservable {
    pub const fn name(self) -> &'static str {
        match self {
            Self::ActiveVertexCount => "ActiveVertexCount",
            Self::ActiveEdgeCount => "ActiveEdgeCount",
            Self::ComponentCount => "ComponentCount",
            Self::LargestSize => "LargestSize",
            Self::GiantFraction => "GiantFraction",
            Self::RawSecondMoment => "RawSecondMoment",
            Self::FiniteClusterSusceptibility => "FiniteClusterSusceptibility",
            Self::Wrapping => "Wrapping",
            Self::WindingRank => "WindingRank",
        }
    }

    /// Whether recording this observable requires a periodic embedding (W2).
    pub(crate) const fn needs_embedding(self) -> bool {
        matches!(self, Self::Wrapping | Self::WindingRank)
    }
}

/// One realized configuration as seen by an observable.
///
/// `total_vertices` is the topology size (occupied or not), which the
/// built-in `GiantFraction` uses as its denominator. `winding` is present
/// only when the adapter computed a winding summary; plans that request
/// winding are rejected at construction without an embedding, so a custom
/// observable declaring [`ScalarObservable::requires_winding`] always sees `Some`
/// when it is evaluated. `labels` is present only when the plan or adapter
/// asked the analyzer to materialize canonical per-vertex labels.
#[derive(Clone, Copy, Debug)]
pub struct SampleView<'a> {
    total_vertices: usize,
    summary: ComponentSummary,
    boundary_outcomes: &'a [bool],
    winding: Option<WindingSummary>,
    labels: Option<&'a [ComponentLabel]>,
}

impl<'a> SampleView<'a> {
    pub(crate) const fn new(
        total_vertices: usize,
        summary: ComponentSummary,
        boundary_outcomes: &'a [bool],
        winding: Option<WindingSummary>,
        labels: Option<&'a [ComponentLabel]>,
    ) -> Self {
        Self {
            total_vertices,
            summary,
            boundary_outcomes,
            winding,
            labels,
        }
    }

    /// Zero-vertex realization for adapter vacuum paths that skip analysis.
    pub(crate) const fn vacuum(boundary_outcomes: &'a [bool]) -> Self {
        Self::new(0, ComponentSummary::EMPTY, boundary_outcomes, None, None)
    }

    #[inline]
    pub const fn total_vertices(&self) -> usize {
        self.total_vertices
    }

    #[inline]
    pub const fn summary(&self) -> ComponentSummary {
        self.summary
    }

    #[inline]
    pub const fn boundary_outcomes(&self) -> &'a [bool] {
        self.boundary_outcomes
    }

    #[inline]
    pub const fn winding(&self) -> Option<WindingSummary> {
        self.winding
    }

    #[inline]
    pub const fn labels(&self) -> Option<&'a [ComponentLabel]> {
        self.labels
    }
}

/// Open extension point for a user-defined scalar observable.
///
/// The built-in [`StaticObservable`] set stays the zero-cost default path;
/// this trait is how downstream crates add a column without editing CMC.rs.
/// Implementations are evaluated once per realized sample from the shared
/// [`SampleView`]. Returning `None` leaves the column absent for that
/// sample, matching the built-in undefined-observable convention (`n_bins`
/// below the sweep count).
///
/// Capability methods are declarations, not hints: the plan aggregates them
/// so adapters materialize exactly what the observables need. Label
/// materialization and winding scans are skipped otherwise.
///
/// Names are validated when the plan is assembled with
/// [`ObservablePlan::with_custom`]: they must be non-empty and must not
/// collide with a planned built-in, its indicator column, a boundary column,
/// or another custom observable.
///
/// ```
/// use std::sync::Arc;
/// use cmc_rs::{ScalarObservable, ObservablePlan, SampleView};
///
/// #[derive(Debug)]
/// struct Occupancy { name: &'static str }
///
/// impl ScalarObservable for Occupancy {
///     fn name(&self) -> &str { self.name }
///     fn evaluate(&self, sample: &SampleView<'_>) -> Option<f64> {
///         let vertices = sample.total_vertices();
///         (vertices != 0).then(|| {
///             sample.summary().active_vertex_count as f64 / vertices as f64
///         })
///     }
/// }
///
/// let plan = ObservablePlan::all(Vec::new())
///     .with_custom([Arc::new(Occupancy { name: "Occupancy" }) as Arc<dyn ScalarObservable>])
///     .expect("unique name");
/// assert!(plan.has_custom());
/// ```
pub trait ScalarObservable: fmt::Debug + Send + Sync {
    /// Stable Results column name.
    fn name(&self) -> &str;

    /// Whether evaluation needs canonical per-vertex component labels.
    fn requires_labels(&self) -> bool {
        false
    }

    /// Whether the plan needs a periodic embedding and winding summary.
    fn requires_winding(&self) -> bool {
        false
    }

    /// Whether evaluation consumes moments that must stay exact in `f64`.
    fn requires_exact_moments(&self) -> bool {
        false
    }

    /// Evaluate one realization; `None` skips this sample for this column.
    fn evaluate(&self, sample: &SampleView<'_>) -> Option<f64>;
}

/// Emission protocol for one built-in plan entry.
///
/// A spec may record more than one Results column (the finite-cluster
/// susceptibility writes its presence indicator alongside the value), which
/// is why this is separate from the single-column [`ScalarObservable`].
/// Built-ins are statically dispatched: the plan stores
/// [`StaticObservable`], whose impl delegates to one zero-sized spec per
/// variant. Third-party observables use [`ScalarObservable`] instead.
pub trait ObservableSpec {
    fn name(&self) -> &'static str;

    fn needs_embedding(&self) -> bool {
        false
    }

    fn needs_exact_moments(&self) -> bool {
        false
    }

    fn indicator_names(&self) -> &'static [&'static str] {
        &[]
    }

    fn record(&self, sample: &SampleView<'_>, ctx: &mut Context<Xoshiro256PlusPlus>);
}

#[derive(Clone, Copy, Debug, Default)]
struct ActiveVertexCountSpec;

impl ObservableSpec for ActiveVertexCountSpec {
    fn name(&self) -> &'static str {
        StaticObservable::ActiveVertexCount.name()
    }

    fn record(&self, sample: &SampleView<'_>, ctx: &mut Context<Xoshiro256PlusPlus>) {
        ctx.measure(self.name(), sample.summary().active_vertex_count as f64);
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct ActiveEdgeCountSpec;

impl ObservableSpec for ActiveEdgeCountSpec {
    fn name(&self) -> &'static str {
        StaticObservable::ActiveEdgeCount.name()
    }

    fn record(&self, sample: &SampleView<'_>, ctx: &mut Context<Xoshiro256PlusPlus>) {
        ctx.measure(self.name(), sample.summary().active_edge_count as f64);
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct ComponentCountSpec;

impl ObservableSpec for ComponentCountSpec {
    fn name(&self) -> &'static str {
        StaticObservable::ComponentCount.name()
    }

    fn record(&self, sample: &SampleView<'_>, ctx: &mut Context<Xoshiro256PlusPlus>) {
        ctx.measure(self.name(), sample.summary().component_count as f64);
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct LargestSizeSpec;

impl ObservableSpec for LargestSizeSpec {
    fn name(&self) -> &'static str {
        StaticObservable::LargestSize.name()
    }

    fn record(&self, sample: &SampleView<'_>, ctx: &mut Context<Xoshiro256PlusPlus>) {
        ctx.measure(self.name(), sample.summary().largest_component_size as f64);
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct GiantFractionSpec;

impl ObservableSpec for GiantFractionSpec {
    fn name(&self) -> &'static str {
        StaticObservable::GiantFraction.name()
    }

    fn record(&self, sample: &SampleView<'_>, ctx: &mut Context<Xoshiro256PlusPlus>) {
        let total = sample.total_vertices();
        if total != 0 {
            ctx.measure(
                self.name(),
                sample.summary().largest_component_size as f64 / total as f64,
            );
        }
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct RawSecondMomentSpec;

impl ObservableSpec for RawSecondMomentSpec {
    fn name(&self) -> &'static str {
        StaticObservable::RawSecondMoment.name()
    }

    fn needs_exact_moments(&self) -> bool {
        true
    }

    fn record(&self, sample: &SampleView<'_>, ctx: &mut Context<Xoshiro256PlusPlus>) {
        ctx.measure(self.name(), sample.summary().raw_second_moment as f64);
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct FiniteClusterSusceptibilitySpec;

const SUSCEPTIBILITY_INDICATORS: &[&str] = &["FiniteClusterSusceptibilityDefined"];

impl ObservableSpec for FiniteClusterSusceptibilitySpec {
    fn name(&self) -> &'static str {
        StaticObservable::FiniteClusterSusceptibility.name()
    }

    fn needs_exact_moments(&self) -> bool {
        true
    }

    fn indicator_names(&self) -> &'static [&'static str] {
        SUSCEPTIBILITY_INDICATORS
    }

    fn record(&self, sample: &SampleView<'_>, ctx: &mut Context<Xoshiro256PlusPlus>) {
        let summary = sample.summary();
        let finite_vertices = summary
            .active_vertex_count
            .saturating_sub(summary.largest_component_size);
        let largest_squared =
            (summary.largest_component_size as u128) * (summary.largest_component_size as u128);
        let finite_second_moment = summary.raw_second_moment.saturating_sub(largest_squared);
        let defined = finite_vertices != 0;
        ctx.measure(SUSCEPTIBILITY_INDICATORS[0], f64::from(u8::from(defined)));
        if defined {
            ctx.measure(
                self.name(),
                finite_second_moment as f64 / finite_vertices as f64,
            );
        }
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct WrappingSpec;

impl ObservableSpec for WrappingSpec {
    fn name(&self) -> &'static str {
        StaticObservable::Wrapping.name()
    }

    fn needs_embedding(&self) -> bool {
        true
    }

    fn record(&self, sample: &SampleView<'_>, ctx: &mut Context<Xoshiro256PlusPlus>) {
        debug_assert!(
            sample.winding().is_some(),
            "plan validation pairs winding observables with an embedding"
        );
        if let Some(winding) = sample.winding() {
            ctx.measure(
                self.name(),
                f64::from(u8::from(winding.wrapping_component_count > 0)),
            );
        }
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct WindingRankSpec;

impl ObservableSpec for WindingRankSpec {
    fn name(&self) -> &'static str {
        StaticObservable::WindingRank.name()
    }

    fn needs_embedding(&self) -> bool {
        true
    }

    fn record(&self, sample: &SampleView<'_>, ctx: &mut Context<Xoshiro256PlusPlus>) {
        debug_assert!(
            sample.winding().is_some(),
            "plan validation pairs winding observables with an embedding"
        );
        if let Some(winding) = sample.winding() {
            ctx.measure(self.name(), winding.largest_winding_rank as f64);
        }
    }
}

impl ObservableSpec for StaticObservable {
    fn name(&self) -> &'static str {
        StaticObservable::name(*self)
    }

    fn needs_embedding(&self) -> bool {
        StaticObservable::needs_embedding(*self)
    }

    fn needs_exact_moments(&self) -> bool {
        matches!(
            self,
            Self::RawSecondMoment | Self::FiniteClusterSusceptibility
        )
    }

    fn indicator_names(&self) -> &'static [&'static str] {
        match self {
            Self::FiniteClusterSusceptibility => SUSCEPTIBILITY_INDICATORS,
            _ => &[],
        }
    }

    fn record(&self, sample: &SampleView<'_>, ctx: &mut Context<Xoshiro256PlusPlus>) {
        match self {
            Self::ActiveVertexCount => ActiveVertexCountSpec.record(sample, ctx),
            Self::ActiveEdgeCount => ActiveEdgeCountSpec.record(sample, ctx),
            Self::ComponentCount => ComponentCountSpec.record(sample, ctx),
            Self::LargestSize => LargestSizeSpec.record(sample, ctx),
            Self::GiantFraction => GiantFractionSpec.record(sample, ctx),
            Self::RawSecondMoment => RawSecondMomentSpec.record(sample, ctx),
            Self::FiniteClusterSusceptibility => {
                FiniteClusterSusceptibilitySpec.record(sample, ctx);
            }
            Self::Wrapping => WrappingSpec.record(sample, ctx),
            Self::WindingRank => WindingRankSpec.record(sample, ctx),
        }
    }
}

/// Explicit requests for scalar component observables and boundary events.
///
/// Duplicate requests are collapsed, so every planned observable records
/// exactly one sample per realization regardless of request order. Boundary
/// outcomes are measured in plan order as `BoundaryCrossing0`,
/// `BoundaryCrossing1`, and so on. Unlike scalars, duplicate boundary queries
/// are kept: each gets its own numbered name (identical duplicate queries
/// therefore yield identically-valued rows with independent `n_bins`). Labels
/// are intentionally not requested:
/// all listed scalars follow from [`ComponentSummary`], whose canonical
/// lowest-ID largest component proves the one-component tie policy.
#[derive(Clone, Debug)]
pub struct ObservablePlan {
    observables: Vec<StaticObservable>,
    custom: Vec<Arc<dyn ScalarObservable>>,
    boundary_queries: Vec<BoundaryQuery>,
    boundary_names: Vec<String>,
}

impl ObservablePlan {
    pub fn new(
        observables: impl IntoIterator<Item = StaticObservable>,
        boundary_queries: Vec<BoundaryQuery>,
    ) -> Self {
        // A full sort+dedup (not just consecutive dedup) guarantees one
        // sample per name for interleaved duplicates such as [A, B, A]:
        // double recording would inflate n_bins and corrupt stderr while
        // leaving the mean seemingly unchanged.
        let mut observables = observables.into_iter().collect::<Vec<_>>();
        observables.sort_unstable();
        observables.dedup();
        let boundary_names = (0..boundary_queries.len())
            .map(|index| format!("BoundaryCrossing{index}"))
            .collect();
        Self {
            observables,
            custom: Vec::new(),
            boundary_queries,
            boundary_names,
        }
    }

    pub fn all(boundary_queries: Vec<BoundaryQuery>) -> Self {
        Self::new(
            [
                StaticObservable::ActiveVertexCount,
                StaticObservable::ActiveEdgeCount,
                StaticObservable::ComponentCount,
                StaticObservable::LargestSize,
                StaticObservable::GiantFraction,
                StaticObservable::RawSecondMoment,
                StaticObservable::FiniteClusterSusceptibility,
            ],
            boundary_queries,
        )
    }

    pub fn observables(&self) -> &[StaticObservable] {
        &self.observables
    }

    /// Add user-defined scalar observables to this plan.
    ///
    /// Names must be non-empty and unique across built-ins, the
    /// `FiniteClusterSusceptibilityDefined` indicator, boundary columns,
    /// and one another. The returned plan is otherwise unchanged, so a
    /// failed assembly leaves the caller's original plan untouched.
    pub fn with_custom(
        mut self,
        observables: impl IntoIterator<Item = Arc<dyn ScalarObservable>>,
    ) -> Result<Self, ObservablePlanError> {
        let custom: Vec<Arc<dyn ScalarObservable>> = observables.into_iter().collect();
        let reserved = self.reserved_names();
        let mut seen: Vec<String> = Vec::with_capacity(custom.len());
        for observable in &custom {
            let name = observable.name();
            if name.is_empty() {
                return Err(ObservablePlanError::EmptyCustomName);
            }
            if seen.iter().any(|existing| existing == name) {
                return Err(ObservablePlanError::DuplicateCustomName {
                    name: name.to_owned(),
                });
            }
            if reserved.contains(&name) {
                return Err(ObservablePlanError::CustomNameCollision {
                    name: name.to_owned(),
                });
            }
            seen.push(name.to_owned());
        }
        self.custom = custom;
        Ok(self)
    }

    /// The user-defined observables attached to this plan.
    pub fn custom_observables(&self) -> &[Arc<dyn ScalarObservable>] {
        &self.custom
    }

    /// Whether any user-defined observable is attached.
    pub fn has_custom(&self) -> bool {
        !self.custom.is_empty()
    }

    fn reserved_names(&self) -> Vec<&str> {
        let mut names: Vec<&str> = self.observables.iter().map(|value| value.name()).collect();
        for observable in &self.observables {
            names.extend(observable.indicator_names());
        }
        names.extend(self.boundary_names.iter().map(String::as_str));
        names
    }

    pub fn boundary_queries(&self) -> &[BoundaryQuery] {
        &self.boundary_queries
    }

    /// Whether any planned observable requires a periodic embedding (W2):
    /// the caller must supply one or reject the plan loudly.
    pub fn requests_winding(&self) -> bool {
        self.observables
            .iter()
            .any(|observable| observable.needs_embedding())
            || self
                .custom
                .iter()
                .any(|observable| observable.requires_winding())
    }

    /// Whether any user-defined observable needs canonical component labels.
    /// Adapters pass this to `analyze_with_labels`, so the label array is
    /// materialized only on request.
    pub fn requests_labels(&self) -> bool {
        self.custom
            .iter()
            .any(|observable| observable.requires_labels())
    }

    /// Reject moment observables whose worst-case `u128 -> f64` conversion
    /// could silently round, i.e. topologies with `V^2 > 2^53`.
    pub fn validate(&self, total_vertices: usize) -> Result<(), ObservablePlanError> {
        // Both moment-based observables convert u128 cluster-size squares of
        // magnitude up to V^2 into f64; beyond 2^53 that conversion would
        // silently round, so the plan is rejected at construction instead.
        let needs_exact_moments = self
            .observables
            .iter()
            .any(|observable| observable.needs_exact_moments())
            || self
                .custom
                .iter()
                .any(|observable| observable.requires_exact_moments());
        if needs_exact_moments {
            let max_moment = (total_vertices as u128) * (total_vertices as u128);
            if max_moment > (1_u128 << f64::MANTISSA_DIGITS) {
                return Err(ObservablePlanError::InexactMomentConversion { total_vertices });
            }
        }
        Ok(())
    }

    /// Record every planned column for one realization.
    ///
    /// Built-ins are emitted first in their canonical enum order, then
    /// user-defined observables in plan order. Plans without customs
    /// therefore keep the exact v2 recording sequence.
    pub(crate) fn record(&self, sample: &SampleView<'_>, ctx: &mut Context<Xoshiro256PlusPlus>) {
        for observable in &self.observables {
            observable.record(sample, ctx);
        }
        for (name, &crosses) in self.boundary_names.iter().zip(sample.boundary_outcomes()) {
            ctx.measure(name, f64::from(u8::from(crosses)));
        }
        for observable in &self.custom {
            if let Some(value) = observable.evaluate(sample) {
                ctx.measure(observable.name(), value);
            }
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ObservablePlanError {
    InexactMomentConversion {
        total_vertices: usize,
    },
    /// A custom observable returned an empty column name.
    EmptyCustomName,
    /// Two custom observables use the same column name.
    DuplicateCustomName {
        name: String,
    },
    /// A custom observable collides with a built-in or boundary column.
    CustomNameCollision {
        name: String,
    },
}

impl fmt::Display for ObservablePlanError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InexactMomentConversion { total_vertices } => write!(
                formatter,
                "moment observables may exceed exact f64 integer range for V={total_vertices}"
            ),
            Self::EmptyCustomName => write!(formatter, "custom observable names must not be empty"),
            Self::DuplicateCustomName { name } => {
                write!(formatter, "duplicate custom observable name {name:?}")
            }
            Self::CustomNameCollision { name } => write!(
                formatter,
                "custom observable name {name:?} collides with a built-in or boundary column"
            ),
        }
    }
}

impl std::error::Error for ObservablePlanError {}

#[cfg(test)]
mod tests {
    use super::*;
    use carlo_rs::Context;
    use rand::SeedableRng;

    use crate::build_chain;

    fn summary(active: usize, largest: usize, second: u128) -> ComponentSummary {
        ComponentSummary {
            active_vertex_count: active,
            active_edge_count: active.saturating_sub(1),
            component_count: usize::from(active != 0),
            largest_component: None,
            largest_component_size: largest,
            raw_second_moment: second,
        }
    }

    fn context() -> Context<Xoshiro256PlusPlus> {
        Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(1), 0, 1)
    }

    fn sample<'a>(
        total: usize,
        summary: ComponentSummary,
        outcomes: &'a [bool],
        winding: Option<WindingSummary>,
    ) -> SampleView<'a> {
        SampleView::new(total, summary, outcomes, winding, None)
    }

    #[test]
    fn records_units_total_vertex_normalization_and_multiple_queries() {
        let graph = build_chain(10, false);
        let plan = ObservablePlan::all(vec![
            BoundaryQuery::new(&graph, &[0], &[9]).unwrap(),
            BoundaryQuery::new(&graph, &[1], &[8]).unwrap(),
        ]);
        let mut ctx = context();
        plan.record(
            &sample(10, summary(6, 3, 14), &[true, false], None),
            &mut ctx,
        );
        let values = ctx.finalize_measurements();
        assert_eq!(values["ActiveVertexCount"].mean, 6.0);
        assert_eq!(values["ActiveEdgeCount"].mean, 5.0);
        assert_eq!(values["LargestSize"].mean, 3.0);
        assert_eq!(values["GiantFraction"].mean, 0.3);
        assert_eq!(values["RawSecondMoment"].mean, 14.0);
        assert_eq!(values["FiniteClusterSusceptibility"].mean, 5.0 / 3.0);
        assert_eq!(values["FiniteClusterSusceptibilityDefined"].mean, 1.0);
        assert_eq!(values["BoundaryCrossing0"].mean, 1.0);
        assert_eq!(values["BoundaryCrossing1"].mean, 0.0);
    }

    #[test]
    fn undefined_values_are_missing_not_nan_and_have_an_indicator() {
        let plan = ObservablePlan::new(
            [
                StaticObservable::GiantFraction,
                StaticObservable::FiniteClusterSusceptibility,
            ],
            Vec::new(),
        );
        let mut empty = context();
        plan.record(&sample(0, summary(0, 0, 0), &[], None), &mut empty);
        let values = empty.finalize_measurements();
        assert!(!values.contains_key("GiantFraction"));
        assert!(!values.contains_key("FiniteClusterSusceptibility"));
        assert_eq!(values["FiniteClusterSusceptibilityDefined"].mean, 0.0);

        let mut giant = context();
        plan.record(&sample(4, summary(4, 4, 16), &[], None), &mut giant);
        let values = giant.finalize_measurements();
        assert_eq!(values["GiantFraction"].mean, 1.0);
        assert!(!values.contains_key("FiniteClusterSusceptibility"));
        assert_eq!(values["FiniteClusterSusceptibilityDefined"].mean, 0.0);
    }

    #[test]
    fn susceptibility_excludes_exactly_one_largest_per_sample() {
        let plan = ObservablePlan::new([StaticObservable::FiniteClusterSusceptibility], Vec::new());
        let mut ctx = context();
        // Components [3, 3, 1]: canonical tie policy removes one size-3 cluster.
        plan.record(&sample(7, summary(7, 3, 19), &[], None), &mut ctx);
        assert_eq!(
            ctx.finalize_measurements()["FiniteClusterSusceptibility"].mean,
            2.5
        );
    }

    #[test]
    fn susceptibility_is_averaged_per_sample_not_as_ratio_of_means() {
        let plan = ObservablePlan::new([StaticObservable::FiniteClusterSusceptibility], Vec::new());
        let mut ctx = context();
        // Remaining clusters are [1] and [2, 2], giving chi 1 and 2.
        plan.record(&sample(4, summary(4, 3, 10), &[], None), &mut ctx);
        plan.record(&sample(7, summary(7, 3, 17), &[], None), &mut ctx);
        let per_sample_mean = ctx.finalize_measurements()["FiniteClusterSusceptibility"].mean;
        assert_eq!(per_sample_mean, 1.5);
        assert_ne!(per_sample_mean, (1.0 + 8.0) / (1.0 + 4.0));
    }

    #[test]
    fn winding_scalars_record_from_the_winding_summary() {
        let plan = ObservablePlan::new(
            [StaticObservable::Wrapping, StaticObservable::WindingRank],
            Vec::new(),
        );
        assert!(plan.requests_winding());
        let winding = WindingSummary {
            active_vertex_count: 4,
            active_edge_count: 8,
            component_count: 1,
            wrapping_component_count: 1,
            largest_component: None,
            largest_component_size: 4,
            largest_winding_rank: 2,
        };
        let contractible = WindingSummary {
            wrapping_component_count: 0,
            largest_winding_rank: 0,
            ..winding
        };
        let mut ctx = context();
        plan.record(&sample(4, summary(4, 4, 16), &[], Some(winding)), &mut ctx);
        plan.record(
            &sample(4, summary(4, 4, 16), &[], Some(contractible)),
            &mut ctx,
        );
        let values = ctx.finalize_measurements();
        assert_eq!(values["Wrapping"].mean, 0.5);
        assert_eq!(values["WindingRank"].mean, 1.0);
    }

    #[test]
    fn plain_plans_do_not_request_winding() {
        let plan = ObservablePlan::all(Vec::new());
        assert!(!plan.requests_winding());
    }

    #[test]
    fn interleaved_duplicate_requests_collapse_to_one_sample_per_name() {
        let plan = ObservablePlan::new(
            [
                StaticObservable::LargestSize,
                StaticObservable::ComponentCount,
                StaticObservable::LargestSize,
            ],
            Vec::new(),
        );
        assert_eq!(
            plan.observables(),
            &[
                StaticObservable::ComponentCount,
                StaticObservable::LargestSize
            ]
        );
        let mut ctx = context();
        plan.record(&sample(4, summary(4, 3, 11), &[], None), &mut ctx);
        let values = ctx.finalize_measurements();
        // One sample per sweep, not two: the mean is unchanged by duplicates
        // but n_bins and stderr would silently double without collapsing.
        assert_eq!(values["LargestSize"].n_bins, 1);
    }

    #[derive(Debug)]
    struct MomentConsumer;

    impl ScalarObservable for MomentConsumer {
        fn name(&self) -> &str {
            "MomentConsumer"
        }

        fn requires_exact_moments(&self) -> bool {
            true
        }

        fn evaluate(&self, sample: &SampleView<'_>) -> Option<f64> {
            Some(sample.summary().raw_second_moment as f64)
        }
    }

    #[test]
    fn custom_capabilities_feed_plan_validation_and_winding_requests() {
        let plan = ObservablePlan::new(Vec::<StaticObservable>::new(), Vec::new())
            .with_custom([
                std::sync::Arc::new(MomentConsumer) as std::sync::Arc<dyn ScalarObservable>
            ])
            .expect("plan");
        const EXACT_LIMIT: usize = 94_906_265;
        assert!(plan.validate(EXACT_LIMIT).is_ok());
        assert_eq!(
            plan.validate(EXACT_LIMIT + 1),
            Err(ObservablePlanError::InexactMomentConversion {
                total_vertices: EXACT_LIMIT + 1,
            }),
        );
    }

    #[test]
    fn moment_observables_require_worst_case_exact_f64_conversion() {
        const EXACT_LIMIT: usize = 94_906_265;
        let plan = ObservablePlan::new([StaticObservable::RawSecondMoment], Vec::new());
        assert!(plan.validate(EXACT_LIMIT).is_ok());
        assert_eq!(
            plan.validate(EXACT_LIMIT + 1),
            Err(ObservablePlanError::InexactMomentConversion {
                total_vertices: EXACT_LIMIT + 1,
            }),
        );
        // The susceptibility squares sizes in its numerator, so the same
        // worst-case bound guards its per-sample conversion.
        let susceptibility =
            ObservablePlan::new([StaticObservable::FiniteClusterSusceptibility], Vec::new());
        assert_eq!(
            susceptibility.validate(EXACT_LIMIT + 1),
            Err(ObservablePlanError::InexactMomentConversion {
                total_vertices: EXACT_LIMIT + 1,
            }),
        );
        assert!(
            ObservablePlan::new([StaticObservable::LargestSize], Vec::new())
                .validate(usize::MAX)
                .is_ok()
        );
    }
}
