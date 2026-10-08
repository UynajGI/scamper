use core::fmt;

use carlo_rs::Context;
use rand_xoshiro::Xoshiro256PlusPlus;

use super::{BoundaryQuery, ComponentSummary};

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
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub enum StaticObservable {
    ActiveVertexCount,
    ActiveEdgeCount,
    ComponentCount,
    LargestSize,
    GiantFraction,
    RawSecondMoment,
    FiniteClusterSusceptibility,
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

    pub fn boundary_queries(&self) -> &[BoundaryQuery] {
        &self.boundary_queries
    }

    /// Reject moment observables whose worst-case `u128 -> f64` conversion
    /// could silently round, i.e. topologies with `V^2 > 2^53`.
    pub fn validate(&self, total_vertices: usize) -> Result<(), ObservablePlanError> {
        // Both moment-based observables convert u128 cluster-size squares of
        // magnitude up to V^2 into f64; beyond 2^53 that conversion would
        // silently round, so the plan is rejected at construction instead.
        let needs_exact_moments = self.observables.iter().any(|observable| {
            matches!(
                observable,
                StaticObservable::RawSecondMoment | StaticObservable::FiniteClusterSusceptibility
            )
        });
        if needs_exact_moments {
            let max_moment = (total_vertices as u128) * (total_vertices as u128);
            if max_moment > (1_u128 << f64::MANTISSA_DIGITS) {
                return Err(ObservablePlanError::InexactMomentConversion { total_vertices });
            }
        }
        Ok(())
    }

    pub(crate) fn record(
        &self,
        total_vertices: usize,
        summary: ComponentSummary,
        boundary_outcomes: &[bool],
        ctx: &mut Context<Xoshiro256PlusPlus>,
    ) {
        for observable in &self.observables {
            match observable {
                StaticObservable::ActiveVertexCount => {
                    ctx.measure(observable.name(), summary.active_vertex_count as f64);
                }
                StaticObservable::ActiveEdgeCount => {
                    ctx.measure(observable.name(), summary.active_edge_count as f64);
                }
                StaticObservable::ComponentCount => {
                    ctx.measure(observable.name(), summary.component_count as f64);
                }
                StaticObservable::LargestSize => {
                    ctx.measure(observable.name(), summary.largest_component_size as f64);
                }
                StaticObservable::GiantFraction if total_vertices != 0 => ctx.measure(
                    observable.name(),
                    summary.largest_component_size as f64 / total_vertices as f64,
                ),
                StaticObservable::GiantFraction => {}
                StaticObservable::RawSecondMoment => {
                    ctx.measure(observable.name(), summary.raw_second_moment as f64);
                }
                StaticObservable::FiniteClusterSusceptibility => {
                    let finite_vertices = summary
                        .active_vertex_count
                        .saturating_sub(summary.largest_component_size);
                    let largest_squared = (summary.largest_component_size as u128)
                        * (summary.largest_component_size as u128);
                    let finite_second_moment =
                        summary.raw_second_moment.saturating_sub(largest_squared);
                    let defined = finite_vertices != 0;
                    ctx.measure(
                        "FiniteClusterSusceptibilityDefined",
                        f64::from(u8::from(defined)),
                    );
                    if defined {
                        ctx.measure(
                            observable.name(),
                            finite_second_moment as f64 / finite_vertices as f64,
                        );
                    }
                }
            }
        }
        for (name, &crosses) in self.boundary_names.iter().zip(boundary_outcomes) {
            ctx.measure(name, f64::from(u8::from(crosses)));
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ObservablePlanError {
    InexactMomentConversion { total_vertices: usize },
}

impl fmt::Display for ObservablePlanError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InexactMomentConversion { total_vertices } => write!(
                formatter,
                "moment observables may exceed exact f64 integer range for V={total_vertices}"
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

    #[test]
    fn records_units_total_vertex_normalization_and_multiple_queries() {
        let graph = build_chain(10, false);
        let plan = ObservablePlan::all(vec![
            BoundaryQuery::new(&graph, &[0], &[9]).unwrap(),
            BoundaryQuery::new(&graph, &[1], &[8]).unwrap(),
        ]);
        let mut ctx = context();
        plan.record(10, summary(6, 3, 14), &[true, false], &mut ctx);
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
        plan.record(0, summary(0, 0, 0), &[], &mut empty);
        let values = empty.finalize_measurements();
        assert!(!values.contains_key("GiantFraction"));
        assert!(!values.contains_key("FiniteClusterSusceptibility"));
        assert_eq!(values["FiniteClusterSusceptibilityDefined"].mean, 0.0);

        let mut giant = context();
        plan.record(4, summary(4, 4, 16), &[], &mut giant);
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
        plan.record(7, summary(7, 3, 19), &[], &mut ctx);
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
        plan.record(4, summary(4, 3, 10), &[], &mut ctx);
        plan.record(7, summary(7, 3, 17), &[], &mut ctx);
        let per_sample_mean = ctx.finalize_measurements()["FiniteClusterSusceptibility"].mean;
        assert_eq!(per_sample_mean, 1.5);
        assert_ne!(per_sample_mean, (1.0 + 8.0) / (1.0 + 4.0));
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
        plan.record(4, summary(4, 3, 11), &[], &mut ctx);
        let values = ctx.finalize_measurements();
        // One sample per sweep, not two: the mean is unchanged by duplicates
        // but n_bins and stderr would silently double without collapsing.
        assert_eq!(values["LargestSize"].n_bins, 1);
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
