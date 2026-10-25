//! A1a: open custom-observable extension.
//!
//! This integration test lives outside the library crate, so implementing a
//! [`CustomObservable`] here proves downstream crates can add a Results column
//! without editing CMC.rs. The built-in v2 plan path is unchanged; a plan with
//! no custom observables must behave exactly as before.

use std::sync::Arc;

use carlo_rs::{Context, Run, RunConfig, RunId, TaskId};
use cmc_rs::{
    build_chain, AnalysisResult, CustomObservable, ObservablePlan, ObservablePlanError, StaticLaw,
    StaticObservable, StaticPercolationMC,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

#[derive(Debug)]
struct DoubleActiveEdges;

impl CustomObservable for DoubleActiveEdges {
    fn name(&self) -> &str {
        "CustomDoubleActiveEdges"
    }

    fn evaluate(&self, analysis: &AnalysisResult<'_>) -> Option<f64> {
        Some(2.0 * analysis.summary().active_edge_count as f64)
    }
}

#[derive(Debug)]
struct NeverDefined;

impl CustomObservable for NeverDefined {
    fn name(&self) -> &str {
        "CustomNeverDefined"
    }

    fn evaluate(&self, _analysis: &AnalysisResult<'_>) -> Option<f64> {
        None
    }
}

#[derive(Debug)]
struct FixedName(&'static str);

impl CustomObservable for FixedName {
    fn name(&self) -> &str {
        self.0
    }

    fn evaluate(&self, _analysis: &AnalysisResult<'_>) -> Option<f64> {
        Some(1.0)
    }
}

fn custom(name: &'static str) -> Arc<dyn CustomObservable> {
    Arc::new(FixedName(name))
}

fn run(plan: ObservablePlan, sweeps: u64, binsize: usize) -> carlo_rs::Results {
    let lattice = build_chain(4, false);
    let adapter = StaticPercolationMC::new(lattice, StaticLaw::site(1.0).expect("law"), plan)
        .expect("adapter");
    let seed = 0x4131_0001;
    let context = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(seed), 0, binsize);
    let mut run = Run::from_parts(
        context,
        adapter,
        TaskId::new(0),
        RunId::new(0),
        RunConfig {
            thermalization_sweeps: 0,
            measurement_sweeps: sweeps,
            binsize,
            base_seed: seed,
            ..Default::default()
        },
    );
    run.run(sweeps);
    run.finalize(seed)
}

#[test]
fn external_custom_observable_is_recorded_through_the_scheduler() {
    let plan = ObservablePlan::new([StaticObservable::ActiveVertexCount], Vec::new())
        .with_custom([
            Arc::new(DoubleActiveEdges) as Arc<dyn CustomObservable>,
            Arc::new(NeverDefined) as Arc<dyn CustomObservable>,
        ])
        .expect("custom plan");
    assert_eq!(plan.custom_observables().len(), 2);
    assert!(plan.has_custom());

    let results = run(plan, 4, 1);
    let built_in = results.get("ActiveVertexCount").expect("built-in column");
    assert_eq!(built_in.mean, 4.0);
    assert_eq!(built_in.n_bins, 4);

    // Chain of 4 at p=1 has 3 physical edges; the custom column doubles that.
    let custom = results
        .get("CustomDoubleActiveEdges")
        .expect("custom column");
    assert_eq!(custom.mean, 6.0);
    assert_eq!(custom.n_bins, 4);
    // A custom observable may leave a sample undefined; no column is written.
    assert!(results.get("CustomNeverDefined").is_none());
}

#[test]
fn custom_names_are_validated_against_builtins_boundaries_and_each_other() {
    let base = || ObservablePlan::new([StaticObservable::LargestSize], Vec::new());
    assert!(matches!(
        base().with_custom([custom("LargestSize")]),
        Err(ObservablePlanError::CustomNameCollision { .. })
    ));
    assert!(matches!(
        base().with_custom([
            Arc::new(DoubleActiveEdges) as Arc<dyn CustomObservable>,
            Arc::new(DoubleActiveEdges) as Arc<dyn CustomObservable>,
        ]),
        Err(ObservablePlanError::DuplicateCustomName { .. })
    ));
    assert!(matches!(
        base().with_custom([custom("")]),
        Err(ObservablePlanError::EmptyCustomName)
    ));

    // The susceptibility indicator and boundary columns are reserved too.
    let susceptibility =
        ObservablePlan::new([StaticObservable::FiniteClusterSusceptibility], Vec::new());
    assert!(matches!(
        susceptibility.with_custom([custom("FiniteClusterSusceptibilityDefined")]),
        Err(ObservablePlanError::CustomNameCollision { .. })
    ));

    let lattice = build_chain(4, false);
    let query = cmc_rs::BoundaryQuery::new(&lattice, &[0], &[3]).expect("query");
    let boundary = ObservablePlan::new([StaticObservable::ActiveVertexCount], vec![query]);
    assert!(matches!(
        boundary.with_custom([custom("BoundaryCrossing0")]),
        Err(ObservablePlanError::CustomNameCollision { .. })
    ));
}

#[test]
fn a_plan_without_custom_observables_is_byte_identical_to_the_v2_path() {
    let baseline = run(
        ObservablePlan::new([StaticObservable::ActiveVertexCount], Vec::new()),
        4,
        1,
    );
    let empty_custom = run(
        ObservablePlan::new([StaticObservable::ActiveVertexCount], Vec::new())
            .with_custom(Vec::new())
            .expect("empty custom list"),
        4,
        1,
    );
    let left = baseline.get("ActiveVertexCount").expect("baseline");
    let right = empty_custom.get("ActiveVertexCount").expect("custom");
    assert_eq!(left.mean.to_bits(), right.mean.to_bits());
    assert_eq!(left.stderr.to_bits(), right.stderr.to_bits());
    assert_eq!(left.n_bins, right.n_bins);
}
