//! A1: open observable extension and capability declarations.
//!
//! This integration test lives outside the library crate, so implementing an
//! [`ScalarObservable`] here proves downstream crates can add a Results column
//! without editing CMC.rs. Plans without customs must keep the v2 recording
//! path bit-for-bit, and capability declarations (labels, winding, exact
//! moments) must control what adapters materialize.

use std::sync::Arc;

use carlo_rs::{Context, MonteCarlo, Results, Run, RunConfig, RunId, TaskId};
use cmc_rs::{
    build_chain, build_square_with_embedding, ClusterQ, ObservablePlan, ObservablePlanError,
    RandomClusterMC, SampleView, ScalarObservable, StaticLaw, StaticObservable,
    StaticPercolationError, StaticPercolationMC,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

#[derive(Debug)]
struct DoubleActiveEdges;

impl ScalarObservable for DoubleActiveEdges {
    fn name(&self) -> &str {
        "CustomDoubleActiveEdges"
    }

    fn evaluate(&self, sample: &SampleView<'_>) -> Option<f64> {
        Some(2.0 * sample.summary().active_edge_count as f64)
    }
}

#[derive(Debug)]
struct NeverDefined;

impl ScalarObservable for NeverDefined {
    fn name(&self) -> &str {
        "CustomNeverDefined"
    }

    fn evaluate(&self, _sample: &SampleView<'_>) -> Option<f64> {
        None
    }
}

#[derive(Debug)]
struct FixedName(&'static str);

impl ScalarObservable for FixedName {
    fn name(&self) -> &str {
        self.0
    }

    fn evaluate(&self, _sample: &SampleView<'_>) -> Option<f64> {
        Some(1.0)
    }
}

/// Reports the label-array length, or stays undefined when labels were not
/// materialized. `declared` toggles the capability so one observable can prove
/// both the requested and the unrequested path.
#[derive(Debug)]
struct LabelCount {
    name: &'static str,
    declared: bool,
}

impl ScalarObservable for LabelCount {
    fn name(&self) -> &str {
        self.name
    }

    fn requires_labels(&self) -> bool {
        self.declared
    }

    fn evaluate(&self, sample: &SampleView<'_>) -> Option<f64> {
        sample.labels().map(|labels| labels.len() as f64)
    }
}

#[derive(Debug)]
struct LargestWindingRank;

impl ScalarObservable for LargestWindingRank {
    fn name(&self) -> &str {
        "CustomLargestWindingRank"
    }

    fn requires_winding(&self) -> bool {
        true
    }

    fn evaluate(&self, sample: &SampleView<'_>) -> Option<f64> {
        sample
            .winding()
            .map(|winding| winding.largest_winding_rank as f64)
    }
}

fn custom(name: &'static str) -> Arc<dyn ScalarObservable> {
    Arc::new(FixedName(name))
}

fn run_adapter<A>(adapter: A, seed: u64, sweeps: u64, binsize: usize) -> Results
where
    A: MonteCarlo<Rng = Xoshiro256PlusPlus>,
{
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

fn run_static(plan: ObservablePlan, sweeps: u64, binsize: usize) -> Results {
    let lattice = build_chain(4, false);
    let adapter = StaticPercolationMC::new(lattice, StaticLaw::site(1.0).expect("law"), plan)
        .expect("adapter");
    run_adapter(adapter, 0x4131_0001, sweeps, binsize)
}

#[test]
fn external_custom_observable_is_recorded_through_the_scheduler() {
    let plan = ObservablePlan::new([StaticObservable::ActiveVertexCount], Vec::new())
        .with_custom([
            Arc::new(DoubleActiveEdges) as Arc<dyn ScalarObservable>,
            Arc::new(NeverDefined) as Arc<dyn ScalarObservable>,
        ])
        .expect("custom plan");
    assert_eq!(plan.custom_observables().len(), 2);
    assert!(plan.has_custom());

    let results = run_static(plan, 4, 1);
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
fn labels_are_materialized_only_when_a_custom_observable_requires_them() {
    let requested = ObservablePlan::new(Vec::<StaticObservable>::new(), Vec::new())
        .with_custom([Arc::new(LabelCount {
            name: "RequestedLabelCount",
            declared: true,
        }) as Arc<dyn ScalarObservable>])
        .expect("plan");
    assert!(requested.requests_labels());
    let results = run_static(requested, 2, 1);
    assert_eq!(
        results.get("RequestedLabelCount").expect("labels").mean,
        4.0
    );

    let undeclared = ObservablePlan::new(Vec::<StaticObservable>::new(), Vec::new())
        .with_custom([Arc::new(LabelCount {
            name: "UndeclaredLabelCount",
            declared: false,
        }) as Arc<dyn ScalarObservable>])
        .expect("plan");
    assert!(!undeclared.requests_labels());
    let results = run_static(undeclared, 2, 1);
    // The analyzer skipped labels, so the observable correctly stays undefined
    // instead of inventing an empty array.
    assert!(results.get("UndeclaredLabelCount").is_none());
}

#[test]
fn bond_adapter_honors_the_label_capability() {
    let plan = ObservablePlan::new(Vec::<StaticObservable>::new(), Vec::new())
        .with_custom([Arc::new(LabelCount {
            name: "BondLabelCount",
            declared: true,
        }) as Arc<dyn ScalarObservable>])
        .expect("plan");
    let lattice = build_chain(4, false);
    let adapter = RandomClusterMC::new(lattice, ClusterQ::q1(), 1.0, plan).expect("adapter");
    let results = run_adapter(adapter, 0x4131_0002, 2, 1);
    assert_eq!(results.get("BondLabelCount").expect("labels").mean, 4.0);
}

#[test]
fn custom_winding_capability_requires_an_embedding_and_sees_the_summary() {
    let plan = ObservablePlan::new(Vec::<StaticObservable>::new(), Vec::new())
        .with_custom([Arc::new(LargestWindingRank) as Arc<dyn ScalarObservable>])
        .expect("plan");
    assert!(plan.requests_winding());

    let lattice = build_chain(4, false);
    assert!(matches!(
        StaticPercolationMC::new(lattice, StaticLaw::site(1.0).expect("law"), plan.clone()),
        Err(StaticPercolationError::WindingWithoutEmbedding)
    ));

    let (lattice, embedding) = build_square_with_embedding(4, 4, true).expect("torus");
    let adapter = StaticPercolationMC::with_embedding(
        lattice,
        embedding,
        StaticLaw::site(1.0).expect("law"),
        plan,
    )
    .expect("embedded adapter");
    let results = run_adapter(adapter, 0x4131_0003, 2, 1);
    // Fully occupied 4x4 torus: the largest component wraps both axes.
    assert_eq!(
        results
            .get("CustomLargestWindingRank")
            .expect("winding column")
            .mean,
        2.0
    );
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
            Arc::new(DoubleActiveEdges) as Arc<dyn ScalarObservable>,
            Arc::new(DoubleActiveEdges) as Arc<dyn ScalarObservable>,
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

/// A1 allocation parity: the custom slot must not allocate per sample. The
/// context binsize exceeds warmup + measured samples, so no completed-bin
/// allocation can fall inside the measured window; each column's first-sample
/// insertion is also warmed up.
#[cfg(feature = "allocation-probe")]
#[test]
fn custom_observable_recording_is_allocation_free_in_steady_state() {
    use allocation_counter::measure;

    const WARMUP: usize = 64;
    const SAMPLES: usize = 256;

    fn steady_state_allocations(plan: ObservablePlan) -> u64 {
        let lattice = build_chain(64, false);
        let mut adapter =
            StaticPercolationMC::new(lattice, StaticLaw::site(0.5).expect("law"), plan)
                .expect("adapter");
        let mut context = Context::new_with_binsize(
            Xoshiro256PlusPlus::seed_from_u64(0x4131_0004),
            0,
            WARMUP + SAMPLES + 1,
        );
        for _ in 0..WARMUP {
            adapter.sweep(&mut context);
            adapter.measure(&mut context);
        }
        let info = measure(|| {
            for _ in 0..SAMPLES {
                adapter.sweep(&mut context);
                adapter.measure(&mut context);
            }
        });
        info.count_total
    }

    let builtin = ObservablePlan::new([StaticObservable::ActiveVertexCount], Vec::new());
    let with_custom = builtin
        .clone()
        .with_custom([Arc::new(DoubleActiveEdges) as Arc<dyn ScalarObservable>])
        .expect("custom plan");
    assert_eq!(steady_state_allocations(builtin), 0);
    assert_eq!(steady_state_allocations(with_custom), 0);
}

#[test]
fn a_plan_without_custom_observables_is_byte_identical_to_the_v2_path() {
    let baseline = run_static(
        ObservablePlan::new([StaticObservable::ActiveVertexCount], Vec::new()),
        4,
        1,
    );
    let empty_custom = run_static(
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
