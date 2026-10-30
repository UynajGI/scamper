//! W2 winding-analyzer cost: per-sample F3-only vs winding vs F3+winding on
//! builder-embedded tori (D = 1, 2, 3), with lift-BFS parity assertions on
//! non-clipped samples. Criterion IDs live under `percolation_winding/...`.

#[allow(dead_code)]
mod percolation_support;

use std::time::Duration;

use carlo_rs::{Context, MonteCarlo};
use cmc_rs::{
    analyze, analyze_winding, BondBernoulli, BoundaryQuery, ComponentWorkspace, GraphView,
    ObservablePlan, PeriodicEmbedding, Probability, StaticConfiguration, StaticLaw,
    StaticObservable, StaticPercolationMC, UndirectedGraphView, WindingWorkspace,
};
use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use percolation_support::{lift_winding_reference, torus_cases, P_BOND};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

const PARITY_SAMPLES: usize = 16;
const TIMED_INITIAL_SAMPLES: usize = 1;

/// Cross-checks the production winding analyzer against the finite-cover
/// lift reference on freshly sampled configurations (skipping cover-clipped
/// samples) and against F3 for shared connectivity facts.
fn assert_parity(torus: &percolation_support::TorusCase) {
    let mut configuration =
        StaticConfiguration::new(torus.lattice.vertex_count(), torus.lattice.edge_count());
    let mut winding_workspace = WindingWorkspace::new();
    let mut component_workspace = ComponentWorkspace::new();
    let query = BoundaryQuery::new(&torus.lattice, &[0], &[torus.lattice.vertex_count() - 1])
        .expect("ring boundary query");
    component_workspace
        .prepare(torus.lattice.vertex_count(), 1, false)
        .unwrap();
    let bond = BondBernoulli::new(Probability::new(P_BOND).unwrap().into());
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x5749_4e32);
    for _ in 0..PARITY_SAMPLES {
        bond.sample(&torus.lattice, &mut configuration, &mut rng)
            .expect("bench dimensions match");
        let result = analyze_winding(
            &torus.lattice,
            &torus.embedding,
            &configuration,
            &configuration,
            &mut winding_workspace,
        )
        .expect("bench winding analysis");
        let f3 = analyze(
            &torus.lattice,
            &configuration,
            &configuration,
            std::slice::from_ref(&query),
            &mut component_workspace,
        )
        .unwrap();
        assert_eq!(
            result.summary().component_count,
            f3.summary().component_count,
            "{}/component count",
            torus.name
        );
        assert_eq!(
            result.summary().largest_component,
            f3.summary().largest_component,
            "{}/largest identity",
            torus.name
        );
        let (reference, clipped) =
            lift_winding_reference(&torus.lattice, &torus.embedding, &configuration);
        if clipped {
            continue;
        }
        assert_eq!(
            result.components().len(),
            reference.len(),
            "{}/reference component count",
            torus.name
        );
        for (info, (canonical, size, rank)) in result.components().iter().zip(reference.iter()) {
            assert_eq!(
                (info.component().index(), info.size(), info.rank()),
                (*canonical, *size, *rank),
                "{}/reference parity",
                torus.name
            );
        }
    }
}

fn bench_percolation_winding(criterion: &mut Criterion) {
    for torus in torus_cases() {
        assert_parity(&torus);
        let vertex_count = torus.lattice.vertex_count();
        let edge_count = torus.lattice.edge_count();
        let dimensions = torus.embedding.periodic_dimensions();

        let bond = BondBernoulli::new(Probability::new(P_BOND).unwrap().into());
        let mut configuration = StaticConfiguration::new(vertex_count, edge_count);
        let mut winding_workspace = WindingWorkspace::new();
        winding_workspace.prepare(vertex_count, dimensions).unwrap();
        let mut component_workspace = ComponentWorkspace::new();
        let query = BoundaryQuery::new(&torus.lattice, &[0], &[vertex_count - 1])
            .expect("ring boundary query");
        let queries = vec![query];
        component_workspace
            .prepare(vertex_count, queries.len(), false)
            .unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x4633_5753_4e44);
        for _ in 0..TIMED_INITIAL_SAMPLES {
            bond.sample(&torus.lattice, &mut configuration, &mut rng)
                .expect("bench dimensions match");
        }

        let mut group = criterion.benchmark_group(format!(
            "percolation_winding/{}/V={}/E={}/D={}",
            torus.name, vertex_count, edge_count, dimensions
        ));
        group.throughput(Throughput::Elements(edge_count as u64));

        group.bench_function("sample+f3-analyze", |bencher| {
            bencher.iter(|| {
                bond.sample(&torus.lattice, &mut configuration, &mut rng)
                    .expect("bench dimensions match");
                black_box(
                    analyze(
                        &torus.lattice,
                        &configuration,
                        &configuration,
                        &queries,
                        &mut component_workspace,
                    )
                    .unwrap(),
                );
            });
        });
        group.bench_function("sample+winding-analyze", |bencher| {
            bencher.iter(|| {
                bond.sample(&torus.lattice, &mut configuration, &mut rng)
                    .expect("bench dimensions match");
                black_box(
                    analyze_winding(
                        &torus.lattice,
                        &torus.embedding,
                        &configuration,
                        &configuration,
                        &mut winding_workspace,
                    )
                    .unwrap(),
                );
            });
        });
        group.bench_function("sample+f3+winding", |bencher| {
            bencher.iter(|| {
                bond.sample(&torus.lattice, &mut configuration, &mut rng)
                    .expect("bench dimensions match");
                black_box(
                    analyze(
                        &torus.lattice,
                        &configuration,
                        &configuration,
                        &queries,
                        &mut component_workspace,
                    )
                    .unwrap(),
                );
                black_box(
                    analyze_winding(
                        &torus.lattice,
                        &torus.embedding,
                        &configuration,
                        &configuration,
                        &mut winding_workspace,
                    )
                    .unwrap(),
                );
            });
        });

        // Full adapter path with winding observables requested: sweep +
        // measure runs the F3 analysis, the winding analysis, and records
        // the scalars through a Carlo.rs Context.
        let plan = ObservablePlan::new(
            [
                StaticObservable::ComponentCount,
                StaticObservable::LargestSize,
                StaticObservable::GiantFraction,
                StaticObservable::Wrapping,
                StaticObservable::WindingRank,
            ],
            vec![
                BoundaryQuery::new(&torus.lattice, &[0], &[vertex_count - 1])
                    .expect("ring boundary query"),
            ],
        );
        let mut adapter = StaticPercolationMC::with_embedding(
            torus.lattice.clone(),
            torus.embedding.clone(),
            StaticLaw::bond(P_BOND).unwrap(),
            plan,
        )
        .unwrap();
        let mut mc_context =
            Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(0x4633_5753_4e44), 0, 1024);
        for _ in 0..TIMED_INITIAL_SAMPLES {
            adapter.sweep(&mut mc_context);
            adapter.measure(&mut mc_context);
        }
        group.bench_function("adapter-sweep+measure-winding", |bencher| {
            bencher.iter(|| {
                adapter.sweep(&mut mc_context);
                adapter.measure(&mut mc_context);
                black_box(&mut mc_context);
            });
        });

        let plain_plan = ObservablePlan::new(
            [
                StaticObservable::ComponentCount,
                StaticObservable::LargestSize,
                StaticObservable::GiantFraction,
            ],
            vec![
                BoundaryQuery::new(&torus.lattice, &[0], &[vertex_count - 1])
                    .expect("ring boundary query"),
            ],
        );
        let mut plain_adapter = StaticPercolationMC::new(
            torus.lattice.clone(),
            StaticLaw::bond(P_BOND).unwrap(),
            plain_plan,
        )
        .unwrap();
        let mut plain_context =
            Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(0x4633_5753_4e44), 0, 1024);
        for _ in 0..TIMED_INITIAL_SAMPLES {
            plain_adapter.sweep(&mut plain_context);
            plain_adapter.measure(&mut plain_context);
        }
        group.bench_function("adapter-sweep+measure-plain", |bencher| {
            bencher.iter(|| {
                plain_adapter.sweep(&mut plain_context);
                plain_adapter.measure(&mut plain_context);
                black_box(&mut plain_context);
            });
        });
        group.finish();
    }
}

criterion_group! {
    name = benches;
    config = Criterion::default()
        .sample_size(10)
        .warm_up_time(Duration::from_millis(100))
        .measurement_time(Duration::from_millis(250));
    targets = bench_percolation_winding
}
criterion_main!(benches);
