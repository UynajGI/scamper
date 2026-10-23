//! R2 q = 1 random-cluster realization benchmark.
//!
//! Realization throughput of the `RandomClusterMC` sweep+measure path
//! (BondBernoulli delegate -> dense-mask mirror -> R1 external-bond
//! analysis) versus the static bond path (BondBernoulli draw ->
//! `StaticConfiguration` analysis) on the five shared topologies at
//! p = 0.5. Both pipelines end in the same F3 analyzer over the same
//! activity; the q = 1 adapter's extra work is exactly the mask mirror
//! that feeds the R1 external path, which is what this bench isolates.
//! Observable-plan recording (`Context::measure`) is excluded from both
//! paths equally. `BondSample` materialization (the R3 contract surface)
//! is benched separately.

#[allow(dead_code)]
mod percolation_support;

use std::time::Duration;

use cmc_rs::{
    analyze, analyze_bonds, BondBernoulli, BondSample, BoundaryQuery, ComponentWorkspace,
    EdgeActivity, GraphView, Probability, StaticConfiguration, UndirectedGraphView,
};
use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use percolation_support::{cases, P_BOND};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

const DRAW_SEED: u64 = 20261022;

/// One q = 1 realization through the adapter's exact sweep+measure work
/// (minus plan recording): draw via the BondBernoulli delegate, mirror the
/// draw into the reused dense mask, analyze through the R1 external path.
fn q1_realization(
    lattice: &impl UndirectedGraphView,
    law: &BondBernoulli,
    configuration: &mut StaticConfiguration,
    mask: &mut [bool],
    workspace: &mut ComponentWorkspace,
    queries: &[BoundaryQuery],
    rng: &mut Xoshiro256PlusPlus,
) {
    law.sample(lattice, configuration, rng)
        .expect("bench dimensions match");
    for (index, open) in mask.iter_mut().enumerate() {
        let edge = lattice.edge_id(index).expect("dense edge index");
        *open = configuration.edge_active(edge);
    }
    black_box(analyze_bonds(lattice, mask, queries, workspace).expect("bench dimensions match"));
}

/// One static-bond realization through `StaticPercolationMC`'s exact
/// sweep+measure work (minus plan recording).
fn static_realization(
    lattice: &impl UndirectedGraphView,
    law: &BondBernoulli,
    configuration: &mut StaticConfiguration,
    workspace: &mut ComponentWorkspace,
    queries: &[BoundaryQuery],
    rng: &mut Xoshiro256PlusPlus,
) {
    law.sample(lattice, configuration, rng)
        .expect("bench dimensions match");
    black_box(
        analyze(lattice, configuration, configuration, queries, workspace)
            .expect("bench dimensions match"),
    );
}

fn bench_random_cluster_q1(criterion: &mut Criterion) {
    for case in cases() {
        let queries = [BoundaryQuery::new(&case.lattice, &case.from, &case.to).unwrap()];
        let edge_count = case.lattice.edge_count();
        let law = BondBernoulli::new(Probability::new(P_BOND).unwrap().into());

        let mut q1_configuration =
            StaticConfiguration::new(case.lattice.vertex_count(), edge_count);
        let mut q1_mask = vec![false; edge_count];
        let mut q1_workspace = ComponentWorkspace::new();
        q1_workspace
            .prepare(case.lattice.vertex_count(), queries.len(), false)
            .unwrap();
        let mut static_configuration =
            StaticConfiguration::new(case.lattice.vertex_count(), edge_count);
        let mut static_workspace = ComponentWorkspace::new();
        static_workspace
            .prepare(case.lattice.vertex_count(), queries.len(), false)
            .unwrap();

        // Parity outside timing: identical draws (same seed, one f64 per
        // edge) analyze to identical summaries through both paths.
        q1_realization(
            &case.lattice,
            &law,
            &mut q1_configuration,
            &mut q1_mask,
            &mut q1_workspace,
            &queries,
            &mut Xoshiro256PlusPlus::seed_from_u64(DRAW_SEED),
        );
        static_realization(
            &case.lattice,
            &law,
            &mut static_configuration,
            &mut static_workspace,
            &queries,
            &mut Xoshiro256PlusPlus::seed_from_u64(DRAW_SEED),
        );
        let q1_summary = analyze_bonds(&case.lattice, &q1_mask, &queries, &mut q1_workspace)
            .unwrap()
            .summary();
        let static_summary = analyze(
            &case.lattice,
            &static_configuration,
            &static_configuration,
            &queries,
            &mut static_workspace,
        )
        .unwrap()
        .summary();
        assert_eq!(
            q1_summary, static_summary,
            "{}: the two realization paths diverge",
            case.name
        );

        let mut group = criterion.benchmark_group(format!(
            "random_cluster_q1/{}/V={}/E={}/p_bond={}",
            case.name,
            case.lattice.vertex_count(),
            edge_count,
            P_BOND
        ));
        group.throughput(Throughput::Elements(case.lattice.vertex_count() as u64));
        group.bench_function("q1 realization (external path)", |bencher| {
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(DRAW_SEED);
            bencher.iter(|| {
                q1_realization(
                    &case.lattice,
                    &law,
                    &mut q1_configuration,
                    &mut q1_mask,
                    &mut q1_workspace,
                    &queries,
                    &mut rng,
                );
            });
        });
        group.bench_function("static bond realization", |bencher| {
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(DRAW_SEED);
            bencher.iter(|| {
                static_realization(
                    &case.lattice,
                    &law,
                    &mut static_configuration,
                    &mut static_workspace,
                    &queries,
                    &mut rng,
                );
            });
        });
        group.bench_function("q1 BondSample materialization", |bencher| {
            bencher.iter(|| {
                black_box(BondSample::from_mask(&case.lattice, &q1_mask).unwrap());
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
    targets = bench_random_cluster_q1
}
criterion_main!(benches);
