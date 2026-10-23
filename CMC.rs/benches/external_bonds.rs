//! R1 external bond-configuration benchmark.
//!
//! Mask-to-analysis round trip on the five shared topologies at p = 0.5,
//! with a reused workspace, compared against the `StaticConfiguration`
//! path holding the SAME bonds (drawn through `BondBernoulli` with an
//! identically seeded RNG — the law's one-f64-per-edge draw contract makes
//! the two masks bit-identical, asserted once below). The two paths execute
//! the same F3 analyzer over the same activity, so the recorded times must
//! be equal to noise; any systematic gap would mean the external-mask
//! wrapper adds work, which is exactly what R1 must not do.

#[allow(dead_code)]
mod percolation_support;

use std::time::Duration;

use cmc_rs::{
    analyze, analyze_bonds, BondBernoulli, BoundaryQuery, ComponentWorkspace, EdgeActivity,
    GraphView, Probability, StaticConfiguration, UndirectedGraphView,
};
use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use percolation_support::{cases, P_BOND};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

const MASK_SEED: u64 = 0x4d41_534b;

fn bench_external_bonds(criterion: &mut Criterion) {
    for case in cases() {
        let queries = [BoundaryQuery::new(&case.lattice, &case.from, &case.to).unwrap()];
        let edge_count = case.lattice.edge_count();

        // One shared bond draw at p = P_BOND: raw mask for the external path,
        // StaticConfiguration for the law path. Identically seeded and the
        // law draws exactly one f64 per edge in order, so the masks match.
        let mut mask_rng = Xoshiro256PlusPlus::seed_from_u64(MASK_SEED);
        let mask: Vec<bool> = (0..edge_count)
            .map(|_| mask_rng.random::<f64>() < P_BOND)
            .collect();
        let mut configuration = StaticConfiguration::new(case.lattice.vertex_count(), edge_count);
        BondBernoulli::new(Probability::new(P_BOND).unwrap().into())
            .sample(
                &case.lattice,
                &mut configuration,
                &mut Xoshiro256PlusPlus::seed_from_u64(MASK_SEED),
            )
            .expect("benchmark dimensions match");
        for (index, &open) in mask.iter().enumerate() {
            let edge = case.lattice.edge_id(index).expect("dense edge ID");
            assert_eq!(
                configuration.edge_active(edge),
                open,
                "{}: the two paths must analyze the same bonds",
                case.name
            );
        }

        let mut external_workspace = ComponentWorkspace::new();
        external_workspace
            .prepare(case.lattice.vertex_count(), queries.len(), false)
            .unwrap();
        let mut static_workspace = ComponentWorkspace::new();
        static_workspace
            .prepare(case.lattice.vertex_count(), queries.len(), false)
            .unwrap();

        // Parity outside timing: both paths produce the identical summary.
        let external_summary =
            analyze_bonds(&case.lattice, &mask, &queries, &mut external_workspace)
                .unwrap()
                .summary();
        let static_summary = analyze(
            &case.lattice,
            &configuration,
            &configuration,
            &queries,
            &mut static_workspace,
        )
        .unwrap()
        .summary();
        assert_eq!(
            external_summary, static_summary,
            "{}: external and static paths diverge",
            case.name
        );

        let mut group = criterion.benchmark_group(format!(
            "external_bonds/{}/V={}/E={}/p_bond={}",
            case.name,
            case.lattice.vertex_count(),
            edge_count,
            P_BOND
        ));
        group.throughput(Throughput::Elements(edge_count as u64));
        group.bench_function("analyze_bonds (external mask)", |bencher| {
            bencher.iter(|| {
                black_box(
                    analyze_bonds(&case.lattice, &mask, &queries, &mut external_workspace).unwrap(),
                );
            });
        });
        group.bench_function("analyze (StaticConfiguration)", |bencher| {
            bencher.iter(|| {
                black_box(
                    analyze(
                        &case.lattice,
                        &configuration,
                        &configuration,
                        &queries,
                        &mut static_workspace,
                    )
                    .unwrap(),
                );
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
    targets = bench_external_bonds
}
criterion_main!(benches);
