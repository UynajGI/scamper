//! G2 benchmark: Achlioptas process steps per second.
//!
//! One iteration = 2 complete realizations on the shared benchmark
//! topologies at V = 4096 (product rule, m = 2, without replacement, full
//! edge universe, `Every` grid with 64 windows). Each realization
//! re-initializes the V vertices, draws 2 candidates per step, scores both
//! through the G1 `component_size` primitive, and accepts one edge, running
//! to complete occupation (E steps). The throughput unit is process steps
//! (candidate draws + scoring + one union each).
//!
//! Run with: cargo bench -p cmc-rs --bench achlioptas_process

// Shared bench-support items are dead in this target.
#![allow(dead_code)]

mod percolation_support;

use cmc_rs::{
    AchlioptasGrid, AchlioptasProcessMC, AchlioptasRule, CandidateSampling, EdgeUniverse,
    GraphView, TiePolicy, UndirectedGraphView,
};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use percolation_support::cases;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

const REALIZATIONS: usize = 2;

fn bench_achlioptas_steps(c: &mut Criterion) {
    let mut group = c.benchmark_group("achlioptas_process");
    for case in cases() {
        let edges = case.lattice.edge_count();
        let mut process = AchlioptasProcessMC::new(
            case.lattice.clone(),
            AchlioptasRule::Product,
            CandidateSampling::WithoutReplacement,
            2,
            TiePolicy::LowestEdgeIndex,
            EdgeUniverse::AllEdges,
            AchlioptasGrid::Every {
                step: (edges / 64).max(1),
            },
        )
        .expect("bench process");
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x00C2_4247);
        let vertices = case.lattice.vertex_count();

        // Every realization runs E steps to complete occupation.
        group.throughput(criterion::Throughput::Elements(
            (REALIZATIONS * edges) as u64,
        ));
        group.bench_with_input(
            BenchmarkId::new("steps", format!("{}_V{vertices}", case.name)),
            &case.name,
            |b, _| {
                b.iter(|| {
                    for _ in 0..REALIZATIONS {
                        process.run_realization(&mut rng);
                        std::hint::black_box(process.largest_curve());
                    }
                });
            },
        );
    }
    group.finish();
}

criterion_group!(benches, bench_achlioptas_steps);
criterion_main!(benches);
