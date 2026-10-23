//! G1 benchmark: incremental unions per second.
//!
//! One iteration = 4 complete realizations on the shared benchmark
//! topologies at V = 4096. Each realization resets the reused workspace,
//! activates every physical edge (waiting), then activates all vertices in a
//! fresh random order — the Z1 site law expressed through the generic
//! substrate, V + E accepted operations per realization. The throughput unit
//! is accepted operations (unions).
//!
//! Run with: cargo bench -p cmc-rs --bench incremental_process

// Shared bench-support items are dead in this target.
#![allow(dead_code)]

mod percolation_support;

use cmc_rs::{GraphView, IncrementalUnionWorkspace, UndirectedGraphView};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use percolation_support::cases;
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

const REALIZATIONS: usize = 4;

fn bench_incremental_unions(c: &mut Criterion) {
    let mut group = c.benchmark_group("incremental_process");
    for case in cases() {
        let vertices = case.lattice.vertex_count();
        let edges = case.lattice.edge_count();
        let mut workspace = IncrementalUnionWorkspace::new();
        workspace.prepare(vertices, edges).expect("bench prepare");
        let mut order: Vec<usize> = (0..vertices).collect();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x00C1_6731);

        // Every realization accepts each vertex and each edge exactly once.
        group.throughput(criterion::Throughput::Elements(
            (REALIZATIONS * (vertices + edges)) as u64,
        ));
        group.bench_with_input(
            BenchmarkId::new("unions", format!("{}_V{vertices}", case.name)),
            &case.name,
            |b, _| {
                b.iter(|| {
                    for _ in 0..REALIZATIONS {
                        workspace.reset();
                        for index in 0..edges {
                            let edge = case.lattice.edge_id(index).expect("dense edge");
                            std::hint::black_box(
                                workspace
                                    .activate_edge(&case.lattice, edge, index)
                                    .expect("bench edge activation"),
                            );
                        }
                        for index in (1..order.len()).rev() {
                            let swap_with = rng.random_range(0..=index);
                            order.swap(index, swap_with);
                        }
                        for (step, &vertex) in order.iter().enumerate() {
                            let vertex_id = case.lattice.vertex_id(vertex).expect("dense vertex");
                            std::hint::black_box(
                                workspace
                                    .activate_vertex(&case.lattice, vertex_id, step)
                                    .expect("bench vertex activation"),
                            );
                        }
                    }
                    std::hint::black_box(workspace.trajectory().summaries());
                });
            },
        );
    }
    group.finish();
}

criterion_group!(benches, bench_incremental_unions);
criterion_main!(benches);
