//! E2 benchmark: bootstrap-percolation closure realizations.
//!
//! One iteration = 8 complete realizations (density seed sampling plus one
//! O(V + E) closure with the E1 queue/counter/stamp buffers) at r = 2 on
//! the shared benchmark topologies, for both recorded clocks, plus the
//! Holroyd-anchor operating point on the periodic 128 x 128 square at its
//! measured half-majority density.
//!
//! Run with: cargo bench -p cmc-rs --bench bootstrap_percolation

// Shared bench support items are dead in this target.
#![allow(dead_code)]

mod percolation_support;

use carlo_rs::{Context, MonteCarlo};
use cmc_rs::{
    build_square, BootstrapDynamics, BootstrapPercolationMC, BootstrapSeeds, BootstrapThreshold,
    GraphView, UndirectedGraphView,
};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use percolation_support::cases;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

const REALIZATIONS: usize = 8;
const SEED_DENSITY: f64 = 0.1;
const STEPS: usize = 16;

fn bench_bootstrap_closure(c: &mut Criterion) {
    let mut group = c.benchmark_group("bootstrap_percolation");
    for case in cases() {
        let vertices = case.lattice.vertex_count();
        let edges = case.lattice.edge_count();
        group.throughput(criterion::Throughput::Elements(
            (REALIZATIONS * (vertices + edges)) as u64,
        ));
        for dynamics in [
            BootstrapDynamics::SynchronousRounds,
            BootstrapDynamics::QueueClosure,
        ] {
            let label = match dynamics {
                BootstrapDynamics::SynchronousRounds => "sync_rounds",
                BootstrapDynamics::QueueClosure => "queue_closure",
            };
            let mut process = BootstrapPercolationMC::new(
                case.lattice.clone(),
                BootstrapThreshold::Homogeneous { value: 2 },
                BootstrapSeeds::Density { p: SEED_DENSITY },
                dynamics,
                STEPS,
            )
            .expect("bench adapter");
            let mut ctx =
                Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(0x00e2), 0, 1 << 12);
            group.bench_with_input(
                BenchmarkId::new(label, format!("{}_V{vertices}", case.name)),
                &case.name,
                |b, _| {
                    b.iter(|| {
                        for _ in 0..REALIZATIONS {
                            process.sweep(&mut ctx);
                            std::hint::black_box(process.last_trajectory());
                        }
                    });
                },
            );
        }
    }
    // The Holroyd-anchor operating point: 2D r = 2 on the periodic
    // 128 x 128 square near the measured half-majority density.
    let width = 128usize;
    let anchor_density = 0.048;
    let vertices = width * width;
    let edges = 2 * vertices;
    group.throughput(criterion::Throughput::Elements(
        (REALIZATIONS * (vertices + edges)) as u64,
    ));
    for dynamics in [
        BootstrapDynamics::SynchronousRounds,
        BootstrapDynamics::QueueClosure,
    ] {
        let label = match dynamics {
            BootstrapDynamics::SynchronousRounds => "sync_rounds",
            BootstrapDynamics::QueueClosure => "queue_closure",
        };
        let mut process = BootstrapPercolationMC::new(
            build_square(width, width, true),
            BootstrapThreshold::Homogeneous { value: 2 },
            BootstrapSeeds::Density { p: anchor_density },
            dynamics,
            STEPS,
        )
        .expect("bench adapter");
        let mut ctx =
            Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(0x00e3), 0, 1 << 12);
        group.bench_with_input(
            BenchmarkId::new(label, "torus128_holroyd"),
            &label,
            |b, _| {
                b.iter(|| {
                    for _ in 0..REALIZATIONS {
                        process.sweep(&mut ctx);
                        std::hint::black_box(process.last_trajectory());
                    }
                });
            },
        );
    }
    group.finish();
}

criterion_group!(benches, bench_bootstrap_closure);
criterion_main!(benches);
