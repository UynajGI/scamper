//! E3 benchmark: k-core peeling realizations.
//!
//! One iteration = 8 complete realizations (site-keep activity draw plus
//! deterministic peeling and summary recording) on the shared benchmark
//! topologies, at k = 1 (light peeling) and k = 3 (full cascade). The
//! recorded-round budget is 128; peeling itself always runs to completion.
//!
//! Run with: cargo bench -p cmc-rs --bench kcore_percolation

// Shared bench support items are dead in this target.
#![allow(dead_code)]

mod percolation_support;

use carlo_rs::{Context, MonteCarlo};
use cmc_rs::{CoreThresholds, GraphView, KCoreMC};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use percolation_support::cases;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

const REALIZATIONS: usize = 8;
const KEEP: f64 = 0.5927;
const RECORDED_ROUNDS: usize = 128;

fn bench_kcore(c: &mut Criterion) {
    let mut group = c.benchmark_group("kcore_percolation");
    for case in cases() {
        let vertices = case.lattice.vertex_count();
        for k in [1usize, 3] {
            let mut mc = KCoreMC::new(
                case.lattice.clone(),
                KEEP,
                CoreThresholds::uniform(k).expect("bench threshold"),
                RECORDED_ROUNDS,
            )
            .expect("bench adapter");
            let mut ctx =
                Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(0xE3B3), 0, 64);
            // One realization peels every vertex at most once: the
            // throughput unit is vertices.
            group.throughput(criterion::Throughput::Elements(vertices as u64));
            group.bench_with_input(
                BenchmarkId::new(
                    format!("realization_k{k}"),
                    format!("{}_V{vertices}", case.name),
                ),
                &case.name,
                |b, _| {
                    b.iter(|| {
                        for _ in 0..REALIZATIONS {
                            mc.sweep(&mut ctx);
                            mc.measure(&mut ctx);
                        }
                    });
                },
            );
        }
    }
    group.finish();
}

criterion_group!(benches, bench_kcore);
criterion_main!(benches);
