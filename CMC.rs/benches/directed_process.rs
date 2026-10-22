//! D3 benchmark: layered directed-percolation realizations.
//!
//! One complete space-time realization per iteration (T layers over a
//! periodic chain, single seed), with the frontier bitset pair as the only
//! process state.
//!
//! Run with: cargo bench -p cmc-rs --bench directed_process

use carlo_rs::{Context, MonteCarlo};
use cmc_rs::{build_chain, DpSeedPolicy, LayeredDirectedPercolationMC};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

fn bench_directed_process(c: &mut Criterion) {
    let mut group = c.benchmark_group("directed_process");
    for (label, width, layers) in [
        ("chain512_T1200", 512_usize, 1200_usize),
        ("chain1024_T2000", 1024, 2000),
    ] {
        let center = width / 2;
        let mut process = LayeredDirectedPercolationMC::new(
            build_chain(width, true),
            0.644_700_1,
            layers,
            DpSeedPolicy::Single { site: center },
        )
        .unwrap();
        let mut ctx = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(0x0d38), 0, 64);
        // Site-layer operations per realization (the throughput unit).
        group.throughput(criterion::Throughput::Elements((width * layers) as u64));
        group.bench_with_input(BenchmarkId::new("realization", label), &label, |b, _| {
            b.iter(|| {
                for _ in 0..8 {
                    process.sweep(&mut ctx);
                    std::hint::black_box(process.layer_active_counts());
                }
            });
        });
    }
    group.finish();
}

criterion_group!(benches, bench_directed_process);
criterion_main!(benches);
