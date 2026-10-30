//! N2 benchmark: random-failure sampling with degree observables.
//!
//! Measures `RandomFailureMC` (site law + labeled component analysis + six
//! degree arrays) against a `StaticPercolationMC` scalar baseline (no
//! labels, no arrays) on the same topologies, isolating the cost of the
//! degree-conditioned observables.
//!
//! Run with: cargo bench -p cmc-rs --bench percolation_network

// Shared bench support items are dead in this target.
#![allow(dead_code)]

mod percolation_support;

use carlo_rs::{Context, MonteCarlo};
use cmc_rs::{ObservablePlan, RandomFailureMC, StaticLaw, StaticObservable, StaticPercolationMC};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use percolation_support::{cases, P_SITE};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

const SWEEPS: usize = 64;

fn bench_random_failure(c: &mut Criterion) {
    let mut group = c.benchmark_group("percolation_network");
    for case in cases() {
        // Degree observables: full sweep+measure pipeline of the facade.
        let mut facade = RandomFailureMC::new(
            case.lattice.clone(),
            P_SITE,
            ObservablePlan::all(Vec::new()),
        )
        .expect("bench facade");
        let mut ctx = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(0x4e32), 0, 1);
        group.bench_with_input(
            BenchmarkId::new("random_failure_degree", case.name),
            &case.name,
            |b, _| {
                b.iter(|| {
                    for _ in 0..SWEEPS {
                        facade.sweep(&mut ctx);
                        facade.measure(&mut ctx);
                    }
                });
            },
        );

        // Scalar baseline: same sampler and analysis without labels/arrays.
        let law = StaticLaw::site(P_SITE).expect("bench probability");
        let plan = ObservablePlan::new(
            [
                StaticObservable::GiantFraction,
                StaticObservable::LargestSize,
            ],
            Vec::new(),
        );
        let mut adapter =
            StaticPercolationMC::new(case.lattice.clone(), law, plan).expect("bench adapter");
        let mut ctx = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(0x4e32), 0, 1);
        group.bench_with_input(
            BenchmarkId::new("scalar_baseline", case.name),
            &case.name,
            |b, _| {
                b.iter(|| {
                    for _ in 0..SWEEPS {
                        adapter.sweep(&mut ctx);
                        adapter.measure(&mut ctx);
                    }
                });
            },
        );
    }
    group.finish();
}

criterion_group!(benches, bench_random_failure);
criterion_main!(benches);
