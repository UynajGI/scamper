//! R3 Edwards-Sokal exposure benchmark.
//!
//! Throughput of the full R3 cycle on the five shared topologies
//! (V = 4096): `es_sweep` is one Edwards-Sokal step — the conditional FK
//! bond draw given the current spins, the union-find cluster construction,
//! and the `BondSample` exposure of that actual draw; `es_sweep+measure`
//! adds the F3 analysis of the exposed sample plus the nine-scalar
//! observable recording through a Carlo.rs `Context` (the same composition
//! the F4 `static-mc-sweep+measure` IDs measure for the static family).
//! `beta = 0.5` with unit couplings gives `p_e = 1 - exp(-1) ~ 0.6327`, a
//! mid-occupation regime; thermalization is irrelevant to throughput but a
//! handful of warm sweeps precede Criterion as in the F4 bench.

#[allow(dead_code)]
mod percolation_support;

use std::time::Duration;

use carlo_rs::{Context, MonteCarlo};
use cmc_rs::{BoundaryQuery, EdwardsSokalMC, GraphView, ObservablePlan, UndirectedGraphView};
use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use percolation_support::cases;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

const BETA: f64 = 0.5;
const TIMED_SAMPLE_SEED: u64 = 0x4553_5345_4550;
const TIMED_INITIAL_SAMPLES: usize = 4;

fn bench_edwards_sokal(criterion: &mut Criterion) {
    for case in cases() {
        let queries = vec![BoundaryQuery::new(&case.lattice, &case.from, &case.to).unwrap()];
        let plan = ObservablePlan::all(queries);
        let mut adapter = EdwardsSokalMC::new(case.lattice.clone(), BETA, plan).unwrap();
        let mut mc_context = Context::new_with_binsize(
            Xoshiro256PlusPlus::seed_from_u64(TIMED_SAMPLE_SEED),
            0,
            1024,
        );
        for _ in 0..TIMED_INITIAL_SAMPLES {
            adapter.sweep(&mut mc_context);
            adapter.measure(&mut mc_context);
        }

        let edge_count = case.lattice.edge_count();
        let mut group = criterion.benchmark_group(format!(
            "edwards_sokal/{}/V={}/E={}/beta={BETA}",
            case.name,
            case.lattice.vertex_count(),
            edge_count,
        ));
        group.throughput(Throughput::Elements(edge_count as u64));
        group.bench_function("es_sweep (FK draw + exposure)", |bencher| {
            bencher.iter(|| {
                adapter.sweep(&mut mc_context);
                black_box(adapter.last_bond_sample());
            });
        });
        group.bench_function("es_sweep+measure (analyze exposed draw)", |bencher| {
            bencher.iter(|| {
                adapter.sweep(&mut mc_context);
                adapter.measure(&mut mc_context);
                black_box(&mut mc_context);
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
    targets = bench_edwards_sokal
}
criterion_main!(benches);
