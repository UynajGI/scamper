//! G3 benchmark: invasion growth steps per second.
//!
//! One iteration = REALIZATIONS complete invasion realizations on the five
//! shared benchmark topologies at V = 4096 (site and bond modes): each
//! realization freezes a uniform threshold field (the adapter's documented
//! dense-order draw contract), grows from an explicit seed through the
//! binary-heap minimum-threshold frontier until frontier exhaustion, and
//! drives the G1 substrate with every accepted operation. The throughput
//! unit is accepted growth steps.
//!
//! The heap-versus-complete-frontier-scan question is a CORRECTNESS
//! comparison (tests/physics/invasion_percolation.rs), not a benchmark
//! axis: the scan is O(frontier) per step by construction and is not
//! production code.
//!
//! Run with: cargo bench -p cmc-rs --bench invasion_percolation

// Shared bench-support items are dead in this target.
#![allow(dead_code)]

mod percolation_support;

use cmc_rs::{GraphView, InvasionMode, InvasionPercolationMC, InvasionStop, InvasionTrapping};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use percolation_support::cases;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

const REALIZATIONS: usize = 4;

fn bench_invasion_steps(c: &mut Criterion) {
    let mut group = c.benchmark_group("invasion_percolation");
    for case in cases() {
        let vertices = case.lattice.vertex_count();
        // A fixed-length C0 shape is required by the adapter; 4096 covers
        // the largest possible GrowAll cluster (V - 1 growth steps).
        for mode in [InvasionMode::Site, InvasionMode::Bond] {
            let mut adapter = InvasionPercolationMC::new(
                case.lattice.clone(),
                mode,
                vec![0],
                InvasionStop::GrowAll,
                InvasionTrapping::NonTrapping,
                vertices.saturating_sub(1).max(1),
            )
            .expect("bench adapter constructs");
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x6333);
            adapter.realize(&mut rng); // warm every buffer once
            let expected = adapter.last_outcome().accepted_total as u64;
            eprintln!(
                "invasion bench: {} {} V={} E={} accepted_steps_per_realization={} \
                 max_accepted_threshold={:.4}",
                case.name,
                mode_label(mode),
                vertices,
                case.lattice.n_edges(),
                expected,
                adapter.last_outcome().max_accepted_threshold
            );
            group.throughput(criterion::Throughput::Elements(
                REALIZATIONS as u64 * expected,
            ));
            group.bench_with_input(
                BenchmarkId::new(
                    "steps",
                    format!("{}_{}_V{vertices}", case.name, mode_label(mode)),
                ),
                &case.name,
                |b, _| {
                    b.iter(|| {
                        for _ in 0..REALIZATIONS {
                            adapter.realize(&mut rng);
                            std::hint::black_box(adapter.last_outcome());
                        }
                    });
                },
            );
        }
    }
    group.finish();
}

fn mode_label(mode: InvasionMode) -> &'static str {
    match mode {
        InvasionMode::Site => "site",
        InvasionMode::Bond => "bond",
    }
}

criterion_group!(benches, bench_invasion_steps);
criterion_main!(benches);
