//! N3 benchmark: degree-ordered attack trajectories.
//!
//! One iteration = 8 complete attack trajectories (order construction plus
//! a full checkpoint walk with a fresh component analysis per checkpoint)
//! on the shared benchmark topologies, for both tie policies. The grid
//! step adapts so every case walks ~17 checkpoints; the Full grid costs
//! O(V * (V + E)) per trajectory and is intentionally not benchmarked.
//!
//! Run with: cargo bench -p cmc-rs --bench percolation_attack

// Shared bench support items are dead in this target.
#![allow(dead_code)]

mod percolation_support;

use carlo_rs::{Context, MonteCarlo};
use cmc_rs::{AttackGrid, AttackPolicy, DegreeAttackMC, GraphView};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use percolation_support::cases;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

const TRAJECTORIES: usize = 8;

fn bench_degree_attack(c: &mut Criterion) {
    let mut group = c.benchmark_group("percolation_attack");
    for case in cases() {
        let vertices = case.lattice.vertex_count();
        let step = (vertices / 16).max(1);
        for policy in [
            AttackPolicy::DegreeDescendingStableId,
            AttackPolicy::DegreeDescendingSeededShuffle,
        ] {
            let label = match policy {
                AttackPolicy::DegreeDescendingStableId => "stable_id",
                AttackPolicy::DegreeDescendingSeededShuffle => "seeded_shuffle",
            };
            let mut mc =
                DegreeAttackMC::new(case.lattice.clone(), policy, AttackGrid::Every { step })
                    .expect("bench adapter");
            let checkpoints = mc.checkpoints().len();
            let mut ctx =
                Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(0xa11a), 0, 64);
            group.bench_with_input(
                BenchmarkId::new(
                    format!("{label}_ck{checkpoints}"),
                    format!("{}_V{vertices}", case.name),
                ),
                &case.name,
                |b, _| {
                    b.iter(|| {
                        for _ in 0..TRAJECTORIES {
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

criterion_group!(benches, bench_degree_attack);
criterion_main!(benches);
