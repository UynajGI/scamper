//! N3 benchmark: degree-ordered attack trajectories, plus the T1 adaptive
//! residual-degree order generation (lazy heap vs naive full rescan).
//!
//! One iteration = 8 complete attack trajectories (order construction plus
//! a full checkpoint walk with a fresh component analysis per checkpoint)
//! on the shared benchmark topologies, for both tie policies. The grid
//! step adapts so every case walks ~17 checkpoints; the Full grid costs
//! O(V * (V + E)) per trajectory and is intentionally not benchmarked.
//!
//! The T1 section isolates ORDER GENERATION on the Holme et al. (2002)
//! Fig. 7(e) configuration (BA, N = 10^4, m0 = 5, m = 3): the production
//! lazy max-heap (`DegreeAttackMC::sweep`, O((V + E) log V)) against the
//! naive full-rescan reference (every step recomputes every active
//! vertex's residual degree from adjacency, O(V * (V + E))). The naive
//! reference is deliberately NOT production code; it duplicates the
//! test-local reference of `tests/physics/percolation_attack.rs` (kept in
//! sync by the exact-parity gate there).
//!
//! Run with: cargo bench -p cmc-rs --bench percolation_attack

// Shared bench support items are dead in this target.
#![allow(dead_code)]

mod percolation_support;

use std::time::Duration;

use carlo_rs::{Context, MonteCarlo};
use cmc_rs::{AttackGrid, AttackPolicy, Bond, BondType, CsrLattice, DegreeAttackMC, GraphView};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use percolation_support::cases;
use rand::{RngExt, SeedableRng};
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
                AttackPolicy::ResidualDegreeDescendingStableId => "residual_stable_id",
                AttackPolicy::ResidualDegreeDescendingSeededKey => "residual_seeded_key",
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

/// Barabási–Albert graph (Holme et al. 2002 Sec. IV E): initial complete
/// graph on `m0` vertices; every later vertex attaches to `m` distinct
/// existing vertices drawn from a degree-proportional endpoint pool.
/// Identical construction to the test builder (documented there).
fn barabasi_albert_graph(n: usize, m0: usize, m: usize, seed: u64) -> CsrLattice {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut edges = Vec::with_capacity(m0 * (m0 - 1) / 2 + m * (n - m0));
    let mut pool: Vec<usize> = Vec::with_capacity(2 * edges.capacity());
    for left in 0..m0 {
        for right in (left + 1)..m0 {
            edges.push(Bond::new(left, right, BondType::Generic, 1.0));
            pool.push(left);
            pool.push(right);
        }
    }
    for vertex in m0..n {
        let mut targets = Vec::with_capacity(m);
        while targets.len() < m.min(vertex) {
            let target = pool[rng.random_range(0..pool.len())];
            if !targets.contains(&target) {
                targets.push(target);
            }
        }
        for target in targets {
            edges.push(Bond::new(vertex, target, BondType::Generic, 1.0));
            pool.push(vertex);
            pool.push(target);
        }
    }
    CsrLattice::from_edges(n, edges)
}

/// Naive full-rescan reference (bench-side duplicate of the test-local
/// reference): at EVERY step recompute every active vertex's residual
/// degree from adjacency and remove the maximum (stable ascending-ID ties).
fn naive_residual_order(adjacency: &[Vec<usize>]) -> Vec<usize> {
    let vertices = adjacency.len();
    let mut removed = vec![false; vertices];
    let mut order = Vec::with_capacity(vertices);
    for _ in 0..vertices {
        let mut victim = usize::MAX;
        let mut best_degree = 0_u32;
        for vertex in 0..vertices {
            if removed[vertex] {
                continue;
            }
            let residual = adjacency[vertex]
                .iter()
                .filter(|&&neighbor| !removed[neighbor])
                .count() as u32;
            if residual > best_degree || (residual == best_degree && vertex < victim) {
                best_degree = residual;
                victim = vertex;
            }
        }
        removed[victim] = true;
        order.push(victim);
    }
    order
}

/// T1 order-generation comparison on the Holme Fig. 7(e) configuration:
/// production lazy heap vs the naive full-rescan reference.
fn bench_adaptive_order_generation(c: &mut Criterion) {
    let lattice = barabasi_albert_graph(10_000, 5, 3, 0x0BA1_BE4C);
    let vertices = lattice.vertex_count();
    let adjacency: Vec<Vec<usize>> = {
        let mut adjacency = vec![Vec::new(); vertices];
        for bond in &lattice.edges {
            adjacency[bond.source].push(bond.target);
            adjacency[bond.target].push(bond.source);
        }
        adjacency
    };

    let mut group = c.benchmark_group("percolation_attack_order");
    // The naive reference costs O(V * (V + E)) per trajectory (~seconds at
    // this size); keep criterion's sample count at the minimum so the
    // comparison stays a measurement, not a marathon.
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(20));

    let mut mc = DegreeAttackMC::new(
        lattice,
        AttackPolicy::ResidualDegreeDescendingSeededKey,
        AttackGrid::Every { step: vertices },
    )
    .expect("bench adapter");
    let mut ctx = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(0x7e51), 0, 8);
    group.bench_with_input(
        BenchmarkId::new("lazy_heap", format!("ba_V{vertices}_m0_5_m_3")),
        &(),
        |b, _| {
            b.iter(|| {
                for _ in 0..TRAJECTORIES {
                    mc.sweep(&mut ctx);
                    std::hint::black_box(mc.removal_order());
                }
            });
        },
    );
    group.bench_with_input(
        BenchmarkId::new("naive_rescan", format!("ba_V{vertices}_m0_5_m_3")),
        &(),
        |b, _| {
            b.iter(|| {
                for _ in 0..TRAJECTORIES {
                    std::hint::black_box(naive_residual_order(&adjacency));
                }
            });
        },
    );
    group.finish();
}

criterion_group!(
    benches,
    bench_degree_attack,
    bench_adaptive_order_generation
);
criterion_main!(benches);
