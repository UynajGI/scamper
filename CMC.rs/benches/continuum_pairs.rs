//! Continuum candidate-pair benchmark (P1).
//!
//! Records the two throughput figures of the per-analysis pipeline at unit
//! number density (`rho = 1`, so the box side is `N^(1/D)`):
//!
//! - `enumerate`: candidate pairs emitted per second — one full transient
//!   `CandidatePairGraph::pairs()` (coordinate copy + cell-list build +
//!   stencil walk) per iteration, the honest per-analysis cost under the
//!   no-persistent-CSR rule;
//! - `analyze`: complete snapshot analyses per second — enumeration plus
//!   the F3 component pass with an all-active site mask.
//!
//! The cutoff is set from a target mean coordination
//! `z = rho * V_D * r_c^D` (`V_D` is the unit-ball volume): at fixed
//! density the cutoff is the density knob of the pair criterion. Systems:
//! 2D and 3D, N = 4096 and 16384, z in {2, 6, 12}.
//!
//! Run with: cargo bench -p cmc-rs --bench continuum_pairs

use cmc_rs::{
    analyze, CandidatePairGraph, ComponentWorkspace, EdgeActivity, EdgeId, GraphView,
    OrthorhombicCell, VertexActivity, VertexId,
};
use criterion::measurement::WallTime;
use criterion::{
    criterion_group, criterion_main, BenchmarkGroup, BenchmarkId, Criterion, Throughput,
};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

struct AllVertices(usize);

impl VertexActivity for AllVertices {
    fn vertex_count(&self) -> usize {
        self.0
    }

    fn vertex_active(&self, _vertex: VertexId) -> bool {
        true
    }
}

struct AllEdges(usize);

impl EdgeActivity for AllEdges {
    fn edge_count(&self) -> usize {
        self.0
    }

    fn edge_active(&self, _edge: EdgeId) -> bool {
        true
    }
}

/// Flat coordinates, periodic cell, and cutoff for one benchmark system.
fn flat_config<const D: usize>(
    n: usize,
    target_coordination: f64,
    seed: u64,
) -> (Vec<f64>, OrthorhombicCell<D>, f64) {
    let side = (n as f64).powf(1.0 / D as f64);
    let unit_ball = match D {
        2 => std::f64::consts::PI,
        3 => (4.0 / 3.0) * std::f64::consts::PI,
        _ => panic!("benchmark supports 2D and 3D"),
    };
    let cutoff = (target_coordination / unit_ball).powf(1.0 / D as f64);
    let mut lengths = [0.0; D];
    lengths.fill(side);
    let cell = OrthorhombicCell::new(lengths).expect("benchmark cell is valid");
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut flat = Vec::with_capacity(n * D);
    for _ in 0..n {
        for _ in 0..D {
            flat.push(rng.random_range(0.0..side));
        }
    }
    (flat, cell, cutoff)
}

fn bench_dimension<const D: usize>(group: &mut BenchmarkGroup<'_, WallTime>, label: &str) {
    for (n, coordination) in [
        (4_096_usize, 2.0_f64),
        (4_096, 6.0),
        (4_096, 12.0),
        (16_384, 2.0),
        (16_384, 6.0),
        (16_384, 12.0),
    ] {
        let (flat, cell, cutoff) = flat_config::<D>(n, coordination, 0xa11ce);
        let adapter = CandidatePairGraph::new(&flat, &cell, cutoff).expect("valid system");
        let pair_count = adapter.pairs().unwrap().pair_count();
        let system = format!("{label}_n{n}_z{}", coordination as usize);

        group.throughput(Throughput::Elements(pair_count as u64));
        group.bench_with_input(
            BenchmarkId::new("enumerate", &system),
            &adapter,
            |b, adapter| {
                b.iter(|| adapter.pairs().unwrap().pair_count());
            },
        );

        let mut workspace = ComponentWorkspace::new();
        group.throughput(Throughput::Elements(1));
        group.bench_with_input(
            BenchmarkId::new("analyze", &system),
            &adapter,
            |b, adapter| {
                b.iter(|| {
                    let graph = adapter.pairs().unwrap();
                    let summary = analyze(
                        &graph,
                        &AllVertices(graph.vertex_count()),
                        &AllEdges(graph.pair_count()),
                        &[],
                        &mut workspace,
                    )
                    .expect("analysis succeeds")
                    .summary();
                    std::hint::black_box(summary);
                });
            },
        );

        eprintln!(
            "CONTINUUM_BENCH dim={label} n={n} z={coordination} cutoff={cutoff:.4} \
             pairs={pair_count} pairs_per_particle={:.2}",
            pair_count as f64 / n as f64,
        );
    }
}

fn bench_continuum_pairs(c: &mut Criterion) {
    let mut group = c.benchmark_group("continuum_pairs");
    bench_dimension::<2>(&mut group, "2d");
    bench_dimension::<3>(&mut group, "3d");
    group.finish();
}

criterion_group!(benches, bench_continuum_pairs);
criterion_main!(benches);
