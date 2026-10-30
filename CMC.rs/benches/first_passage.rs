//! F1 benchmark: single-source weighted shortest paths — production
//! Dijkstra versus an independent test-side Bellman-Ford, plus Dijkstra
//! throughput and the fixed-average-degree size sweep (doubling V should
//! roughly double the time — the O((V+E) log V) complexity sanity check).
//!
//! The Bellman-Ford here is a from-scratch edge-list reference (the parity
//! gate lives in tests/physics/first_passage.rs); it reuses one distance
//! buffer across calls so the comparison measures the algorithm, not the
//! allocator, and its early exit (no change in a round) is the honest
//! sparse-graph Bellman-Ford — O(rounds * E) with rounds the hop-diameter
//! of the reached ball, never the textbook V*E worst case on these graphs.
//!
//! Run with: cargo bench -p cmc-rs --bench first_passage

use cmc_rs::{
    build_square, first_passage_run, Bond, BondType, CsrLattice, FirstPassageWorkspace,
    GeodesicTracking, GraphView, UndirectedGraphView,
};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

/// Sparse ER multigraph: `m` uniform random pairs (self-loops excluded,
/// duplicates kept), weights uniform [0,1). Returns the production graph,
/// the `Bond.weight` weight slice, and the plain edge triples the
/// Bellman-Ford reference consumes.
fn sparse_er(n: usize, m: usize, seed: u64) -> (CsrLattice, Vec<f64>, Vec<(usize, usize, f64)>) {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut edges = Vec::with_capacity(m);
    while edges.len() < m {
        let u = rng.random_range(0..n);
        let v = rng.random_range(0..n);
        if u != v {
            edges.push((u, v, rng.random::<f64>()));
        }
    }
    let lattice = CsrLattice::from_edges(
        n,
        edges
            .iter()
            .map(|&(source, target, weight)| Bond::new(source, target, BondType::Generic, weight))
            .collect(),
    );
    let weights: Vec<f64> = lattice.edges.iter().map(|bond| bond.weight).collect();
    (lattice, weights, edges)
}

/// Independent Bellman-Ford with a reusable distance buffer (zero
/// allocation in steady state; `clear` + `resize` writes only).
fn bellman_ford_reused(
    vertices: usize,
    edges: &[(usize, usize, f64)],
    dist: &mut Vec<f64>,
    source: usize,
) {
    dist.clear();
    dist.resize(vertices, f64::INFINITY);
    dist[source] = 0.0;
    for _ in 0..vertices {
        let mut changed = false;
        for &(u, v, w) in edges {
            if dist[u] + w < dist[v] {
                dist[v] = dist[u] + w;
                changed = true;
            }
            if dist[v] + w < dist[u] {
                dist[u] = dist[v] + w;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
}

/// One comparison case: the production graph, its `Bond.weight` slice, and
/// the plain edge triples the Bellman-Ford reference consumes.
struct CompareCase {
    label: &'static str,
    lattice: CsrLattice,
    weights: Vec<f64>,
    edges: Vec<(usize, usize, f64)>,
}

/// Dijkstra vs Bellman-Ford on two large sparse graphs where their
/// asymptotics point opposite ways:
/// - a small-diameter ER multigraph (V = 100k, E = 300k, hop diameter ~7):
///   early-exit Bellman-Ford converges in a handful of sequential rounds and
///   is expected to WIN here — recorded honestly, like the T2 sub-1
///   speedups;
/// - a large-diameter open 200x200 grid with random weights (V = 40k,
///   E ~ 79.6k, hop diameter 398): Bellman-Ford needs ~diameter rounds, so
///   Dijkstra is expected to win — the percolation-relevant regime (F2 runs
///   corridors on Z^2).
fn bench_dijkstra_vs_bellman_ford(c: &mut Criterion) {
    let mut group = c.benchmark_group("first_passage_dijkstra_vs_bellman_ford");

    let er = {
        let (lattice, weights, edges) = sparse_er(100_000, 300_000, 0x00F1_5350);
        CompareCase {
            label: "er_multigraph",
            lattice,
            weights,
            edges,
        }
    };
    let grid = {
        let base = build_square(200, 200, false);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x00F1_6721_4400);
        let bonds = base
            .edges
            .iter()
            .map(|bond| Bond::new(bond.source, bond.target, bond.kind, rng.random::<f64>()))
            .collect();
        let lattice = CsrLattice::from_edges(base.n_sites, bonds);
        let weights: Vec<f64> = lattice.edges.iter().map(|bond| bond.weight).collect();
        let edges: Vec<(usize, usize, f64)> = lattice
            .edges
            .iter()
            .map(|bond| (bond.source, bond.target, bond.weight))
            .collect();
        CompareCase {
            label: "open_grid_random_weights",
            lattice,
            weights,
            edges,
        }
    };

    for case in [&er, &grid] {
        let CompareCase {
            label,
            lattice,
            weights,
            edges,
        } = case;
        let n = lattice.vertex_count();
        let m = lattice.edge_count();
        let source = n / 2;
        let mut workspace = FirstPassageWorkspace::new();
        let outcome = first_passage_run(
            lattice,
            weights,
            source,
            GeodesicTracking::Predecessors,
            &mut workspace,
        )
        .expect("bench production pass");
        let mut bf_dist = Vec::with_capacity(n);
        bellman_ford_reused(n, edges, &mut bf_dist, source);
        let mismatches = workspace
            .distances()
            .iter()
            .zip(bf_dist.iter())
            .filter(|(&produced, &expected)| {
                produced.is_finite() != expected.is_finite()
                    || (expected.is_finite()
                        && (produced - expected).abs() > 1e-12 * expected.abs().max(1.0))
            })
            .count();
        eprintln!(
            "first passage bench: {label} V={n} E={m} reached={} bf_mismatches={mismatches}",
            outcome.reached
        );

        group.bench_function(BenchmarkId::new("dijkstra_production", label), |b| {
            b.iter(|| {
                std::hint::black_box(
                    first_passage_run(
                        lattice,
                        weights,
                        source,
                        GeodesicTracking::Predecessors,
                        &mut workspace,
                    )
                    .expect("bench production pass"),
                )
            })
        });
        group.bench_function(BenchmarkId::new("bellman_ford_reference", label), |b| {
            b.iter(|| {
                bellman_ford_reused(n, edges, &mut bf_dist, source);
                std::hint::black_box(&bf_dist);
            })
        });
    }
    group.finish();
}

/// Dijkstra throughput (vertices/s) across doubling sizes at fixed average
/// degree 3 (E = 3V): the complexity sanity check — doubling V at fixed
/// degree should roughly double the mean time.
fn bench_dijkstra_throughput(c: &mut Criterion) {
    let mut group = c.benchmark_group("first_passage_dijkstra");
    for n in [25_000_usize, 50_000, 100_000, 200_000, 400_000] {
        let (lattice, weights, _edges) = sparse_er(n, 3 * n, 0x00F1_5448 + n as u64);
        let source = n / 2;
        let mut workspace = FirstPassageWorkspace::new();
        first_passage_run(
            &lattice,
            &weights,
            source,
            GeodesicTracking::Predecessors,
            &mut workspace,
        )
        .expect("bench warmup pass");
        group.throughput(criterion::Throughput::Elements(n as u64));
        group.bench_with_input(
            BenchmarkId::new("single_source", format!("V{n}")),
            &n,
            |b, &_| {
                b.iter(|| {
                    std::hint::black_box(
                        first_passage_run(
                            &lattice,
                            &weights,
                            source,
                            GeodesicTracking::Predecessors,
                            &mut workspace,
                        )
                        .expect("bench pass"),
                    )
                })
            },
        );
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_dijkstra_vs_bellman_ford,
    bench_dijkstra_throughput
);
criterion_main!(benches);
