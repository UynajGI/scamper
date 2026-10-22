//! D2 benchmark: weak components, SCC, and reachability throughput.
//!
//! Measures one full analysis of each kind on a directed grid (V=16384,
//! A=40576 with back arcs) and a random digraph and a random digraph, with workspaces
//! prepared once so the timed region is the analysis itself.
//!
//! Run with: cargo bench -p cmc-rs --bench directed_analysis

use cmc_rs::{
    analyze_strong_components, analyze_weak_components, reachable, Arc, DirectedCsr,
    DirectedGraphView, GraphView, ReachabilityWorkspace, StrongComponentWorkspace,
    WeakComponentWorkspace,
};
use criterion::{criterion_group, criterion_main, Criterion};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

fn directed_grid(width: usize, height: usize) -> DirectedCsr {
    let index = |x: usize, y: usize| y * width + x;
    let mut arcs = Vec::with_capacity(3 * width * height);
    for y in 0..height {
        for x in 0..width {
            if x + 1 < width {
                arcs.push(Arc::new(index(x, y), index(x + 1, y)));
            }
            if y + 1 < height {
                arcs.push(Arc::new(index(x, y), index(x, y + 1)));
                if x % 2 == 0 && y > 0 {
                    // Occasional backward arcs create non-trivial SCCs.
                    arcs.push(Arc::new(index(x, y), index(x, y - 1)));
                }
            }
        }
    }
    DirectedCsr::from_arcs(width * height, arcs).expect("valid arcs")
}

fn random_digraph(vertices: usize, arc_count: usize, seed: u64) -> DirectedCsr {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let arcs = (0..arc_count)
        .map(|_| Arc::new(rng.random_range(0..vertices), rng.random_range(0..vertices)))
        .collect();
    DirectedCsr::from_arcs(vertices, arcs).expect("valid arcs")
}

fn bench_directed_analysis(c: &mut Criterion) {
    let cases: Vec<(&str, DirectedCsr)> = vec![
        ("grid", directed_grid(128, 128)),
        ("random", random_digraph(16_384, 49_152, 0xd2_b2)),
    ];
    for (name, graph) in &cases {
        let vertices = graph.vertex_count();
        let arcs = graph.arc_count();

        let mut weak = WeakComponentWorkspace::default();
        let mut strong = StrongComponentWorkspace::default();
        let mut reach = ReachabilityWorkspace::default();
        analyze_weak_components(graph, &mut weak).unwrap();
        analyze_strong_components(graph, &mut strong).unwrap();
        reachable(
            graph,
            graph.vertex_id(0).unwrap(),
            graph.vertex_id(vertices - 1).unwrap(),
            &mut reach,
        )
        .unwrap();

        let mut group =
            c.benchmark_group(format!("directed_analysis/{name}/V={vertices}/A={arcs}"));
        group.throughput(criterion::Throughput::Elements((vertices + arcs) as u64));
        group.bench_function("weak_components", |b| {
            b.iter(|| {
                let result = analyze_weak_components(graph, &mut weak).unwrap();
                std::hint::black_box(result.summary());
            });
        });
        group.bench_function("strong_components", |b| {
            b.iter(|| {
                let result = analyze_strong_components(graph, &mut strong).unwrap();
                std::hint::black_box(result.summary());
            });
        });
        // Far-corner queries visit (nearly) the whole reachable subgraph.
        group.bench_function("reachability_far_corner", |b| {
            b.iter(|| {
                std::hint::black_box(
                    reachable(
                        graph,
                        graph.vertex_id(0).unwrap(),
                        graph.vertex_id(vertices - 1).unwrap(),
                        &mut reach,
                    )
                    .unwrap(),
                )
            });
        });
        group.finish();
    }
}

criterion_group!(benches, bench_directed_analysis);
criterion_main!(benches);
