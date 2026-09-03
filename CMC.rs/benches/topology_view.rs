use std::time::Duration;

use cmc_rs::{build_square, BorrowedUndirectedCsr, CsrLattice, UndirectedGraphView};
use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};

const WIDTH: usize = 128;
const HEIGHT: usize = 128;

#[inline]
fn direct_edge_scan(graph: &CsrLattice) -> usize {
    graph.edges.iter().fold(0usize, |sum, edge| {
        sum.wrapping_add(edge.source).wrapping_add(edge.target)
    })
}

#[inline]
fn generic_edge_scan<G: UndirectedGraphView>(graph: &G) -> usize {
    (0..graph.edge_count()).fold(0usize, |sum, edge| {
        let [left, right] = graph.edge_endpoints(graph.edge_id(edge).expect("edge in dense range"));
        sum.wrapping_add(left.index()).wrapping_add(right.index())
    })
}

#[inline]
fn direct_incidence_scan(graph: &CsrLattice) -> usize {
    (0..graph.n_sites).fold(0usize, |sum, vertex| {
        graph.incidences(vertex).fold(sum, |sum, (neighbor, edge)| {
            sum.wrapping_add(neighbor).wrapping_add(edge)
        })
    })
}

#[inline]
fn generic_incidence_scan<G: UndirectedGraphView>(graph: &G) -> usize {
    graph.vertex_ids().fold(0usize, |sum, vertex| {
        graph.incidences(vertex).fold(sum, |sum, incidence| {
            sum.wrapping_add(incidence.neighbor.index())
                .wrapping_add(incidence.edge.index())
        })
    })
}

fn bench_topology_view(criterion: &mut Criterion) {
    let graph = build_square(WIDTH, HEIGHT, false);
    let endpoints: Vec<_> = graph
        .edges
        .iter()
        .map(|edge| [edge.source, edge.target])
        .collect();
    let borrowed = BorrowedUndirectedCsr::new(
        &graph.offsets,
        &graph.neighbors,
        &graph.edge_ids,
        &endpoints,
    )
    .expect("builder output must validate as borrowed CSR");

    let direct_edges = direct_edge_scan(&graph);
    assert_eq!(direct_edges, generic_edge_scan(&graph));
    assert_eq!(direct_edges, generic_edge_scan(&borrowed));
    let direct_incidences = direct_incidence_scan(&graph);
    assert_eq!(direct_incidences, generic_incidence_scan(&graph));
    assert_eq!(direct_incidences, generic_incidence_scan(&borrowed));

    let mut edges = criterion.benchmark_group(format!(
        "topology_view/edge_scan/V={}/E={}",
        graph.n_sites,
        graph.n_edges()
    ));
    edges.throughput(Throughput::Elements(graph.n_edges() as u64));
    edges.bench_function("direct_owned", |bencher| {
        bencher.iter(|| black_box(direct_edge_scan(black_box(&graph))));
    });
    edges.bench_function("generic_owned", |bencher| {
        bencher.iter(|| black_box(generic_edge_scan(black_box(&graph))));
    });
    edges.bench_function("generic_borrowed", |bencher| {
        bencher.iter(|| black_box(generic_edge_scan(black_box(&borrowed))));
    });
    edges.finish();

    let mut incidences = criterion.benchmark_group(format!(
        "topology_view/incidence_scan/V={}/I={}",
        graph.n_sites, graph.n_bonds
    ));
    incidences.throughput(Throughput::Elements(graph.n_bonds as u64));
    incidences.bench_function("direct_owned", |bencher| {
        bencher.iter(|| black_box(direct_incidence_scan(black_box(&graph))));
    });
    incidences.bench_function("generic_owned", |bencher| {
        bencher.iter(|| black_box(generic_incidence_scan(black_box(&graph))));
    });
    incidences.bench_function("generic_borrowed", |bencher| {
        bencher.iter(|| black_box(generic_incidence_scan(black_box(&borrowed))));
    });
    incidences.finish();
}

criterion_group! {
    name = benches;
    config = Criterion::default()
        .sample_size(20)
        .warm_up_time(Duration::from_millis(250))
        .measurement_time(Duration::from_secs(1));
    targets = bench_topology_view
}
criterion_main!(benches);
