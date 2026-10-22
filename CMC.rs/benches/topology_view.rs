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

/// D1 gate: generic directed adjacency scans must stay within the same
/// noise band as direct field access on the owned representation.
fn direct_directed_arc_scan(graph: &cmc_rs::DirectedCsr) -> usize {
    graph.arcs.iter().fold(0usize, |sum, arc| {
        sum.wrapping_add(arc.source).wrapping_add(arc.target)
    })
}

fn generic_directed_arc_scan<G: cmc_rs::DirectedGraphView>(graph: &G) -> usize {
    graph.arc_ids().fold(0usize, |sum, arc| {
        let [source, target] = graph.arc_endpoints(arc);
        sum.wrapping_add(source.index())
            .wrapping_add(target.index())
    })
}

fn direct_directed_out_scan(graph: &cmc_rs::DirectedCsr) -> usize {
    (0..graph.vertex_count).fold(0usize, |sum, vertex| {
        let (start, end) = (graph.out_offsets[vertex], graph.out_offsets[vertex + 1]);
        (start..end).fold(sum, |sum, slot| {
            sum.wrapping_add(graph.out_neighbors[slot])
                .wrapping_add(graph.out_arc_ids[slot])
        })
    })
}

fn generic_directed_out_scan<G: cmc_rs::DirectedOutGraphView>(graph: &G) -> usize {
    graph.vertex_ids().fold(0usize, |sum, vertex| {
        graph.out_arcs(vertex).fold(sum, |sum, incidence| {
            sum.wrapping_add(incidence.neighbor.index())
                .wrapping_add(incidence.arc.index())
        })
    })
}

fn direct_directed_in_scan(graph: &cmc_rs::DirectedCsr) -> usize {
    (0..graph.vertex_count).fold(0usize, |sum, vertex| {
        let (start, end) = (graph.in_offsets[vertex], graph.in_offsets[vertex + 1]);
        (start..end).fold(sum, |sum, slot| {
            sum.wrapping_add(graph.in_neighbors[slot])
                .wrapping_add(graph.in_arc_ids[slot])
        })
    })
}

fn generic_directed_in_scan<G: cmc_rs::DirectedInGraphView>(graph: &G) -> usize {
    graph.vertex_ids().fold(0usize, |sum, vertex| {
        graph.in_arcs(vertex).fold(sum, |sum, incidence| {
            sum.wrapping_add(incidence.neighbor.index())
                .wrapping_add(incidence.arc.index())
        })
    })
}

fn bench_directed_topology_view(criterion: &mut Criterion) {
    // Orient every square bond +x/+y, plus reverse arcs on every other
    // column to exercise mixed in/out degrees.
    let square = build_square(WIDTH, HEIGHT, false);
    let mut arcs = Vec::with_capacity(square.n_edges() * 2);
    for edge in &square.edges {
        let (source, target) = if edge.source < edge.target {
            (edge.source, edge.target)
        } else {
            (edge.target, edge.source)
        };
        arcs.push(cmc_rs::Arc::new(source, target));
        if source % 2 == 0 {
            arcs.push(cmc_rs::Arc::new(target, source));
        }
    }
    let directed = cmc_rs::DirectedCsr::from_arcs(square.n_sites, arcs).expect("valid arcs");
    let endpoints: Vec<_> = directed
        .arcs
        .iter()
        .map(|arc| [arc.source, arc.target])
        .collect();
    let borrowed = cmc_rs::BorrowedDirectedCsr::new(
        &directed.out_offsets,
        &directed.out_neighbors,
        &directed.out_arc_ids,
        &endpoints,
        &directed.in_offsets,
        &directed.in_neighbors,
        &directed.in_arc_ids,
    )
    .expect("builder output must validate");

    let arc_count = cmc_rs::DirectedGraphView::arc_count(&directed);
    assert_eq!(
        direct_directed_arc_scan(&directed),
        generic_directed_arc_scan(&directed)
    );
    assert_eq!(
        direct_directed_out_scan(&directed),
        generic_directed_out_scan(&directed)
    );
    assert_eq!(
        direct_directed_out_scan(&directed),
        generic_directed_out_scan(&borrowed)
    );
    assert_eq!(
        direct_directed_in_scan(&directed),
        generic_directed_in_scan(&directed)
    );
    assert_eq!(
        direct_directed_in_scan(&directed),
        generic_directed_in_scan(&borrowed)
    );

    for (label, elements) in [
        ("arc_scan", arc_count),
        ("out_scan", arc_count),
        ("in_scan", arc_count),
    ] {
        let mut group = criterion.benchmark_group(format!(
            "topology_view/directed_{label}/V={}/A={}",
            directed.vertex_count, arc_count
        ));
        group.throughput(Throughput::Elements(elements as u64));
        match label {
            "arc_scan" => {
                group.bench_function("direct_owned", |bencher| {
                    bencher.iter(|| black_box(direct_directed_arc_scan(black_box(&directed))));
                });
                group.bench_function("generic_owned", |bencher| {
                    bencher.iter(|| black_box(generic_directed_arc_scan(black_box(&directed))));
                });
                group.bench_function("generic_borrowed", |bencher| {
                    bencher.iter(|| black_box(generic_directed_arc_scan(black_box(&borrowed))));
                });
            }
            "out_scan" => {
                group.bench_function("direct_owned", |bencher| {
                    bencher.iter(|| black_box(direct_directed_out_scan(black_box(&directed))));
                });
                group.bench_function("generic_owned", |bencher| {
                    bencher.iter(|| black_box(generic_directed_out_scan(black_box(&directed))));
                });
                group.bench_function("generic_borrowed", |bencher| {
                    bencher.iter(|| black_box(generic_directed_out_scan(black_box(&borrowed))));
                });
            }
            _ => {
                group.bench_function("direct_owned", |bencher| {
                    bencher.iter(|| black_box(direct_directed_in_scan(black_box(&directed))));
                });
                group.bench_function("generic_owned", |bencher| {
                    bencher.iter(|| black_box(generic_directed_in_scan(black_box(&directed))));
                });
                group.bench_function("generic_borrowed", |bencher| {
                    bencher.iter(|| black_box(generic_directed_in_scan(black_box(&borrowed))));
                });
            }
        }
        group.finish();
    }
}

criterion_group! {
    name = benches;
    config = Criterion::default()
        .sample_size(20)
        .warm_up_time(Duration::from_millis(250))
        .measurement_time(Duration::from_secs(1));
    targets = bench_topology_view, bench_directed_topology_view
}
criterion_main!(benches);
