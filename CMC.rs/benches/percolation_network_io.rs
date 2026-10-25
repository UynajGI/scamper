//! X1 network-document I/O benchmark.
//!
//! Measures both versioned formats (edge-list v1 text and multiplex JSON v1)
//! and the layer-to-CSR projection at the shared 16k-edge scale.
//!
//! Run with: cargo bench -p cmc-rs --features io --bench percolation_network_io

use cmc_rs::{EdgeListSpec, NetworkDocument, NetworkEdge, UndirectedGraphView};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;
use std::collections::HashSet;

const VERTICES: usize = 4_096;
const EDGES: usize = 16_384;

fn document() -> NetworkDocument {
    let mut document = NetworkDocument::new(VERTICES, ["main"]).expect("document");
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x1e57_1042);
    let mut seen: HashSet<(usize, usize)> = HashSet::with_capacity(EDGES);
    while seen.len() < EDGES {
        let source = rng.random_range(0..VERTICES);
        let target = rng.random_range(0..VERTICES);
        if source == target {
            continue;
        }
        let pair = (source.min(target), source.max(target));
        if seen.insert(pair) {
            document
                .push_edge(NetworkEdge::new(pair.0, pair.1, 1.0, [0]).expect("edge"))
                .expect("edge");
        }
    }
    document
}

fn bench_network_io(c: &mut Criterion) {
    let document = document();
    let edge_list = document.to_edge_list_v1().expect("serialize text");
    let json = document.to_json().expect("serialize json");
    eprintln!(
        "NETWORK_IO vertices={VERTICES} edges={EDGES} text_bytes={} json_bytes={}",
        edge_list.len(),
        json.len()
    );

    let mut group = c.benchmark_group("percolation_network_io");
    group.throughput(Throughput::Elements(EDGES as u64));
    group.bench_function(BenchmarkId::new("edge_list_v1", "parse"), |b| {
        b.iter(|| {
            let parsed = NetworkDocument::from_edge_list_v1(&edge_list).expect("parse");
            std::hint::black_box(parsed.edges().len());
        });
    });
    group.bench_function(BenchmarkId::new("edge_list_v1", "write"), |b| {
        b.iter(|| std::hint::black_box(document.to_edge_list_v1().expect("write")));
    });
    group.bench_function(BenchmarkId::new("json_v1", "parse"), |b| {
        b.iter(|| {
            let parsed = NetworkDocument::from_json(&json).expect("parse");
            std::hint::black_box(parsed.edges().len());
        });
    });
    group.bench_function(BenchmarkId::new("json_v1", "write"), |b| {
        b.iter(|| std::hint::black_box(document.to_json().expect("write")));
    });
    group.bench_function(BenchmarkId::new("layer_to_csr", "main"), |b| {
        b.iter(|| {
            let lattice = document
                .to_csr_layer(0, EdgeListSpec::simple(VERTICES))
                .expect("project");
            std::hint::black_box(lattice.edge_count());
        });
    });
    group.finish();
}

criterion_group!(benches, bench_network_io);
criterion_main!(benches);
