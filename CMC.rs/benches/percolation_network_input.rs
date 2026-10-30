//! Network-input benchmark: text edge-list ingestion into owned CSR.
//!
//! The roadmap defers any binary/memory-mapped format until million-edge
//! ingestion evidence justifies one; this bench records that evidence —
//! parse+build throughput at 16k edges (the shared benchmark scale) and at
//! 1M edges.
//!
//! Run with: cargo bench -p cmc-rs --bench percolation_network_input

use cmc_rs::{csr_from_edge_list, EdgeListSpec, UndirectedGraphView};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;
use std::collections::HashSet;

/// Deterministic G(n, m) edge-list text with distinct non-loop edges.
fn edge_list_text(vertices: usize, edge_count: usize, seed: u64) -> String {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut seen: HashSet<(usize, usize)> = HashSet::with_capacity(edge_count);
    let mut text = String::with_capacity(edge_count * 14);
    while seen.len() < edge_count {
        let left = rng.random_range(0..vertices);
        let right = rng.random_range(0..vertices);
        if left != right {
            let pair = if left < right {
                (left, right)
            } else {
                (right, left)
            };
            if seen.insert(pair) {
                text.push_str(&format!("{} {}\n", pair.0, pair.1));
            }
        }
    }
    text
}

fn bench_network_input(c: &mut Criterion) {
    let mut group = c.benchmark_group("percolation_network_input");
    for (label, vertices, edges) in [
        ("er16k", 4_096_usize, 16_384_usize),
        ("er1m", 250_000, 1_000_000),
    ] {
        let text = edge_list_text(vertices, edges, 0x1e57);
        let approx_bytes = text.len();
        group.throughput(criterion::Throughput::Elements(edges as u64));
        group.bench_with_input(
            BenchmarkId::new("parse_and_build", label),
            &(text, vertices),
            |b, (text, vertices)| {
                b.iter(|| {
                    let lattice = csr_from_edge_list(EdgeListSpec::simple(*vertices), text)
                        .expect("generated edge list parses");
                    std::hint::black_box(lattice.edge_count());
                });
            },
        );
        eprintln!(
            "NETWORK_INPUT text={label} edges={edges} vertices={vertices} \
             text_bytes={approx_bytes} bytes_per_edge={:.1}",
            approx_bytes as f64 / edges as f64
        );
    }
    group.finish();
}

criterion_group!(benches, bench_network_input);
criterion_main!(benches);
