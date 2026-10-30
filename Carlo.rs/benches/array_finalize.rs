//! C0 benchmark: element-preserving array finalize.
//!
//! Compares `Accumulator::finalize_array` (one pass over stored bins,
//! per-element statistics) against the cost of the legacy
//! element-averaging `finalize`, for array observables that downstream
//! adapters (percolation robustness curves, survival profiles) record.
//!
//! Run with: cargo bench -p carlo-rs --bench array_finalize

use carlo_rs::Accumulator;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};

fn populated(bin_capacity: usize, shape: &[usize], n_samples: usize) -> Accumulator {
    let flat: usize = shape.iter().product();
    let mut acc = Accumulator::with_shape(bin_capacity, shape);
    let mut sample = vec![0.0f64; flat];
    for s in 0..n_samples {
        for (i, value) in sample.iter_mut().enumerate() {
            *value = (s % 97) as f64 + i as f64 * 1e-3;
        }
        acc.add_array(&sample).expect("fixed shape");
    }
    acc
}

fn bench_array_finalize(c: &mut Criterion) {
    let mut group = c.benchmark_group("array_finalize");
    group.sample_size(30);

    for (elements, n_bins) in [(64usize, 128usize), (256, 128), (1024, 64)] {
        // Bins are pre-populated so both variants measure pure finalize cost.
        let acc = populated(8, &[elements], n_bins * 8);
        assert_eq!(acc.num_bins(), n_bins);

        group.bench_with_input(
            BenchmarkId::new("finalize_array", format!("E={elements}_B={n_bins}")),
            &acc,
            |b, acc| {
                b.iter(|| {
                    let out = acc.finalize_array().expect("array accumulator");
                    criterion::black_box(out);
                });
            },
        );

        group.bench_with_input(
            BenchmarkId::new(
                "finalize_legacy_collapse",
                format!("E={elements}_B={n_bins}"),
            ),
            &acc,
            |b, acc| {
                b.iter(|| {
                    let out = acc.finalize();
                    criterion::black_box(out);
                });
            },
        );
    }

    // Full pipeline check: sampling plus per-element finalize is the cost a
    // curve-producing adapter actually pays per replica.
    for elements in [64usize, 256] {
        let n_samples = 1024;
        group.bench_with_input(
            BenchmarkId::new("sample_and_finalize", format!("E={elements}_N={n_samples}")),
            &elements,
            |b, &elements| {
                b.iter_batched(
                    || (0.0f64, vec![0.0; elements]),
                    |(mut phase, mut sample)| {
                        let mut acc = Accumulator::with_shape(8, &[elements]);
                        for s in 0..n_samples {
                            phase += 0.37;
                            for (i, value) in sample.iter_mut().enumerate() {
                                *value = phase + i as f64 * 1e-3 + (s % 13) as f64;
                            }
                            acc.add_array(&sample).expect("fixed shape");
                        }
                        criterion::black_box(acc.finalize_array())
                    },
                    criterion::BatchSize::SmallInput,
                );
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_array_finalize);
criterion_main!(benches);
