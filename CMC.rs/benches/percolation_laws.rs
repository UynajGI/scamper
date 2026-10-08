use std::time::Duration;

use cmc_rs::{
    build_square, BondBernoulli, GraphView, MixedBernoulli, Probability, ProbabilityField,
    SiteBernoulli, StaticConfiguration, UndirectedGraphView,
};
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

const WIDTH: usize = 64;
const HEIGHT: usize = 64;

fn probability(value: f64) -> Probability {
    Probability::new(value).expect("benchmark probability")
}

fn bench_percolation_laws(criterion: &mut Criterion) {
    let graph = build_square(WIDTH, HEIGHT, false);
    let vertex_count = graph.vertex_count();
    let edge_count = graph.edge_count();

    for p in [0.01, 0.5, 0.99] {
        let site = SiteBernoulli::new(probability(p).into());
        let bond = BondBernoulli::new(probability(p).into());
        let mut site_configuration = StaticConfiguration::new(vertex_count, edge_count);
        let mut bond_configuration = StaticConfiguration::new(vertex_count, edge_count);
        let mut site_rng = Xoshiro256PlusPlus::seed_from_u64(0x0000_5349_5445);
        let mut bond_rng = Xoshiro256PlusPlus::seed_from_u64(0x0000_424f_4e44);

        let mut site_group = criterion.benchmark_group(format!(
            "percolation_laws/site/V={vertex_count}/E={edge_count}/p={p}"
        ));
        site_group.throughput(Throughput::Elements(vertex_count as u64));
        site_group.bench_function(BenchmarkId::new("sample", "uniform"), |bencher| {
            bencher.iter(|| {
                site.sample(
                    black_box(&graph),
                    black_box(&mut site_configuration),
                    black_box(&mut site_rng),
                )
                .expect("benchmark domains match");
            });
        });
        site_group.finish();

        let mut bond_group = criterion.benchmark_group(format!(
            "percolation_laws/bond/V={vertex_count}/E={edge_count}/p={p}"
        ));
        bond_group.throughput(Throughput::Elements(edge_count as u64));
        bond_group.bench_function(BenchmarkId::new("sample", "uniform"), |bencher| {
            bencher.iter(|| {
                bond.sample(
                    black_box(&graph),
                    black_box(&mut bond_configuration),
                    black_box(&mut bond_rng),
                )
                .expect("benchmark domains match");
            });
        });
        bond_group.finish();
    }

    for (index, (p_vertex, p_edge)) in [
        (0.01, 0.01),
        (0.5, 0.5),
        (0.99, 0.99),
        (0.01, 0.5),
        (0.99, 0.5),
        (0.5, 0.01),
        (0.5, 0.99),
    ]
    .into_iter()
    .enumerate()
    {
        let mixed = MixedBernoulli::new(probability(p_vertex).into(), probability(p_edge).into());
        let mut mixed_configuration = StaticConfiguration::new(vertex_count, edge_count);
        let mut mixed_rng = Xoshiro256PlusPlus::seed_from_u64(0x4d49_5845_4400 + index as u64);
        let mut mixed_group = criterion.benchmark_group(format!(
            "percolation_laws/mixed/V={vertex_count}/E={edge_count}/p_vertex={p_vertex}/p_edge={p_edge}"
        ));
        mixed_group.throughput(Throughput::Elements((vertex_count + edge_count) as u64));
        mixed_group.bench_function("sample", |bencher| {
            bencher.iter(|| {
                mixed
                    .sample(
                        black_box(&graph),
                        black_box(&mut mixed_configuration),
                        black_box(&mut mixed_rng),
                    )
                    .expect("benchmark domains match");
            });
        });
        mixed_group.finish();
    }

    let vertex_probabilities = (0..vertex_count)
        .map(|index| probability(if index % 2 == 0 { 0.25 } else { 0.75 }))
        .collect::<Vec<_>>();
    let edge_probabilities = (0..edge_count)
        .map(|index| probability(if index % 3 == 0 { 0.2 } else { 0.6 }))
        .collect::<Vec<_>>();

    // N1 uniform-versus-heterogeneous sampling comparison, one heterogeneous
    // field per family next to its uniform rows above.
    let heterogeneous_site = SiteBernoulli::new(ProbabilityField::Borrowed(&vertex_probabilities));
    let mut heterogeneous_site_configuration = StaticConfiguration::new(vertex_count, edge_count);
    let mut heterogeneous_site_rng = Xoshiro256PlusPlus::seed_from_u64(0x5349_5445_4845_5452);
    let mut heterogeneous_site_group = criterion.benchmark_group(format!(
        "percolation_laws/heterogeneous-site/V={vertex_count}/E={edge_count}"
    ));
    heterogeneous_site_group.throughput(Throughput::Elements(vertex_count as u64));
    heterogeneous_site_group.bench_function(
        BenchmarkId::new("sample", "heterogeneous"),
        |bencher| {
            bencher.iter(|| {
                heterogeneous_site
                    .sample(
                        black_box(&graph),
                        black_box(&mut heterogeneous_site_configuration),
                        black_box(&mut heterogeneous_site_rng),
                    )
                    .expect("benchmark domains match");
            });
        },
    );
    heterogeneous_site_group.finish();

    let heterogeneous_bond = BondBernoulli::new(ProbabilityField::Borrowed(&edge_probabilities));
    let mut heterogeneous_bond_configuration = StaticConfiguration::new(vertex_count, edge_count);
    let mut heterogeneous_bond_rng = Xoshiro256PlusPlus::seed_from_u64(0x424f_4e44_4845_5452);
    let mut heterogeneous_bond_group = criterion.benchmark_group(format!(
        "percolation_laws/heterogeneous-bond/V={vertex_count}/E={edge_count}"
    ));
    heterogeneous_bond_group.throughput(Throughput::Elements(edge_count as u64));
    heterogeneous_bond_group.bench_function(
        BenchmarkId::new("sample", "heterogeneous"),
        |bencher| {
            bencher.iter(|| {
                heterogeneous_bond
                    .sample(
                        black_box(&graph),
                        black_box(&mut heterogeneous_bond_configuration),
                        black_box(&mut heterogeneous_bond_rng),
                    )
                    .expect("benchmark domains match");
            });
        },
    );
    heterogeneous_bond_group.finish();

    let heterogeneous = MixedBernoulli::new(
        ProbabilityField::Borrowed(&vertex_probabilities),
        ProbabilityField::Borrowed(&edge_probabilities),
    );
    let mut heterogeneous_configuration = StaticConfiguration::new(vertex_count, edge_count);
    let mut heterogeneous_rng = Xoshiro256PlusPlus::seed_from_u64(0x4845_5445_524f);
    let mut heterogeneous_group = criterion.benchmark_group(format!(
        "percolation_laws/heterogeneous-mixed/V={vertex_count}/E={edge_count}"
    ));
    heterogeneous_group.throughput(Throughput::Elements((vertex_count + edge_count) as u64));
    heterogeneous_group.bench_function("sample", |bencher| {
        bencher.iter(|| {
            heterogeneous
                .sample(
                    black_box(&graph),
                    black_box(&mut heterogeneous_configuration),
                    black_box(&mut heterogeneous_rng),
                )
                .expect("benchmark domains match");
        });
    });
    heterogeneous_group.finish();
}

criterion_group! {
    name = benches;
    config = Criterion::default()
        .sample_size(10)
        .warm_up_time(Duration::from_millis(100))
        .measurement_time(Duration::from_millis(250));
    targets = bench_percolation_laws
}
criterion_main!(benches);
