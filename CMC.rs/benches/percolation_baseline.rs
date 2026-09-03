mod percolation_support;

use std::time::Duration;

use cmc_rs::{OccupancyState, PercolationMode};
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use percolation_support::{cases, sample_and_analyze, P_BOND, P_SITE};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

fn bench_percolation_baseline(criterion: &mut Criterion) {
    for case in cases() {
        for mode in [
            PercolationMode::Site,
            PercolationMode::Bond,
            PercolationMode::SiteBond,
        ] {
            let mut occupancy = OccupancyState::new(&case.lattice, mode);
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x0042_454e_4348);
            let mut group = criterion.benchmark_group(format!(
                "percolation_reference/{}/V={}/E={}/p_site={}/p_bond={}",
                case.name,
                case.lattice.n_sites,
                case.lattice.n_edges(),
                P_SITE,
                P_BOND
            ));
            group.throughput(Throughput::Elements(1));
            group.bench_function(
                BenchmarkId::new("sample+analyze", mode.as_label()),
                |bencher| {
                    bencher.iter(|| black_box(sample_and_analyze(&case, &mut occupancy, &mut rng)));
                },
            );
            group.finish();
        }
    }
}

criterion_group! {
    name = benches;
    config = Criterion::default()
        .sample_size(10)
        .warm_up_time(Duration::from_millis(100))
        .measurement_time(Duration::from_millis(250));
    targets = bench_percolation_baseline
}
criterion_main!(benches);
