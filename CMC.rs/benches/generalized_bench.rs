use cmc_rs::{
    build_chain, canonical_reweight, enumerate_ising_density_of_states, Algorithm, BinnedAxis,
    IsingModel, LogDensityOfStates, MacrostateAxis, SimulationPhase, System, WangLandauConfig,
    WangLandauCore,
};
use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use rand::SeedableRng;

const ISING_SITES: usize = 12;
const REWEIGHT_BINS: usize = 4096;

fn bench_wang_landau_sweep(criterion: &mut Criterion) {
    let lattice = build_chain(ISING_SITES, true);
    let model = IsingModel::new(1.0);
    let axis = enumerate_ising_density_of_states(&lattice, &model)
        .unwrap()
        .axis()
        .unwrap();
    let config = WangLandauConfig {
        flatness_check_interval: u64::MAX,
        max_adaptation_sweeps: 0,
        ..Default::default()
    };
    let mut kernel = WangLandauCore::new(axis, config).unwrap();
    let mut system = System::new(lattice, 1, 1.0, 0.0);
    system.recompute_energy(&model);
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(0x0057_4C44_4F53);

    let mut group = criterion.benchmark_group("wang_landau");
    group.throughput(Throughput::Elements(ISING_SITES as u64));
    group.bench_function("local_attempts", |bencher| {
        bencher.iter(|| {
            kernel.sweep_with_phase(
                &mut system,
                &model,
                &mut rng,
                SimulationPhase::Thermalization,
            );
            black_box(kernel.estimator().log_f());
        });
    });
    group.finish();
}

fn bench_canonical_reweighting(criterion: &mut Criterion) {
    let axis = BinnedAxis::new(-2048.0, 2048.0, REWEIGHT_BINS).unwrap();
    let values = (0..axis.bins())
        .map(|bin| 0.0001 * axis.center(bin).powi(2))
        .collect();
    let density = LogDensityOfStates::from_values(values, vec![true; REWEIGHT_BINS]).unwrap();

    let mut group = criterion.benchmark_group("generalized_reweight");
    group.throughput(Throughput::Elements(REWEIGHT_BINS as u64));
    group.bench_function("canonical_log_sum_exp", |bencher| {
        bencher.iter(|| {
            black_box(canonical_reweight(&axis, &density, black_box(0.7)).unwrap());
        });
    });
    group.finish();
}

criterion_group!(
    benches,
    bench_wang_landau_sweep,
    bench_canonical_reweighting
);
criterion_main!(benches);
