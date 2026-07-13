use cmc_rs::{
    build_square, statistical_efficiency, IsingGraphWormMC, IsingGraphWormModel, WormConfig,
    WormKernel,
};
use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use std::time::Instant;

const LINEAR_SIZE: usize = 8;
const SITES: usize = LINEAR_SIZE * LINEAR_SIZE;

fn worm_kernel(track_endpoint_pairs: bool) -> WormKernel<IsingGraphWormModel> {
    let model =
        IsingGraphWormModel::new(build_square(LINEAR_SIZE, LINEAR_SIZE, true), 0.44, 1.0).unwrap();
    let configuration = model.empty_configuration();
    let local_updates_per_sweep = model.lattice().n_edges();
    WormKernel::new(
        model,
        configuration,
        WormConfig {
            local_updates_per_sweep,
            close_probability: 0.25,
            log_worm_fugacity: (0.1 / SITES as f64).ln(),
            track_endpoint_pairs,
            cache_audit_interval: 0,
        },
    )
    .unwrap()
}

fn bench_worm_transitions(criterion: &mut Criterion) {
    let mut kernel = worm_kernel(false);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x51a9e);
    let mut group = criterion.benchmark_group("classical_worm");
    group.throughput(Throughput::Elements(1));
    group.bench_function("local_transition", |bencher| {
        bencher.iter(|| black_box(kernel.local_update(&mut rng).unwrap()));
    });
    group.finish();

    let mut kernel = worm_kernel(false);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x51a9f);
    let updates = kernel.config().local_updates_per_sweep;
    let mut group = criterion.benchmark_group("classical_worm_sweep");
    group.throughput(Throughput::Elements(updates as u64));
    group.bench_function("local_updates", |bencher| {
        bencher.iter(|| {
            kernel.sweep(&mut rng).unwrap();
            black_box(kernel.state().configuration().occupied_edges());
        });
    });
    group.finish();
}

fn bench_endpoint_observation(criterion: &mut Criterion) {
    let mut kernel = worm_kernel(true);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x51aa0);
    let mut group = criterion.benchmark_group("classical_worm_endpoint_pairs");
    group.throughput(Throughput::Elements(1));
    group.bench_function("tracked_local_transition", |bencher| {
        bencher.iter(|| black_box(kernel.local_update(&mut rng).unwrap()));
    });
    group.finish();
}

fn bench_checkpoint_serialization(criterion: &mut Criterion) {
    let model =
        IsingGraphWormModel::new(build_square(LINEAR_SIZE, LINEAR_SIZE, true), 0.44, 1.0).unwrap();
    let config = WormConfig {
        local_updates_per_sweep: model.lattice().n_edges(),
        close_probability: 0.25,
        log_worm_fugacity: (0.1 / SITES as f64).ln(),
        track_endpoint_pairs: true,
        cache_audit_interval: 0,
    };
    let mut mc = IsingGraphWormMC::new(model, config).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x51aa1);
    for _ in 0..2_048 {
        mc.kernel_mut().local_update(&mut rng).unwrap();
    }
    let mut group = criterion.benchmark_group("classical_worm_checkpoint");
    group.bench_function("json_serialization", |bencher| {
        bencher.iter(|| black_box(serde_json::to_vec(&mc.save_snapshot()).unwrap()));
    });
    group.finish();
}

fn report_statistical_efficiency() {
    let mut kernel = worm_kernel(false);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x51aa2);
    for _ in 0..2_000 {
        kernel.sweep(&mut rng).unwrap();
    }
    let start = Instant::now();
    let mut samples = Vec::with_capacity(8_192);
    while samples.len() < samples.capacity() {
        kernel.sweep(&mut rng).unwrap();
        if kernel.state().is_physical() {
            samples.push(kernel.state().configuration().occupied_edges() as f64);
        }
    }
    let efficiency = statistical_efficiency(&samples, start.elapsed().as_secs_f64());
    println!(
        "STAT_EFF classical_worm tau_int={:.6} ess={:.3} ess_per_second={:.3}",
        efficiency.integrated_autocorrelation_time,
        efficiency.effective_samples,
        efficiency.effective_samples_per_second,
    );
}

fn worm_benchmarks(criterion: &mut Criterion) {
    report_statistical_efficiency();
    bench_worm_transitions(criterion);
    bench_endpoint_observation(criterion);
    bench_checkpoint_serialization(criterion);
}

criterion_group!(benches, worm_benchmarks);
criterion_main!(benches);
