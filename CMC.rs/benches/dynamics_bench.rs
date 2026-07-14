use cmc_rs::{
    build_square, statistical_efficiency, Algorithm, BklIsingKernel, GillespieKernel,
    HardSphereEventChain, IsingModel, KawasakiCore, KineticIsingModel, KineticRateLaw,
    OrthorhombicCell, ParticleConfiguration, System,
};
use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use std::sync::Once;
use std::time::Instant;

const ISING_LENGTH: usize = 16;
const ISING_SITES: usize = ISING_LENGTH * ISING_LENGTH;
const HARD_SPHERE_AXIS: usize = 8;
const HARD_SPHERES: usize = HARD_SPHERE_AXIS * HARD_SPHERE_AXIS;
const EVENT_CHAIN_LENGTH: f64 = 8.0;

fn ising_state(beta: f64) -> System {
    let lattice = build_square(ISING_LENGTH, ISING_LENGTH, true);
    let mut state = System::new(lattice, 1, 1.0, beta);
    for (site, spin) in state.spins.iter_mut().enumerate() {
        let x = site % ISING_LENGTH;
        let y = site / ISING_LENGTH;
        *spin = if (x + y).is_multiple_of(2) { 1.0 } else { -1.0 };
    }
    state.recompute_energy(&IsingModel::new(1.0));
    state
}

fn kinetic_model() -> KineticIsingModel {
    KineticIsingModel::new(1.0, KineticRateLaw::glauber(1.0).unwrap()).unwrap()
}

fn bkl_kernel() -> BklIsingKernel {
    BklIsingKernel::new(kinetic_model(), ising_state(0.44), 0).unwrap()
}

fn hard_sphere_kernel() -> HardSphereEventChain<2> {
    let box_length = 16.0;
    let spacing = box_length / HARD_SPHERE_AXIS as f64;
    let mut positions = Vec::with_capacity(HARD_SPHERES);
    for y in 0..HARD_SPHERE_AXIS {
        for x in 0..HARD_SPHERE_AXIS {
            positions.push([(x as f64 + 0.5) * spacing, (y as f64 + 0.5) * spacing]);
        }
    }
    let cell = OrthorhombicCell::new([box_length; 2]).unwrap();
    let configuration = ParticleConfiguration::new(positions, vec![0; HARD_SPHERES], cell).unwrap();
    HardSphereEventChain::new(configuration, 1.0, EVENT_CHAIN_LENGTH, 0).unwrap()
}

fn report_statistical_efficiency_once() {
    static ONCE: Once = Once::new();
    ONCE.call_once(|| {
        let mut kernel = bkl_kernel();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x424B_4C45_5353);
        let mut samples = Vec::with_capacity(4096);
        let started = Instant::now();
        for _ in 0..4096 {
            kernel.advance_by(0.25, &mut rng).unwrap();
            samples.push(kernel.state().energy);
        }
        let efficiency = statistical_efficiency(&samples, started.elapsed().as_secs_f64());
        eprintln!(
            "STAT_EFF bkl_ising tau_int={:.6} ess={:.3} ess_per_second={:.3}",
            efficiency.integrated_autocorrelation_time,
            efficiency.effective_samples,
            efficiency.effective_samples_per_second
        );
    });
}

fn bench_kawasaki(criterion: &mut Criterion) {
    report_statistical_efficiency_once();
    let model = IsingModel::new(1.0);
    let mut state = ising_state(0.44);
    let mut kernel = KawasakiCore::new(ISING_SITES);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(1);
    let mut group = criterion.benchmark_group("classical_dynamics_kawasaki");
    group.throughput(Throughput::Elements(ISING_SITES as u64));
    group.bench_function("attempted_exchanges", |bencher| {
        bencher.iter(|| {
            kernel.sweep(&mut state, &model, &mut rng);
            black_box((state.energy, kernel.last_accepts()));
        });
    });
    group.finish();
}

fn bench_rejection_free(criterion: &mut Criterion) {
    let model = kinetic_model();
    let state = ising_state(0.44);
    let mut direct = GillespieKernel::new(model, state).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(2);
    let mut group = criterion.benchmark_group("classical_dynamics_rejection_free");
    group.throughput(Throughput::Elements(1));
    group.bench_function("direct_gillespie_event", |bencher| {
        bencher.iter(|| black_box(direct.step(&mut rng).unwrap()));
    });

    let mut bkl = bkl_kernel();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(3);
    group.bench_function("bkl_fenwick_event", |bencher| {
        bencher.iter(|| black_box(bkl.step(&mut rng).unwrap()));
    });
    group.finish();
}

fn bench_event_chain(criterion: &mut Criterion) {
    let mut kernel = hard_sphere_kernel();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(4);
    let mut group = criterion.benchmark_group("hard_sphere_event_chain");
    group.throughput(Throughput::Elements(EVENT_CHAIN_LENGTH as u64));
    group.bench_function("lifted_distance", |bencher| {
        bencher.iter(|| black_box(kernel.step(&mut rng).unwrap()));
    });
    group.finish();
}

fn bench_dynamics_checkpoint(criterion: &mut Criterion) {
    let mut kernel = bkl_kernel();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(5);
    for _ in 0..64 {
        kernel.step(&mut rng).unwrap();
    }
    let mut group = criterion.benchmark_group("classical_dynamics_checkpoint");
    group.bench_function("bkl_json_serialization", |bencher| {
        bencher.iter(|| black_box(serde_json::to_vec(&kernel.save_snapshot()).unwrap()));
    });
    group.finish();
}

criterion_group!(
    benches,
    bench_kawasaki,
    bench_rejection_free,
    bench_event_chain,
    bench_dynamics_checkpoint
);
criterion_main!(benches);
