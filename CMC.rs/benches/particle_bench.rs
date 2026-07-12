use cmc_rs::{
    compute_particle_energy, CutoffTreatment, LennardJones, OrthorhombicCell, ParticleAlgorithm,
    ParticleConfiguration, ParticleMetropolisCore, ParticleSystem,
};
use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use rand::SeedableRng;

const PARTICLES_PER_AXIS: usize = 6;
const N_PARTICLES: usize = PARTICLES_PER_AXIS * PARTICLES_PER_AXIS * PARTICLES_PER_AXIS;

fn dense_lj_system() -> (ParticleSystem<3>, LennardJones) {
    let density = 0.7;
    let length = (N_PARTICLES as f64 / density).cbrt();
    let cell = OrthorhombicCell::new([length; 3]).unwrap();
    let spacing = length / PARTICLES_PER_AXIS as f64;
    let mut positions = Vec::with_capacity(N_PARTICLES);
    for x in 0..PARTICLES_PER_AXIS {
        for y in 0..PARTICLES_PER_AXIS {
            for z in 0..PARTICLES_PER_AXIS {
                positions.push([
                    (x as f64 + 0.5) * spacing,
                    (y as f64 + 0.5) * spacing,
                    (z as f64 + 0.5) * spacing,
                ]);
            }
        }
    }
    let configuration = ParticleConfiguration::new(positions, vec![0; N_PARTICLES], cell).unwrap();
    let potential =
        LennardJones::with_treatment(1.0, 1.0, 2.5, CutoffTreatment::ShiftedPotential).unwrap();
    let system = ParticleSystem::new(configuration, &potential, 1.0).unwrap();
    (system, potential)
}

fn bench_lj_translation_sweep(criterion: &mut Criterion) {
    let (mut system, potential) = dense_lj_system();
    let mut kernel = ParticleMetropolisCore::new(0.1);
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(0xC311_1157);
    let mut group = criterion.benchmark_group("lj_nvt");
    group.throughput(Throughput::Elements(N_PARTICLES as u64));
    group.bench_function("cell_list_translation_attempts", |bencher| {
        bencher.iter(|| {
            kernel.sweep(&mut system, &potential, &mut rng);
            black_box(system.energy);
        });
    });
    group.finish();
}

fn bench_lj_full_energy(criterion: &mut Criterion) {
    let (system, potential) = dense_lj_system();
    let mut group = criterion.benchmark_group("lj_reference");
    group.throughput(Throughput::Elements(N_PARTICLES as u64));
    group.bench_function("quadratic_full_energy", |bencher| {
        bencher.iter(|| {
            black_box(compute_particle_energy(system.configuration(), &potential));
        });
    });
    group.finish();
}

criterion_group!(benches, bench_lj_translation_sweep, bench_lj_full_energy);
criterion_main!(benches);
