use cmc_rs::{
    build_square, enumerate_ising_density_of_states, statistical_efficiency, Algorithm,
    BatchEnergyPatch, BatchSpinMove, CellList, CutoffTreatment, IsingModel, LennardJones,
    MetropolisCore, OrthorhombicCell, ParticleAlgorithm, ParticleConfiguration,
    ParticleMetropolisCore, ParticleSystem, SWCore, SimulationPhase, System, TrialEvaluator,
    WangLandauConfig, WangLandauCore, WolffCore,
};
use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use std::sync::Once;
use std::time::Instant;

const ISING_LENGTH: usize = 24;
const ISING_SITES: usize = ISING_LENGTH * ISING_LENGTH;
const LJ_AXIS: usize = 5;
const LJ_PARTICLES: usize = LJ_AXIS * LJ_AXIS * LJ_AXIS;

fn ising_system() -> (System, IsingModel) {
    let lattice = build_square(ISING_LENGTH, ISING_LENGTH, true);
    let model = IsingModel::new(1.0);
    let mut system = System::new(lattice, 1, 1.0, 0.44);
    system.recompute_energy(&model);
    (system, model)
}

fn lj_system() -> (ParticleSystem<3>, LennardJones) {
    let density = 0.55;
    let length = (LJ_PARTICLES as f64 / density).cbrt();
    let cell = OrthorhombicCell::new([length; 3]).unwrap();
    let spacing = length / LJ_AXIS as f64;
    let mut positions = Vec::with_capacity(LJ_PARTICLES);
    for x in 0..LJ_AXIS {
        for y in 0..LJ_AXIS {
            for z in 0..LJ_AXIS {
                positions.push([
                    (x as f64 + 0.5) * spacing,
                    (y as f64 + 0.5) * spacing,
                    (z as f64 + 0.5) * spacing,
                ]);
            }
        }
    }
    let configuration = ParticleConfiguration::new(positions, vec![0; LJ_PARTICLES], cell).unwrap();
    let potential =
        LennardJones::with_treatment(1.0, 1.0, 2.5, CutoffTreatment::ShiftedPotential).unwrap();
    let system = ParticleSystem::new(configuration, &potential, 1.0).unwrap();
    (system, potential)
}

fn print_efficiency(name: &str, samples: &[f64], elapsed: f64) {
    let result = statistical_efficiency(samples, elapsed);
    eprintln!(
        "STAT_EFF {name} tau_int={:.6} ess={:.3} ess_per_second={:.3}",
        result.integrated_autocorrelation_time,
        result.effective_samples,
        result.effective_samples_per_second
    );
}

fn report_statistical_efficiency_once() {
    static ONCE: Once = Once::new();
    ONCE.call_once(|| {
        let (mut system, model) = ising_system();
        let mut kernel = MetropolisCore::new();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x4D45_5452_4F50);
        let mut samples = Vec::with_capacity(1024);
        let started = Instant::now();
        for _ in 0..1024 {
            kernel.sweep(&mut system, &model, &mut rng);
            samples.push(system.energy);
        }
        print_efficiency(
            "ising_metropolis",
            &samples,
            started.elapsed().as_secs_f64(),
        );

        let (mut system, model) = ising_system();
        let mut kernel = WolffCore::new();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x574F_4C46);
        samples.clear();
        let started = Instant::now();
        for _ in 0..1024 {
            kernel.sweep(&mut system, &model, &mut rng);
            samples.push(system.energy);
        }
        print_efficiency("ising_wolff", &samples, started.elapsed().as_secs_f64());

        let (mut system, model) = ising_system();
        let mut kernel = SWCore::new();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x5357_4544);
        samples.clear();
        let started = Instant::now();
        for _ in 0..512 {
            kernel.sweep(&mut system, &model, &mut rng);
            samples.push(system.energy);
        }
        print_efficiency(
            "ising_swendsen_wang",
            &samples,
            started.elapsed().as_secs_f64(),
        );

        let (mut system, potential) = lj_system();
        let mut kernel = ParticleMetropolisCore::new(0.1);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x004C_4A4E_5654);
        samples.clear();
        let started = Instant::now();
        for _ in 0..128 {
            kernel.sweep(&mut system, &potential, &mut rng);
            samples.push(system.energy);
        }
        print_efficiency("lj_translation", &samples, started.elapsed().as_secs_f64());
    });
}

fn bench_lattice_updates(criterion: &mut Criterion) {
    report_statistical_efficiency_once();

    let (mut system, model) = ising_system();
    let mut metropolis = MetropolisCore::new();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(1);
    let mut group = criterion.benchmark_group("lattice_updates");
    group.throughput(Throughput::Elements(ISING_SITES as u64));
    group.bench_function("ising_metropolis_attempted_updates", |bencher| {
        bencher.iter(|| {
            metropolis.sweep(&mut system, &model, &mut rng);
            black_box(system.energy);
        });
    });
    group.finish();

    let (mut system, model) = ising_system();
    let mut wolff = WolffCore::new();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(2);
    let mut cluster_total = 0usize;
    for _ in 0..128 {
        wolff.sweep(&mut system, &model, &mut rng);
        cluster_total += wolff.last_cluster_size();
    }
    let mean_cluster = (cluster_total / 128).max(1);
    let mut group = criterion.benchmark_group("wolff_cluster");
    group.throughput(Throughput::Elements(mean_cluster as u64));
    group.bench_function("cluster_sites", |bencher| {
        bencher.iter(|| {
            wolff.sweep(&mut system, &model, &mut rng);
            black_box(wolff.last_cluster_size());
        });
    });
    group.finish();

    let (mut system, model) = ising_system();
    let edges = system.lattice.n_edges();
    let mut sw = SWCore::new();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(3);
    let mut group = criterion.benchmark_group("swendsen_wang");
    group.throughput(Throughput::Elements(edges as u64));
    group.bench_function("physical_edges", |bencher| {
        bencher.iter(|| {
            sw.sweep(&mut system, &model, &mut rng);
            black_box(system.energy);
        });
    });
    group.finish();
}

fn bench_particle_paths(criterion: &mut Criterion) {
    let (mut system, potential) = lj_system();
    let mut kernel = ParticleMetropolisCore::new(0.1);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(4);
    let mut group = criterion.benchmark_group("particle_updates");
    group.throughput(Throughput::Elements(LJ_PARTICLES as u64));
    group.bench_function("lj_trial_translations", |bencher| {
        bencher.iter(|| {
            kernel.sweep(&mut system, &potential, &mut rng);
            black_box(system.energy);
        });
    });
    group.finish();

    let cell_list: &CellList<3> = system.cell_list();
    let cell_index = cell_list.particle_cell(0);
    let mut candidates = Vec::new();
    let mut group = criterion.benchmark_group("cell_list");
    group.throughput(Throughput::Elements(1));
    group.bench_function("neighbor_query", |bencher| {
        bencher.iter(|| {
            cell_list.fill_candidates(cell_index, &mut candidates);
            black_box(candidates.len());
        });
    });
    group.finish();
}

fn bench_batch_delta_energy(criterion: &mut Criterion) {
    let (system, model) = ising_system();
    let mut movement = BatchSpinMove::with_capacity(1, 32);
    for site in (0..ISING_SITES).step_by(17).take(32) {
        movement.push(site, &[-system.spins[site]]);
    }
    let mut patch = BatchEnergyPatch::default();
    let mut group = criterion.benchmark_group("batch_move");
    group.throughput(Throughput::Elements(movement.len() as u64));
    group.bench_function("delta_energy", |bencher| {
        bencher.iter(|| {
            black_box(system.evaluate_trial(&model, &movement, &mut patch));
        });
    });
    group.finish();
}

fn bench_checkpoint_serialization(criterion: &mut Criterion) {
    let lattice = cmc_rs::build_chain(12, true);
    let model = IsingModel::new(1.0);
    let axis = enumerate_ising_density_of_states(&lattice, &model)
        .unwrap()
        .axis()
        .unwrap();
    let mut kernel = WangLandauCore::new(axis, WangLandauConfig::default()).unwrap();
    let mut system = System::new(lattice, 1, 1.0, 0.0);
    system.recompute_energy(&model);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(5);
    for _ in 0..64 {
        kernel.sweep_with_phase(
            &mut system,
            &model,
            &mut rng,
            SimulationPhase::Thermalization,
        );
    }
    let mut group = criterion.benchmark_group("checkpoint");
    group.bench_function("wang_landau_json_serialization", |bencher| {
        bencher.iter(|| black_box(serde_json::to_vec(&kernel.save_snapshot()).unwrap()));
    });
    group.finish();
}

criterion_group!(
    benches,
    bench_lattice_updates,
    bench_particle_paths,
    bench_batch_delta_energy,
    bench_checkpoint_serialization
);
criterion_main!(benches);
