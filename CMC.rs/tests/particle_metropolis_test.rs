use cmc_rs::{
    CutoffTreatment, LennardJones, OrthorhombicCell, PairPotential, ParticleAlgorithm,
    ParticleConfiguration, ParticleMetropolisCore, ParticleSystem, SimulationPhase,
    TranslateParticle,
};
use rand::SeedableRng;

fn grid_system() -> (ParticleSystem<2>, LennardJones) {
    let cell = OrthorhombicCell::new([8.0, 8.0]).unwrap();
    let mut positions = Vec::new();
    for x in 0..4 {
        for y in 0..4 {
            positions.push([1.0 + 2.0 * x as f64, 1.0 + 2.0 * y as f64]);
        }
    }
    let configuration = ParticleConfiguration::new(positions, vec![0; 16], cell).unwrap();
    let potential =
        LennardJones::with_treatment(1.0, 1.0, 2.5, CutoffTreatment::ShiftedForce).unwrap();
    let system = ParticleSystem::new(configuration, &potential, 1.0).unwrap();
    (system, potential)
}

#[test]
fn nvt_translation_kernel_preserves_all_caches() {
    let (mut system, potential) = grid_system();
    let translation = TranslateParticle::new(0.25).with_adaptation(0.5, 5, 0.4, 1e-4, 2.0);
    let mut kernel = ParticleMetropolisCore::new(0.25)
        .with_translation(translation)
        .with_energy_check_interval(7);
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(1234);

    for _ in 0..30 {
        kernel.sweep_with_phase(
            &mut system,
            &potential,
            &mut rng,
            SimulationPhase::Thermalization,
        );
    }
    let frozen = kernel.translation().max_displacement();
    for _ in 0..40 {
        kernel.sweep_with_phase(
            &mut system,
            &potential,
            &mut rng,
            SimulationPhase::Measurement,
        );
        assert_eq!(kernel.translation().max_displacement(), frozen);
    }

    assert_eq!(kernel.translation().attempted(), 40 * system.len() as u64);
    assert!((0.0..=1.0).contains(&kernel.translation().acceptance_rate()));
    assert!(system.energy_error(&potential).abs() < 1e-10);
    system.validate(&potential).unwrap();
}

#[test]
fn fixed_seed_reproduces_particle_trajectory() {
    let (mut left, potential_left) = grid_system();
    let (mut right, potential_right) = grid_system();
    let mut kernel_left = ParticleMetropolisCore::new(0.2);
    let mut kernel_right = ParticleMetropolisCore::new(0.2);
    let mut rng_left = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(991);
    let mut rng_right = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(991);

    for _ in 0..25 {
        kernel_left.sweep(&mut left, &potential_left, &mut rng_left);
        kernel_right.sweep(&mut right, &potential_right, &mut rng_right);
    }
    assert_eq!(left.configuration(), right.configuration());
    assert_eq!(left.energy.to_bits(), right.energy.to_bits());
    assert_eq!(
        kernel_left.translation().accepted(),
        kernel_right.translation().accepted()
    );
}

#[test]
fn two_particle_energy_distribution_matches_quadrature() {
    const LENGTH: f64 = 5.0;
    const BETA: f64 = 1.0;
    const GRID: usize = 320;

    let cell = OrthorhombicCell::new([LENGTH, LENGTH]).unwrap();
    let configuration =
        ParticleConfiguration::new(vec![[1.25, 2.5], [3.75, 2.5]], vec![0, 0], cell).unwrap();
    let potential =
        LennardJones::with_treatment(1.0, 1.0, 2.0, CutoffTreatment::ShiftedPotential).unwrap();
    let mut system = ParticleSystem::new(configuration, &potential, BETA).unwrap();
    let mut kernel = ParticleMetropolisCore::new(0.8);
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(0x005E_ED2D);

    for _ in 0..20_000 {
        kernel.sweep(&mut system, &potential, &mut rng);
    }
    let mut sampled_sum = 0.0;
    let samples = 120_000usize;
    for _ in 0..samples {
        kernel.sweep(&mut system, &potential, &mut rng);
        sampled_sum += system.energy;
    }
    let sampled_mean = sampled_sum / samples as f64;

    // Translational invariance reduces the exact two-particle integral to the
    // relative displacement over one periodic square. Midpoint quadrature is
    // deterministic and resolves the repulsive core without evaluating r=0.
    let spacing = LENGTH / GRID as f64;
    let mut weighted_energy = 0.0;
    let mut partition = 0.0;
    for ix in 0..GRID {
        let x = (ix as f64 + 0.5) * spacing - 0.5 * LENGTH;
        for idx_y in 0..GRID {
            let y = (idx_y as f64 + 0.5) * spacing - 0.5 * LENGTH;
            let energy = potential.energy(0, 0, x * x + y * y);
            let weight = (-BETA * energy).exp();
            weighted_energy += energy * weight;
            partition += weight;
        }
    }
    let reference_mean = weighted_energy / partition;
    assert!(
        (sampled_mean - reference_mean).abs() < 0.02,
        "sampled mean {sampled_mean} differs from quadrature {reference_mean}"
    );
}
