use approx::assert_abs_diff_eq;
use cmc_rs::{
    compute_particle_energy, metropolis_hastings_step, CanonicalEnsemble, CellList,
    CutoffTreatment, LennardJones, MetropolisHastingsAcceptance, OrthorhombicCell, PairPotential,
    ParticleConfiguration, ParticleEnergyPatch, ParticleSystem, ParticleTranslation, ProposedMove,
    SimulationCell, TrialEvaluator,
};
use rand::SeedableRng;

fn simple_system() -> (ParticleSystem<2>, LennardJones) {
    let cell = OrthorhombicCell::new([8.0, 8.0]).unwrap();
    let configuration = ParticleConfiguration::new(
        vec![[0.5, 0.5], [2.0, 0.5], [4.0, 4.0], [7.5, 7.5]],
        vec![0; 4],
        cell,
    )
    .unwrap();
    let potential =
        LennardJones::with_treatment(1.0, 1.0, 3.0, CutoffTreatment::ShiftedPotential).unwrap();
    let system = ParticleSystem::new(configuration, &potential, 1.0).unwrap();
    (system, potential)
}

#[test]
fn two_particle_energy_matches_analytic_pair_value() {
    let cell = OrthorhombicCell::new([10.0, 10.0, 10.0]).unwrap();
    let configuration =
        ParticleConfiguration::new(vec![[1.0, 1.0, 1.0], [2.5, 1.0, 1.0]], vec![0, 0], cell)
            .unwrap();
    let potential =
        LennardJones::with_treatment(1.0, 2.0, 4.0, CutoffTreatment::Truncated).unwrap();
    let expected = 8.0 * ((1.0 / 1.5f64).powi(12) - (1.0 / 1.5f64).powi(6));
    assert_abs_diff_eq!(
        compute_particle_energy(&configuration, &potential),
        expected,
        epsilon = 1e-14
    );
}

#[test]
fn minimum_image_interaction_crosses_periodic_boundary() {
    let cell = OrthorhombicCell::new([10.0, 10.0]).unwrap();
    let configuration =
        ParticleConfiguration::new(vec![[0.2, 5.0], [9.0, 5.0]], vec![0, 0], cell).unwrap();
    let potential =
        LennardJones::with_treatment(1.0, 1.0, 2.0, CutoffTreatment::Truncated).unwrap();
    let distance_squared = cell.distance_squared(&[0.2, 5.0], &[9.0, 5.0]);
    assert_abs_diff_eq!(distance_squared, 1.44, epsilon = 1e-14);
    assert_abs_diff_eq!(
        compute_particle_energy(&configuration, &potential),
        potential.energy(0, 0, 1.44),
        epsilon = 1e-8
    );
}

#[test]
fn rejected_trial_keeps_position_energy_and_cell_membership_identical() {
    let (mut system, potential) = simple_system();
    let before_position = *system.configuration().position(0);
    let before_energy = system.energy;
    let before_cell = system.cell_list().particle_cell(0);
    let proposal = ProposedMove::symmetric(ParticleTranslation::new(
        0,
        *system.configuration().position(1),
    ));
    let mut patch = ParticleEnergyPatch::default();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(9);
    let outcome = metropolis_hastings_step(
        &mut system,
        &potential,
        &proposal,
        &CanonicalEnsemble::new(1.0),
        &MetropolisHastingsAcceptance,
        &mut patch,
        &mut rng,
    );
    assert!(!outcome.accepted);
    assert_eq!(*system.configuration().position(0), before_position);
    assert_eq!(system.energy, before_energy);
    assert_eq!(system.cell_list().particle_cell(0), before_cell);
    system.cell_list().validate(system.configuration()).unwrap();
}

#[test]
fn crossing_boundary_commit_patches_position_cell_list_and_energy() {
    let (mut system, potential) = simple_system();
    let movement = ParticleTranslation::new(0, [7.9, 0.5]);
    let mut patch = ParticleEnergyPatch::default();
    let before_energy = system.energy;
    let delta = system.evaluate_trial(&potential, &movement, &mut patch);
    assert_eq!(*system.configuration().position(0), [0.5, 0.5]);
    <ParticleSystem<2> as TrialEvaluator<LennardJones, ParticleTranslation<2>>>::commit_trial(
        &mut system,
        &movement,
        &patch,
    );
    assert_abs_diff_eq!(system.energy_error(&potential), 0.0, epsilon = 1e-12);
    assert_abs_diff_eq!(system.energy, before_energy + delta.energy, epsilon = 1e-14);
    system.cell_list().validate(system.configuration()).unwrap();
}

#[test]
fn cell_list_candidate_energy_matches_brute_force_after_many_commits() {
    let cell = OrthorhombicCell::new([12.0, 12.0, 12.0]).unwrap();
    let mut positions = Vec::new();
    for x in 0..3 {
        for y in 0..3 {
            for z in 0..3 {
                positions.push([
                    1.0 + 3.5 * x as f64,
                    1.0 + 3.5 * y as f64,
                    1.0 + 3.5 * z as f64,
                ]);
            }
        }
    }
    let configuration = ParticleConfiguration::new(positions, vec![0; 27], cell).unwrap();
    let potential = LennardJones::new(1.0, 1.0, 2.5).unwrap();
    let mut system = ParticleSystem::new(configuration, &potential, 0.0).unwrap();
    let mut patch = ParticleEnergyPatch::default();

    for step in 0..80 {
        let particle = step % system.len();
        let old = *system.configuration().position(particle);
        let movement = ParticleTranslation::new(
            particle,
            [
                (old[0] + 0.37).rem_euclid(12.0),
                (old[1] + 0.19).rem_euclid(12.0),
                (old[2] + 0.11).rem_euclid(12.0),
            ],
        );
        let delta = system.evaluate_trial(&potential, &movement, &mut patch);
        if delta.energy.is_finite() {
            <ParticleSystem<3> as TrialEvaluator<LennardJones, ParticleTranslation<3>>>::commit_trial(
                &mut system,
                &movement,
                &patch,
            );
        }
        assert_abs_diff_eq!(system.energy_error(&potential), 0.0, epsilon = 2e-11);
        system.cell_list().validate(system.configuration()).unwrap();
    }
}

#[test]
fn cell_list_neighbor_query_contains_every_pair_inside_cutoff() {
    let cell = OrthorhombicCell::new([9.0, 9.0]).unwrap();
    let configuration = ParticleConfiguration::new(
        vec![[0.1, 0.1], [8.8, 0.2], [4.2, 4.3], [5.9, 4.3], [2.1, 7.8]],
        vec![0; 5],
        cell,
    )
    .unwrap();
    let cutoff = 2.0;
    let list = CellList::new(&configuration, cutoff * cutoff).unwrap();
    let mut candidates = Vec::new();
    for particle in 0..configuration.len() {
        let cell_index = list.particle_cell(particle);
        list.fill_candidates(cell_index, &mut candidates);
        for other in 0..configuration.len() {
            if other == particle {
                continue;
            }
            let r2 = cell.distance_squared(
                configuration.position(particle),
                configuration.position(other),
            );
            if r2 < cutoff * cutoff {
                assert!(candidates.contains(&other));
            }
        }
    }
}

#[test]
fn two_and_three_dimensional_cells_construct() {
    let cell_2d = OrthorhombicCell::new([5.0, 6.0]).unwrap();
    let cell_3d = OrthorhombicCell::new([5.0, 6.0, 7.0]).unwrap();
    assert_eq!(cell_2d.volume(), 30.0);
    assert_eq!(cell_3d.volume(), 210.0);
}
