use cmc_rs::{
    compute_particle_energy, metropolis_hastings_step, GrandCanonical, InsertDeleteParticle,
    IsothermalIsobaric, IsotropicVolumeChange, LennardJones, MetropolisHastingsAcceptance,
    MoleculeTopology, PairPotential, ParticleBatchMove, ParticleBatchPatch, ParticleConfiguration,
    ParticleDeletion, ParticleInsertion, ParticleSystem, ProposedMove, RigidMoleculeRotation,
    SimulationCell, TorsionDefinition, TorsionRotation, TrialEvaluator, VolumeChangePatch,
};
use rand::SeedableRng;

#[derive(Debug, Clone, Copy)]
struct ZeroPotential {
    cutoff_squared: f64,
}

impl PairPotential for ZeroPotential {
    fn cutoff_squared(&self) -> f64 {
        self.cutoff_squared
    }

    fn energy(&self, _species_i: u16, _species_j: u16, _distance_squared: f64) -> f64 {
        0.0
    }
}

fn lj_system() -> (ParticleSystem<2>, LennardJones) {
    let cell = cmc_rs::OrthorhombicCell::new([8.0, 8.0]).unwrap();
    let configuration = ParticleConfiguration::new(
        vec![[1.0, 1.0], [2.1, 1.2], [5.0, 5.0], [6.0, 5.4]],
        vec![0; 4],
        cell,
    )
    .unwrap();
    let potential = LennardJones::new(1.0, 1.0, 2.5).unwrap();
    let system = ParticleSystem::new(configuration, &potential, 1.0).unwrap();
    (system, potential)
}

#[test]
fn batch_trial_is_transactional_and_commit_matches_exact_energy() {
    let (mut system, potential) = lj_system();
    let old_positions = system.configuration().positions().to_vec();
    let old_energy = system.energy;
    let old_cells: Vec<_> = (0..system.len())
        .map(|particle| system.cell_list().particle_cell(particle))
        .collect();

    let movement = ParticleBatchMove::new(vec![0, 1], vec![[7.7, 1.1], [0.8, 1.3]]).unwrap();
    let mut patch = ParticleBatchPatch::default();
    let delta = system.evaluate_trial(&potential, &movement, &mut patch);

    assert_eq!(system.configuration().positions(), old_positions);
    assert_eq!(system.energy, old_energy);
    for (particle, old_cell) in old_cells.into_iter().enumerate() {
        assert_eq!(system.cell_list().particle_cell(particle), old_cell);
    }

    let mut expected_positions = old_positions;
    expected_positions[0] = [7.7, 1.1];
    expected_positions[1] = [0.8, 1.3];
    let expected_configuration = ParticleConfiguration::new(
        expected_positions,
        vec![0; 4],
        *system.configuration().cell(),
    )
    .unwrap();
    let exact_delta = compute_particle_energy(&expected_configuration, &potential) - old_energy;
    assert!((delta.energy - exact_delta).abs() < 1e-12);

    <ParticleSystem<2> as TrialEvaluator<LennardJones, ParticleBatchMove<2>>>::commit_trial(
        &mut system,
        &movement,
        &patch,
    );
    system.validate(&potential).unwrap();
}

#[test]
fn rigid_rotation_preserves_minimum_image_bond_lengths_across_boundary() {
    let cell = cmc_rs::OrthorhombicCell::new([10.0, 10.0]).unwrap();
    let configuration =
        ParticleConfiguration::new(vec![[9.8, 5.0], [0.2, 5.0], [0.0, 5.3]], vec![0; 3], cell)
            .unwrap();
    let topology = MoleculeTopology::new(3, vec![vec![0, 1, 2]]).unwrap();
    let rotation = RigidMoleculeRotation::new(0.8).unwrap();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(17);
    let proposal = rotation.propose(&configuration, &topology, 0, &mut rng);

    let before = configuration
        .cell()
        .distance_squared(configuration.position(0), configuration.position(1));
    let after = configuration.cell().distance_squared(
        &proposal.movement.positions()[0],
        &proposal.movement.positions()[1],
    );
    assert!((before - after).abs() < 1e-12);
}

#[test]
fn torsion_rotation_preserves_axial_and_radial_coordinates() {
    let cell = cmc_rs::OrthorhombicCell::new([10.0, 10.0, 10.0]).unwrap();
    let configuration = ParticleConfiguration::new(
        vec![[9.5, 5.0, 5.0], [0.5, 5.0, 5.0], [0.0, 6.0, 5.0]],
        vec![0; 3],
        cell,
    )
    .unwrap();
    let definition = TorsionDefinition::new(3, 0, 1, vec![2]).unwrap();
    let torsion = TorsionRotation::new(definition, 1.0).unwrap();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(29);
    let proposal = torsion.propose(&configuration, &mut rng).unwrap();

    let anchor = configuration.position(0);
    let axis = configuration
        .cell()
        .displacement(anchor, configuration.position(1));
    let before = configuration
        .cell()
        .displacement(anchor, configuration.position(2));
    let after = configuration
        .cell()
        .displacement(anchor, &proposal.movement.positions()[0]);
    let norm = axis.iter().map(|value| value * value).sum::<f64>().sqrt();
    let unit = [axis[0] / norm, axis[1] / norm, axis[2] / norm];
    let axial_before = before.iter().zip(unit).map(|(x, u)| x * u).sum::<f64>();
    let axial_after = after.iter().zip(unit).map(|(x, u)| x * u).sum::<f64>();
    let radial_before = before.iter().map(|x| x * x).sum::<f64>() - axial_before.powi(2);
    let radial_after = after.iter().map(|x| x * x).sum::<f64>() - axial_after.powi(2);
    assert!((axial_before - axial_after).abs() < 1e-12);
    assert!((radial_before - radial_after).abs() < 1e-12);
}

#[test]
fn npt_volume_delta_contains_pressure_coordinate_and_proposal_terms() {
    let potential = ZeroPotential {
        cutoff_squared: 1.0,
    };
    let cell = cmc_rs::OrthorhombicCell::new([10.0, 10.0, 10.0]).unwrap();
    let configuration =
        ParticleConfiguration::new(vec![[1.0; 3], [4.0; 3]], vec![0; 2], cell).unwrap();
    let mut system = ParticleSystem::new(configuration, &potential, 0.5).unwrap();
    let movement = IsotropicVolumeChange::new(8.0f64.ln());
    let proposal = ProposedMove::new(movement, movement.log_volume_ratio);
    let mut patch = VolumeChangePatch::default();
    let delta = system.evaluate_trial(&potential, &movement, &mut patch);

    assert!((delta.volume - 7000.0).abs() < 1e-9);
    assert!((delta.log_jacobian - 2.0 * 8.0f64.ln()).abs() < 1e-12);
    let target = IsothermalIsobaric::new(0.5, 0.002);
    let log_acceptance = cmc_rs::AcceptanceRule::log_acceptance(
        &MetropolisHastingsAcceptance,
        &target,
        &delta,
        proposal.log_reverse_over_forward,
    );
    let expected = -0.5 * 0.002 * 7000.0 + 3.0 * 8.0f64.ln();
    assert!((log_acceptance - expected).abs() < 1e-12);

    <ParticleSystem<3> as TrialEvaluator<ZeroPotential, IsotropicVolumeChange>>::commit_trial(
        &mut system,
        &movement,
        &patch,
    );
    assert!((system.configuration().cell().volume() - 8000.0).abs() < 1e-9);
    assert_eq!(system.configuration().position(0), &[2.0; 3]);
    system.validate(&potential).unwrap();
}

#[test]
fn invalid_volume_contraction_is_rejected_without_mutation() {
    let (mut system, potential) = lj_system();
    let positions = system.configuration().positions().to_vec();
    let volume = system.configuration().cell().volume();
    let energy = system.energy;
    let proposal = ProposedMove::new(IsotropicVolumeChange::new(-10.0), -10.0);
    let target = IsothermalIsobaric::new(system.beta, 1.0);
    let mut patch = VolumeChangePatch::default();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(3);
    let outcome = metropolis_hastings_step(
        &mut system,
        &potential,
        &proposal,
        &target,
        &MetropolisHastingsAcceptance,
        &mut patch,
        &mut rng,
    );
    assert!(!outcome.accepted);
    assert_eq!(system.configuration().positions(), positions);
    assert_eq!(system.configuration().cell().volume(), volume);
    assert_eq!(system.energy, energy);
}

#[test]
fn insertion_and_swap_remove_deletion_preserve_all_caches() {
    let potential = ZeroPotential {
        cutoff_squared: 1.0,
    };
    let cell = cmc_rs::OrthorhombicCell::new([6.0, 6.0]).unwrap();
    let configuration =
        ParticleConfiguration::new(vec![[0.5, 0.5], [2.5, 2.5], [5.5, 5.5]], vec![0; 3], cell)
            .unwrap();
    let mut system = ParticleSystem::new(configuration, &potential, 1.0).unwrap();
    let mut patch = cmc_rs::GrandCanonicalPatch::default();

    let insertion = cmc_rs::GrandCanonicalMove::Insert(ParticleInsertion {
        species: 0,
        position: [4.0, 1.0],
    });
    system.evaluate_trial(&potential, &insertion, &mut patch);
    <ParticleSystem<2> as TrialEvaluator<ZeroPotential, cmc_rs::GrandCanonicalMove<2>>>::commit_trial(
        &mut system,
        &insertion,
        &patch,
    );
    assert_eq!(system.len(), 4);
    system.validate(&potential).unwrap();

    let deletion = cmc_rs::GrandCanonicalMove::Delete(ParticleDeletion { particle: 1 });
    system.evaluate_trial(&potential, &deletion, &mut patch);
    <ParticleSystem<2> as TrialEvaluator<ZeroPotential, cmc_rs::GrandCanonicalMove<2>>>::commit_trial(
        &mut system,
        &deletion,
        &patch,
    );
    assert_eq!(system.len(), 3);
    assert_eq!(system.configuration().position(1), &[4.0, 1.0]);
    system.validate(&potential).unwrap();
}

#[test]
fn ideal_gas_grand_canonical_number_mean_is_poisson() {
    let potential = ZeroPotential {
        cutoff_squared: 0.25,
    };
    let cell = cmc_rs::OrthorhombicCell::new([5.0]).unwrap();
    let configuration = ParticleConfiguration::new(Vec::new(), Vec::new(), cell).unwrap();
    let mut system = ParticleSystem::new(configuration, &potential, 1.0).unwrap();
    let proposal_kernel = InsertDeleteParticle::new(0);
    let expected_mean: f64 = 3.0;
    let target = GrandCanonical::new(1.0, (expected_mean / 5.0).ln());
    let acceptance = MetropolisHastingsAcceptance;
    let mut patch = cmc_rs::GrandCanonicalPatch::default();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(1234);

    let mut sum = 0.0;
    let mut samples = 0usize;
    for step in 0..160_000 {
        let proposal = proposal_kernel.propose(&system, &mut rng);
        metropolis_hastings_step(
            &mut system,
            &potential,
            &proposal,
            &target,
            &acceptance,
            &mut patch,
            &mut rng,
        );
        if step >= 20_000 {
            sum += system.len() as f64;
            samples += 1;
        }
    }
    let mean = sum / samples as f64;
    assert!((mean - expected_mean).abs() < 0.08, "mean N was {mean}");
    system.validate(&potential).unwrap();
}

#[test]
fn stage3_kernels_run_and_preserve_invariants() {
    use cmc_rs::{
        MolecularMetropolisCore, ParticleAlgorithm, ParticleGrandCanonicalCore,
        ParticleNptMetropolisCore, SimulationPhase,
    };

    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(77);

    let (mut molecular_system, molecular_potential) = lj_system();
    let topology = MoleculeTopology::new(4, vec![vec![0, 1], vec![2, 3]]).unwrap();
    let mut molecular = MolecularMetropolisCore::<2>::new(topology, 0.1, 0.2)
        .unwrap()
        .with_energy_check_interval(1);
    for _ in 0..20 {
        molecular.sweep_with_phase(
            &mut molecular_system,
            &molecular_potential,
            &mut rng,
            SimulationPhase::Thermalization,
        );
    }
    molecular.sweep_with_phase(
        &mut molecular_system,
        &molecular_potential,
        &mut rng,
        SimulationPhase::Measurement,
    );
    assert!(molecular.move_mixture().is_frozen());
    molecular_system.validate(&molecular_potential).unwrap();

    let zero = ZeroPotential {
        cutoff_squared: 0.25,
    };
    let npt_cell = cmc_rs::OrthorhombicCell::new([5.0, 5.0]).unwrap();
    let npt_configuration =
        ParticleConfiguration::new(vec![[1.0, 1.0], [3.0, 3.0]], vec![0; 2], npt_cell).unwrap();
    let mut npt_system = ParticleSystem::new(npt_configuration, &zero, 1.0).unwrap();
    let mut npt = ParticleNptMetropolisCore::<2>::new(0.1, 0.01, 1.0).with_energy_check_interval(1);
    for _ in 0..30 {
        npt.sweep_with_phase(
            &mut npt_system,
            &zero,
            &mut rng,
            SimulationPhase::Measurement,
        );
    }
    npt_system.validate(&zero).unwrap();

    let gcmc_cell = cmc_rs::OrthorhombicCell::new([5.0]).unwrap();
    let gcmc_configuration = ParticleConfiguration::new(Vec::new(), Vec::new(), gcmc_cell).unwrap();
    let mut gcmc_system = ParticleSystem::new(gcmc_configuration, &zero, 1.0).unwrap();
    let mut gcmc = ParticleGrandCanonicalCore::<1>::new(
        0.1,
        InsertDeleteParticle::new(0),
        (2.0f64 / 5.0).ln(),
    )
    .with_exchange_attempts(2)
    .with_energy_check_interval(1);
    for _ in 0..100 {
        gcmc.sweep_with_phase(
            &mut gcmc_system,
            &zero,
            &mut rng,
            SimulationPhase::Measurement,
        );
    }
    gcmc_system.validate(&zero).unwrap();
}

#[test]
fn periodic_energy_audit_detects_cache_corruption_instead_of_masking_it() {
    use cmc_rs::{ParticleAlgorithm, ParticleMetropolisCore, SimulationPhase};

    let (mut system, potential) = lj_system();
    system.energy += 1.0;
    let mut algorithm = ParticleMetropolisCore::<2>::new(0.01).with_energy_check_interval(1);
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(88);
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        algorithm.sweep_with_phase(
            &mut system,
            &potential,
            &mut rng,
            SimulationPhase::Measurement,
        );
    }));
    assert!(result.is_err());
}
