use super::common::assert_close;
use cmc_rs::{
    compute_particle_energy, AcceptanceRule, CutoffTreatment, IsothermalIsobaric,
    IsotropicVolumeChange, LennardJones, MetropolisHastingsAcceptance, OrthorhombicCell,
    PairPotential, ParticleBatchMove, ParticleBatchPatch, ParticleConfiguration, ParticleSystem,
    ParticleTranslation, ProposedMove, SimulationCell, TrialEvaluator, VolumeChangePatch,
};

#[derive(Debug, Clone, Copy)]
struct ZeroPotential(f64);
impl PairPotential for ZeroPotential {
    fn cutoff_squared(&self) -> f64 {
        self.0
    }
    fn energy(&self, _: u16, _: u16, _: f64) -> f64 {
        0.0
    }
}

#[test]
fn periodic_minimum_image_and_wrap_are_exact_for_large_offsets() {
    let cell = OrthorhombicCell::new([10.0, 8.0, 6.0]).unwrap();
    let mut position = [31.25, -17.5, 18.25];
    cell.wrap(&mut position);
    assert_eq!(position, [1.25, 6.5, 0.25]);
    let displacement = cell.displacement(&[0.2, 7.7, 0.1], &[9.8, 0.3, 5.9]);
    assert_close(displacement[0], -0.4, 2e-15);
    assert_close(displacement[1], 0.6, 2e-15);
    assert_close(displacement[2], -0.2, 2e-15);
}

#[test]
fn lennard_jones_minimum_and_cutoff_conditions_match_analytic_values() {
    let minimum = 2.0f64.powf(1.0 / 6.0);
    let raw = LennardJones::with_treatment(1.0, 2.5, 3.0, CutoffTreatment::Truncated).unwrap();
    assert_close(raw.energy(0, 0, minimum * minimum), -2.5, 2e-14);

    let shifted =
        LennardJones::with_treatment(1.0, 1.0, 2.5, CutoffTreatment::ShiftedPotential).unwrap();
    assert_eq!(shifted.energy(0, 0, 2.5f64.powi(2)), 0.0);
    assert!(shifted.energy(0, 0, (2.5_f64 - 1e-8).powi(2)).abs() < 1e-8);

    let shifted_force =
        LennardJones::with_treatment(1.0, 1.0, 2.5, CutoffTreatment::ShiftedForce).unwrap();
    let h = 1e-5;
    let derivative = (shifted_force.energy(0, 0, (2.5_f64 - h).powi(2))
        - shifted_force.energy(0, 0, (2.5_f64 - 2.0 * h).powi(2)))
        / h;
    assert!(derivative.abs() < 2e-4, "cutoff derivative={derivative}");
}

#[test]
fn particle_translation_and_batch_deltas_equal_full_n_squared_recomputation() {
    let cell = OrthorhombicCell::new([8.0, 8.0]).unwrap();
    let configuration = ParticleConfiguration::new(
        vec![[0.2, 0.3], [1.3, 0.7], [7.7, 0.4], [4.2, 4.7], [5.0, 4.8]],
        vec![0; 5],
        cell,
    )
    .unwrap();
    let potential = LennardJones::new(1.0, 0.7, 2.5).unwrap();
    let mut system = ParticleSystem::new(configuration, &potential, 0.9).unwrap();

    let movement = ParticleTranslation {
        particle: 2,
        position: [0.8, 0.55],
    };
    let mut patch = cmc_rs::ParticleEnergyPatch::default();
    let delta = system.evaluate_trial(&potential, &movement, &mut patch);
    let mut positions = system.configuration().positions().to_vec();
    positions[2] = movement.position;
    let expected = ParticleConfiguration::new(
        positions,
        system.configuration().species().to_vec(),
        *system.configuration().cell(),
    )
    .unwrap();
    assert_close(
        delta.energy,
        compute_particle_energy(&expected, &potential) - system.energy,
        2e-12,
    );
    let original = system.configuration().positions().to_vec();
    assert_eq!(system.configuration().positions(), original);
    <ParticleSystem<2> as TrialEvaluator<LennardJones, ParticleTranslation<2>>>::commit_trial(
        &mut system,
        &movement,
        &patch,
    );
    system.validate(&potential).unwrap();

    let batch = ParticleBatchMove::new(vec![0, 4], vec![[7.85, 0.25], [4.8, 5.2]]).unwrap();
    let mut batch_patch = ParticleBatchPatch::default();
    let before = system.energy;
    let delta = system.evaluate_trial(&potential, &batch, &mut batch_patch);
    let mut positions = system.configuration().positions().to_vec();
    positions[0] = [7.85, 0.25];
    positions[4] = [4.8, 5.2];
    let expected = ParticleConfiguration::new(
        positions,
        system.configuration().species().to_vec(),
        *system.configuration().cell(),
    )
    .unwrap();
    assert_close(
        delta.energy,
        compute_particle_energy(&expected, &potential) - before,
        2e-12,
    );
}

#[test]
fn npt_log_volume_move_contains_state_jacobian_pressure_and_hastings_terms() {
    let potential = ZeroPotential(1.0);
    let cell = OrthorhombicCell::new([10.0, 10.0, 10.0]).unwrap();
    let configuration =
        ParticleConfiguration::new(vec![[1.0; 3], [4.0; 3]], vec![0; 2], cell).unwrap();
    let system = ParticleSystem::new(configuration, &potential, 0.5).unwrap();
    let movement = IsotropicVolumeChange::new(8.0f64.ln());
    let proposal = ProposedMove::new(movement, movement.log_volume_ratio);
    let mut patch = VolumeChangePatch::default();
    let delta = system.evaluate_trial(&potential, &movement, &mut patch);
    assert_close(delta.volume, 7_000.0, 1e-9);
    assert_close(delta.log_jacobian, 2.0 * 8.0f64.ln(), 1e-13);
    let target = IsothermalIsobaric::new(0.5, 0.002);
    let log_acceptance = MetropolisHastingsAcceptance.log_acceptance(
        &target,
        &delta,
        proposal.log_reverse_over_forward,
    );
    assert_close(
        log_acceptance,
        -0.5 * 0.002 * 7_000.0 + 3.0 * 8.0f64.ln(),
        2e-13,
    );
}

#[test]
fn ideal_gas_insertion_and_deletion_proposal_ratios_are_exact_reciprocals() {
    use cmc_rs::{GrandCanonicalMove, InsertDeleteParticle};
    use rand::SeedableRng;
    let potential = ZeroPotential(0.25);
    let cell = OrthorhombicCell::new([5.0]).unwrap();
    let empty = ParticleConfiguration::new(Vec::new(), Vec::new(), cell).unwrap();
    let system0 = ParticleSystem::new(empty, &potential, 1.0).unwrap();
    let proposal_kernel = InsertDeleteParticle::new(0)
        .with_particle_bounds(0, Some(1))
        .unwrap();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(0x1DEA1);
    let insertion = proposal_kernel.propose(&system0, &mut rng);
    assert!(matches!(insertion.movement, GrandCanonicalMove::Insert(_)));
    assert_close(insertion.log_reverse_over_forward, 5.0f64.ln(), 1e-15);

    let occupied = ParticleConfiguration::new(vec![[1.0]], vec![0], cell).unwrap();
    let system1 = ParticleSystem::new(occupied, &potential, 1.0).unwrap();
    let deletion = proposal_kernel.propose(&system1, &mut rng);
    assert!(matches!(deletion.movement, GrandCanonicalMove::Delete(_)));
    assert_close(deletion.log_reverse_over_forward, -5.0f64.ln(), 1e-15);
    assert_close(
        insertion.log_reverse_over_forward + deletion.log_reverse_over_forward,
        0.0,
        1e-15,
    );
}

#[test]
fn rigid_molecular_moves_preserve_internal_geometry_in_minimum_image_convention() {
    use cmc_rs::{MoleculeTopology, RigidMoleculeRotation, RigidMoleculeTranslation};
    use rand::SeedableRng;
    let cell = OrthorhombicCell::new([10.0, 10.0]).unwrap();
    let configuration =
        ParticleConfiguration::new(vec![[9.8, 5.0], [0.2, 5.0], [0.0, 5.3]], vec![0; 3], cell)
            .unwrap();
    let topology = MoleculeTopology::new(3, vec![vec![0, 1, 2]]).unwrap();
    let before = |a: usize, b: usize| {
        configuration
            .cell()
            .distance_squared(configuration.position(a), configuration.position(b))
    };
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(0xB00D);
    let translation =
        RigidMoleculeTranslation::new(1.0)
            .unwrap()
            .propose(&configuration, &topology, 0, &mut rng);
    let rotation =
        RigidMoleculeRotation::new(1.0)
            .unwrap()
            .propose(&configuration, &topology, 0, &mut rng);
    for (local_a, &a) in [0usize, 1, 2].iter().enumerate() {
        for (local_b, &b) in [0usize, 1, 2].iter().enumerate().skip(local_a + 1) {
            let _ = local_b;
            assert_close(
                configuration.cell().distance_squared(
                    &translation.movement.positions()[local_a],
                    &translation.movement.positions()[local_b],
                ),
                before(a, b),
                3e-13,
            );
            assert_close(
                configuration.cell().distance_squared(
                    &rotation.movement.positions()[local_a],
                    &rotation.movement.positions()[local_b],
                ),
                before(a, b),
                3e-13,
            );
        }
    }
}
