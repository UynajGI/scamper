use cmc_rs::{
    build_chain, BatchEnergyPatch, BatchSpinMove, CanonicalEnsemble, Ensemble, Hamiltonian,
    IsingModel, SiteSpinMove, Spin, System, ThermodynamicDelta, TrialEvaluator,
};

#[test]
fn site_trial_is_transactional_until_commit() {
    let model = IsingModel::new(1.0);
    let mut system = System::new(build_chain(4, true), 1, 1.0, 0.7);
    system.recompute_energy(&model);
    let before_spins = system.spins.clone();
    let before_energy = system.energy;
    let movement = SiteSpinMove::new(0, Spin::from_slice(&[-1.0]));
    let mut patch = cmc_rs::EnergyPatch::default();

    let delta = system.evaluate_trial(&model, &movement, &mut patch);
    assert_eq!(system.spins, before_spins);
    assert_eq!(system.energy, before_energy);
    assert_eq!(delta.energy, patch.delta_energy);

    <System as TrialEvaluator<IsingModel, SiteSpinMove>>::commit_trial(
        &mut system,
        &movement,
        &patch,
    );
    assert!(system.energy_error(&model).abs() < 1e-12);
}

#[test]
fn batch_trial_matches_full_recomputation() {
    let model = IsingModel::new(1.0);
    let mut system = System::new(build_chain(12, true), 1, 1.0, 0.4);
    system.recompute_energy(&model);
    let mut movement = BatchSpinMove::with_capacity(1, 3);
    movement.push(1, &[-1.0]);
    movement.push(4, &[-1.0]);
    movement.push(10, &[-1.0]);
    let mut patch = BatchEnergyPatch::default();

    let delta = system.evaluate_trial(&model, &movement, &mut patch);
    let mut expected_spins = system.spins.clone();
    expected_spins[1] = -1.0;
    expected_spins[4] = -1.0;
    expected_spins[10] = -1.0;
    let expected =
        model.compute_total_energy(&expected_spins, &system.lattice, 1.0) - system.energy;
    assert!((delta.energy - expected).abs() < 1e-12);

    <System as TrialEvaluator<IsingModel, BatchSpinMove>>::commit_trial(
        &mut system,
        &movement,
        &patch,
    );
    assert!(system.energy_error(&model).abs() < 1e-12);
}

#[test]
fn canonical_ensemble_is_independent_of_move_representation() {
    let target = CanonicalEnsemble::new(2.0);
    let delta = ThermodynamicDelta {
        energy: 1.5,
        log_jacobian: 0.25,
        ..Default::default()
    };
    assert!((target.log_weight_ratio(&delta) + 2.75).abs() < 1e-14);
}

struct ThreeBodyProduct;

impl Hamiltonian for ThreeBodyProduct {
    fn spin_dim(&self) -> usize {
        1
    }

    fn coupling(&self) -> f64 {
        1.0
    }

    fn local_energy(
        &self,
        spins: &[f64],
        _lattice: &cmc_rs::CsrLattice,
        site: usize,
        _beta: f64,
        proposed: &[f64],
    ) -> f64 {
        -spins
            .iter()
            .enumerate()
            .map(|(index, spin)| if index == site { proposed[0] } else { *spin })
            .product::<f64>()
    }

    fn compute_total_energy(
        &self,
        spins: &[f64],
        _lattice: &cmc_rs::CsrLattice,
        _beta: f64,
    ) -> f64 {
        -spins.iter().product::<f64>()
    }
}

#[test]
fn direct_multibody_hamiltonian_uses_exact_scratch_batch_path() {
    let model = ThreeBodyProduct;
    let mut system = System::new(build_chain(3, true), 1, 1.0, 0.6);
    system.recompute_energy(&model);
    let mut movement = BatchSpinMove::with_capacity(1, 2);
    movement.push(0, &[-1.0]);
    movement.push(1, &[0.5]);
    let mut patch = BatchEnergyPatch::default();

    let delta = system.evaluate_trial(&model, &movement, &mut patch);
    assert!((delta.energy - 1.5).abs() < 1e-12);
    <System as TrialEvaluator<ThreeBodyProduct, BatchSpinMove>>::commit_trial(
        &mut system,
        &movement,
        &patch,
    );
    assert!(system.energy_error(&model).abs() < 1e-12);
}
