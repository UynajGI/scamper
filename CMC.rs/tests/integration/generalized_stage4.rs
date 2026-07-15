use cmc_rs::{
    build_chain, canonical_reweight, enumerate_ising_density_of_states, Algorithm, BinnedAxis,
    DiscreteAxis, EnergyBiasCore, FixedBias, HarmonicUmbrellaBias, Histogram, IsingModel, LogBias,
    LogDensityOfStates, MacrostateAxis, MulticanonicalBias, PairPotential, SimulationPhase, System,
    WangLandauConfig, WangLandauCore, WangLandauPhase, WangLandauState, WangLandauTermination,
};
use rand::SeedableRng;

#[derive(Debug, Clone, Copy)]
struct SpeciesZeroOnly;

impl PairPotential for SpeciesZeroOnly {
    fn cutoff_squared(&self) -> f64 {
        1.0
    }

    fn energy(&self, _species_i: u16, _species_j: u16, _distance_squared: f64) -> f64 {
        0.0
    }

    fn supports_species(&self, species: u16) -> bool {
        species == 0
    }
}

#[test]
fn continuous_and_discrete_axes_handle_boundaries_explicitly() {
    let binned = BinnedAxis::new(-2.0, 2.0, 4).unwrap();
    assert_eq!(binned.bin(-2.0), Some(0));
    assert_eq!(binned.bin(-1.0), Some(1));
    assert_eq!(binned.bin(2.0), Some(3));
    let just_below_maximum = f64::from_bits(2.0_f64.to_bits() - 1);
    assert_eq!(binned.bin(just_below_maximum), Some(3));
    assert_eq!(binned.bin(2.0 + 1e-12), None);
    assert_eq!(binned.center(0), -1.5);

    let discrete = DiscreteAxis::new(vec![4.0, -4.0, 0.0]).unwrap();
    assert_eq!(discrete.values(), &[-4.0, 0.0, 4.0]);
    assert_eq!(discrete.bin(4.0 + 1e-11), Some(2));
    assert_eq!(discrete.bin(1.0), None);
}

#[test]
fn histogram_flatness_requires_the_requested_coverage() {
    let mut histogram = Histogram::new(4).unwrap();
    for bin in 0..3 {
        for _ in 0..10 {
            histogram.record(bin);
        }
    }
    assert!(histogram.is_flat(0.8, 0.75));
    assert!(!histogram.is_flat(0.8, 1.0));
    for _ in 0..10 {
        histogram.record(3);
    }
    assert!(histogram.is_flat(1.0, 1.0));
}

#[test]
fn exact_four_site_ising_dos_has_known_levels_and_degeneracies() {
    let lattice = build_chain(4, true);
    let exact = enumerate_ising_density_of_states(&lattice, &IsingModel::new(1.0)).unwrap();
    assert_eq!(exact.energies(), &[-4.0, 0.0, 4.0]);
    assert_eq!(exact.degeneracies(), &[2, 12, 2]);
    assert_eq!(exact.states(), 16);
}

#[test]
fn canonical_reweighting_matches_direct_exact_sum() {
    let lattice = build_chain(4, true);
    let exact = enumerate_ising_density_of_states(&lattice, &IsingModel::new(1.0)).unwrap();
    let axis = exact.axis().unwrap();
    let dos = exact.log_density().unwrap();
    let beta = 0.7;
    let estimate = canonical_reweight(&axis, &dos, beta).unwrap();

    let weights: Vec<_> = exact
        .energies()
        .iter()
        .zip(exact.degeneracies())
        .map(|(&energy, &degeneracy)| degeneracy as f64 * (-beta * energy).exp())
        .collect();
    let partition: f64 = weights.iter().sum();
    let direct_mean: f64 = exact
        .energies()
        .iter()
        .zip(&weights)
        .map(|(&energy, &weight)| energy * weight)
        .sum::<f64>()
        / partition;
    assert!((estimate.mean_energy() - direct_mean).abs() < 1e-12);
    assert!((estimate.probabilities().iter().sum::<f64>() - 1.0).abs() < 1e-12);
}

#[test]
fn frozen_biases_use_complete_log_target_weights_and_preserve_unvisited_bins() {
    let axis = DiscreteAxis::new(vec![-1.0, 0.0, 1.0]).unwrap();
    let umbrella = HarmonicUmbrellaBias::new(&axis, 2.0, 0.0, 4.0).unwrap();
    assert!((umbrella.log_weight(0) - 0.0).abs() < 1e-14);
    assert!((umbrella.log_weight(1) - 0.0).abs() < 1e-14);
    assert!((umbrella.log_weight(2) + 4.0).abs() < 1e-14);

    let dos =
        LogDensityOfStates::from_values(vec![0.0, 3.0, 1.0], vec![true, false, true]).unwrap();
    let multicanonical = MulticanonicalBias::from_log_density(&dos).unwrap();
    assert_eq!(multicanonical.log_weight(0), 0.0);
    assert_eq!(multicanonical.log_weight(1), f64::NEG_INFINITY);
    assert_eq!(multicanonical.log_weight(2), -1.0);
}

#[test]
fn wang_landau_state_refines_freezes_and_round_trips_checkpoint() {
    let config = WangLandauConfig {
        initial_log_f: 1.0,
        final_log_f: 0.2,
        flatness: 1.0,
        flatness_check_interval: 1,
        discovery_sweeps: 0,
        one_over_t_threshold: 0.0,
        max_adaptation_sweeps: 100,
        minimum_visited_fraction: 1.0,
    };
    let mut state = WangLandauState::new(3, config).unwrap();
    for _ in 0..3 {
        for bin in 0..3 {
            state.record_visit(bin);
        }
        state.finish_sweep();
    }
    assert_eq!(state.phase(), WangLandauPhase::FrozenProduction);
    assert_eq!(state.termination(), Some(WangLandauTermination::Converged));
    assert_eq!(state.flatness_passes(), 3);
    assert!(state
        .log_density()
        .values()
        .iter()
        .all(|&value| value == 0.0));

    let snapshot = state.save_snapshot();
    let restored = WangLandauState::load_snapshot(&snapshot).unwrap();
    assert_eq!(restored, state);
}

#[test]
fn frozen_bias_rejects_moves_outside_the_represented_axis_transactionally() {
    let lattice = build_chain(2, true);
    let model = IsingModel::new(1.0);
    let mut system = System::new(lattice, 1, 1.0, 0.0);
    system.recompute_energy(&model);
    let original_spins = system.spins.clone();
    let original_energy = system.energy;

    let axis = DiscreteAxis::new(vec![original_energy]).unwrap();
    let bias = FixedBias::new(vec![0.0]).unwrap();
    let mut kernel = EnergyBiasCore::new(axis, bias);
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(9);
    kernel.sweep_with_phase(&mut system, &model, &mut rng, SimulationPhase::Measurement);

    assert_eq!(system.spins, original_spins);
    assert_eq!(system.energy, original_energy);
    assert_eq!(kernel.histogram().total(), 2);
    assert_eq!(kernel.out_of_range_proposals(), 2);
}

#[test]
fn two_site_wang_landau_estimates_equal_degeneracies_then_runs_frozen_production() {
    let lattice = build_chain(2, true);
    let model = IsingModel::new(1.0);
    let exact = enumerate_ising_density_of_states(&lattice, &model).unwrap();
    assert_eq!(exact.degeneracies(), &[2, 2]);
    let axis = exact.axis().unwrap();
    let config = WangLandauConfig {
        initial_log_f: 1.0,
        final_log_f: 1.0 / 64.0,
        flatness: 0.7,
        flatness_check_interval: 2,
        discovery_sweeps: 0,
        one_over_t_threshold: 0.0,
        max_adaptation_sweeps: 10_000,
        minimum_visited_fraction: 1.0,
    };
    let mut kernel = WangLandauCore::new(axis, config).unwrap();
    let mut system = System::new(lattice, 1, 1.0, 0.0);
    system.recompute_energy(&model);
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(42);

    while kernel.estimator().is_adaptive() {
        kernel.sweep_with_phase(
            &mut system,
            &model,
            &mut rng,
            SimulationPhase::Thermalization,
        );
    }
    assert_eq!(
        kernel.estimator().phase(),
        WangLandauPhase::FrozenProduction
    );
    assert_eq!(kernel.estimator().log_density().visited_bins(), 2);
    let shifted = kernel.estimator().log_density().shifted_to_max_zero();
    let difference = (shifted[0].unwrap() - shifted[1].unwrap()).abs();
    assert!(
        difference < 0.75,
        "equal exact degeneracies differed by {difference}"
    );

    for _ in 0..10 {
        kernel.sweep_with_phase(&mut system, &model, &mut rng, SimulationPhase::Measurement);
    }
    assert_eq!(kernel.estimator().production_sweeps(), 10);
    assert_eq!(kernel.estimator().production_histogram().total(), 20);
    assert!(system.energy_error(&model).abs() < 1e-12);
}

#[test]
fn four_site_wang_landau_dos_reweights_close_to_the_exact_reference() {
    let lattice = build_chain(4, true);
    let model = IsingModel::new(1.0);
    let exact = enumerate_ising_density_of_states(&lattice, &model).unwrap();
    let axis = exact.axis().unwrap();
    let config = WangLandauConfig {
        initial_log_f: 1.0,
        final_log_f: 1.0 / 4096.0,
        flatness: 0.8,
        flatness_check_interval: 50,
        discovery_sweeps: 0,
        one_over_t_threshold: 0.0,
        max_adaptation_sweeps: 100_000,
        minimum_visited_fraction: 1.0,
    };
    let mut kernel = WangLandauCore::new(axis.clone(), config).unwrap();
    let mut system = System::new(lattice, 1, 1.0, 0.0);
    system.recompute_energy(&model);
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(77);
    while kernel.estimator().is_adaptive() {
        kernel.sweep_with_phase(
            &mut system,
            &model,
            &mut rng,
            SimulationPhase::Thermalization,
        );
    }
    assert_eq!(
        kernel.estimator().termination(),
        Some(WangLandauTermination::Converged)
    );

    let estimated = kernel.estimator().log_density().shifted_to_max_zero();
    let exact_shifted = exact.log_density().unwrap().shifted_to_max_zero();
    for (estimate, reference) in estimated.iter().zip(&exact_shifted) {
        assert!((estimate.unwrap() - reference.unwrap()).abs() < 0.4);
    }

    let beta = 0.7;
    let wl_thermodynamics =
        canonical_reweight(&axis, kernel.estimator().log_density(), beta).unwrap();
    let exact_dos = exact.log_density().unwrap();
    let exact_thermodynamics = canonical_reweight(&axis, &exact_dos, beta).unwrap();
    assert!((wl_thermodynamics.mean_energy() - exact_thermodynamics.mean_energy()).abs() < 0.35);
    assert!((wl_thermodynamics.heat_capacity() - exact_thermodynamics.heat_capacity()).abs() < 0.3);
}

#[test]
fn dos_normalization_is_an_additive_gauge_transformation() {
    let mut dos = LogDensityOfStates::from_values(vec![8.0, 10.0, 9.0], vec![true; 3]).unwrap();
    let before = dos.value(0) - dos.value(2);
    dos.normalize_max_zero();
    assert_eq!(dos.value(1), 0.0);
    assert_eq!(dos.value(0) - dos.value(2), before);
}

#[test]
fn stage_three_review_rejects_unsupported_insertion_species_up_front() {
    let proposal = cmc_rs::InsertDeleteParticle::from_species(vec![(0, 1.0), (1, 1.0)]).unwrap();
    let error = proposal.validate_potential(&SpeciesZeroOnly).unwrap_err();
    assert!(error.to_string().contains("species 1"));
}

#[test]
fn scheduler_ready_ising_wang_landau_uses_adaptive_run_control() {
    use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
    use cmc_rs::{IsingWangLandau, WangLandauRunControl};

    let mut params = Params::new();
    params.set("L", 2usize);
    params.set("initial_state", "cold");
    params.set("wl_initial_log_f", 1.0);
    params.set("wl_final_log_f", 0.25);
    params.set("wl_flatness", 0.5);
    params.set("wl_flatness_check_interval", 1u64);
    params.set("wl_one_over_t_threshold", 0.0);
    params.set("wl_max_adaptation_sweeps", 1_000u64);
    let scheduler = Scheduler::new(
        RayonBackend::new(1),
        RunConfig {
            base_seed: 12,
            binsize: 2,
            ..Default::default()
        },
    );
    let (simulation, results) = scheduler
        .run_controlled_with_state::<IsingWangLandau, _>(&params, WangLandauRunControl::new(5))
        .unwrap();
    assert!(results.get("Energy").is_some());
    assert_eq!(simulation.estimator().phase(), WangLandauPhase::Finished);
    assert_eq!(simulation.estimator().production_sweeps(), 5);
    assert_eq!(simulation.estimator().log_density().visited_bins(), 2);
}

#[test]
fn discovery_delays_flatness_and_resets_the_adaptation_window() {
    let config = WangLandauConfig {
        initial_log_f: 1.0,
        final_log_f: 0.1,
        flatness: 1.0,
        flatness_check_interval: 1,
        discovery_sweeps: 2,
        one_over_t_threshold: 0.0,
        max_adaptation_sweeps: 100,
        minimum_visited_fraction: 1.0,
    };
    let mut state = WangLandauState::new(1, config).unwrap();
    state.record_visit(0);
    state.finish_sweep();
    assert_eq!(state.phase(), WangLandauPhase::Discovery);
    assert_eq!(state.flatness_checks(), 0);

    state.record_visit(0);
    state.finish_sweep();
    assert_eq!(state.phase(), WangLandauPhase::Adaptation);
    assert_eq!(state.adaptation_histogram().total(), 0);
    assert_eq!(state.refinement_visits(), 0);
    assert_eq!(state.flatness_checks(), 0);

    state.record_visit(0);
    state.finish_sweep();
    assert_eq!(state.flatness_passes(), 1);
}

#[test]
fn maximum_sweep_guard_finishes_an_unconverged_estimate() {
    let config = WangLandauConfig {
        initial_log_f: 1.0,
        final_log_f: 0.1,
        flatness: 1.0,
        flatness_check_interval: 1,
        discovery_sweeps: 0,
        one_over_t_threshold: 0.0,
        max_adaptation_sweeps: 2,
        minimum_visited_fraction: 1.0,
    };
    let mut state = WangLandauState::new(2, config).unwrap();
    for _ in 0..2 {
        state.record_visit(0);
        state.finish_sweep();
    }
    assert_eq!(state.phase(), WangLandauPhase::Finished);
    assert_eq!(
        state.termination(),
        Some(WangLandauTermination::MaximumSweeps)
    );
    assert_eq!(
        WangLandauState::load_snapshot(&state.save_snapshot()).unwrap(),
        state
    );
}

#[test]
fn one_over_t_refinement_converges_without_another_flatness_reset() {
    let config = WangLandauConfig {
        initial_log_f: 1.0,
        final_log_f: 0.1,
        flatness: 1.0,
        flatness_check_interval: 1,
        discovery_sweeps: 0,
        one_over_t_threshold: 0.75,
        max_adaptation_sweeps: 100,
        minimum_visited_fraction: 1.0,
    };
    let mut state = WangLandauState::new(2, config).unwrap();
    state.record_visit(0);
    state.record_visit(1);
    state.finish_sweep();
    assert_eq!(state.refinement(), cmc_rs::WangLandauRefinement::OneOverT);
    for visit in 0..20 {
        state.record_visit(visit % 2);
    }
    state.finish_sweep();
    assert_eq!(state.phase(), WangLandauPhase::FrozenProduction);
    assert_eq!(state.termination(), Some(WangLandauTermination::Converged));
}

#[test]
fn kernel_checkpoint_validates_physical_axis_not_only_bin_count() {
    let axis = DiscreteAxis::new(vec![-2.0, 2.0]).unwrap();
    let config = WangLandauConfig {
        initial_log_f: 1.0,
        final_log_f: 0.1,
        flatness: 0.8,
        flatness_check_interval: 10,
        discovery_sweeps: 0,
        one_over_t_threshold: 0.0,
        max_adaptation_sweeps: 100,
        minimum_visited_fraction: 1.0,
    };
    let kernel = WangLandauCore::new(axis.clone(), config).unwrap();
    let snapshot = kernel.save_snapshot();
    let restored =
        WangLandauCore::from_snapshot(axis, cmc_rs::StandardStrategy::new(), &snapshot).unwrap();
    assert_eq!(restored.estimator(), kernel.estimator());

    let wrong_axis = DiscreteAxis::new(vec![-3.0, 3.0]).unwrap();
    assert!(
        WangLandauCore::from_snapshot(wrong_axis, cmc_rs::StandardStrategy::new(), &snapshot,)
            .is_err()
    );
}

#[test]
fn checkpoint_rejects_inconsistent_lifecycle_metadata() {
    let config = WangLandauConfig {
        initial_log_f: 1.0,
        final_log_f: 0.1,
        flatness: 0.8,
        flatness_check_interval: 10,
        discovery_sweeps: 0,
        one_over_t_threshold: 0.0,
        max_adaptation_sweeps: 100,
        minimum_visited_fraction: 1.0,
    };
    let state = WangLandauState::new(2, config).unwrap();
    let mut snapshot = state.save_snapshot();
    snapshot["termination"] = serde_json::json!("converged");
    assert!(WangLandauState::load_snapshot(&snapshot).is_err());

    let mut inconsistent_counts = state.save_snapshot();
    inconsistent_counts["refinement_visits"] = serde_json::json!(1_u64);
    assert!(WangLandauState::load_snapshot(&inconsistent_counts).is_err());
}

#[test]
fn kernel_checkpoint_restores_schedule_and_counters() {
    let axis = DiscreteAxis::new(vec![-2.0, 2.0]).unwrap();
    let config = WangLandauConfig {
        initial_log_f: 1.0,
        final_log_f: 0.1,
        flatness: 0.8,
        flatness_check_interval: 10,
        discovery_sweeps: 0,
        one_over_t_threshold: 0.0,
        max_adaptation_sweeps: 100,
        minimum_visited_fraction: 1.0,
    };
    let kernel = WangLandauCore::new(axis.clone(), config)
        .unwrap()
        .with_visit_schedule(cmc_rs::VisitSchedule::Sequential)
        .with_energy_check_interval(17);
    let snapshot = kernel.save_snapshot();
    let restored =
        WangLandauCore::from_snapshot(axis, cmc_rs::StandardStrategy::new(), &snapshot).unwrap();
    assert_eq!(restored.sweeps(), 0);
    assert_eq!(restored.out_of_range_proposals(), 0);
    let restored_snapshot = restored.save_snapshot();
    assert_eq!(
        restored_snapshot["visit_schedule"].as_str(),
        Some("sequential")
    );
    assert_eq!(
        restored_snapshot["energy_check_interval"].as_u64(),
        Some(17)
    );
}
