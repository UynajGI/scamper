use carlo_rs::{Context, FromParams, MonteCarlo, Params};
use cmc_rs::{
    build_chain, enumerate_ising_graph_expansion, Bond, BondType, CsrLattice, IsingGraphWormMC,
    IsingGraphWormModel, IsingWormStep, WormConfig, WormKernel, WormModel, WormSector, WormState,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

fn ring_model(beta: f64) -> IsingGraphWormModel {
    IsingGraphWormModel::new(build_chain(4, true), beta, 1.0).unwrap()
}

fn kernel(beta: f64, track_pairs: bool) -> WormKernel<IsingGraphWormModel> {
    let model = ring_model(beta);
    let config = WormConfig {
        local_updates_per_sweep: 4,
        close_probability: 0.25,
        log_worm_fugacity: (1.0f64 / 4.0).ln(),
        track_endpoint_pairs: track_pairs,
        cache_audit_interval: 17,
    };
    let configuration = model.empty_configuration();
    WormKernel::new(model, configuration, config).unwrap()
}

#[test]
fn generic_worm_state_enforces_sector_structure() {
    let mut state = WormState::new(vec![false; 2]);
    assert_eq!(state.sector(), WormSector::Physical);
    state.open(3usize).unwrap();
    assert_eq!(state.head(), Some(&3));
    assert!(state.close().is_ok());
    assert!(WormState::from_parts(vec![false], WormSector::Physical, Some(0), None).is_err());
}

#[test]
fn ising_worm_rejects_sign_problem_and_self_loops() {
    assert!(IsingGraphWormModel::new(build_chain(4, true), 0.5, -1.0).is_err());
    let lattice = CsrLattice::from_edges(2, vec![Bond::new(0, 0, BondType::Generic, 1.0)]);
    assert!(IsingGraphWormModel::new(lattice, 0.5, 1.0).is_err());
}

#[test]
fn local_step_is_transactional_and_reversible() {
    let model = ring_model(0.4);
    let mut state = WormState::new(model.empty_configuration());
    state.open(0).unwrap();
    let edge_id = model.lattice().edge_ids[model.lattice().offsets[0]];
    let step = IsingWormStep { edge_id };
    let mut patch = Default::default();
    let forward = model.evaluate_step(&state, &step, &mut patch).unwrap();
    assert!(!state.configuration().occupied()[edge_id]);
    model.commit_step(&mut state, &step, &patch);
    state.move_head(forward.new_head).unwrap();
    model.validate_state(&state).unwrap();
    assert!(state.configuration().occupied()[edge_id]);

    let mut reverse_patch = Default::default();
    let reverse = model
        .evaluate_step(&state, &step, &mut reverse_patch)
        .unwrap();
    assert!((forward.log_weight_ratio + reverse.log_weight_ratio).abs() < 1e-14);
    model.commit_step(&mut state, &step, &reverse_patch);
    state.move_head(reverse.new_head).unwrap();
    model.validate_state(&state).unwrap();
    assert!(state.head() == state.tail());
    assert!(state
        .configuration()
        .occupied()
        .iter()
        .all(|occupied| !occupied));
}

#[test]
fn open_close_and_step_branch_hastings_terms_are_reciprocal() {
    let n_sites = 7.0f64;
    let close_probability = 0.3f64;
    let log_eta = -1.7f64;
    let open = log_eta + close_probability.ln() + n_sites.ln();
    let close = -log_eta - close_probability.ln() - n_sites.ln();
    assert!((open + close).abs() < 1e-14);

    let coincident_to_open = 0.0 - (1.0 - close_probability).ln();
    let open_to_coincident = (1.0 - close_probability).ln() - 0.0;
    assert!((coincident_to_open + open_to_coincident).abs() < 1e-14);
}

#[test]
fn irregular_weighted_parallel_graph_runs_with_exact_cache_audits() {
    let lattice = CsrLattice::from_edges(
        4,
        vec![
            Bond::new(0, 1, BondType::Generic, 1.0),
            Bond::new(0, 1, BondType::Generic, 0.6),
            Bond::new(1, 2, BondType::Generic, 0.8),
            Bond::new(2, 3, BondType::Generic, 1.2),
            Bond::new(3, 0, BondType::Generic, 0.7),
        ],
    );
    let model = IsingGraphWormModel::new(lattice, 0.5, 1.0).unwrap();
    let configuration = model.empty_configuration();
    let mut kernel = WormKernel::new(
        model,
        configuration,
        WormConfig {
            local_updates_per_sweep: 9,
            close_probability: 0.2,
            log_worm_fugacity: -1.0,
            track_endpoint_pairs: true,
            cache_audit_interval: 1,
        },
    )
    .unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0xabc1);
    for _ in 0..2_000 {
        kernel.sweep(&mut rng).unwrap();
    }
    kernel.validate().unwrap();
    assert!(kernel.statistics().completed_worms > 0);
}

#[test]
fn zero_temperature_parameter_beta_zero_remains_finite_and_empty() {
    let mut kernel = kernel(0.0, false);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0xabc2);
    for _ in 0..10_000 {
        kernel.local_update(&mut rng).unwrap();
    }
    kernel.validate().unwrap();
    assert_eq!(kernel.state().configuration().occupied_edges(), 0);
    assert_eq!(
        kernel
            .model()
            .energy_estimator(kernel.state().configuration()),
        0.0
    );
}

#[test]
fn exact_ring_graph_expansion_has_only_empty_and_full_cycles() {
    let model = ring_model(0.45);
    let exact = enumerate_ising_graph_expansion(&model).unwrap();
    assert_eq!(exact.physical_configurations, 2);
    let t = 0.45f64.tanh();
    let full_probability = t.powi(4) / (1.0 + t.powi(4));
    assert!((exact.mean_occupied_edges - 4.0 * full_probability).abs() < 1e-13);
    for probability in exact.edge_occupation_probabilities {
        assert!((probability - full_probability).abs() < 1e-13);
    }
}

#[test]
fn persistent_worm_matches_exact_graph_energy_and_endpoint_correlation() {
    let beta = 0.45;
    let mut kernel = kernel(beta, true);
    let exact = enumerate_ising_graph_expansion(kernel.model()).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x5eed);

    for _ in 0..50_000 {
        kernel.local_update(&mut rng).unwrap();
    }
    let mut physical_samples = 0u64;
    let mut occupied_sum = 0.0;
    let mut energy_sum = 0.0;
    for _ in 0..600_000 {
        kernel.local_update(&mut rng).unwrap();
        if kernel.state().is_physical() {
            physical_samples += 1;
            occupied_sum += kernel.state().configuration().occupied_edges() as f64;
            energy_sum += kernel
                .model()
                .energy_estimator(kernel.state().configuration());
        }
    }
    assert!(physical_samples > 50_000);
    let mean_occupied = occupied_sum / physical_samples as f64;
    let mean_energy = energy_sum / physical_samples as f64;
    assert!((mean_occupied - exact.mean_occupied_edges).abs() < 0.035);
    assert!((mean_energy - exact.mean_energy).abs() < 0.06);

    let histogram = kernel.endpoint_pairs().unwrap();
    let measured = histogram.correlation_ratio(0, 1).unwrap();
    let t = beta.tanh();
    let expected = (t + t.powi(3)) / (1.0 + t.powi(4));
    assert!(
        (measured - expected).abs() < 0.06,
        "{measured} vs {expected}"
    );
}

#[test]
fn json_snapshot_restores_exact_future_kernel_trajectory() {
    let model = ring_model(0.37);
    let config = WormConfig {
        local_updates_per_sweep: 7,
        close_probability: 0.31,
        log_worm_fugacity: -1.2,
        track_endpoint_pairs: true,
        cache_audit_interval: 11,
    };
    let mut original = IsingGraphWormMC::new(model.clone(), config.clone()).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(91);
    for _ in 0..300 {
        original.kernel_mut().local_update(&mut rng).unwrap();
    }
    let snapshot = original.save_snapshot();
    let mut restored = IsingGraphWormMC::new(model, config).unwrap();
    restored.load_snapshot(&snapshot).unwrap();
    assert_eq!(restored.save_snapshot(), snapshot);

    let mut restored_rng = rng.clone();
    for _ in 0..1_000 {
        original.kernel_mut().local_update(&mut rng).unwrap();
        restored
            .kernel_mut()
            .local_update(&mut restored_rng)
            .unwrap();
    }
    assert_eq!(original.save_snapshot(), restored.save_snapshot());
}

#[test]
fn checkpoint_rejects_inconsistent_transition_counters() {
    let model = ring_model(0.37);
    let config = WormConfig {
        local_updates_per_sweep: 4,
        close_probability: 0.25,
        log_worm_fugacity: -1.0,
        track_endpoint_pairs: false,
        cache_audit_interval: 0,
    };
    let mut mc = IsingGraphWormMC::new(model.clone(), config.clone()).unwrap();
    let mut snapshot = mc.save_snapshot();
    snapshot["runtime"]["statistics"]["close_attempts"] = serde_json::json!(1);
    snapshot["runtime"]["statistics"]["close_accepts"] = serde_json::json!(1);
    snapshot["runtime"]["statistics"]["completed_worms"] = serde_json::json!(1);
    snapshot["runtime"]["statistics"]["physical_visits"] = serde_json::json!(1);
    let error = mc.load_snapshot(&snapshot).unwrap_err();
    assert!(error
        .to_string()
        .contains("closures exceed accepted openings"));
}

#[test]
fn scheduler_ready_adapter_constructs_and_measures() {
    let mut params = Params::new();
    params.set("lattice_type", "chain");
    params.set("L", 6);
    params.set("pbc", true);
    params.set("beta", 0.4);
    params.set("J", 1.0);
    params.set("worm_updates_per_sweep", 12);
    params.set("worm_track_endpoint_pairs", true);

    let mut rng = Xoshiro256PlusPlus::seed_from_u64(17);
    IsingGraphWormMC::validate_params(&params).unwrap();
    let mut mc = IsingGraphWormMC::from_params(&params, &mut rng).unwrap();
    let mut context = Context::new(rng, 0);
    for _ in 0..200 {
        mc.sweep(&mut context);
        mc.measure(&mut context);
        context.advance_sweep();
    }
    mc.kernel().validate().unwrap();
    let results = context.finalize_measurements();
    assert!(results.contains_key("WormSector"));
    assert!(results.contains_key("PhysicalSector"));
    assert!(results.contains_key("Energy"));
}
