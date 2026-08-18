//! Convergence robustness of the Wang–Landau flatness gate.
//!
//! MATURITY_ASSESSMENT.md item 19: "`minimum_visited_fraction` needs
//! hand-tuning to the reachable-bin count (silent non-convergence risk
//! otherwise)". The estimator now auto-derives the reachable set as the
//! discovery plateau of the walk: when the configured fraction demands more
//! visited bins than the plateau provides, the estimate terminates loudly
//! with [`WangLandauTermination::UnreachableBins`] instead of silently
//! burning sweeps up to the maximum-sweep guard with an unconverged DOS.
//!
//! The system is the validated weighted 6-ring on a 14-bin `BinnedAxis`
//! where exactly 6 bins are physically reachable (`ceil(0.4 · 14) = 6`):
//! the historical hand-tuned configuration. Demanding all 14 bins is the
//! historical silent-failure configuration.

use cmc_rs::{
    Algorithm, BinnedAxis, Bond, BondType, CsrLattice, IsingModel, MacrostateAxis, SimulationPhase,
    System, WangLandauConfig, WangLandauCore, WangLandauState, WangLandauTermination,
};
use rand::SeedableRng;

/// Weighted periodic 6-ring with alternating bond weights (spectrum split
/// into eight levels, only 6 of the 14 unit-width bins occupied).
fn weighted_ring() -> CsrLattice {
    let edges: Vec<Bond> = (0..6)
        .map(|site| {
            let weight = if site % 2 == 0 { 1.0 } else { 1.1 };
            let target = (site + 1) % 6;
            Bond::new(site, target, BondType::Generic, weight)
        })
        .collect();
    CsrLattice::from_edges(6, edges)
}

fn test_config(minimum_visited_fraction: f64) -> WangLandauConfig {
    WangLandauConfig {
        initial_log_f: 1.0,
        final_log_f: 1.0 / 4096.0,
        flatness: 0.8,
        // Small interval so the 500-check stall limit fires after ~2 500
        // sweeps; the default interval of 100 would take ~50 000.
        flatness_check_interval: 5,
        discovery_sweeps: 0,
        one_over_t_threshold: 0.0,
        max_adaptation_sweeps: 2_000_000,
        minimum_visited_fraction,
    }
}

/// Adapt until the estimator leaves the adaptive phases; returns the kernel
/// with a sweep cap so a regression to silent non-convergence cannot hang.
fn run_until_terminal(
    axis: &BinnedAxis,
    lattice: &CsrLattice,
    model: &IsingModel,
    config: &WangLandauConfig,
    seed: u64,
    sweep_cap: u64,
) -> WangLandauCore<cmc_rs::BinnedAxis> {
    let mut system = System::new(lattice.clone(), 1, 1.0, 0.0);
    system.recompute_energy(model);
    assert!(
        axis.bin(system.energy).is_some(),
        "cold start energy must lie on the axis"
    );
    let mut kernel = WangLandauCore::new(*axis, config.clone()).unwrap();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut sweeps = 0_u64;
    while kernel.estimator().is_adaptive() {
        kernel.sweep_with_phase(
            &mut system,
            model,
            &mut rng,
            SimulationPhase::Thermalization,
        );
        sweeps += 1;
        assert!(
            sweeps < sweep_cap,
            "Wang-Landau did not terminate within {sweep_cap} sweeps \
             (silent non-convergence regression)"
        );
    }
    kernel
}

#[test]
fn wang_landau_unattainable_visited_fraction_fails_loudly() {
    let lattice = weighted_ring();
    let model = IsingModel::new(1.0);
    let axis = BinnedAxis::new(-7.0, 7.0, 14).unwrap();
    // The historical silent-failure configuration: every axis bin demanded,
    // but 8 of the 14 bins contain no physical state.
    let kernel = run_until_terminal(&axis, &lattice, &model, &test_config(1.0), 7, 20_000);
    let state = kernel.estimator();

    assert_eq!(
        state.termination(),
        Some(WangLandauTermination::UnreachableBins),
        "an unattainable visited fraction must terminate loudly, got {:?}",
        state.termination()
    );
    assert!(!state.is_frozen(), "no production may follow a failed gate");
    assert!(
        state.adaptation_sweeps() < 20_000,
        "loud termination must fire far below the 2M sweep guard, took {}",
        state.adaptation_sweeps()
    );

    // The auto-derived reachable set is the true one: exactly the 6 bins
    // that carry enumerated states.
    let visited: Vec<usize> = (0..axis.bins())
        .filter(|&bin| state.log_density().is_visited(bin))
        .collect();
    assert_eq!(visited, vec![0, 4, 5, 8, 9, 13]);

    // Versioned checkpoints round-trip the failure and its evidence.
    let snapshot = state.save_snapshot();
    let restored = WangLandauState::load_snapshot(&snapshot).unwrap();
    assert_eq!(
        restored.termination(),
        Some(WangLandauTermination::UnreachableBins)
    );
    assert_eq!(restored.phase(), state.phase());
}

#[test]
fn wang_landau_sane_fraction_on_the_same_axis_still_converges() {
    // The hand-tuned historical value ceil(0.4·14) = 6 must keep converging
    // on the same system: the stall detector may not fire when the gate is
    // attainable.
    let lattice = weighted_ring();
    let model = IsingModel::new(1.0);
    let axis = BinnedAxis::new(-7.0, 7.0, 14).unwrap();
    let kernel = run_until_terminal(&axis, &lattice, &model, &test_config(0.4), 11, 500_000);
    assert_eq!(
        kernel.estimator().termination(),
        Some(WangLandauTermination::Converged)
    );
    assert!(kernel.estimator().is_frozen());
}

#[test]
fn wang_landau_checkpoint_rejects_fabricated_unreachable_termination() {
    let lattice = weighted_ring();
    let model = IsingModel::new(1.0);
    let axis = BinnedAxis::new(-7.0, 7.0, 14).unwrap();
    let kernel = run_until_terminal(&axis, &lattice, &model, &test_config(1.0), 13, 20_000);
    let mut snapshot = kernel.estimator().save_snapshot();

    // Strip the discovery-stall evidence: a checkpoint claiming the loud
    // termination without the plateau proof must be rejected.
    snapshot["discovery_stall_checks"] = serde_json::json!(0);
    let error = WangLandauState::load_snapshot(&snapshot);
    assert!(
        error.is_err(),
        "unreachable-bins checkpoints require the stall evidence"
    );

    // Claiming enough visited bins for the gate is equally inconsistent.
    let mut snapshot = kernel.estimator().save_snapshot();
    snapshot["last_check_visited_bins"] = serde_json::json!(14);
    assert!(WangLandauState::load_snapshot(&snapshot).is_err());
}

#[test]
fn wang_landau_version_one_checkpoint_without_stall_fields_still_loads() {
    // Version-1 checkpoints predate the discovery-stall fields; loading must
    // keep working (the fields default to a fresh plateau detection).
    let lattice = weighted_ring();
    let model = IsingModel::new(1.0);
    let axis = BinnedAxis::new(-7.0, 7.0, 14).unwrap();
    let kernel = run_until_terminal(&axis, &lattice, &model, &test_config(0.4), 17, 500_000);
    let mut snapshot = kernel.estimator().save_snapshot();
    snapshot
        .as_object_mut()
        .unwrap()
        .remove("last_check_visited_bins");
    snapshot
        .as_object_mut()
        .unwrap()
        .remove("discovery_stall_checks");
    let restored = WangLandauState::load_snapshot(&snapshot)
        .expect("legacy checkpoint without stall fields must load");
    assert_eq!(
        restored.termination(),
        Some(WangLandauTermination::Converged)
    );
}
