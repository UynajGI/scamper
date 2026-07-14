use super::common::{assert_close, exact_ising_moments};
use cmc_rs::{
    build_chain, canonical_reweight, enumerate_ising_density_of_states, LogDensityOfStates,
};

#[test]
fn exact_ising_dos_counts_every_microstate_once() {
    for n in 2..=8 {
        let exact =
            enumerate_ising_density_of_states(&build_chain(n, true), &cmc_rs::IsingModel::new(1.0))
                .unwrap();
        assert_eq!(exact.states(), 1u64 << n);
        assert_eq!(exact.degeneracies().iter().sum::<u64>(), 1u64 << n);
    }
}

#[test]
fn discrete_dos_reweighting_matches_direct_microstate_enumeration() {
    let lattice = build_chain(6, true);
    let exact = enumerate_ising_density_of_states(&lattice, &cmc_rs::IsingModel::new(1.0)).unwrap();
    let axis = exact.axis().unwrap();
    let dos = exact.log_density().unwrap();
    for beta in [0.0, 0.17, 0.63, 1.4] {
        let reweighted = canonical_reweight(&axis, &dos, beta).unwrap();
        let (_, mean_e, mean_e2, _) = exact_ising_moments(&lattice, 1.0, beta);
        assert_close(reweighted.mean_energy(), mean_e, 3e-13);
        assert_close(reweighted.mean_energy_squared(), mean_e2, 3e-12);
        assert_close(reweighted.probabilities().iter().sum(), 1.0, 3e-15);
        assert_close(
            reweighted.heat_capacity(),
            beta * beta * (mean_e2 - mean_e * mean_e),
            4e-12,
        );
    }
}

#[test]
fn dos_additive_gauge_does_not_change_canonical_probabilities() {
    let lattice = build_chain(4, true);
    let exact = enumerate_ising_density_of_states(&lattice, &cmc_rs::IsingModel::new(1.0)).unwrap();
    let axis = exact.axis().unwrap();
    let original = exact.log_density().unwrap();
    let shifted = LogDensityOfStates::from_values(
        original
            .values()
            .iter()
            .map(|value| value + 1_000.0)
            .collect(),
        vec![true; original.bins()],
    )
    .unwrap();
    let left = canonical_reweight(&axis, &original, 0.71).unwrap();
    let right = canonical_reweight(&axis, &shifted, 0.71).unwrap();
    for (&a, &b) in left.probabilities().iter().zip(right.probabilities()) {
        assert_close(a, b, 2e-14);
    }
    assert_close(left.mean_energy(), right.mean_energy(), 2e-13);
}

#[test]
fn wang_landau_lifecycle_and_checkpoint_are_exact_state_machine_properties() {
    use cmc_rs::{WangLandauConfig, WangLandauPhase, WangLandauState, WangLandauTermination};
    let config = WangLandauConfig {
        initial_log_f: 1.0,
        final_log_f: 0.25,
        flatness: 1.0,
        flatness_check_interval: 1,
        discovery_sweeps: 0,
        one_over_t_threshold: 0.0,
        max_adaptation_sweeps: 100,
        minimum_visited_fraction: 1.0,
    };
    let mut state = WangLandauState::new(2, config).unwrap();
    assert_eq!(state.phase(), WangLandauPhase::Adaptation);
    for _ in 0..2 {
        state.record_visit(0);
        state.record_visit(1);
        state.finish_sweep();
    }
    assert_eq!(state.phase(), WangLandauPhase::FrozenProduction);
    assert_eq!(state.termination(), Some(WangLandauTermination::Converged));
    let snapshot = state.save_snapshot();
    let restored = WangLandauState::load_snapshot(&snapshot).unwrap();
    assert_eq!(restored, state);
}

#[test]
fn exact_dos_preserves_the_four_physically_distinct_weighted_levels() {
    use cmc_rs::{Bond, BondType, CsrLattice, IsingModel};
    let epsilon = 1e-11;
    let lattice = CsrLattice::from_edges(
        3,
        vec![
            Bond::new(0, 1, BondType::Generic, 1.0),
            Bond::new(1, 2, BondType::Generic, 1.0 + epsilon),
        ],
    );
    let exact = enumerate_ising_density_of_states(&lattice, &IsingModel::new(1.0)).unwrap();

    // E = -(s0*s1 + (1+epsilon)*s1*s2), so the exact physical
    // levels are ±(2+epsilon) and ±epsilon.  The small middle levels
    // are real splittings caused by the unequal bond weights, not roundoff
    // that may be merged into E=0.
    let expected = [-2.0 - epsilon, -epsilon, epsilon, 2.0 + epsilon];
    assert_eq!(
        exact.energies().len(),
        expected.len(),
        "levels={:?}",
        exact.energies()
    );
    for (&actual, &target) in exact.energies().iter().zip(&expected) {
        assert_close(actual, target, 2e-15);
    }
    assert_eq!(exact.degeneracies(), &[2, 2, 2, 2]);
    assert_eq!(exact.states(), 8);
}
