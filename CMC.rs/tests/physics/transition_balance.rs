use super::common::{assert_close, enumerate_ising};
use cmc_rs::{
    Bond, BondType, CanonicalEnsemble, CsrLattice, Ensemble, Hamiltonian, IsingModel,
    MetropolisHastingsAcceptance, ThermodynamicDelta,
};

#[test]
fn exact_single_spin_metropolis_matrix_obeys_detailed_balance() {
    let beta = 0.73;
    let lattice = CsrLattice::from_edges(
        4,
        vec![
            Bond::new(0, 1, BondType::Generic, 0.8),
            Bond::new(0, 1, BondType::Generic, -0.15),
            Bond::new(1, 2, BondType::Generic, 1.1),
            Bond::new(2, 3, BondType::Generic, -0.35),
            Bond::new(3, 0, BondType::Generic, 0.6),
            Bond::new(2, 2, BondType::Generic, 7.0),
        ],
    );
    let model = IsingModel::new(0.9);
    let states = enumerate_ising(4);
    let energies: Vec<_> = states
        .iter()
        .map(|state| model.compute_total_energy(state, &lattice, 1.0))
        .collect();
    let weights: Vec<_> = energies
        .iter()
        .map(|energy| (-beta * energy).exp())
        .collect();
    let partition: f64 = weights.iter().sum();
    let probabilities: Vec<_> = weights.iter().map(|weight| weight / partition).collect();

    let mut transition = vec![vec![0.0; states.len()]; states.len()];
    for (x, state) in states.iter().enumerate() {
        for site in 0..4 {
            let mut target = state.clone();
            target[site] = -target[site];
            let y = states
                .iter()
                .position(|candidate| candidate == &target)
                .unwrap();
            let delta = energies[y] - energies[x];
            let acceptance = (-beta * delta).min(0.0).exp();
            transition[x][y] += acceptance / 4.0;
            transition[x][x] += (1.0 - acceptance) / 4.0;
        }
        assert_close(transition[x].iter().sum(), 1.0, 2e-15);
    }

    for x in 0..states.len() {
        for y in 0..states.len() {
            assert_close(
                probabilities[x] * transition[x][y],
                probabilities[y] * transition[y][x],
                3e-15,
            );
        }
    }
}

#[test]
fn asymmetric_hastings_correction_restores_exact_balance() {
    let beta = 0.4;
    let energy_minus = 0.7;
    let energy_plus = -1.2;
    let q_plus_given_minus: f64 = 0.8;
    let q_minus_given_plus: f64 = 0.2;
    let forward_log_ratio = (q_minus_given_plus / q_plus_given_minus).ln();
    let reverse_log_ratio = -forward_log_ratio;
    let ensemble = CanonicalEnsemble::new(beta);
    let rule = MetropolisHastingsAcceptance;

    let forward_delta = ThermodynamicDelta::energy(energy_plus - energy_minus);
    let reverse_delta = ThermodynamicDelta::energy(energy_minus - energy_plus);
    let forward_log =
        cmc_rs::AcceptanceRule::log_acceptance(&rule, &ensemble, &forward_delta, forward_log_ratio);
    let reverse_log =
        cmc_rs::AcceptanceRule::log_acceptance(&rule, &ensemble, &reverse_delta, reverse_log_ratio);
    let p_minus = (-beta * energy_minus).exp();
    let p_plus = (-beta * energy_plus).exp();
    assert_close(
        p_minus * q_plus_given_minus * forward_log.min(0.0).exp(),
        p_plus * q_minus_given_plus * reverse_log.min(0.0).exp(),
        2e-15,
    );
}

#[test]
fn canonical_target_contains_beta_exactly_once() {
    let ensemble = CanonicalEnsemble::new(2.3);
    let delta = ThermodynamicDelta {
        energy: -1.7,
        log_jacobian: 0.4,
        ..Default::default()
    };
    assert_close(ensemble.log_weight_ratio(&delta), 2.3 * 1.7 + 0.4, 1e-15);
}

#[test]
fn cluster_activation_probabilities_match_fortuin_kasteleyn_weights() {
    use cmc_rs::{Bond, BondType, ClusterAuxiliary, ClusterModel, ONModel, PottsModel};
    use smallvec::smallvec;

    let bond = Bond::new(0, 1, BondType::Generic, 0.7);
    let beta = 0.8;
    let ising = IsingModel::new(1.3);
    let p = ising.cluster_bond_probability(&[1.0], &[1.0], &bond, &ClusterAuxiliary::None, beta);
    assert_close(p, 1.0 - (-2.0 * beta * 1.3 * 0.7).exp(), 1e-15);
    assert_eq!(
        ising.cluster_bond_probability(&[1.0], &[-1.0], &bond, &ClusterAuxiliary::None, beta,),
        0.0
    );

    let potts = PottsModel::new(1.3, 3);
    let p = potts.cluster_bond_probability(&[2.0], &[2.0], &bond, &ClusterAuxiliary::None, beta);
    assert_close(p, 1.0 - (-beta * 1.3 * 0.7).exp(), 1e-15);

    let on = ONModel::<3>::new(1.3);
    let normal = smallvec![1.0, 0.0, 0.0];
    let p = on.cluster_bond_probability(
        &[0.6, 0.8, 0.0],
        &[0.5, 0.0, 0.8660254037844386],
        &bond,
        &ClusterAuxiliary::Reflection(normal),
        beta,
    );
    assert_close(p, 1.0 - (-2.0 * beta * 1.3 * 0.7 * 0.6 * 0.5).exp(), 2e-15);
}
