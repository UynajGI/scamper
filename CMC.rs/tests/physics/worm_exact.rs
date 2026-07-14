use super::common::{assert_close, enumerate_ising};
use cmc_rs::{
    build_chain, enumerate_ising_graph_expansion, Hamiltonian, IsingGraphWormModel, IsingModel,
};

#[test]
fn high_temperature_graph_partition_identity_matches_spin_enumeration() {
    let beta = 0.47;
    let coupling = 1.2;
    let lattice = build_chain(5, true);
    let worm = IsingGraphWormModel::new(lattice.clone(), beta, coupling).unwrap();
    let graph = enumerate_ising_graph_expansion(&worm).unwrap();

    let spin_partition: f64 = enumerate_ising(lattice.n_sites)
        .iter()
        .map(|spins| {
            let energy = IsingModel::new(coupling).compute_total_energy(spins, &lattice, 1.0);
            (-beta * energy).exp()
        })
        .sum();
    let prefactor = 2.0f64.powi(lattice.n_sites as i32)
        * worm
            .edge_couplings()
            .iter()
            .map(|&edge_coupling| (beta * edge_coupling).cosh())
            .product::<f64>();
    assert_close(
        spin_partition,
        prefactor * graph.log_reduced_partition.exp(),
        2e-12,
    );
}

#[test]
fn graph_energy_estimator_average_matches_exact_spin_energy() {
    let beta = 0.39;
    let coupling = 0.8;
    let lattice = build_chain(6, true);
    let worm = IsingGraphWormModel::new(lattice.clone(), beta, coupling).unwrap();
    let graph = enumerate_ising_graph_expansion(&worm).unwrap();
    let model = IsingModel::new(coupling);
    let mut z = 0.0;
    let mut e = 0.0;
    for spins in enumerate_ising(lattice.n_sites) {
        let energy = model.compute_total_energy(&spins, &lattice, 1.0);
        let weight = (-beta * energy).exp();
        z += weight;
        e += weight * energy;
    }
    assert_close(graph.mean_energy, e / z, 4e-13);
}

#[test]
fn worm_edge_toggle_is_transactional_and_log_weight_reversible() {
    use cmc_rs::{IsingWormStep, WormModel, WormState};
    let model = IsingGraphWormModel::new(build_chain(4, true), 0.4, 1.0).unwrap();
    let mut state = WormState::new(model.empty_configuration());
    state.open(0).unwrap();
    let edge_id = model.lattice().edge_ids[model.lattice().offsets[0]];
    let step = IsingWormStep { edge_id };
    let before = state.configuration().clone();
    let mut patch = Default::default();
    let forward = model.evaluate_step(&state, &step, &mut patch).unwrap();
    assert_eq!(state.configuration(), &before);
    model.commit_step(&mut state, &step, &patch);
    state.move_head(forward.new_head).unwrap();
    let mut reverse_patch = Default::default();
    let reverse = model
        .evaluate_step(&state, &step, &mut reverse_patch)
        .unwrap();
    assert_close(
        forward.log_weight_ratio + reverse.log_weight_ratio,
        0.0,
        2e-15,
    );
    model.commit_step(&mut state, &step, &reverse_patch);
    state.move_head(reverse.new_head).unwrap();
    model.validate_state(&state).unwrap();
    assert_eq!(state.configuration(), &before);
}
