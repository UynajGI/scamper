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

#[test]
fn worm_rejects_multi_component_lattices_loudly() {
    use cmc_rs::{Bond, BondType, CsrLattice};
    // Two disjoint bonds: sites {0,1} and {2,3} form separate components.
    let disconnected = CsrLattice::from_edges(
        4,
        vec![
            Bond::new(0, 1, BondType::Generic, 1.0),
            Bond::new(2, 3, BondType::Generic, 1.0),
        ],
    );
    let error = IsingGraphWormModel::new(disconnected, 0.4, 1.0)
        .expect_err("a multi-component lattice must be rejected at input");
    assert!(
        error.to_string().contains("connected"),
        "the rejection must name the connectivity requirement: {error}"
    );

    // An isolated site is its own component with the same silent-freeze
    // failure mode (a defect opened there can never step).
    let isolated = CsrLattice::from_edges(3, vec![Bond::new(0, 1, BondType::Generic, 1.0)]);
    assert!(IsingGraphWormModel::new(isolated, 0.4, 1.0).is_err());

    // The connected counterpart is accepted.
    let connected = CsrLattice::from_edges(
        4,
        vec![
            Bond::new(0, 1, BondType::Generic, 1.0),
            Bond::new(1, 2, BondType::Generic, 1.0),
            Bond::new(2, 3, BondType::Generic, 1.0),
        ],
    );
    assert!(IsingGraphWormModel::new(connected, 0.4, 1.0).is_ok());
}

#[test]
fn worm_energy_agrees_with_spin_metropolis_cross_solver() {
    // Criterion F for the worm: two independent solvers of the same
    // canonical ensemble — the high-temperature graph worm and local
    // Metropolis on spins — must agree on ⟨E⟩ within pooled errors.
    use super::common::zscore_seed_count;
    use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
    use cmc_rs::{ClassicalMC, MetropolisCore};

    let beta = 0.44;
    let n_seeds = zscore_seed_count(8);
    let run = |solver: &str, seed: u64| -> f64 {
        let mut params = Params::new();
        params.set("lattice_type", "square");
        params.set("Lx", 4usize);
        params.set("Ly", 4usize);
        params.set("beta", beta);
        params.set("J", 1.0);
        let config = RunConfig {
            thermalization_sweeps: 2_000,
            measurement_sweeps: 20_000,
            binsize: 500,
            base_seed: seed,
            ..Default::default()
        };
        let results = match solver {
            "worm" => Scheduler::new(RayonBackend::new(1), config)
                .run_one::<cmc_rs::IsingGraphWormMC>(&params),
            "metropolis" => Scheduler::new(RayonBackend::new(1), config)
                .run_one::<ClassicalMC<cmc_rs::IsingModel, MetropolisCore>>(&params),
            other => panic!("unknown solver {other}"),
        };
        results.get("Energy").expect("Energy observable").mean
    };
    let pool = |solver: &str, base: u64| {
        let means: Vec<f64> = (0..n_seeds as u64)
            .map(|seed| run(solver, base + seed))
            .collect();
        let count = means.len() as f64;
        let mean = means.iter().sum::<f64>() / count;
        // Seed-spread stderr of the pooled mean.
        let variance = means
            .iter()
            .map(|value| (value - mean) * (value - mean))
            .sum::<f64>()
            / (count - 1.0);
        (mean, (variance / count).sqrt())
    };
    let (worm_mean, worm_stderr) = pool("worm", 0x0E51);
    let (metro_mean, metro_stderr) = pool("metropolis", 0x0E52);
    let z =
        (worm_mean - metro_mean) / (worm_stderr * worm_stderr + metro_stderr * metro_stderr).sqrt();
    eprintln!(
        "[worm-cross] 4x4 beta={beta}: worm ⟨E⟩ = {worm_mean:.4} ± {worm_stderr:.4}, \
         metropolis ⟨E⟩ = {metro_mean:.4} ± {metro_stderr:.4}, z = {z:+.2}"
    );
    assert!(
        z.abs() < 4.0,
        "worm vs Metropolis ⟨E⟩ disagree: z = {z:.2} ({worm_mean:.4} vs {metro_mean:.4})"
    );
}
