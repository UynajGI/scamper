use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use qmc_rs::lattice::ContinuousLatticeEngine;
use qmc_rs::QmcKernel;
use qmc_rs::{
    CsrGraph, EdgeCoupling, LatticeConfiguration, LatticeSpinQmc, SpinLatticeModel,
    SpinModelBuilder, SpinSpace, UpdateSchedule,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

#[test]
fn smoke_arbitrary_graph_and_spin_are_first_class() {
    let graph = CsrGraph::from_adjacency(&[
        vec![1, 3],
        vec![0, 2, 4],
        vec![1, 5],
        vec![0, 4],
        vec![1, 3, 5],
        vec![2, 4],
    ])
    .expect("graph");
    let space = SpinSpace::site_resolved(vec![1, 2, 3, 1, 2, 3]).expect("space");
    let model = SpinModelBuilder::new(graph, space)
        .uniform_edge(EdgeCoupling::xxz(-0.6, 0.3))
        .build()
        .expect("model");
    let mut configuration =
        LatticeConfiguration::new(4.0, vec![0, 1, 1, 0, 1, 2], &model).expect("configuration");
    let mut engine = ContinuousLatticeEngine::new(model, UpdateSchedule::new(8, 4, 64));
    engine.set_validate_each_sweep(true);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(90210);
    for _ in 0..200 {
        engine.sweep(&mut configuration, &mut rng).expect("sweep");
    }
    configuration.validate(engine.model()).expect("valid");
}

#[test]
fn smoke_carlo_adapter_runs_on_square_and_edge_list_graphs() {
    for (topology, extra) in [("square", None), ("edges", Some("0-1,1-2,2-3,3-0"))] {
        let mut params = Params::new();
        params.set("beta", 2.0);
        params.set("model", "xxz");
        params.set("topology", topology);
        params.set("two_s", 1);
        params.set("J_xy", -0.5);
        params.set("J_z", 0.2);
        params.set("validate_each_sweep", true);
        if topology == "square" {
            params.set("L", 2);
            params.set("pbc", false);
        } else if let Some(edges) = extra {
            params.set("n_sites", 4);
            params.set("edges", edges);
        }
        let run = RunConfig {
            thermalization_sweeps: 200,
            measurement_sweeps: 400,
            binsize: 20,
            base_seed: 42,
            ..Default::default()
        };
        let results = Scheduler::new(RayonBackend::new(1), run).run_one::<LatticeSpinQmc>(&params);
        assert!(results.get("EnergyPerSite").is_some());
        assert!(results.get("AverageSign").is_some());
    }
}

#[test]
fn antiferromagnetic_heisenberg_uses_marshall_gauge() {
    let graph = CsrGraph::square(3, 2, false).expect("graph");
    let model = SpinLatticeModel::heisenberg(graph, 1, 1.0).expect("bipartite model");
    assert!(model.gauge().contains(&-1));
    assert!(model.gauge().contains(&1));
}

#[test]
fn continuous_time_dimer_energy_matches_exact_partition_function() {
    let beta = 3.0;
    let j_xy = -0.8;
    let j_z = 0.3;
    let graph = CsrGraph::chain(2, false).expect("graph");
    let model = SpinLatticeModel::xxz(graph, 1, j_xy, j_z).expect("model");
    let mut configuration =
        LatticeConfiguration::new(beta, vec![0, 1], &model).expect("configuration");
    let mut engine = ContinuousLatticeEngine::new(model, UpdateSchedule::new(4, 4, 64));
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x5EED_D1A3);
    let mut energy_sum = 0.0;
    let mut samples = 0_u64;
    for sweep in 0..60_000 {
        engine.sweep(&mut configuration, &mut rng).expect("sweep");
        if sweep >= 10_000 && sweep % 5 == 0 {
            energy_sum +=
                engine.model().constant_shift() - configuration.expansion_order() as f64 / beta;
            samples += 1;
        }
    }
    let measured = energy_sum / samples as f64;
    let levels = [
        0.25 * j_z,
        0.25 * j_z,
        -0.25 * j_z + 0.5 * j_xy,
        -0.25 * j_z - 0.5 * j_xy,
    ];
    let partition = levels
        .iter()
        .map(|energy| (-beta * *energy).exp())
        .sum::<f64>();
    let exact = levels
        .iter()
        .map(|energy| *energy * (-beta * *energy).exp())
        .sum::<f64>()
        / partition;
    assert!(
        (measured - exact).abs() < 0.04,
        "measured {measured}, exact {exact}"
    );
}
