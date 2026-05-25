//! Smoke tests: each model + Metropolis runs through Scheduler.

use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use cmc_rs::{
    ClassicalMC, ClusterModel, FromHamiltonianParams, Hamiltonian, HeisenbergModel, Measurable,
    MetropolisCore, Proposable, PottsModel, XYModel,
};

fn run_model<M>(
    extra_params: &[(&str, &str)],
    n_sites_approx: usize,
) -> (f64, f64)
where
    M: Hamiltonian + Measurable + Proposable + ClusterModel + FromHamiltonianParams,
{
    let l = (n_sites_approx as f64).sqrt().round() as usize;
    let mut params = Params::new();
    params.set("Lx", l);
    params.set("Ly", l);
    params.set("J", 1.0);
    params.set("beta", 1.0);
    for &(k, v) in extra_params {
        params.set(k, v);
    }

    let config = RunConfig {
        thermalization_sweeps: 300,
        measurement_sweeps: 500,
        binsize: 100,
        base_seed: 123,
        ..Default::default()
    };

    let backend = RayonBackend::new(1);
    let scheduler = Scheduler::new(backend, config);
    let results = scheduler.run_one::<ClassicalMC<M, MetropolisCore>>(&params);

    let e = results.get("Energy").expect("Energy missing");
    let m = results.get("Magnetization").expect("Magnetization missing");
    (e.mean, m.mean)
}

#[test]
fn test_potts_end_to_end() {
    let (e, m) = run_model::<PottsModel>(&[("q", "4")], 16);
    // At beta=1, ferromagnetic Potts should have negative energy
    assert!(e < 0.0, "Energy should be negative, got {}", e);
    assert!(
        (0.0..=1.0).contains(&m),
        "Magnetization in [0,1], got {}",
        m
    );
}

#[test]
fn test_xy_end_to_end() {
    let (_e, m) = run_model::<XYModel>(&[], 16);
    assert!(
        (0.0..=1.0).contains(&m),
        "Magnetization in [0,1], got {}",
        m
    );
}

#[test]
fn test_heisenberg_end_to_end() {
    let (_e, m) = run_model::<HeisenbergModel>(&[], 16);
    assert!(
        (0.0..=1.0).contains(&m),
        "Magnetization in [0,1], got {}",
        m
    );
}
