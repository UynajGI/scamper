//! Smoke tests: each model + Metropolis runs through Scheduler.

use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use cmc_rs::{
    ClassicalMC, ClusterModel, FromHamiltonianParams, Hamiltonian, HeisenbergModel, Initializable,
    IsingModel, Measurable, MetropolisCore, PottsModel, Proposable, XYModel,
};

fn run_model<M>(extra_params: &[(&str, &str)], n_sites_approx: usize) -> (f64, f64)
where
    M: Hamiltonian + Initializable + Measurable + Proposable + ClusterModel + FromHamiltonianParams,
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

// ── Non-square lattice validation ──────────────────────────

#[test]
fn test_triangular_ising_ferro_cools() {
    let mut params = Params::new();
    params.set("Lx", 4usize);
    params.set("Ly", 4usize);
    params.set("lattice_type", "triangular");
    params.set("J", 1.0);
    params.set("beta", 2.0); // well below Tc (βc ≈ 0.27), ordered

    let config = RunConfig {
        thermalization_sweeps: 500,
        measurement_sweeps: 500,
        binsize: 100,
        base_seed: 42,
        ..Default::default()
    };

    let backend = RayonBackend::new(1);
    let scheduler = Scheduler::new(backend, config);
    let results = scheduler.run_one::<ClassicalMC<IsingModel, MetropolisCore>>(&params);

    let e = results.get("Energy").expect("Energy missing");
    let m = results.get("Magnetization").expect("Magnetization missing");

    // At low T, ferromagnetic triangular Ising should order
    assert!(
        e.mean < -30.0,
        "Expected strongly negative energy, got {}",
        e.mean
    );
    assert!(m.mean > 0.7, "Expected high magnetization, got {}", m.mean);
}

#[test]
fn test_kagome_ising_af_low_t() {
    let mut params = Params::new();
    params.set("Lx", 3usize);
    params.set("Ly", 3usize);
    params.set("lattice_type", "kagome");
    params.set("J", -1.0); // antiferromagnetic
    params.set("beta", 5.0);

    let config = RunConfig {
        thermalization_sweeps: 500,
        measurement_sweeps: 500,
        binsize: 100,
        base_seed: 42,
        ..Default::default()
    };

    let backend = RayonBackend::new(1);
    let scheduler = Scheduler::new(backend, config);
    let results = scheduler.run_one::<ClassicalMC<IsingModel, MetropolisCore>>(&params);

    let e = results.get("Energy").expect("Energy missing");
    let m = results.get("Magnetization").expect("Magnetization missing");

    // AF kagome: energy should be negative, magnetization near 0 (no FM order)
    assert!(e.mean < 0.0, "AF energy should be negative, got {}", e.mean);
    assert!(
        m.mean < 0.5,
        "AF magnetization should be low, got {}",
        m.mean
    );
}
