//! Slow stochastic checks. These are secondary regression tests, not the
//! primary proof of correctness. Every target comes from exact enumeration.

use super::common::exact_ising_moments;
use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use cmc_rs::{build_chain, ClassicalMC, IsingModel, MetropolisCore};

#[test]
#[ignore = "long statistical regression against exact finite-system moments"]
fn metropolis_energy_matches_exact_four_site_distribution_with_error_budget() {
    let beta = 0.61;
    let lattice = build_chain(4, true);
    let (_, exact_energy, _, _) = exact_ising_moments(&lattice, 1.0, beta);
    let mut params = Params::new();
    params.set("lattice_type", "chain");
    params.set("L", 4usize);
    params.set("pbc", true);
    params.set("beta", beta);
    params.set("J", 1.0);
    let config = RunConfig {
        thermalization_sweeps: 20_000,
        measurement_sweeps: 400_000,
        binsize: 200,
        base_seed: 0x5059_5349_4353,
        ..Default::default()
    };
    let results = Scheduler::new(RayonBackend::new(1), config)
        .run_one::<ClassicalMC<IsingModel, MetropolisCore>>(&params);
    let estimate = results.get("Energy").unwrap();
    let budget = (6.0 * estimate.stderr).max(0.02);
    assert!(
        (estimate.mean - exact_energy).abs() <= budget,
        "measured={} exact={} stderr={} budget={}",
        estimate.mean,
        exact_energy,
        estimate.stderr,
        budget
    );
}
