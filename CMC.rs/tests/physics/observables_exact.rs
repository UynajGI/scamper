use super::common::{assert_close, exact_ising_moments};
use carlo_rs::{Estimate, Results};
use cmc_rs::{
    binder_cumulant, build_chain, connected_order_parameter_fluctuation, specific_heat,
    susceptibility, zero_field_ising_susceptibility,
};

fn estimate(mean: f64) -> Estimate {
    Estimate {
        mean,
        stderr: 0.0,
        autocorr_time: 1.0,
        n_bins: 1,
    }
}

#[test]
fn thermodynamic_postprocessing_matches_exact_finite_ising_moments() {
    let beta = 0.58;
    let lattice = build_chain(6, true);
    let (_, e, e2, m2) = exact_ising_moments(&lattice, 1.0, beta);
    let mut results = Results::new();
    results.add("Energy", estimate(e));
    results.add("E2", estimate(e2));
    results.add("M2", estimate(m2));

    assert_close(
        specific_heat(&results, beta, lattice.n_sites).unwrap(),
        beta * beta * (e2 - e * e) / lattice.n_sites as f64,
        2e-13,
    );
    assert_close(
        zero_field_ising_susceptibility(&results, beta, lattice.n_sites).unwrap(),
        beta * lattice.n_sites as f64 * m2,
        2e-14,
    );
}

#[test]
fn magnitude_fluctuation_is_not_mislabeled_as_zero_field_ising_response() {
    let mut results = Results::new();
    results.add("Magnetization", estimate(0.6));
    results.add("M2", estimate(0.5));
    let connected = connected_order_parameter_fluctuation(&results, 2.0, 10).unwrap();
    let compatibility = susceptibility(&results, 2.0, 10).unwrap();
    let response = zero_field_ising_susceptibility(&results, 2.0, 10).unwrap();
    assert_close(connected, 2.8, 1e-15);
    assert_eq!(compatibility, connected);
    assert_close(response, 10.0, 1e-15);
    assert_ne!(connected, response);
}

#[test]
fn binder_cumulant_uses_even_moments_exactly() {
    let mut results = Results::new();
    results.add("M2", estimate(0.25));
    results.add("M4", estimate(0.125));
    assert_close(binder_cumulant(&results).unwrap(), 1.0 / 3.0, 1e-15);
}
