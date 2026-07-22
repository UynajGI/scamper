//! Ergodicity tests for the occupation worldline solver.
//!
//! Verifies that independent runs from different RNG seeds converge to
//! the same thermal expectation values for ⟨σz⟩ and ⟨n⟩, confirming
//! that the occupation-basis sampler is ergodic across the full Hilbert
//! space (spin × boson occupations).

use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use qmc_rs::OccupationWorldlineQmc;

fn run_occupation_rabi(seed: u64) -> (f64, f64, f64, f64) {
    let mut params = Params::new();
    params.set("beta", 4.0);
    params.set("kind", "rabi");
    params.set("spin_splitting", 1.0);
    params.set("g", 0.3);
    params.set("omega0", 1.0);
    params.set("cutoff", 10);

    let run = RunConfig {
        thermalization_sweeps: 2000,
        measurement_sweeps: 10000,
        binsize: 100,
        base_seed: seed,
        ..Default::default()
    };
    let results =
        Scheduler::new(RayonBackend::new(1), run).run_one::<OccupationWorldlineQmc>(&params);
    let sigma_z = results.get("OccupationSigmaZ").expect("OccupationSigmaZ");
    let boson_n = results
        .get("OccupationBosonNumber")
        .expect("OccupationBosonNumber");
    (sigma_z.mean, sigma_z.stderr, boson_n.mean, boson_n.stderr)
}

#[test]
fn occupation_ergodicity_multi_seed_convergence() {
    let seeds = [42u64, 123, 456, 789];
    let results: Vec<(f64, f64, f64, f64)> =
        seeds.iter().map(|&s| run_occupation_rabi(s)).collect();

    // ⟨σz⟩ consistency: max−min < 4 × max(stderr)
    let sz_values: Vec<f64> = results.iter().map(|r| r.0).collect();
    let sz_stderrs: Vec<f64> = results.iter().map(|r| r.1).collect();
    let sz_spread = sz_values.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
        - sz_values.iter().cloned().fold(f64::INFINITY, f64::min);
    let sz_max_stderr = sz_stderrs.iter().cloned().fold(0.0_f64, f64::max);
    assert!(
        sz_spread < 4.0 * sz_max_stderr.max(0.01),
        "⟨σz⟩ spread={sz_spread:.6} exceeds 4σ={:.6} (values: {sz_values:?})",
        4.0 * sz_max_stderr.max(0.01)
    );

    // ⟨n⟩ consistency: max−min < 4 × max(stderr)
    let n_values: Vec<f64> = results.iter().map(|r| r.2).collect();
    let n_stderrs: Vec<f64> = results.iter().map(|r| r.3).collect();
    let n_spread = n_values.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
        - n_values.iter().cloned().fold(f64::INFINITY, f64::min);
    let n_max_stderr = n_stderrs.iter().cloned().fold(0.0_f64, f64::max);
    assert!(
        n_spread < 4.0 * n_max_stderr.max(0.01),
        "⟨n⟩ spread={n_spread:.6} exceeds 4σ={:.6} (values: {n_values:?})",
        4.0 * n_max_stderr.max(0.01)
    );
}
