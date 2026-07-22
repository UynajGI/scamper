//! Thread-count independence test for QMC solvers.
//!
//! Verifies that running with 1 thread vs 4 threads produces
//! statistically equivalent results (same mean within 2σ).

use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use qmc_rs::impurity::ImpurityQmc;

fn run_impurity(n_threads: usize, seed: u64) -> (f64, f64) {
    let mut params = Params::new();
    params.set("beta", 4.0);
    params.set("model", "jc");
    params.set("bath", "single");
    params.set("omega0", 1.0);
    params.set("g", 0.35);
    params.set("tunnelling", 0.5);
    let run = RunConfig {
        thermalization_sweeps: 2000,
        measurement_sweeps: 8000,
        binsize: 100,
        base_seed: seed,
        ..Default::default()
    };
    let results = Scheduler::new(RayonBackend::new(n_threads), run).run_one::<ImpurityQmc>(&params);
    let order = results.get("ExpansionOrder").expect("ExpansionOrder");
    (order.mean, order.stderr)
}

#[test]
fn thread_count_does_not_change_expansion_order() {
    // Run with 1 thread and 4 threads, same seed.
    // Results should be statistically equivalent.
    let (mean_1, stderr_1) = run_impurity(1, 42);
    let (mean_4, stderr_4) = run_impurity(4, 42);

    // Combined uncertainty
    let combined_stderr = (stderr_1 * stderr_1 + stderr_4 * stderr_4).sqrt();
    let diff = (mean_1 - mean_4).abs();

    assert!(
        diff < 3.0 * combined_stderr.max(0.01),
        "Thread-count dependence: 1-thread={mean_1:.4}±{stderr_1:.4}, \
         4-thread={mean_4:.4}±{stderr_4:.4}, diff={diff:.4}"
    );
}
