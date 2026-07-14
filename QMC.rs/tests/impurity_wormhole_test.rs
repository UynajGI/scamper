use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use qmc_rs::impurity::ImpurityQmc;

fn run_model(model: &str) -> carlo_rs::Results {
    let mut params = Params::new();
    params.set("beta", 4.0);
    params.set("model", model);
    params.set("bath", "single");
    params.set("omega0", 1.0);
    params.set("g", 0.35);
    params.set("g_xy", 0.35);
    params.set("g_x", 0.35);
    params.set("g_y", 0.20);
    params.set("g_z", 0.10);
    params.set("h_z", 0.2);
    params.set("crw_ratio", 0.2);
    params.set("tunnelling", 0.2);
    params.set("validate_each_sweep", true);
    let run = RunConfig {
        thermalization_sweeps: 300,
        measurement_sweeps: 600,
        binsize: 30,
        base_seed: 2026,
        ..Default::default()
    };
    Scheduler::new(RayonBackend::new(1), run).run_one::<ImpurityQmc>(&params)
}

#[test]
fn all_impurity_catalogs_run_through_carlo() {
    for model in ["jc", "rw_crw", "xxz", "xyz", "rabi"] {
        let results = run_model(model);
        let order = results
            .get("ExpansionOrder")
            .unwrap_or_else(|| panic!("missing ExpansionOrder for {model}"));
        assert!(order.mean >= 0.0);
        assert!(results.get("WormholeFraction").is_some());
        assert!(results.get("LoopAbortFraction").is_some());
    }
}

#[test]
fn free_limit_has_bounded_magnetization() {
    let mut params = Params::new();
    params.set("beta", 2.0);
    params.set("model", "xxz");
    params.set("bath", "single");
    params.set("omega0", 1.0);
    params.set("lambda_xy", 0.0);
    params.set("lambda_z", 0.0);
    params.set("h_z", 0.0);
    let run = RunConfig {
        thermalization_sweeps: 100,
        measurement_sweeps: 500,
        binsize: 25,
        base_seed: 77,
        ..Default::default()
    };
    let results = Scheduler::new(RayonBackend::new(1), run).run_one::<ImpurityQmc>(&params);
    let magnetization = results
        .get("MagnetizationSigmaZ")
        .expect("MagnetizationSigmaZ missing");
    assert!(magnetization.mean.abs() < 0.4);
}
