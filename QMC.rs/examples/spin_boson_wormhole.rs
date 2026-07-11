use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use qmc_rs::SpinBosonQmc;

fn main() {
    let mut params = Params::new();
    params.set("model", "jc");
    params.set("bath", "single");
    params.set("beta", 8.0);
    params.set("omega0", 1.0);
    params.set("g", 0.35);
    params.set("h_z", 0.4);
    params.set("adaptive_schedule", true);

    let run = RunConfig {
        thermalization_sweeps: 5_000,
        measurement_sweeps: 20_000,
        binsize: 100,
        base_seed: 2_026,
        ..Default::default()
    };
    let results = Scheduler::new(RayonBackend::new(1), run).run_one::<SpinBosonQmc>(&params);

    for name in [
        "MagnetizationSz",
        "ChiZ",
        "ExpansionOrder",
        "DiagonalAcceptance",
        "WormholeFraction",
    ] {
        if let Some(estimate) = results.get(name) {
            println!("{name}: {} +/- {}", estimate.mean, estimate.stderr);
        }
    }
}
