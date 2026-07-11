use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use qmc_rs::LatticeSpinQmc;

fn main() {
    let mut params = Params::new();
    params.set("beta", 8.0);
    params.set("model", "heisenberg");
    params.set("topology", "square");
    params.set("Lx", 4);
    params.set("Ly", 4);
    params.set("pbc", true);
    params.set("spin", 0.5);
    params.set("J", 1.0);
    params.set("gauge", "auto");
    params.set("adaptive_schedule", true);

    let run = RunConfig {
        thermalization_sweeps: 2_000,
        measurement_sweeps: 10_000,
        binsize: 100,
        base_seed: 2026,
        ..Default::default()
    };
    let results = Scheduler::new(RayonBackend::new(1), run).run_one::<LatticeSpinQmc>(&params);
    for (name, estimate) in results.estimates() {
        println!("{name}: {}", estimate.format());
    }
}
