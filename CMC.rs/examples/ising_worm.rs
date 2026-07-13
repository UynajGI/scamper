use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use cmc_rs::IsingGraphWormMC;

fn main() {
    let mut params = Params::new();
    params.set("lattice_type", "square");
    params.set("Lx", 8);
    params.set("Ly", 8);
    params.set("pbc", true);
    params.set("beta", 0.44);
    params.set("J", 1.0);
    params.set("worm_updates_per_sweep", 128);
    params.set("worm_track_endpoint_pairs", false);

    let config = RunConfig {
        thermalization_sweeps: 1_000,
        measurement_sweeps: 5_000,
        ..RunConfig::default()
    };
    let _results =
        Scheduler::new(RayonBackend::new(1), config).run_one::<IsingGraphWormMC>(&params);
}
