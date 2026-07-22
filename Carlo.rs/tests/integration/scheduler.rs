use carlo_rs::{
    CarloError, Context, FromParams, MonteCarlo, Params, RayonBackend, RunConfig, Scheduler,
};
use rand_xoshiro::Xoshiro256PlusPlus;

/// A simple counting MC for testing.
struct CountingMC {
    sweep_count: u64,
}

impl MonteCarlo for CountingMC {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        self.sweep_count += 1;
        // Record energy = sweep_count (for testing)
        if ctx.is_thermalized() {
            ctx.measure("SweepCount", self.sweep_count as f64);
        }
    }
}

impl FromParams for CountingMC {
    fn from_params(_params: &Params, _rng: &mut Self::Rng) -> Result<Self, CarloError> {
        Ok(Self { sweep_count: 0 })
    }
}

#[test]
fn test_scheduler_single_task() {
    let config = RunConfig {
        thermalization_sweeps: 100,
        measurement_sweeps: 1000,
        binsize: 100,
        base_seed: 42,
        ..Default::default()
    };

    let backend = RayonBackend::default();
    let scheduler = Scheduler::new(backend, config);
    let params = Params::new();

    let results = scheduler.run_one::<CountingMC>(&params);

    // Should have SweepCount observable
    assert!(results.get("SweepCount").is_some());

    // Mean of sweep counts 101..=1100 = (101+1100)/2 = 600.5
    let est = results.get("SweepCount").unwrap();
    assert!(
        (est.mean - 600.5).abs() < 50.0,
        "expected mean ~600.5, got {}",
        est.mean
    );
}

#[test]
fn test_scheduler_parallel_tasks() {
    let config = RunConfig {
        thermalization_sweeps: 10,
        measurement_sweeps: 100,
        binsize: 10,
        base_seed: 12345,
        ..Default::default()
    };

    let backend = RayonBackend::new(2);
    let scheduler = Scheduler::new(backend, config);
    let params = Params::new();

    let results = scheduler.run_parallel::<CountingMC>(4, &params);

    assert_eq!(results.len(), 4);
    for r in &results {
        assert!(r.get("SweepCount").is_some());
    }

    // Mean of sweep counts 11..=110 = (11+110)/2 = 60.5
    let est = results[0].get("SweepCount").unwrap();
    assert!(
        (est.mean - 60.5).abs() < 20.0,
        "expected mean ~60.5, got {}",
        est.mean
    );
}
