use carlo_rs::{
    CarloError, Context, FromParams, MonteCarlo, Params, RayonBackend, RunConfig, Scheduler,
};
use rand_core::Rng;
use rand_xoshiro::Xoshiro256PlusPlus;

/// Simple MC that generates deterministic values from RNG.
struct DeterministicMC;

impl MonteCarlo for DeterministicMC {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        // Generate deterministic value from RNG
        let raw_val = ctx.rng.next_u64();
        let val = (raw_val as f64) / (u64::MAX as f64);
        if ctx.is_thermalized() {
            ctx.measure("RandomValue", val);
        }
    }
}

impl FromParams for DeterministicMC {
    fn from_params(_params: &Params, _rng: &mut Self::Rng) -> Result<Self, CarloError> {
        Ok(Self)
    }
}

#[test]
fn test_deterministic_same_seed_same_result() {
    let config = RunConfig {
        thermalization_sweeps: 10,
        measurement_sweeps: 100,
        binsize: 10,
        base_seed: 12345,
        ..Default::default()
    };

    let backend = RayonBackend::new(1); // Single thread for reproducibility
    let params = Params::new();

    // Run twice with same configuration
    let results1 =
        Scheduler::new(backend.clone(), config.clone()).run_one::<DeterministicMC>(&params);
    let results2 =
        Scheduler::new(backend.clone(), config.clone()).run_one::<DeterministicMC>(&params);

    // Same seed should produce exactly the same mean
    let est1 = results1.get("RandomValue").unwrap();
    let est2 = results2.get("RandomValue").unwrap();

    // The mean should be bit-for-bit identical
    assert!(
        (est1.mean - est2.mean).abs() < 1e-15,
        "Means should be identical: {} vs {}",
        est1.mean,
        est2.mean
    );
}

#[test]
fn test_different_seed_different_result() {
    let backend = RayonBackend::new(1);
    let params = Params::new();

    let config1 = RunConfig {
        base_seed: 111,
        thermalization_sweeps: 10,
        measurement_sweeps: 100,
        binsize: 10,
        ..Default::default()
    };

    let config2 = RunConfig {
        base_seed: 222,
        ..config1.clone()
    };

    let results1 = Scheduler::new(backend.clone(), config1).run_one::<DeterministicMC>(&params);
    let results2 = Scheduler::new(backend.clone(), config2).run_one::<DeterministicMC>(&params);

    let est1 = results1.get("RandomValue").unwrap();
    let est2 = results2.get("RandomValue").unwrap();

    // Different seeds should produce different results
    assert!(
        (est1.mean - est2.mean).abs() > 1e-10,
        "Different seeds should produce different results: {} vs {}",
        est1.mean,
        est2.mean
    );
}
