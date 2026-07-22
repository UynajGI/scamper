use carlo_rs::{Context, MonteCarlo};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

struct DummyMC {
    sweep_count: u64,
}

impl MonteCarlo for DummyMC {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Xoshiro256PlusPlus>) {
        self.sweep_count += 1;
        ctx.measure("sweeps", 1.0);
    }
}

#[test]
fn test_monte_carlo_sweep() {
    let mut mc = DummyMC { sweep_count: 0 };
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut ctx = Context::new(rng, 10);

    mc.sweep(&mut ctx);

    assert_eq!(mc.sweep_count, 1);

    // Verify the measurement was recorded
    let estimates = ctx.finalize_measurements();
    assert!(
        estimates.contains_key("sweeps"),
        "sweeps observable should exist after sweep"
    );
    let est = &estimates["sweeps"];
    assert!(
        (est.mean - 1.0).abs() < 1e-10,
        "expected mean 1.0, got {}",
        est.mean
    );
}
