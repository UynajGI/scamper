use carlo_rs::{
    CarloError, Context, FromParams, MonteCarlo, Params, Run, RunConfig, RunId, TaskId,
};
use rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

struct TestMC {
    value: f64,
}

impl MonteCarlo for TestMC {
    type Rng = Xoshiro256PlusPlus;
    fn sweep(&mut self, _ctx: &mut Context<Self::Rng>) {}
}

impl FromParams for TestMC {
    fn from_params(_params: &Params, _rng: &mut Self::Rng) -> Result<Self, CarloError> {
        Ok(TestMC { value: 0.0 })
    }
}

#[test]
fn test_run_creation() {
    let params = Params::new();
    let config = RunConfig {
        measurement_sweeps: 100,
        thermalization_sweeps: 10,
        binsize: 10,
        base_seed: 42,
        progress_interval: 100,
        checkpoint_interval: 0,
    };
    let run: Run<TestMC, Xoshiro256PlusPlus> =
        Run::new(&params, TaskId::new(0), RunId::new(1), &config, 42)
            .expect("Failed to create run");
    assert_eq!(run.sweep_count(), 0);
    assert_eq!(run.sweeps_done(), 0);
}

#[test]
fn test_run_step() {
    let params = Params::new();
    let config = RunConfig {
        measurement_sweeps: 100,
        thermalization_sweeps: 10,
        binsize: 10,
        base_seed: 42,
        progress_interval: 100,
        checkpoint_interval: 0,
    };
    let mut run: Run<TestMC, Xoshiro256PlusPlus> =
        Run::new(&params, TaskId::new(0), RunId::new(1), &config, 42)
            .expect("Failed to create run");

    // Run until thermalized
    for _ in 0..15 {
        run.step();
    }

    assert!(run.is_thermalized());
    assert!(run.sweeps_done() > 0);
}

#[test]
fn test_run_from_context() {
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let context = Context::new(rng, 100);
    let mc = TestMC { value: 0.0 };
    let run: Run<TestMC, Xoshiro256PlusPlus> = Run::from_context(context, mc);
    assert_eq!(run.sweep_count(), 0);
}

#[test]
fn test_run_task_id() {
    let params = Params::new();
    let config = RunConfig::default();
    let run: Run<TestMC, Xoshiro256PlusPlus> =
        Run::new(&params, TaskId::new(5), RunId::new(3), &config, 42)
            .expect("Failed to create run");
    assert_eq!(run.task_id().as_usize(), 5);
    assert_eq!(run.run_id().as_u64(), 3);
}

#[test]
fn test_run_remaining_sweeps() {
    let params = Params::new();
    let config = RunConfig {
        measurement_sweeps: 100,
        thermalization_sweeps: 10,
        binsize: 10,
        base_seed: 42,
        progress_interval: 100,
        checkpoint_interval: 0,
    };
    let run: Run<TestMC, Xoshiro256PlusPlus> =
        Run::new(&params, TaskId::new(0), RunId::new(1), &config, 42)
            .expect("Failed to create run");
    assert_eq!(run.remaining_sweeps(), 100);
}

#[test]
fn test_run_is_complete() {
    let params = Params::new();
    let config = RunConfig {
        measurement_sweeps: 5,
        thermalization_sweeps: 0,
        binsize: 10,
        base_seed: 42,
        progress_interval: 100,
        checkpoint_interval: 0,
    };
    let mut run: Run<TestMC, Xoshiro256PlusPlus> =
        Run::new(&params, TaskId::new(0), RunId::new(1), &config, 42)
            .expect("Failed to create run");

    assert!(!run.is_complete());

    for _ in 0..5 {
        run.step();
    }

    assert!(run.is_complete());
}

#[test]
fn test_run_run_method() {
    let params = Params::new();
    let config = RunConfig {
        measurement_sweeps: 100,
        thermalization_sweeps: 0,
        binsize: 10,
        base_seed: 42,
        progress_interval: 100,
        checkpoint_interval: 0,
    };
    let mut run: Run<TestMC, Xoshiro256PlusPlus> =
        Run::new(&params, TaskId::new(0), RunId::new(1), &config, 42)
            .expect("Failed to create run");

    run.run(50);
    assert!(run.sweeps_done() >= 50);
}

#[test]
fn test_run_context_access() {
    let params = Params::new();
    let config = RunConfig::default();
    let mut run: Run<TestMC, Xoshiro256PlusPlus> =
        Run::new(&params, TaskId::new(0), RunId::new(1), &config, 42)
            .expect("Failed to create run");

    // Test context access
    let ctx = run.context();
    assert_eq!(ctx.sweep_count(), 0);

    // Test mutable context access
    run.context_mut().advance_sweep();
    assert_eq!(run.context().sweep_count(), 1);
}

#[test]
fn test_run_mc_access() {
    let params = Params::new();
    let config = RunConfig::default();
    let run: Run<TestMC, Xoshiro256PlusPlus> =
        Run::new(&params, TaskId::new(0), RunId::new(1), &config, 42)
            .expect("Failed to create run");

    // Test MC access
    let _mc = run.mc();
}

// ── Additional run lifecycle tests ────────────────────────────────────────

#[test]
fn test_run_finalize() {
    let params = Params::new();
    let config = RunConfig {
        measurement_sweeps: 10,
        thermalization_sweeps: 0,
        binsize: 5,
        base_seed: 42,
        ..Default::default()
    };
    let mut run: Run<TestMC, Xoshiro256PlusPlus> =
        Run::new(&params, TaskId::new(0), RunId::new(1), &config, 42)
            .expect("Failed to create run");

    run.run(10);
    let results = run.finalize(42);
    assert!(results.metadata().measurement_sweeps == 10);
    assert_eq!(results.metadata().base_seed, 42);
}

#[test]
fn test_run_finalize_with_mc() {
    let params = Params::new();
    let config = RunConfig {
        measurement_sweeps: 5,
        thermalization_sweeps: 0,
        binsize: 5,
        ..Default::default()
    };
    let mut run: Run<TestMC, Xoshiro256PlusPlus> =
        Run::new(&params, TaskId::new(0), RunId::new(1), &config, 42)
            .expect("Failed to create run");

    run.run(5);
    let (results, mc) = run.finalize_with_mc(42);
    assert!(results.metadata().measurement_sweeps == 5);
    assert!((mc.value - 0.0).abs() < 1e-10);
}

#[test]
fn test_run_from_parts() {
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let context = Context::new_with_binsize(rng, 0, 10);
    let mc = TestMC { value: 1.0 };
    let config = RunConfig {
        measurement_sweeps: 50,
        thermalization_sweeps: 0,
        binsize: 10,
        ..Default::default()
    };
    let run: Run<TestMC, Xoshiro256PlusPlus> =
        Run::from_parts(context, mc, TaskId::new(2), RunId::new(3), config);
    assert_eq!(run.task_id().as_usize(), 2);
    assert_eq!(run.run_id().as_u64(), 3);
    assert_eq!(run.target_sweeps(), 50);
}

#[test]
fn test_run_target_sweeps() {
    let params = Params::new();
    let config = RunConfig {
        measurement_sweeps: 777,
        thermalization_sweeps: 10,
        binsize: 10,
        ..Default::default()
    };
    let run: Run<TestMC, Xoshiro256PlusPlus> =
        Run::new(&params, TaskId::new(0), RunId::new(0), &config, 42)
            .expect("Failed to create run");
    assert_eq!(run.target_sweeps(), 777);
}

#[test]
fn test_run_mc_mut() {
    let params = Params::new();
    let config = RunConfig::default();
    let mut run: Run<TestMC, Xoshiro256PlusPlus> =
        Run::new(&params, TaskId::new(0), RunId::new(1), &config, 42)
            .expect("Failed to create run");

    run.mc_mut().value = 42.0;
    assert!((run.mc().value - 42.0).abs() < 1e-10);
}

#[test]
fn test_run_timing_observables_recorded() {
    let params = Params::new();
    let config = RunConfig {
        measurement_sweeps: 10,
        thermalization_sweeps: 0,
        binsize: 5,
        ..Default::default()
    };
    let mut run: Run<TestMC, Xoshiro256PlusPlus> =
        Run::new(&params, TaskId::new(0), RunId::new(1), &config, 42)
            .expect("Failed to create run");

    run.run(10);
    let results = run.finalize(42);

    assert!(results.get("_ll_sweep_time").is_some());
    assert!(results.get("_ll_measure_time").is_some());
}

#[test]
fn test_run_zero_measurement_sweeps() {
    let params = Params::new();
    let config = RunConfig {
        measurement_sweeps: 0,
        thermalization_sweeps: 0,
        binsize: 10,
        ..Default::default()
    };
    let run: Run<TestMC, Xoshiro256PlusPlus> =
        Run::new(&params, TaskId::new(0), RunId::new(1), &config, 42)
            .expect("Failed to create run");
    assert!(run.is_complete());
}

#[test]
fn test_runid_taskid() {
    let rid = RunId::new(42);
    assert_eq!(rid.as_u64(), 42);
    let rid2 = RunId::new(42);
    assert_eq!(rid, rid2);

    let tid = TaskId::new(7);
    assert_eq!(tid.as_usize(), 7);
    let tid2 = TaskId::new(7);
    assert_eq!(tid, tid2);
    assert_ne!(tid, TaskId::new(8));
}
