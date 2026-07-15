use carlo_rs::{
    CarloError, Context, FromParams, MonteCarlo, Params, RayonBackend, RunConfig, RunPhase,
    Scheduler,
};
use rand_xoshiro::Xoshiro256PlusPlus;
use std::sync::Mutex;

static EVENTS: Mutex<Vec<(bool, RunPhase)>> = Mutex::new(Vec::new());

struct LifecycleMC;

impl MonteCarlo for LifecycleMC {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        assert!(matches!(
            ctx.phase(),
            RunPhase::Thermalization | RunPhase::Measurement
        ));
    }

    fn on_phase_start(&mut self, phase: RunPhase, _ctx: &mut Context<Self::Rng>) {
        EVENTS
            .lock()
            .expect("events mutex poisoned")
            .push((true, phase));
    }

    fn on_phase_end(&mut self, phase: RunPhase, _ctx: &mut Context<Self::Rng>) {
        EVENTS
            .lock()
            .expect("events mutex poisoned")
            .push((false, phase));
    }
}

impl FromParams for LifecycleMC {
    fn from_params(_params: &Params, _rng: &mut Self::Rng) -> Result<Self, CarloError> {
        Ok(Self)
    }
}

#[test]
fn scheduler_emits_explicit_lifecycle_boundaries() {
    EVENTS.lock().expect("events mutex poisoned").clear();
    let scheduler = Scheduler::new(
        RayonBackend::new(1),
        RunConfig {
            thermalization_sweeps: 2,
            measurement_sweeps: 3,
            binsize: 1,
            progress_interval: 0,
            ..Default::default()
        },
    );
    let _ = scheduler.run_one::<LifecycleMC>(&Params::new());
    let events = EVENTS.lock().expect("events mutex poisoned").clone();
    assert_eq!(
        events,
        vec![
            (true, RunPhase::Thermalization),
            (false, RunPhase::Thermalization),
            (true, RunPhase::Measurement),
            (false, RunPhase::Measurement),
            (true, RunPhase::Finished),
        ]
    );

    EVENTS.lock().expect("events mutex poisoned").clear();
    let scheduler = Scheduler::new(
        RayonBackend::new(1),
        RunConfig {
            thermalization_sweeps: 0,
            measurement_sweeps: 1,
            binsize: 1,
            progress_interval: 0,
            ..Default::default()
        },
    );
    let _ = scheduler.run_one::<LifecycleMC>(&Params::new());
    let events = EVENTS.lock().expect("events mutex poisoned").clone();
    assert_eq!(
        events,
        vec![
            (true, RunPhase::Measurement),
            (false, RunPhase::Measurement),
            (true, RunPhase::Finished),
        ]
    );
}
