use carlo_rs::{
    AdaptiveRunControl, CarloError, Context, FromParams, MonteCarlo, Params, RayonBackend,
    RunConfig, RunDecision, RunPhase, Scheduler,
};
use rand_xoshiro::Xoshiro256PlusPlus;

struct ControlledMC {
    sweep_count: u64,
}

impl MonteCarlo for ControlledMC {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        assert!(matches!(
            ctx.phase(),
            RunPhase::Thermalization | RunPhase::Measurement
        ));
        self.sweep_count += 1;
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        ctx.measure("Value", 1.0);
    }
}

impl FromParams for ControlledMC {
    fn from_params(_params: &Params, _rng: &mut Self::Rng) -> Result<Self, CarloError> {
        Ok(Self { sweep_count: 0 })
    }
}

struct SweepBudget {
    adaptation_seen: u64,
    production_seen: u64,
}

impl AdaptiveRunControl<ControlledMC> for SweepBudget {
    fn after_sweep(
        &mut self,
        _mc: &ControlledMC,
        ctx: &Context<Xoshiro256PlusPlus>,
    ) -> RunDecision {
        match ctx.phase() {
            RunPhase::Thermalization => {
                self.adaptation_seen += 1;
                if self.adaptation_seen == 3 {
                    RunDecision::BeginProduction
                } else {
                    RunDecision::ContinueAdaptation
                }
            }
            RunPhase::Measurement => {
                self.production_seen += 1;
                if self.production_seen == 4 {
                    RunDecision::Stop
                } else {
                    RunDecision::ContinueProduction
                }
            }
            RunPhase::Initialization | RunPhase::Finished => unreachable!(),
        }
    }
}

#[test]
fn controlled_run_uses_algorithm_driven_phase_lengths() {
    let scheduler = Scheduler::new(
        RayonBackend::new(1),
        RunConfig {
            thermalization_sweeps: 1,
            measurement_sweeps: 999,
            binsize: 1,
            progress_interval: 0,
            ..Default::default()
        },
    );
    let results = scheduler
        .run_controlled::<ControlledMC, _>(
            &Params::new(),
            SweepBudget {
                adaptation_seen: 0,
                production_seen: 0,
            },
        )
        .expect("controlled run should succeed");

    assert_eq!(results.metadata().thermalization_sweeps, 3);
    assert_eq!(results.metadata().measurement_sweeps, 4);
    let value = results.get("Value").expect("measurement should exist");
    assert!((value.mean - 1.0).abs() < 1e-14);
}

// ── run_controlled_with_state ─────────────────────────────────────────────

#[test]
fn run_controlled_with_state_returns_mc_and_results() {
    let scheduler = Scheduler::new(
        RayonBackend::new(1),
        RunConfig {
            thermalization_sweeps: 1,
            measurement_sweeps: 999,
            binsize: 1,
            progress_interval: 0,
            ..Default::default()
        },
    );

    let (mc, results) = scheduler
        .run_controlled_with_state::<ControlledMC, _>(
            &Params::new(),
            SweepBudget {
                adaptation_seen: 0,
                production_seen: 0,
            },
        )
        .expect("controlled run with state should succeed");

    // MC state is preserved — total sweeps = 3 adaptation + 4 production
    assert_eq!(mc.sweep_count, 7);

    assert_eq!(results.metadata().thermalization_sweeps, 3);
    assert_eq!(results.metadata().measurement_sweeps, 4);
    assert!(results.get("Value").is_some());
}

#[test]
fn run_controlled_with_state_immediate_stop() {
    struct StopImmediately;
    impl AdaptiveRunControl<ControlledMC> for StopImmediately {
        fn initial_phase(&self) -> RunPhase {
            RunPhase::Measurement
        }
        fn after_sweep(
            &mut self,
            _mc: &ControlledMC,
            _ctx: &Context<Xoshiro256PlusPlus>,
        ) -> RunDecision {
            RunDecision::Stop
        }
    }

    let scheduler = Scheduler::new(RayonBackend::new(1), RunConfig::default());

    let (mc, results) = scheduler
        .run_controlled_with_state::<ControlledMC, _>(&Params::new(), StopImmediately)
        .expect("immediate stop should succeed");

    assert_eq!(mc.sweep_count, 1);
    assert_eq!(results.metadata().measurement_sweeps, 1);
    assert_eq!(results.metadata().thermalization_sweeps, 0);
}

// ── Error paths ───────────────────────────────────────────────────────────

#[test]
fn run_controlled_rejects_initialization_phase() {
    struct BadInitial;
    impl AdaptiveRunControl<ControlledMC> for BadInitial {
        fn initial_phase(&self) -> RunPhase {
            RunPhase::Initialization
        }
        fn after_sweep(
            &mut self,
            _mc: &ControlledMC,
            _ctx: &Context<Xoshiro256PlusPlus>,
        ) -> RunDecision {
            RunDecision::Stop
        }
    }

    let scheduler = Scheduler::new(RayonBackend::new(1), RunConfig::default());
    let err = scheduler
        .run_controlled::<ControlledMC, _>(&Params::new(), BadInitial)
        .unwrap_err();
    assert!(
        matches!(err, CarloError::InvalidConfig { ref field, .. } if field == "run_control.initial_phase")
    );
}

#[test]
fn run_controlled_rejects_finished_phase() {
    struct BadInitial;
    impl AdaptiveRunControl<ControlledMC> for BadInitial {
        fn initial_phase(&self) -> RunPhase {
            RunPhase::Finished
        }
        fn after_sweep(
            &mut self,
            _mc: &ControlledMC,
            _ctx: &Context<Xoshiro256PlusPlus>,
        ) -> RunDecision {
            RunDecision::Stop
        }
    }

    let scheduler = Scheduler::new(RayonBackend::new(1), RunConfig::default());
    let err = scheduler
        .run_controlled::<ControlledMC, _>(&Params::new(), BadInitial)
        .unwrap_err();
    assert!(
        matches!(err, CarloError::InvalidConfig { ref field, .. } if field == "run_control.initial_phase")
    );
}

#[test]
fn run_controlled_rejects_production_during_thermalization() {
    struct MismatchDecision;
    impl AdaptiveRunControl<ControlledMC> for MismatchDecision {
        fn initial_phase(&self) -> RunPhase {
            RunPhase::Thermalization
        }
        fn after_sweep(
            &mut self,
            _mc: &ControlledMC,
            _ctx: &Context<Xoshiro256PlusPlus>,
        ) -> RunDecision {
            RunDecision::ContinueProduction // Wrong: should be ContinueAdaptation or BeginProduction
        }
    }

    let scheduler = Scheduler::new(RayonBackend::new(1), RunConfig::default());
    let err = scheduler
        .run_controlled::<ControlledMC, _>(&Params::new(), MismatchDecision)
        .unwrap_err();
    assert!(
        matches!(err, CarloError::InvalidConfig { ref field, .. } if field == "run_control.decision")
    );
}

#[test]
fn run_controlled_rejects_adaptation_during_measurement() {
    struct StartInMeasurement;
    impl AdaptiveRunControl<ControlledMC> for StartInMeasurement {
        fn initial_phase(&self) -> RunPhase {
            RunPhase::Measurement
        }
        fn after_sweep(
            &mut self,
            _mc: &ControlledMC,
            _ctx: &Context<Xoshiro256PlusPlus>,
        ) -> RunDecision {
            RunDecision::ContinueAdaptation // Wrong: measurement phase can't adapt
        }
    }

    let scheduler = Scheduler::new(RayonBackend::new(1), RunConfig::default());
    let err = scheduler
        .run_controlled::<ControlledMC, _>(&Params::new(), StartInMeasurement)
        .unwrap_err();
    assert!(
        matches!(err, CarloError::InvalidConfig { ref field, .. } if field == "run_control.decision")
    );
}

#[test]
fn run_controlled_rejects_begin_production_during_measurement() {
    struct StartInMeasurement;
    impl AdaptiveRunControl<ControlledMC> for StartInMeasurement {
        fn initial_phase(&self) -> RunPhase {
            RunPhase::Measurement
        }
        fn after_sweep(
            &mut self,
            _mc: &ControlledMC,
            _ctx: &Context<Xoshiro256PlusPlus>,
        ) -> RunDecision {
            RunDecision::BeginProduction // Already in production
        }
    }

    let scheduler = Scheduler::new(RayonBackend::new(1), RunConfig::default());
    let err = scheduler
        .run_controlled::<ControlledMC, _>(&Params::new(), StartInMeasurement)
        .unwrap_err();
    assert!(
        matches!(err, CarloError::InvalidConfig { ref field, .. } if field == "run_control.decision")
    );
}
