use carlo_rs::{
    AdaptiveRunControl, CarloError, Context, FromParams, MonteCarlo, Params, RayonBackend,
    RunConfig, RunDecision, RunPhase, Scheduler,
};
use rand_xoshiro::Xoshiro256PlusPlus;

struct ControlledMC;

impl MonteCarlo for ControlledMC {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        assert!(matches!(
            ctx.phase(),
            RunPhase::Thermalization | RunPhase::Measurement
        ));
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        ctx.measure("Value", 1.0);
    }
}

impl FromParams for ControlledMC {
    fn from_params(_params: &Params, _rng: &mut Self::Rng) -> Result<Self, CarloError> {
        Ok(Self)
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
            // This fixed warmup setting is intentionally shorter than the
            // controller's convergence condition. Explicit phase must win.
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
