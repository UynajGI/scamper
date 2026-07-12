use carlo_rs::{Context, MonteCarlo, RunPhase};
use rand_xoshiro::Xoshiro256PlusPlus;

use crate::target::LogDensity;
use crate::{
    EuclideanState, McmcError, SamplingPhase, TraceStore, TransitionKernel, TransitionReport,
};

/// Adapter that runs an MCMC target/kernel pair inside Carlo.rs lifecycle hooks.
pub struct McmcSampler<T, K, Tr> {
    target: T,
    kernel: K,
    state: EuclideanState,
    trace: Tr,
    chain_id: usize,
    last_report: TransitionReport,
}

impl<T, K, Tr> McmcSampler<T, K, Tr>
where
    T: LogDensity<[f64]>,
{
    pub fn new(
        mut target: T,
        kernel: K,
        initial_position: Vec<f64>,
        trace: Tr,
        chain_id: usize,
    ) -> Result<Self, McmcError> {
        let state = EuclideanState::initialize(&mut target, initial_position)?;
        Ok(Self {
            target,
            kernel,
            state,
            trace,
            chain_id,
            last_report: TransitionReport::default(),
        })
    }

    pub const fn state(&self) -> &EuclideanState {
        &self.state
    }

    pub const fn trace(&self) -> &Tr {
        &self.trace
    }

    pub fn into_trace(self) -> Tr {
        self.trace
    }
}

impl<T, K, Tr> MonteCarlo for McmcSampler<T, K, Tr>
where
    T: LogDensity<[f64]>,
    K: TransitionKernel,
    Tr: TraceStore,
{
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, context: &mut Context<Self::Rng>) {
        let phase = map_phase(context.phase());
        self.last_report = self
            .kernel
            .transition(&mut self.target, &mut self.state, &mut context.rng, phase)
            .unwrap_or_else(|error| panic!("MCMC transition failed: {error}"));
    }

    fn measure(&mut self, context: &mut Context<Self::Rng>) {
        self.trace
            .record(self.chain_id, &self.state, &self.last_report)
            .unwrap_or_else(|error| panic!("MCMC trace write failed: {error}"));
        context.measure("LogDensity", self.state.log_density());
        if let Some(accepted) = self.last_report.accepted {
            context.measure("Acceptance", if accepted { 1.0 } else { 0.0 });
        } else if let Some(rate) = self.last_report.acceptance_rate() {
            context.measure("Acceptance", rate);
        }
        context.measure(
            "TargetEvaluations",
            f64::from(self.last_report.target_evaluations),
        );
        context.measure(
            "Divergent",
            if self.last_report.divergent { 1.0 } else { 0.0 },
        );
        if let Some(scale) = self.last_report.proposal_scale {
            context.measure("ProposalScale", scale);
        }
    }

    fn on_phase_start(&mut self, phase: RunPhase, _context: &mut Context<Self::Rng>) {
        if let Some(sampling_phase) = optional_map_phase(phase) {
            self.kernel
                .on_phase_start(sampling_phase, &self.state)
                .unwrap_or_else(|error| panic!("MCMC phase start failed: {error}"));
        }
    }

    fn on_phase_end(&mut self, phase: RunPhase, _context: &mut Context<Self::Rng>) {
        if let Some(sampling_phase) = optional_map_phase(phase) {
            self.kernel
                .on_phase_end(sampling_phase, &self.state)
                .unwrap_or_else(|error| panic!("MCMC phase end failed: {error}"));
        }
    }

    fn name(&self) -> &'static str {
        self.kernel.name()
    }
}

const fn map_phase(phase: RunPhase) -> SamplingPhase {
    match phase {
        RunPhase::Initialization | RunPhase::Thermalization => SamplingPhase::Warmup,
        RunPhase::Measurement | RunPhase::Finished => SamplingPhase::Sampling,
    }
}

const fn optional_map_phase(phase: RunPhase) -> Option<SamplingPhase> {
    match phase {
        RunPhase::Initialization | RunPhase::Finished => None,
        RunPhase::Thermalization => Some(SamplingPhase::Warmup),
        RunPhase::Measurement => Some(SamplingPhase::Sampling),
    }
}
