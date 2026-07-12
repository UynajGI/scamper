use rand::Rng;

use crate::target::LogDensity;
use crate::{
    EuclideanState, McmcError, SamplingPhase, TraceStore, TransitionKernel, TransitionReport,
};

/// Stateful single-chain executor usable without Carlo.rs.
pub struct ChainRunner<T, K, Tr> {
    target: T,
    kernel: K,
    state: EuclideanState,
    trace: Tr,
    chain_id: usize,
    last_report: TransitionReport,
}

impl<T, K, Tr> ChainRunner<T, K, Tr> {
    pub const fn new(
        target: T,
        kernel: K,
        state: EuclideanState,
        trace: Tr,
        chain_id: usize,
    ) -> Self {
        Self {
            target,
            kernel,
            state,
            trace,
            chain_id,
            last_report: TransitionReport {
                accepted: None,
                log_acceptance_probability: None,
                proposals: 0,
                acceptances: 0,
                target_evaluations: 0,
                gradient_evaluations: 0,
                divergent: false,
                energy_error: None,
                leapfrog_steps: 0,
                tree_depth: None,
                proposal_scale: None,
                subtransitions: 0,
            },
        }
    }

    pub const fn state(&self) -> &EuclideanState {
        &self.state
    }

    pub const fn kernel(&self) -> &K {
        &self.kernel
    }

    pub fn kernel_mut(&mut self) -> &mut K {
        &mut self.kernel
    }

    pub const fn trace(&self) -> &Tr {
        &self.trace
    }

    pub fn into_parts(self) -> (T, K, EuclideanState, Tr) {
        (self.target, self.kernel, self.state, self.trace)
    }
}

impl<T, K, Tr> ChainRunner<T, K, Tr>
where
    T: LogDensity<[f64]>,
    K: TransitionKernel<T>,
    Tr: TraceStore,
{
    pub fn start_phase(&mut self, phase: SamplingPhase) -> Result<(), McmcError> {
        self.kernel
            .on_phase_start(&mut self.target, phase, &self.state)
    }

    pub fn end_phase(&mut self, phase: SamplingPhase) -> Result<(), McmcError> {
        self.kernel
            .on_phase_end(&mut self.target, phase, &self.state)
    }

    pub fn step<R>(&mut self, rng: &mut R, phase: SamplingPhase) -> Result<(), McmcError>
    where
        R: Rng + ?Sized,
        T: LogDensity<[f64]>,
    {
        self.last_report = self
            .kernel
            .transition(&mut self.target, &mut self.state, rng, phase)?;
        self.last_report.validate()?;
        if phase == SamplingPhase::Sampling {
            let _retained = self
                .trace
                .record(self.chain_id, &self.state, &self.last_report)?;
        }
        Ok(())
    }

    pub fn run<R>(
        &mut self,
        rng: &mut R,
        phase: SamplingPhase,
        transitions: u64,
    ) -> Result<(), McmcError>
    where
        R: Rng + ?Sized,
        T: LogDensity<[f64]>,
    {
        self.start_phase(phase)?;
        for _ in 0..transitions {
            self.step(rng, phase)?;
        }
        self.end_phase(phase)
    }
}
