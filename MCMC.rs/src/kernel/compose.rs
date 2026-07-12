use rand::{Rng, RngExt};
use serde::{Deserialize, Serialize};

use crate::target::LogDensity;
use crate::{EuclideanState, McmcError, SamplingPhase, TransitionKernel, TransitionReport};

/// Apply two kernels sequentially.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Then<A, B> {
    first: A,
    second: B,
}

impl<A, B> Then<A, B> {
    pub const fn new(first: A, second: B) -> Self {
        Self { first, second }
    }

    pub const fn first(&self) -> &A {
        &self.first
    }

    pub const fn second(&self) -> &B {
        &self.second
    }

    pub fn into_inner(self) -> (A, B) {
        (self.first, self.second)
    }
}

impl<T, A, B> TransitionKernel<T> for Then<A, B>
where
    T: LogDensity<[f64]> + ?Sized,
    A: TransitionKernel<T>,
    B: TransitionKernel<T>,
{
    fn transition<R>(
        &mut self,
        target: &mut T,
        state: &mut EuclideanState,
        rng: &mut R,
        phase: SamplingPhase,
    ) -> Result<TransitionReport, McmcError>
    where
        R: Rng + ?Sized,
    {
        let mut report = self.first.transition(target, state, rng, phase)?;
        report.validate()?;
        report.normalize_subtransitions();
        let mut second = self.second.transition(target, state, rng, phase)?;
        second.validate()?;
        second.normalize_subtransitions();
        report.merge(second);
        report.validate()?;
        Ok(report)
    }

    fn on_phase_start(
        &mut self,
        target: &mut T,
        phase: SamplingPhase,
        state: &EuclideanState,
    ) -> Result<(), McmcError> {
        self.first.on_phase_start(target, phase, state)?;
        self.second.on_phase_start(target, phase, state)
    }

    fn on_phase_end(
        &mut self,
        target: &mut T,
        phase: SamplingPhase,
        state: &EuclideanState,
    ) -> Result<(), McmcError> {
        self.first.on_phase_end(target, phase, state)?;
        self.second.on_phase_end(target, phase, state)
    }

    fn name(&self, _target: &T) -> &'static str {
        "Then"
    }
}

/// Apply one kernel a fixed positive number of times per outer transition.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Repeat<K> {
    kernel: K,
    repetitions: usize,
}

impl<K> Repeat<K> {
    pub fn new(kernel: K, repetitions: usize) -> Result<Self, McmcError> {
        if repetitions == 0 {
            return Err(McmcError::InvalidConfig(
                "kernel repetition count must be positive".to_string(),
            ));
        }
        Ok(Self {
            kernel,
            repetitions,
        })
    }

    pub const fn repetitions(&self) -> usize {
        self.repetitions
    }

    pub const fn kernel(&self) -> &K {
        &self.kernel
    }

    pub fn into_inner(self) -> K {
        self.kernel
    }
}

impl<T, K> TransitionKernel<T> for Repeat<K>
where
    T: LogDensity<[f64]> + ?Sized,
    K: TransitionKernel<T>,
{
    fn transition<R>(
        &mut self,
        target: &mut T,
        state: &mut EuclideanState,
        rng: &mut R,
        phase: SamplingPhase,
    ) -> Result<TransitionReport, McmcError>
    where
        R: Rng + ?Sized,
    {
        let mut combined = self.kernel.transition(target, state, rng, phase)?;
        combined.validate()?;
        combined.normalize_subtransitions();
        for _ in 1..self.repetitions {
            let mut report = self.kernel.transition(target, state, rng, phase)?;
            report.validate()?;
            report.normalize_subtransitions();
            combined.merge(report);
        }
        combined.validate()?;
        Ok(combined)
    }

    fn on_phase_start(
        &mut self,
        target: &mut T,
        phase: SamplingPhase,
        state: &EuclideanState,
    ) -> Result<(), McmcError> {
        self.kernel.on_phase_start(target, phase, state)
    }

    fn on_phase_end(
        &mut self,
        target: &mut T,
        phase: SamplingPhase,
        state: &EuclideanState,
    ) -> Result<(), McmcError> {
        self.kernel.on_phase_end(target, phase, state)
    }

    fn name(&self, _target: &T) -> &'static str {
        "Repeat"
    }
}

/// Randomly choose one of two kernels for each transition.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Mixture<A, B> {
    first: A,
    second: B,
    first_probability: f64,
}

impl<A, B> Mixture<A, B> {
    pub fn new(first: A, second: B, first_probability: f64) -> Result<Self, McmcError> {
        if !first_probability.is_finite() || !(0.0..=1.0).contains(&first_probability) {
            return Err(McmcError::InvalidConfig(
                "mixture probability must lie between zero and one".to_string(),
            ));
        }
        Ok(Self {
            first,
            second,
            first_probability,
        })
    }

    pub const fn first_probability(&self) -> f64 {
        self.first_probability
    }

    pub fn into_inner(self) -> (A, B) {
        (self.first, self.second)
    }
}

impl<T, A, B> TransitionKernel<T> for Mixture<A, B>
where
    T: LogDensity<[f64]> + ?Sized,
    A: TransitionKernel<T>,
    B: TransitionKernel<T>,
{
    fn transition<R>(
        &mut self,
        target: &mut T,
        state: &mut EuclideanState,
        rng: &mut R,
        phase: SamplingPhase,
    ) -> Result<TransitionReport, McmcError>
    where
        R: Rng + ?Sized,
    {
        let mut report = if rng.random::<f64>() < self.first_probability {
            self.first.transition(target, state, rng, phase)?
        } else {
            self.second.transition(target, state, rng, phase)?
        };
        report.validate()?;
        report.normalize_subtransitions();
        Ok(report)
    }

    fn on_phase_start(
        &mut self,
        target: &mut T,
        phase: SamplingPhase,
        state: &EuclideanState,
    ) -> Result<(), McmcError> {
        self.first.on_phase_start(target, phase, state)?;
        self.second.on_phase_start(target, phase, state)
    }

    fn on_phase_end(
        &mut self,
        target: &mut T,
        phase: SamplingPhase,
        state: &EuclideanState,
    ) -> Result<(), McmcError> {
        self.first.on_phase_end(target, phase, state)?;
        self.second.on_phase_end(target, phase, state)
    }

    fn name(&self, _target: &T) -> &'static str {
        "Mixture"
    }
}
