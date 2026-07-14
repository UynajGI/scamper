//! Generic direct Gillespie/rejection-free event selection.

use super::DynamicsError;
use rand::{Rng, RngExt};

/// Model contract for a finite catalog of continuous-time events.
pub trait RejectionFreeModel: Send {
    type State;
    type Patch: Default;

    fn event_count(&self, state: &Self::State) -> usize;

    fn event_rate(&self, state: &Self::State, event: usize) -> Result<f64, DynamicsError>;

    /// Prepare the selected event without modifying the accepted state.
    fn prepare_event(
        &self,
        state: &Self::State,
        event: usize,
        patch: &mut Self::Patch,
    ) -> Result<(), DynamicsError>;

    /// Commit exactly one selected event.
    fn commit_event(&self, state: &mut Self::State, event: usize, patch: &Self::Patch);

    fn validate_state(&self, state: &Self::State) -> Result<(), DynamicsError>;
}

/// Outcome of one direct-method event selection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GillespieEvent {
    pub event: usize,
    pub delta_time: f64,
    pub total_rate: f64,
}

/// Generic direct Gillespie kernel.  Rates are recomputed before every event.
pub struct GillespieKernel<M>
where
    M: RejectionFreeModel,
{
    model: M,
    state: M::State,
    patch: M::Patch,
    rates: Vec<f64>,
    event_time: f64,
    events: u64,
}

impl<M> GillespieKernel<M>
where
    M: RejectionFreeModel,
{
    pub fn new(model: M, state: M::State) -> Result<Self, DynamicsError> {
        model.validate_state(&state)?;
        Ok(Self {
            model,
            state,
            patch: M::Patch::default(),
            rates: Vec::new(),
            event_time: 0.0,
            events: 0,
        })
    }

    #[inline]
    pub const fn model(&self) -> &M {
        &self.model
    }

    #[inline]
    pub const fn state(&self) -> &M::State {
        &self.state
    }

    #[inline]
    pub const fn event_time(&self) -> f64 {
        self.event_time
    }

    #[inline]
    pub const fn events(&self) -> u64 {
        self.events
    }

    pub fn validate(&self) -> Result<(), DynamicsError> {
        self.model.validate_state(&self.state)?;
        if !self.event_time.is_finite() || self.event_time < 0.0 {
            return Err(DynamicsError::new("Gillespie event time is invalid"));
        }
        Ok(())
    }

    /// Execute one event.  `Ok(None)` denotes an absorbing zero-rate state.
    pub fn step(&mut self, rng: &mut impl Rng) -> Result<Option<GillespieEvent>, DynamicsError> {
        let count = self.model.event_count(&self.state);
        self.rates.clear();
        self.rates
            .reserve(count.saturating_sub(self.rates.capacity()));
        let mut total_rate = 0.0;
        for event in 0..count {
            let rate = self.model.event_rate(&self.state, event)?;
            if !rate.is_finite() || rate < 0.0 {
                return Err(DynamicsError::new(format!(
                    "event {event} has a non-finite or negative rate"
                )));
            }
            self.rates.push(rate);
            total_rate += rate;
        }
        if total_rate == 0.0 {
            return Ok(None);
        }
        if !total_rate.is_finite() {
            return Err(DynamicsError::new("total event rate overflowed"));
        }

        let event = select_weighted(&self.rates, total_rate, rng)?;
        self.model
            .prepare_event(&self.state, event, &mut self.patch)?;
        let delta_time = exponential_wait(total_rate, rng);
        self.model.commit_event(&mut self.state, event, &self.patch);
        self.event_time += delta_time;
        self.events = self.events.saturating_add(1);
        Ok(Some(GillespieEvent {
            event,
            delta_time,
            total_rate,
        }))
    }

    /// Advance exactly `duration` of event time, executing every event inside the window.
    pub fn advance_by(&mut self, duration: f64, rng: &mut impl Rng) -> Result<u64, DynamicsError> {
        if !duration.is_finite() || duration <= 0.0 {
            return Err(DynamicsError::new(
                "Gillespie advance duration must be finite and positive",
            ));
        }
        let target = self.event_time + duration;
        let before = self.events;
        while self.event_time < target {
            let count = self.model.event_count(&self.state);
            self.rates.clear();
            let mut total_rate = 0.0;
            for event in 0..count {
                let rate = self.model.event_rate(&self.state, event)?;
                if !rate.is_finite() || rate < 0.0 {
                    return Err(DynamicsError::new(format!(
                        "event {event} has a non-finite or negative rate"
                    )));
                }
                self.rates.push(rate);
                total_rate += rate;
            }
            if total_rate == 0.0 {
                self.event_time = target;
                break;
            }
            if !total_rate.is_finite() {
                return Err(DynamicsError::new("total event rate overflowed"));
            }
            let event = select_weighted(&self.rates, total_rate, rng)?;
            let delta_time = exponential_wait(total_rate, rng);
            let remaining = target - self.event_time;
            if delta_time > remaining {
                // Exponential waiting times are memoryless.  No event occurs in
                // this observation window; the unchanged catalog may be sampled
                // afresh in the next window without bias.
                self.event_time = target;
                break;
            }
            self.model
                .prepare_event(&self.state, event, &mut self.patch)?;
            self.model.commit_event(&mut self.state, event, &self.patch);
            self.event_time += delta_time;
            self.events = self.events.saturating_add(1);
        }
        Ok(self.events.saturating_sub(before))
    }
}

#[inline]
pub(crate) fn exponential_wait(total_rate: f64, rng: &mut impl Rng) -> f64 {
    let uniform = rng.random::<f64>().max(f64::MIN_POSITIVE);
    -uniform.ln() / total_rate
}

pub(crate) fn select_weighted(
    rates: &[f64],
    total_rate: f64,
    rng: &mut impl Rng,
) -> Result<usize, DynamicsError> {
    if rates.is_empty() || !total_rate.is_finite() || total_rate <= 0.0 {
        return Err(DynamicsError::new(
            "cannot select from an empty event catalog",
        ));
    }
    let mut threshold = rng.random::<f64>() * total_rate;
    let mut fallback = None;
    for (index, &rate) in rates.iter().enumerate() {
        if rate > 0.0 {
            fallback = Some(index);
        }
        if threshold < rate {
            return Ok(index);
        }
        threshold -= rate;
    }
    fallback.ok_or_else(|| DynamicsError::new("event catalog has zero total rate"))
}
