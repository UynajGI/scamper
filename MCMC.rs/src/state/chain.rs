use serde::{Deserialize, Serialize};

/// Accepted Markov-chain state with synchronized target-density cache.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChainState<S, C = ()> {
    position: S,
    log_density: f64,
    cache: C,
    iteration: u64,
}

impl<S, C> ChainState<S, C> {
    pub const fn new(position: S, log_density: f64, cache: C) -> Self {
        Self {
            position,
            log_density,
            cache,
            iteration: 0,
        }
    }

    pub const fn position(&self) -> &S {
        &self.position
    }

    #[allow(dead_code)]
    pub(crate) fn position_mut_for_cache_rebuild(&mut self) -> &mut S {
        &mut self.position
    }

    pub const fn log_density(&self) -> f64 {
        self.log_density
    }

    pub const fn iteration(&self) -> u64 {
        self.iteration
    }

    pub const fn cache(&self) -> &C {
        &self.cache
    }

    pub fn cache_mut(&mut self) -> &mut C {
        &mut self.cache
    }

    pub fn replace(&mut self, position: S, log_density: f64) {
        self.position = position;
        self.log_density = log_density;
        self.iteration = self.iteration.saturating_add(1);
    }

    pub fn swap_position(&mut self, proposal: &mut S, proposed_log_density: f64) {
        std::mem::swap(&mut self.position, proposal);
        self.log_density = proposed_log_density;
        self.iteration = self.iteration.saturating_add(1);
    }

    pub fn mark_rejected_transition(&mut self) {
        self.iteration = self.iteration.saturating_add(1);
    }
}
