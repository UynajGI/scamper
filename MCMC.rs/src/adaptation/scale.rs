use serde::{Deserialize, Serialize};

use crate::McmcError;

/// Diminishing Robbins-Monro adaptation of a global log proposal scale.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RobbinsMonroScale {
    log_multiplier: f64,
    target_acceptance: f64,
    iteration: u64,
    frozen: bool,
}

impl RobbinsMonroScale {
    pub fn new(target_acceptance: f64) -> Result<Self, McmcError> {
        if !target_acceptance.is_finite() || !(0.0..1.0).contains(&target_acceptance) {
            return Err(McmcError::InvalidConfig(
                "target acceptance must lie strictly between zero and one".to_string(),
            ));
        }
        Ok(Self {
            log_multiplier: 0.0,
            target_acceptance,
            iteration: 0,
            frozen: false,
        })
    }

    pub fn observe(&mut self, acceptance_probability: f64) -> Result<(), McmcError> {
        if self.frozen {
            return Err(McmcError::AdaptationFrozen);
        }
        self.iteration = self.iteration.saturating_add(1);
        let gain = (self.iteration as f64 + 10.0).powf(-0.6);
        self.log_multiplier += gain * (acceptance_probability - self.target_acceptance);
        self.log_multiplier = self.log_multiplier.clamp(-20.0, 20.0);
        Ok(())
    }

    pub fn multiplier(&self) -> f64 {
        self.log_multiplier.exp()
    }

    pub fn freeze(&mut self) {
        self.frozen = true;
    }

    pub const fn is_frozen(&self) -> bool {
        self.frozen
    }
}
