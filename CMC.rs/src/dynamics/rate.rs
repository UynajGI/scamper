//! Detailed-balance single-spin kinetic rate laws.

use super::DynamicsError;

/// Continuous-time rate law for an elementary move with physical energy change `delta_energy`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum KineticRateLaw {
    /// Glauber/heat-bath rate `nu / (1 + exp(beta * delta_energy))`.
    Glauber { attempt_frequency: f64 },
    /// Metropolis kinetic rate `nu * min(1, exp(-beta * delta_energy))`.
    Metropolis { attempt_frequency: f64 },
}

impl Default for KineticRateLaw {
    fn default() -> Self {
        Self::Glauber {
            attempt_frequency: 1.0,
        }
    }
}

impl KineticRateLaw {
    pub fn glauber(attempt_frequency: f64) -> Result<Self, DynamicsError> {
        validate_frequency(attempt_frequency)?;
        Ok(Self::Glauber { attempt_frequency })
    }

    pub fn metropolis(attempt_frequency: f64) -> Result<Self, DynamicsError> {
        validate_frequency(attempt_frequency)?;
        Ok(Self::Metropolis { attempt_frequency })
    }

    pub fn validate(self) -> Result<(), DynamicsError> {
        validate_frequency(self.attempt_frequency())
    }

    #[inline]
    pub const fn attempt_frequency(self) -> f64 {
        match self {
            Self::Glauber { attempt_frequency } | Self::Metropolis { attempt_frequency } => {
                attempt_frequency
            }
        }
    }

    /// Evaluate the rate without forming an unstable positive exponential.
    pub fn rate(self, beta: f64, delta_energy: f64) -> Result<f64, DynamicsError> {
        self.validate()?;
        if !beta.is_finite() || beta < 0.0 {
            return Err(DynamicsError::new(
                "kinetic inverse temperature must be finite and non-negative",
            ));
        }
        if !delta_energy.is_finite() {
            return Err(DynamicsError::new(
                "kinetic event energy change must be finite",
            ));
        }
        let log_ratio = -beta * delta_energy;
        let frequency = self.attempt_frequency();
        let rate = match self {
            Self::Metropolis { .. } => {
                if log_ratio >= 0.0 {
                    frequency
                } else {
                    frequency * log_ratio.exp()
                }
            }
            Self::Glauber { .. } => {
                // sigmoid(log_ratio), split by sign to avoid exp overflow.
                let probability = if log_ratio >= 0.0 {
                    1.0 / (1.0 + (-log_ratio).exp())
                } else {
                    let exp_ratio = log_ratio.exp();
                    exp_ratio / (1.0 + exp_ratio)
                };
                frequency * probability
            }
        };
        if rate.is_finite() && rate >= 0.0 {
            Ok(rate)
        } else {
            Err(DynamicsError::new("kinetic rate is non-finite or negative"))
        }
    }
}

fn validate_frequency(value: f64) -> Result<(), DynamicsError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(DynamicsError::new(
            "attempt frequency must be finite and strictly positive",
        ))
    }
}
