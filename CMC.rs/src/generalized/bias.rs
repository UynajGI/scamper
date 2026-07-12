//! Frozen log-weight functions for umbrella and multicanonical production.

use crate::generalized::{GeneralizedError, LogDensityOfStates, MacrostateAxis};

/// Total logarithmic target weight associated with a macrostate bin.
///
/// This is the complete bin-dependent term used by a biased transition kernel.
/// For example, a canonical umbrella includes both `-βE` and the umbrella;
/// a multicanonical bias contains `-ln g(E)`.
pub trait LogBias: Clone + Send + Sync {
    fn bins(&self) -> usize;

    fn log_weight(&self, bin: usize) -> f64;

    fn log_weight_ratio(&self, old_bin: usize, new_bin: usize) -> f64 {
        let old = self.log_weight(old_bin);
        let new = self.log_weight(new_bin);
        if new == f64::NEG_INFINITY {
            f64::NEG_INFINITY
        } else if old == f64::NEG_INFINITY {
            f64::INFINITY
        } else {
            new - old
        }
    }
}

/// Arbitrary fixed log weights supplied in bin order.
#[derive(Debug, Clone, PartialEq)]
pub struct FixedBias {
    log_weights: Vec<f64>,
}

impl FixedBias {
    pub fn new(log_weights: Vec<f64>) -> Result<Self, GeneralizedError> {
        if log_weights.is_empty() {
            return Err(GeneralizedError::new(
                "fixed bias must contain at least one bin",
            ));
        }
        if log_weights
            .iter()
            .any(|value| value.is_nan() || *value == f64::INFINITY)
        {
            return Err(GeneralizedError::new(
                "fixed bias may contain finite values or negative infinity, but not NaN/+infinity",
            ));
        }
        Ok(Self { log_weights })
    }

    #[inline]
    pub fn log_weights(&self) -> &[f64] {
        &self.log_weights
    }
}

impl LogBias for FixedBias {
    #[inline]
    fn bins(&self) -> usize {
        self.log_weights.len()
    }

    #[inline]
    fn log_weight(&self, bin: usize) -> f64 {
        self.log_weights[bin]
    }
}

/// Canonical harmonic umbrella,
/// `ln π(E) = -βE - 1/2 κ(E-E0)^2`.
#[derive(Debug, Clone, PartialEq)]
pub struct HarmonicUmbrellaBias {
    fixed: FixedBias,
    beta: f64,
    center: f64,
    strength: f64,
}

impl HarmonicUmbrellaBias {
    pub fn new<A: MacrostateAxis>(
        axis: &A,
        beta: f64,
        center: f64,
        strength: f64,
    ) -> Result<Self, GeneralizedError> {
        if !beta.is_finite() || beta < 0.0 {
            return Err(GeneralizedError::new(
                "umbrella beta must be finite and non-negative",
            ));
        }
        if !center.is_finite() || !strength.is_finite() || strength < 0.0 {
            return Err(GeneralizedError::new(
                "umbrella center and non-negative strength must be finite",
            ));
        }
        let log_weights = (0..axis.bins())
            .map(|bin| {
                let value = axis.center(bin);
                -beta * value - 0.5 * strength * (value - center).powi(2)
            })
            .collect();
        Ok(Self {
            fixed: FixedBias::new(log_weights)?,
            beta,
            center,
            strength,
        })
    }

    #[inline]
    pub const fn beta(&self) -> f64 {
        self.beta
    }

    #[inline]
    pub const fn center(&self) -> f64 {
        self.center
    }

    #[inline]
    pub const fn strength(&self) -> f64 {
        self.strength
    }
}

impl LogBias for HarmonicUmbrellaBias {
    #[inline]
    fn bins(&self) -> usize {
        self.fixed.bins()
    }

    #[inline]
    fn log_weight(&self, bin: usize) -> f64 {
        self.fixed.log_weight(bin)
    }
}

/// Frozen multicanonical bias `ln π(E) = -ln g(E)`.
#[derive(Debug, Clone, PartialEq)]
pub struct MulticanonicalBias {
    fixed: FixedBias,
}

impl MulticanonicalBias {
    pub fn from_log_density(log_density: &LogDensityOfStates) -> Result<Self, GeneralizedError> {
        let log_weights = log_density
            .values()
            .iter()
            .zip(log_density.visited())
            .map(
                |(&value, &visited)| {
                    if visited {
                        -value
                    } else {
                        f64::NEG_INFINITY
                    }
                },
            )
            .collect();
        Ok(Self {
            fixed: FixedBias::new(log_weights)?,
        })
    }

    #[inline]
    pub fn log_weights(&self) -> &[f64] {
        self.fixed.log_weights()
    }
}

impl LogBias for MulticanonicalBias {
    #[inline]
    fn bins(&self) -> usize {
        self.fixed.bins()
    }

    #[inline]
    fn log_weight(&self, bin: usize) -> f64 {
        self.fixed.log_weight(bin)
    }
}
