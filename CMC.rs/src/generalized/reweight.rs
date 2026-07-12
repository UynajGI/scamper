//! Stable canonical reweighting from a logarithmic density of states.

use crate::generalized::{GeneralizedError, LogDensityOfStates, MacrostateAxis};

/// Canonical thermodynamic estimate reconstructed from `ln g(E)`.
#[derive(Debug, Clone, PartialEq)]
pub struct CanonicalReweighting {
    beta: f64,
    log_partition_function: f64,
    mean_energy: f64,
    mean_energy_squared: f64,
    heat_capacity: f64,
    probabilities: Vec<f64>,
}

impl CanonicalReweighting {
    #[inline]
    pub const fn beta(&self) -> f64 {
        self.beta
    }

    #[inline]
    pub const fn log_partition_function(&self) -> f64 {
        self.log_partition_function
    }

    #[inline]
    pub const fn mean_energy(&self) -> f64 {
        self.mean_energy
    }

    #[inline]
    pub const fn mean_energy_squared(&self) -> f64 {
        self.mean_energy_squared
    }

    #[inline]
    pub const fn heat_capacity(&self) -> f64 {
        self.heat_capacity
    }

    #[inline]
    pub fn probabilities(&self) -> &[f64] {
        &self.probabilities
    }

    /// Reweight a scalar observable tabulated in the same bins as the DOS.
    pub fn mean_observable(&self, values: &[f64]) -> Result<f64, GeneralizedError> {
        if values.len() != self.probabilities.len() {
            return Err(GeneralizedError::new(
                "observable and reweighting buffers have different lengths",
            ));
        }
        if values.iter().any(|value| !value.is_finite()) {
            return Err(GeneralizedError::new(
                "observable buffer contains a non-finite value",
            ));
        }
        Ok(values
            .iter()
            .zip(&self.probabilities)
            .map(|(&value, &probability)| value * probability)
            .sum())
    }
}

/// Canonically reweight a one-dimensional energy DOS at inverse temperature `beta`.
pub fn canonical_reweight<A: MacrostateAxis>(
    axis: &A,
    log_density: &LogDensityOfStates,
    beta: f64,
) -> Result<CanonicalReweighting, GeneralizedError> {
    if !beta.is_finite() || beta < 0.0 {
        return Err(GeneralizedError::new(
            "reweighting beta must be finite and non-negative",
        ));
    }
    if axis.bins() != log_density.bins() {
        return Err(GeneralizedError::new(
            "axis and density of states have different bin counts",
        ));
    }
    let mut log_terms = vec![f64::NEG_INFINITY; axis.bins()];
    let mut maximum = f64::NEG_INFINITY;
    for (bin, term) in log_terms.iter_mut().enumerate() {
        if log_density.is_visited(bin) {
            *term = log_density.value(bin) - beta * axis.center(bin);
            maximum = maximum.max(*term);
        }
    }
    if maximum == f64::NEG_INFINITY {
        return Err(GeneralizedError::new(
            "cannot reweight a density of states with no visited bins",
        ));
    }

    let mut probabilities = vec![0.0; axis.bins()];
    let normalization: f64 = log_terms
        .iter()
        .filter(|term| term.is_finite())
        .map(|term| (*term - maximum).exp())
        .sum();
    if !normalization.is_finite() || normalization <= 0.0 {
        return Err(GeneralizedError::new(
            "canonical reweighting normalization is invalid",
        ));
    }
    for (probability, &term) in probabilities.iter_mut().zip(&log_terms) {
        if term.is_finite() {
            *probability = (term - maximum).exp() / normalization;
        }
    }

    let mean_energy: f64 = probabilities
        .iter()
        .enumerate()
        .map(|(bin, &probability)| probability * axis.center(bin))
        .sum();
    let mean_energy_squared: f64 = probabilities
        .iter()
        .enumerate()
        .map(|(bin, &probability)| probability * axis.center(bin).powi(2))
        .sum();
    let variance = (mean_energy_squared - mean_energy * mean_energy).max(0.0);
    Ok(CanonicalReweighting {
        beta,
        log_partition_function: maximum + normalization.ln(),
        mean_energy,
        mean_energy_squared,
        heat_capacity: beta * beta * variance,
        probabilities,
    })
}
