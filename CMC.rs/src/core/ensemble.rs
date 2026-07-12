//! Target-distribution abstractions for classical Monte Carlo.
//!
//! Update kernels produce a thermodynamic change.  Ensembles convert that
//! change into a log target-weight ratio.  This keeps proposal mechanics,
//! physical energy evaluation and statistical ensemble policy independent.

/// Changes in thermodynamic extensive variables caused by one trial move.
///
/// Lattice-spin kernels modify only `energy`. Particle kernels also use
/// `log_jacobian`. `particle_count` and `volume` are reserved for NPT/grand
/// canonical extensions.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct ThermodynamicDelta {
    pub energy: f64,
    pub particle_count: i64,
    pub volume: f64,
    /// Coordinate/proposal Jacobian contribution owned by the state backend.
    pub log_jacobian: f64,
}

impl ThermodynamicDelta {
    #[inline]
    pub const fn energy(delta_energy: f64) -> Self {
        Self {
            energy: delta_energy,
            particle_count: 0,
            volume: 0.0,
            log_jacobian: 0.0,
        }
    }

    #[inline]
    pub fn is_finite(self) -> bool {
        self.energy.is_finite() && self.volume.is_finite() && self.log_jacobian.is_finite()
    }
}

/// Converts a state change into `ln π(new) - ln π(old)`.
pub trait Ensemble<D>: Send + Sync {
    fn log_weight_ratio(&self, delta: &D) -> f64;
}

/// Canonical NVT ensemble, `π ∝ exp(-βE)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CanonicalEnsemble {
    beta: f64,
}

impl CanonicalEnsemble {
    pub fn new(beta: f64) -> Self {
        assert!(
            beta.is_finite() && beta >= 0.0,
            "beta must be finite and non-negative"
        );
        Self { beta }
    }

    #[inline]
    pub const fn beta(self) -> f64 {
        self.beta
    }
}

impl Ensemble<ThermodynamicDelta> for CanonicalEnsemble {
    #[inline]
    fn log_weight_ratio(&self, delta: &ThermodynamicDelta) -> f64 {
        let energy_term = if delta.energy == f64::INFINITY {
            f64::NEG_INFINITY
        } else if delta.energy == f64::NEG_INFINITY {
            f64::INFINITY
        } else {
            -self.beta * delta.energy
        };
        energy_term + delta.log_jacobian
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_applies_beta_exactly_once() {
        let ensemble = CanonicalEnsemble::new(2.5);
        let delta = ThermodynamicDelta::energy(1.2);
        assert!((ensemble.log_weight_ratio(&delta) + 3.0).abs() < 1e-14);
    }

    #[test]
    fn canonical_rejects_an_infinite_energy_barrier_even_at_zero_beta() {
        let ensemble = CanonicalEnsemble::new(0.0);
        let delta = ThermodynamicDelta::energy(f64::INFINITY);
        assert_eq!(ensemble.log_weight_ratio(&delta), f64::NEG_INFINITY);
    }
}
