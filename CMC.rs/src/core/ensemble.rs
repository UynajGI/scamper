//! Target-distribution abstractions for classical Monte Carlo.
//!
//! Update kernels produce a thermodynamic change.  Ensembles convert that
//! change into a log target-weight ratio.  This keeps proposal mechanics,
//! physical energy evaluation and statistical ensemble policy independent.

/// Changes in thermodynamic extensive variables caused by one trial move.
///
/// Lattice-spin kernels modify only `energy`. Particle kernels additionally use
/// `particle_count`, `volume` and `log_jacobian` for NPT and grand-canonical
/// transitions.
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

    /// Build an isotropic volume-change delta.
    #[inline]
    pub const fn volume(
        delta_energy: f64,
        delta_volume: f64,
        log_coordinate_jacobian: f64,
    ) -> Self {
        Self {
            energy: delta_energy,
            particle_count: 0,
            volume: delta_volume,
            log_jacobian: log_coordinate_jacobian,
        }
    }

    /// Build a particle-number-changing delta.
    #[inline]
    pub const fn particle(delta_energy: f64, delta_particles: i64) -> Self {
        Self {
            energy: delta_energy,
            particle_count: delta_particles,
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
        canonical_energy_term(self.beta, delta.energy) + delta.log_jacobian
    }
}

/// Isothermal-isobaric NPT ensemble,
/// `π ∝ exp[-β(E + PV)]` in scaled coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IsothermalIsobaric {
    beta: f64,
    pressure: f64,
}

impl IsothermalIsobaric {
    pub fn new(beta: f64, pressure: f64) -> Self {
        assert!(
            beta.is_finite() && beta >= 0.0,
            "beta must be finite and non-negative"
        );
        assert!(pressure.is_finite(), "pressure must be finite");
        Self { beta, pressure }
    }

    #[inline]
    pub const fn beta(self) -> f64 {
        self.beta
    }

    #[inline]
    pub const fn pressure(self) -> f64 {
        self.pressure
    }
}

impl Ensemble<ThermodynamicDelta> for IsothermalIsobaric {
    #[inline]
    fn log_weight_ratio(&self, delta: &ThermodynamicDelta) -> f64 {
        canonical_energy_term(self.beta, delta.energy) - self.beta * self.pressure * delta.volume
            + delta.log_jacobian
    }
}

/// Grand-canonical μVT ensemble parameterized by the logarithmic activity.
///
/// For one species, `log_activity = βμ - D ln Λ`. Proposal-density factors
/// such as `V/(N+1)` belong to the Metropolis-Hastings proposal ratio.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GrandCanonical {
    beta: f64,
    log_activity: f64,
}

impl GrandCanonical {
    pub fn new(beta: f64, log_activity: f64) -> Self {
        assert!(
            beta.is_finite() && beta >= 0.0,
            "beta must be finite and non-negative"
        );
        assert!(log_activity.is_finite(), "log activity must be finite");
        Self { beta, log_activity }
    }

    pub fn from_chemical_potential(
        beta: f64,
        chemical_potential: f64,
        dimension: usize,
        thermal_wavelength: f64,
    ) -> Self {
        assert!(
            chemical_potential.is_finite(),
            "chemical potential must be finite"
        );
        assert!(dimension > 0, "dimension must be positive");
        assert!(
            thermal_wavelength.is_finite() && thermal_wavelength > 0.0,
            "thermal wavelength must be finite and positive"
        );
        Self::new(
            beta,
            beta * chemical_potential - dimension as f64 * thermal_wavelength.ln(),
        )
    }

    #[inline]
    pub const fn beta(self) -> f64 {
        self.beta
    }

    #[inline]
    pub const fn log_activity(self) -> f64 {
        self.log_activity
    }
}

impl Ensemble<ThermodynamicDelta> for GrandCanonical {
    #[inline]
    fn log_weight_ratio(&self, delta: &ThermodynamicDelta) -> f64 {
        canonical_energy_term(self.beta, delta.energy)
            + self.log_activity * delta.particle_count as f64
            + delta.log_jacobian
    }
}

#[inline]
fn canonical_energy_term(beta: f64, delta_energy: f64) -> f64 {
    if delta_energy == f64::INFINITY {
        f64::NEG_INFINITY
    } else if delta_energy == f64::NEG_INFINITY {
        f64::INFINITY
    } else {
        -beta * delta_energy
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

    #[test]
    fn npt_includes_pressure_and_coordinate_jacobian() {
        let ensemble = IsothermalIsobaric::new(2.0, 3.0);
        let delta = ThermodynamicDelta::volume(1.5, 0.25, 0.7);
        assert!((ensemble.log_weight_ratio(&delta) + 3.8).abs() < 1e-14);
    }

    #[test]
    fn grand_canonical_includes_activity_once_per_particle() {
        let ensemble = GrandCanonical::new(1.5, 0.8);
        let delta = ThermodynamicDelta::particle(2.0, 1);
        assert!((ensemble.log_weight_ratio(&delta) + 2.2).abs() < 1e-14);
    }
}
