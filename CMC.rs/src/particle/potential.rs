//! Pair-potential abstraction and Lennard-Jones implementations.

use crate::particle::ParticleError;

/// Isotropic pair potential with a finite radial cutoff.
pub trait PairPotential: Send + Sync {
    /// Squared radial cutoff.
    fn cutoff_squared(&self) -> f64;

    /// Pair energy for two species at squared separation.
    fn energy(&self, species_i: u16, species_j: u16, distance_squared: f64) -> f64;

    /// Whether this potential defines parameters for a species label.
    fn supports_species(&self, _species: u16) -> bool {
        true
    }
}

/// Lennard-Jones parameters for one species.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LennardJonesSpecies {
    /// Length scale σ.
    pub sigma: f64,
    /// Well depth ε.
    pub epsilon: f64,
}

impl LennardJonesSpecies {
    /// Validate one species parameter set.
    pub fn new(sigma: f64, epsilon: f64) -> Result<Self, ParticleError> {
        if !sigma.is_finite() || sigma <= 0.0 {
            return Err(ParticleError::InvalidPotential(
                "sigma must be finite and positive".to_string(),
            ));
        }
        if !epsilon.is_finite() || epsilon < 0.0 {
            return Err(ParticleError::InvalidPotential(
                "epsilon must be finite and non-negative".to_string(),
            ));
        }
        Ok(Self { sigma, epsilon })
    }
}

/// Treatment of the Lennard-Jones potential at the cutoff.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CutoffTreatment {
    /// Set the raw potential to zero beyond the cutoff.
    Truncated,
    /// Subtract `U(r_c)` so the energy is continuous.
    ShiftedPotential,
    /// Also subtract the cutoff slope so energy and force are continuous.
    ShiftedForce,
}

#[derive(Debug, Clone, Copy)]
struct PairParameters {
    sigma_squared: f64,
    four_epsilon: f64,
    energy_at_cutoff: f64,
    derivative_at_cutoff: f64,
}

/// Lennard-Jones 12-6 potential with Lorentz-Berthelot mixing.
#[derive(Debug, Clone)]
pub struct LennardJones {
    species: Vec<LennardJonesSpecies>,
    pairs: Vec<PairParameters>,
    cutoff: f64,
    cutoff_squared: f64,
    treatment: CutoffTreatment,
}

impl LennardJones {
    /// Monatomic Lennard-Jones potential with shifted-potential cutoff.
    pub fn new(sigma: f64, epsilon: f64, cutoff: f64) -> Result<Self, ParticleError> {
        Self::with_treatment(sigma, epsilon, cutoff, CutoffTreatment::ShiftedPotential)
    }

    /// Monatomic Lennard-Jones potential with explicit cutoff treatment.
    pub fn with_treatment(
        sigma: f64,
        epsilon: f64,
        cutoff: f64,
        treatment: CutoffTreatment,
    ) -> Result<Self, ParticleError> {
        let species = LennardJonesSpecies::new(sigma, epsilon)?;
        Self::from_species(vec![species], cutoff, treatment)
    }

    /// Multi-species potential using Lorentz-Berthelot mixed pair parameters.
    pub fn from_species(
        species: Vec<LennardJonesSpecies>,
        cutoff: f64,
        treatment: CutoffTreatment,
    ) -> Result<Self, ParticleError> {
        if species.is_empty() {
            return Err(ParticleError::InvalidPotential(
                "at least one species is required".to_string(),
            ));
        }
        for parameter in &species {
            LennardJonesSpecies::new(parameter.sigma, parameter.epsilon)?;
        }
        if !cutoff.is_finite() || cutoff <= 0.0 {
            return Err(ParticleError::InvalidPotential(
                "cutoff must be finite and positive".to_string(),
            ));
        }

        let count = species.len();
        if count > usize::from(u16::MAX) + 1 {
            return Err(ParticleError::InvalidPotential(
                "species table exceeds u16 label capacity".to_string(),
            ));
        }
        let cutoff_squared = cutoff * cutoff;
        if !cutoff_squared.is_finite() {
            return Err(ParticleError::InvalidPotential(
                "squared cutoff must remain finite".to_string(),
            ));
        }
        let pair_count = count.checked_mul(count).ok_or_else(|| {
            ParticleError::InvalidPotential("pair-parameter table size overflow".to_string())
        })?;
        let mut pairs = Vec::with_capacity(pair_count);
        for left in &species {
            for right in &species {
                let sigma = 0.5 * (left.sigma + right.sigma);
                let epsilon = left.epsilon.sqrt() * right.epsilon.sqrt();
                let sigma_squared = sigma * sigma;
                let four_epsilon = 4.0 * epsilon;
                if !sigma_squared.is_finite() || !four_epsilon.is_finite() {
                    return Err(ParticleError::InvalidPotential(
                        "mixed Lennard-Jones parameters overflow".to_string(),
                    ));
                }
                let reduced_squared = sigma_squared / cutoff_squared;
                let reduced_sixth = reduced_squared * reduced_squared * reduced_squared;
                let reduced_twelfth = reduced_sixth * reduced_sixth;
                let energy_at_cutoff = four_epsilon * (reduced_twelfth - reduced_sixth);
                let derivative_at_cutoff =
                    24.0 * epsilon * (reduced_sixth - 2.0 * reduced_twelfth) / cutoff;
                if !energy_at_cutoff.is_finite() || !derivative_at_cutoff.is_finite() {
                    return Err(ParticleError::InvalidPotential(
                        "Lennard-Jones cutoff shift is non-finite".to_string(),
                    ));
                }
                pairs.push(PairParameters {
                    sigma_squared,
                    four_epsilon,
                    energy_at_cutoff,
                    derivative_at_cutoff,
                });
            }
        }

        Ok(Self {
            species,
            pairs,
            cutoff,
            cutoff_squared,
            treatment,
        })
    }

    /// Radial cutoff.
    #[inline]
    pub const fn cutoff(&self) -> f64 {
        self.cutoff
    }

    /// Selected cutoff treatment.
    #[inline]
    pub const fn treatment(&self) -> CutoffTreatment {
        self.treatment
    }

    /// Number of configured species.
    #[inline]
    pub fn species_count(&self) -> usize {
        self.species.len()
    }

    #[inline]
    fn pair(&self, species_i: u16, species_j: u16) -> PairParameters {
        let count = self.species.len();
        let left = usize::from(species_i);
        let right = usize::from(species_j);
        assert!(
            left < count && right < count,
            "Lennard-Jones species index out of range"
        );
        self.pairs[left * count + right]
    }
}

impl PairPotential for LennardJones {
    #[inline]
    fn cutoff_squared(&self) -> f64 {
        self.cutoff_squared
    }

    #[inline]
    fn energy(&self, species_i: u16, species_j: u16, distance_squared: f64) -> f64 {
        if !distance_squared.is_finite() || distance_squared <= 0.0 {
            return f64::INFINITY;
        }
        if distance_squared >= self.cutoff_squared {
            return 0.0;
        }

        let pair = self.pair(species_i, species_j);
        let reduced_squared = pair.sigma_squared / distance_squared;
        let reduced_sixth = reduced_squared * reduced_squared * reduced_squared;
        let reduced_twelfth = reduced_sixth * reduced_sixth;
        let raw = pair.four_epsilon * (reduced_twelfth - reduced_sixth);
        match self.treatment {
            CutoffTreatment::Truncated => raw,
            CutoffTreatment::ShiftedPotential => raw - pair.energy_at_cutoff,
            CutoffTreatment::ShiftedForce => {
                let distance = distance_squared.sqrt();
                raw - pair.energy_at_cutoff - (distance - self.cutoff) * pair.derivative_at_cutoff
            }
        }
    }

    #[inline]
    fn supports_species(&self, species: u16) -> bool {
        usize::from(species) < self.species.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn raw_lj_has_expected_minimum() {
        let potential =
            LennardJones::with_treatment(1.0, 2.0, 3.0, CutoffTreatment::Truncated).unwrap();
        let r = 2.0f64.powf(1.0 / 6.0);
        assert!((potential.energy(0, 0, r * r) + 2.0).abs() < 1e-12);
    }

    #[test]
    fn shifted_variants_vanish_at_cutoff_from_below() {
        for treatment in [
            CutoffTreatment::ShiftedPotential,
            CutoffTreatment::ShiftedForce,
        ] {
            let potential = LennardJones::with_treatment(1.0, 1.0, 2.5, treatment).unwrap();
            let r = 2.5 - 1e-9;
            assert!(potential.energy(0, 0, r * r).abs() < 1e-8);
            assert_eq!(potential.energy(0, 0, 2.5f64.powi(2)), 0.0);
        }
    }

    #[test]
    fn shifted_force_has_zero_cutoff_slope() {
        let potential =
            LennardJones::with_treatment(1.0, 1.0, 2.5, CutoffTreatment::ShiftedForce).unwrap();
        let h = 1e-5_f64;
        let left = potential.energy(0, 0, (2.5 - h).powi(2));
        let middle = potential.energy(0, 0, (2.5 - 2.0 * h).powi(2));
        let slope = (left - middle) / h;
        assert!(slope.abs() < 2e-4, "cutoff slope was {slope}");
    }

    #[test]
    fn lorentz_berthelot_mixing_matches_pair_minimum() {
        let potential = LennardJones::from_species(
            vec![
                LennardJonesSpecies::new(1.0, 1.0).unwrap(),
                LennardJonesSpecies::new(3.0, 4.0).unwrap(),
            ],
            6.0,
            CutoffTreatment::Truncated,
        )
        .unwrap();
        let mixed_sigma = 2.0;
        let mixed_epsilon = 2.0;
        let minimum = mixed_sigma * 2.0f64.powf(1.0 / 6.0);
        assert!((potential.energy(0, 1, minimum * minimum) + mixed_epsilon).abs() < 1e-12);
    }

    #[test]
    fn exact_overlap_is_an_infinite_barrier() {
        let potential = LennardJones::new(1.0, 1.0, 2.5).unwrap();
        assert_eq!(potential.energy(0, 0, 0.0), f64::INFINITY);
        assert_eq!(potential.energy(0, 0, f64::NAN), f64::INFINITY);
    }
}
