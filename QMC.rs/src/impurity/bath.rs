//! Normalized bath samplers for retarded vertices.
//!
//! The diagonal update samples the factorized probability density
//!
//! `calJ(omega) * P(omega, delta_tau)`,
//!
//! where `calJ` is proportional to `J(omega)/omega` and
//! `P(omega,tau) = omega D(omega,tau)`.  Because this density is included in
//! the proposal, it cancels from the Metropolis acceptance ratio.

use rand::Rng;
use rand::RngExt;

use super::error::ImpurityError;

/// Whether a retarded operator uses the directed propagator `D` or the
/// symmetrized propagator `D_+`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KernelDirection {
    /// Keep the orientation of imaginary-time propagation.  Used by JC.
    Directed,
    /// Sample the two orientations with equal probability.  Used by Hermitian
    /// coordinate couplings such as XXZ, XYZ, and rotated impurity models.
    Symmetric,
}

/// A sampled bath frequency and directed time difference.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BathSample {
    /// Positive bosonic frequency.
    pub omega: f64,
    /// Directed time difference in `[0, beta)`.
    pub delta_tau: f64,
}

/// Discrete single-frequency bath.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SingleModeBath {
    omega: f64,
}

impl SingleModeBath {
    /// Construct a single mode of positive frequency.
    pub fn new(omega: f64) -> Result<Self, ImpurityError> {
        if !omega.is_finite() || omega <= 0.0 {
            return Err(ImpurityError::parameter(
                "omega",
                format!("must be finite and positive, got {omega}"),
            ));
        }
        Ok(Self { omega })
    }

    /// Mode frequency.
    pub fn omega(&self) -> f64 {
        self.omega
    }
}

/// Power-law bath shape with `J(omega) proportional to omega^s` on
/// `(0, omega_c)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PowerLawBath {
    exponent: f64,
    cutoff: f64,
}

impl PowerLawBath {
    /// Construct a normalized power-law proposal.
    pub fn new(exponent: f64, cutoff: f64) -> Result<Self, ImpurityError> {
        if !exponent.is_finite() || exponent <= 0.0 {
            return Err(ImpurityError::parameter(
                "s",
                format!("must be finite and positive, got {exponent}"),
            ));
        }
        if !cutoff.is_finite() || cutoff <= 0.0 {
            return Err(ImpurityError::parameter(
                "omega_c",
                format!("must be finite and positive, got {cutoff}"),
            ));
        }
        Ok(Self { exponent, cutoff })
    }

    /// Bath exponent `s`.
    pub fn exponent(&self) -> f64 {
        self.exponent
    }

    /// Ultraviolet cutoff.
    pub fn cutoff(&self) -> f64 {
        self.cutoff
    }
}

/// Positive discrete spectral measure.
///
/// `weights[k]` is a mass proportional to `J(omega_k)/omega_k`; the
/// constructor normalizes the masses internally.  This representation covers
/// arbitrary finite multimode baths and quadrature discretizations of a
/// continuous spectrum.
#[derive(Debug, Clone, PartialEq)]
pub struct TabulatedBath {
    frequencies: Vec<f64>,
    cumulative: Vec<f64>,
}

impl TabulatedBath {
    /// Construct a tabulated bath from frequencies and positive masses.
    pub fn new(frequencies: Vec<f64>, weights: Vec<f64>) -> Result<Self, ImpurityError> {
        if frequencies.is_empty() {
            return Err(ImpurityError::InvalidBathTable(
                "at least one frequency is required".into(),
            ));
        }
        if frequencies.len() != weights.len() {
            return Err(ImpurityError::InvalidBathTable(format!(
                "frequency count {} differs from weight count {}",
                frequencies.len(),
                weights.len()
            )));
        }
        if frequencies
            .iter()
            .any(|omega| !omega.is_finite() || *omega <= 0.0)
        {
            return Err(ImpurityError::InvalidBathTable(
                "all frequencies must be finite and positive".into(),
            ));
        }
        if weights
            .iter()
            .any(|weight| !weight.is_finite() || *weight < 0.0)
        {
            return Err(ImpurityError::InvalidBathTable(
                "all weights must be finite and non-negative".into(),
            ));
        }
        let total: f64 = weights.iter().sum();
        if total <= 0.0 || !total.is_finite() {
            return Err(ImpurityError::InvalidBathTable(
                "the total spectral mass must be positive and finite".into(),
            ));
        }
        let mut running = 0.0;
        let mut cumulative = Vec::with_capacity(weights.len());
        for weight in weights {
            running += weight / total;
            cumulative.push(running.min(1.0));
        }
        if let Some(last) = cumulative.last_mut() {
            *last = 1.0;
        }
        Ok(Self {
            frequencies,
            cumulative,
        })
    }

    /// Frequencies in the table.
    pub fn frequencies(&self) -> &[f64] {
        &self.frequencies
    }

    /// Normalized cumulative masses.
    pub fn cumulative(&self) -> &[f64] {
        &self.cumulative
    }
}

/// Bath shapes supported by the generic impurity engine.
#[derive(Debug, Clone, PartialEq)]
pub enum Bath {
    /// One oscillator.
    SingleMode(SingleModeBath),
    /// Sharp-cutoff power law.
    PowerLaw(PowerLawBath),
    /// Arbitrary positive discrete spectral measure.
    Tabulated(TabulatedBath),
}

impl Bath {
    /// Sample `(omega, delta_tau)` from the normalized retarded proposal.
    pub fn sample<R: Rng + ?Sized>(
        &self,
        beta: f64,
        direction: KernelDirection,
        rng: &mut R,
    ) -> Result<BathSample, ImpurityError> {
        if !beta.is_finite() || beta <= 0.0 {
            return Err(ImpurityError::parameter(
                "beta",
                format!("must be finite and positive, got {beta}"),
            ));
        }

        let omega = match self {
            Self::SingleMode(bath) => bath.omega,
            Self::PowerLaw(bath) => {
                let u = positive_uniform(rng);
                bath.cutoff * u.powf(1.0 / bath.exponent)
            }
            Self::Tabulated(bath) => {
                let u = rng.random::<f64>();
                let index = bath.cumulative.partition_point(|value| *value < u);
                bath.frequencies[index.min(bath.frequencies.len() - 1)]
            }
        };

        let u = positive_uniform(rng);
        let one_minus_exp = -(-beta * omega).exp_m1();
        let mut delta_tau = -(-u * one_minus_exp).ln_1p() / omega;
        if direction == KernelDirection::Symmetric && rng.random::<bool>() {
            delta_tau = (beta - delta_tau) % beta;
        }
        Ok(BathSample { omega, delta_tau })
    }
}

fn positive_uniform<R: Rng + ?Sized>(rng: &mut R) -> f64 {
    let value = rng.random::<f64>();
    value.max(f64::MIN_POSITIVE)
}

#[cfg(test)]
mod tests {
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256PlusPlus;

    use super::*;

    #[test]
    fn single_mode_samples_valid_time() {
        let bath = Bath::SingleMode(SingleModeBath::new(2.0).expect("valid mode"));
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(7);
        for _ in 0..1_000 {
            let sample = bath
                .sample(5.0, KernelDirection::Directed, &mut rng)
                .expect("sample");
            assert!((sample.omega - 2.0).abs() < f64::EPSILON);
            assert!((0.0..5.0).contains(&sample.delta_tau));
        }
    }

    #[test]
    fn power_law_samples_below_cutoff() {
        let bath = Bath::PowerLaw(PowerLawBath::new(0.8, 3.0).expect("valid bath"));
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(11);
        for _ in 0..1_000 {
            let sample = bath
                .sample(4.0, KernelDirection::Symmetric, &mut rng)
                .expect("sample");
            assert!(sample.omega > 0.0 && sample.omega <= 3.0);
            assert!((0.0..4.0).contains(&sample.delta_tau));
        }
    }

    #[test]
    fn tabulated_validates_lengths() {
        let error = TabulatedBath::new(vec![1.0], vec![1.0, 2.0]);
        assert!(error.is_err());
    }
}
