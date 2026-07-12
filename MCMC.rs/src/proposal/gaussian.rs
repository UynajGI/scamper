use rand::{Rng, RngExt};
use serde::{Deserialize, Serialize};

use crate::McmcError;

/// Scale geometry for a Gaussian random-walk proposal.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum GaussianScale {
    Isotropic(f64),
    Diagonal(Vec<f64>),
}

impl GaussianScale {
    pub fn validate(&self, dimension: usize) -> Result<(), McmcError> {
        match self {
            Self::Isotropic(scale) if scale.is_finite() && *scale > 0.0 => Ok(()),
            Self::Isotropic(_) => Err(McmcError::InvalidConfig(
                "isotropic proposal scale must be finite and positive".to_string(),
            )),
            Self::Diagonal(scales) if scales.len() != dimension => {
                Err(McmcError::DimensionMismatch {
                    expected: dimension,
                    actual: scales.len(),
                })
            }
            Self::Diagonal(scales)
                if scales.iter().all(|scale| scale.is_finite() && *scale > 0.0) =>
            {
                Ok(())
            }
            Self::Diagonal(_) => Err(McmcError::InvalidConfig(
                "diagonal proposal scales must be finite and positive".to_string(),
            )),
        }
    }

    pub fn component(&self, index: usize) -> f64 {
        match self {
            Self::Isotropic(scale) => *scale,
            Self::Diagonal(scales) => scales[index],
        }
    }

    pub fn set_diagonal(&mut self, scales: Vec<f64>) {
        *self = Self::Diagonal(scales);
    }
}

/// Draw one standard normal variate using the polar Box-Muller transform.
pub fn standard_normal<R>(rng: &mut R) -> f64
where
    R: Rng + ?Sized,
{
    loop {
        let x = 2.0 * rng.random::<f64>() - 1.0;
        let y = 2.0 * rng.random::<f64>() - 1.0;
        let radius_squared = x.mul_add(x, y * y);
        if (0.0..1.0).contains(&radius_squared) {
            return x * (-2.0 * radius_squared.ln() / radius_squared).sqrt();
        }
    }
}
