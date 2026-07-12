use rand::{Rng, RngExt};
use serde::{Deserialize, Serialize};

use crate::adaptation::regularized_cholesky;
use crate::McmcError;

/// Scale geometry for a Gaussian random-walk proposal.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum GaussianScale {
    Isotropic(f64),
    Diagonal(Vec<f64>),
    /// Row-major lower Cholesky factor; entries above the diagonal are zero.
    Dense {
        dimension: usize,
        cholesky: Vec<f64>,
    },
}

impl GaussianScale {
    pub fn dense_cholesky(dimension: usize, cholesky: Vec<f64>) -> Result<Self, McmcError> {
        let scale = Self::Dense {
            dimension,
            cholesky,
        };
        scale.validate(dimension)?;
        Ok(scale)
    }

    pub fn dense_from_covariance(
        dimension: usize,
        covariance: &[f64],
        jitter: f64,
    ) -> Result<Self, McmcError> {
        if dimension == 0 {
            return Err(McmcError::InvalidConfig(
                "dense proposal dimension must be positive".to_string(),
            ));
        }
        let expected = dimension.saturating_mul(dimension);
        if covariance.len() != expected {
            return Err(McmcError::DimensionMismatch {
                expected,
                actual: covariance.len(),
            });
        }
        if !jitter.is_finite() || jitter <= 0.0 {
            return Err(McmcError::InvalidConfig(
                "dense proposal jitter must be finite and positive".to_string(),
            ));
        }
        let cholesky = regularized_cholesky(covariance, dimension, jitter).ok_or_else(|| {
            McmcError::InvalidConfig(
                "dense proposal covariance is not positive definite after regularization"
                    .to_string(),
            )
        })?;
        Self::dense_cholesky(dimension, cholesky)
    }

    pub fn validate(&self, dimension: usize) -> Result<(), McmcError> {
        if dimension == 0 {
            return Err(McmcError::InvalidConfig(
                "Gaussian proposal dimension must be positive".to_string(),
            ));
        }
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
            Self::Dense {
                dimension: stored,
                cholesky: _,
            } if *stored != dimension => Err(McmcError::DimensionMismatch {
                expected: dimension,
                actual: *stored,
            }),
            Self::Dense {
                dimension: stored,
                cholesky,
            } if cholesky.len() != (*stored).saturating_mul(*stored) => {
                Err(McmcError::DimensionMismatch {
                    expected: (*stored).saturating_mul(*stored),
                    actual: cholesky.len(),
                })
            }
            Self::Dense {
                dimension: stored,
                cholesky,
            } => {
                for row in 0..*stored {
                    for column in 0..*stored {
                        let value = cholesky[row * *stored + column];
                        if !value.is_finite()
                            || (column > row && value != 0.0)
                            || (column == row && value <= 0.0)
                        {
                            return Err(McmcError::InvalidConfig(
                                "dense proposal Cholesky factor is invalid".to_string(),
                            ));
                        }
                    }
                }
                Ok(())
            }
        }
    }

    /// Marginal standard deviation for one coordinate.
    pub fn component(&self, index: usize) -> f64 {
        match self {
            Self::Isotropic(scale) => *scale,
            Self::Diagonal(scales) => scales[index],
            Self::Dense {
                dimension,
                cholesky,
            } => (0..=index)
                .map(|column| {
                    let value = cholesky[index * *dimension + column];
                    value * value
                })
                .sum::<f64>()
                .sqrt(),
        }
    }

    pub fn fill_displacement<R>(
        &self,
        rng: &mut R,
        normal_buffer: &mut [f64],
        displacement: &mut [f64],
    ) -> Result<(), McmcError>
    where
        R: Rng + ?Sized,
    {
        let dimension = displacement.len();
        if normal_buffer.len() != dimension {
            return Err(McmcError::DimensionMismatch {
                expected: dimension,
                actual: normal_buffer.len(),
            });
        }
        self.validate(dimension)?;
        for normal in normal_buffer.iter_mut() {
            *normal = standard_normal(rng);
        }
        match self {
            Self::Isotropic(scale) => {
                for (output, normal) in displacement.iter_mut().zip(normal_buffer.iter().copied()) {
                    *output = *scale * normal;
                }
            }
            Self::Diagonal(scales) => {
                for ((output, normal), scale) in displacement
                    .iter_mut()
                    .zip(normal_buffer.iter().copied())
                    .zip(scales.iter().copied())
                {
                    *output = scale * normal;
                }
            }
            Self::Dense {
                dimension,
                cholesky,
            } => {
                for (row, output) in displacement.iter_mut().enumerate() {
                    *output = (0..=row)
                        .map(|column| cholesky[row * *dimension + column] * normal_buffer[column])
                        .sum();
                }
            }
        }
        Ok(())
    }

    pub fn set_diagonal(&mut self, scales: Vec<f64>) {
        *self = Self::Diagonal(scales);
    }

    pub fn set_dense_cholesky(&mut self, dimension: usize, cholesky: Vec<f64>) {
        *self = Self::Dense {
            dimension,
            cholesky,
        };
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
        if radius_squared > 0.0 && radius_squared < 1.0 {
            return x * (-2.0 * radius_squared.ln() / radius_squared).sqrt();
        }
    }
}
