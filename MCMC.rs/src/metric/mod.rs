//! Euclidean inverse-mass geometries used by Hamiltonian kernels.
//!
//! A metric stores the inverse mass matrix `G = M^-1`. Momentum is sampled
//! from `N(0, M)`, kinetic energy is `0.5 * p^T G p`, and Hamiltonian velocity
//! is `G p`. Position-covariance adaptation therefore installs the estimated
//! covariance directly as the inverse mass geometry.

use rand::Rng;
use serde::{Deserialize, Serialize};

use crate::adaptation::regularized_cholesky;
use crate::proposal::standard_normal;
use crate::McmcError;

/// Built-in geometry family used to validate warmup adaptation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MetricKind {
    Unit,
    Diagonal,
    Dense,
}

/// Allocation-free Euclidean metric contract for HMC integration.
pub trait Metric: Send {
    fn dimension(&self) -> usize;

    fn kind(&self) -> MetricKind;

    fn sample_momentum<R>(&self, momentum: &mut [f64], rng: &mut R) -> Result<(), McmcError>
    where
        R: Rng + ?Sized;

    fn velocity(&self, momentum: &[f64], output: &mut [f64]) -> Result<(), McmcError>;

    /// Compute `(right_position - left_position)^T G momentum` without
    /// requiring the caller to allocate a velocity or displacement vector.
    fn displacement_dot_velocity(
        &self,
        left_position: &[f64],
        right_position: &[f64],
        momentum: &[f64],
    ) -> Result<f64, McmcError> {
        let dimension = self.dimension();
        check_vector(left_position, dimension)?;
        check_vector(right_position, dimension)?;
        check_vector(momentum, dimension)?;
        let mut velocity = vec![0.0; dimension];
        self.velocity(momentum, &mut velocity)?;
        Ok(left_position
            .iter()
            .zip(right_position.iter())
            .zip(velocity.iter())
            .map(|((left, right), velocity)| (right - left) * velocity)
            .sum())
    }

    /// Compute `(G momentum) · other_momentum` without requiring the caller
    /// to allocate a velocity vector.
    fn velocity_dot_momentum(
        &self,
        momentum: &[f64],
        other_momentum: &[f64],
    ) -> Result<f64, McmcError> {
        let dimension = self.dimension();
        check_vector(momentum, dimension)?;
        check_vector(other_momentum, dimension)?;
        let mut velocity = vec![0.0; dimension];
        self.velocity(momentum, &mut velocity)?;
        Ok(velocity
            .iter()
            .zip(other_momentum.iter())
            .map(|(velocity, other)| velocity * other)
            .sum())
    }

    /// Compute `(G momentum) · (first + second)`. Built-in metrics override
    /// this method so generalized NUTS subtree checks remain allocation-free.
    fn velocity_dot_momentum_sum(
        &self,
        momentum: &[f64],
        first: &[f64],
        second: &[f64],
    ) -> Result<f64, McmcError> {
        let dimension = self.dimension();
        check_vector(first, dimension)?;
        check_vector(second, dimension)?;
        let mut sum = Vec::with_capacity(dimension);
        sum.extend(
            first
                .iter()
                .zip(second.iter())
                .map(|(left, right)| left + right),
        );
        self.velocity_dot_momentum(momentum, &sum)
    }

    fn kinetic_energy(&self, momentum: &[f64]) -> Result<f64, McmcError>;

    /// Install a diagonal position-covariance estimate as the inverse mass.
    fn set_diagonal_inverse_mass(&mut self, _diagonal: &[f64]) -> Result<(), McmcError> {
        Err(McmcError::InvalidConfig(
            "this metric does not support diagonal mass-matrix adaptation".to_string(),
        ))
    }

    /// Install a dense row-major position covariance as the inverse mass.
    fn set_dense_inverse_mass(
        &mut self,
        _dimension: usize,
        _inverse_mass: &[f64],
        _jitter: f64,
    ) -> Result<(), McmcError> {
        Err(McmcError::InvalidConfig(
            "this metric does not support dense mass-matrix adaptation".to_string(),
        ))
    }

    fn name(&self) -> &'static str;
}

/// Identity inverse mass matrix.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
pub struct UnitMetric {
    dimension: usize,
}

impl UnitMetric {
    pub fn new(dimension: usize) -> Result<Self, McmcError> {
        validate_dimension(dimension)?;
        Ok(Self { dimension })
    }
}

impl Metric for UnitMetric {
    fn dimension(&self) -> usize {
        self.dimension
    }

    fn kind(&self) -> MetricKind {
        MetricKind::Unit
    }

    fn sample_momentum<R>(&self, momentum: &mut [f64], rng: &mut R) -> Result<(), McmcError>
    where
        R: Rng + ?Sized,
    {
        check_vector(momentum, self.dimension)?;
        for value in momentum {
            *value = standard_normal(rng);
        }
        Ok(())
    }

    fn velocity(&self, momentum: &[f64], output: &mut [f64]) -> Result<(), McmcError> {
        check_vector(momentum, self.dimension)?;
        check_vector(output, self.dimension)?;
        output.copy_from_slice(momentum);
        Ok(())
    }

    fn displacement_dot_velocity(
        &self,
        left_position: &[f64],
        right_position: &[f64],
        momentum: &[f64],
    ) -> Result<f64, McmcError> {
        check_vector(left_position, self.dimension)?;
        check_vector(right_position, self.dimension)?;
        check_vector(momentum, self.dimension)?;
        Ok(left_position
            .iter()
            .zip(right_position.iter())
            .zip(momentum.iter())
            .map(|((left, right), momentum)| (right - left) * momentum)
            .sum())
    }

    fn velocity_dot_momentum(
        &self,
        momentum: &[f64],
        other_momentum: &[f64],
    ) -> Result<f64, McmcError> {
        check_vector(momentum, self.dimension)?;
        check_vector(other_momentum, self.dimension)?;
        Ok(momentum
            .iter()
            .zip(other_momentum.iter())
            .map(|(left, right)| left * right)
            .sum())
    }

    fn velocity_dot_momentum_sum(
        &self,
        momentum: &[f64],
        first: &[f64],
        second: &[f64],
    ) -> Result<f64, McmcError> {
        check_vector(momentum, self.dimension)?;
        check_vector(first, self.dimension)?;
        check_vector(second, self.dimension)?;
        Ok(momentum
            .iter()
            .zip(first.iter())
            .zip(second.iter())
            .map(|((momentum, first), second)| momentum * (first + second))
            .sum())
    }

    fn kinetic_energy(&self, momentum: &[f64]) -> Result<f64, McmcError> {
        check_vector(momentum, self.dimension)?;
        Ok(0.5 * momentum.iter().map(|value| value * value).sum::<f64>())
    }

    fn name(&self) -> &'static str {
        "UnitMetric"
    }
}

/// Diagonal inverse mass matrix.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct DiagonalMetric {
    inverse_mass: Vec<f64>,
}

impl DiagonalMetric {
    pub fn new(inverse_mass: Vec<f64>) -> Result<Self, McmcError> {
        validate_positive_diagonal(&inverse_mass)?;
        Ok(Self { inverse_mass })
    }

    pub fn unit(dimension: usize) -> Result<Self, McmcError> {
        validate_dimension(dimension)?;
        Self::new(vec![1.0; dimension])
    }

    pub fn inverse_mass(&self) -> &[f64] {
        &self.inverse_mass
    }
}

impl Metric for DiagonalMetric {
    fn dimension(&self) -> usize {
        self.inverse_mass.len()
    }

    fn kind(&self) -> MetricKind {
        MetricKind::Diagonal
    }

    fn sample_momentum<R>(&self, momentum: &mut [f64], rng: &mut R) -> Result<(), McmcError>
    where
        R: Rng + ?Sized,
    {
        check_vector(momentum, self.dimension())?;
        for (output, inverse_mass) in momentum.iter_mut().zip(self.inverse_mass.iter().copied()) {
            *output = standard_normal(rng) / inverse_mass.sqrt();
        }
        Ok(())
    }

    fn velocity(&self, momentum: &[f64], output: &mut [f64]) -> Result<(), McmcError> {
        check_vector(momentum, self.dimension())?;
        check_vector(output, self.dimension())?;
        for ((output, momentum), inverse_mass) in output
            .iter_mut()
            .zip(momentum.iter().copied())
            .zip(self.inverse_mass.iter().copied())
        {
            *output = inverse_mass * momentum;
        }
        Ok(())
    }

    fn displacement_dot_velocity(
        &self,
        left_position: &[f64],
        right_position: &[f64],
        momentum: &[f64],
    ) -> Result<f64, McmcError> {
        check_vector(left_position, self.dimension())?;
        check_vector(right_position, self.dimension())?;
        check_vector(momentum, self.dimension())?;
        Ok(left_position
            .iter()
            .zip(right_position.iter())
            .zip(momentum.iter())
            .zip(self.inverse_mass.iter())
            .map(|(((left, right), momentum), inverse_mass)| {
                (right - left) * inverse_mass * momentum
            })
            .sum())
    }

    fn velocity_dot_momentum(
        &self,
        momentum: &[f64],
        other_momentum: &[f64],
    ) -> Result<f64, McmcError> {
        check_vector(momentum, self.dimension())?;
        check_vector(other_momentum, self.dimension())?;
        Ok(momentum
            .iter()
            .zip(other_momentum.iter())
            .zip(self.inverse_mass.iter())
            .map(|((momentum, other), inverse_mass)| inverse_mass * momentum * other)
            .sum())
    }

    fn velocity_dot_momentum_sum(
        &self,
        momentum: &[f64],
        first: &[f64],
        second: &[f64],
    ) -> Result<f64, McmcError> {
        check_vector(momentum, self.dimension())?;
        check_vector(first, self.dimension())?;
        check_vector(second, self.dimension())?;
        Ok(momentum
            .iter()
            .zip(first.iter())
            .zip(second.iter())
            .zip(self.inverse_mass.iter())
            .map(|(((momentum, first), second), inverse_mass)| {
                inverse_mass * momentum * (first + second)
            })
            .sum())
    }

    fn kinetic_energy(&self, momentum: &[f64]) -> Result<f64, McmcError> {
        check_vector(momentum, self.dimension())?;
        Ok(0.5
            * momentum
                .iter()
                .zip(self.inverse_mass.iter())
                .map(|(momentum, inverse_mass)| momentum * momentum * inverse_mass)
                .sum::<f64>())
    }

    fn set_diagonal_inverse_mass(&mut self, diagonal: &[f64]) -> Result<(), McmcError> {
        if diagonal.len() != self.dimension() {
            return Err(McmcError::DimensionMismatch {
                expected: self.dimension(),
                actual: diagonal.len(),
            });
        }
        validate_positive_diagonal(diagonal)?;
        self.inverse_mass.copy_from_slice(diagonal);
        Ok(())
    }

    fn name(&self) -> &'static str {
        "DiagonalMetric"
    }
}

/// Dense inverse mass matrix with a cached lower Cholesky factor.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct DenseMetric {
    dimension: usize,
    inverse_mass: Vec<f64>,
    cholesky: Vec<f64>,
}

impl DenseMetric {
    pub fn unit(dimension: usize) -> Result<Self, McmcError> {
        validate_dimension(dimension)?;
        let mut inverse_mass = vec![0.0; dimension.saturating_mul(dimension)];
        for index in 0..dimension {
            inverse_mass[index * dimension + index] = 1.0;
        }
        Self::from_inverse_mass(dimension, &inverse_mass, 1.0e-12)
    }

    pub fn from_inverse_mass(
        dimension: usize,
        inverse_mass: &[f64],
        jitter: f64,
    ) -> Result<Self, McmcError> {
        validate_dense_input(dimension, inverse_mass, jitter)?;
        let mut symmetric = symmetrize(inverse_mass, dimension);
        let cholesky = regularized_cholesky(&symmetric, dimension, jitter).ok_or_else(|| {
            McmcError::InvalidConfig(
                "dense inverse mass is not positive definite after regularization".to_string(),
            )
        })?;
        // Store the exact matrix represented by the accepted Cholesky factor so
        // velocity, kinetic energy and momentum sampling remain synchronized.
        symmetric = lower_times_transpose(&cholesky, dimension);
        Ok(Self {
            dimension,
            inverse_mass: symmetric,
            cholesky,
        })
    }

    pub fn inverse_mass(&self) -> &[f64] {
        &self.inverse_mass
    }

    pub fn cholesky(&self) -> &[f64] {
        &self.cholesky
    }
}

impl Metric for DenseMetric {
    fn dimension(&self) -> usize {
        self.dimension
    }

    fn kind(&self) -> MetricKind {
        MetricKind::Dense
    }

    fn sample_momentum<R>(&self, momentum: &mut [f64], rng: &mut R) -> Result<(), McmcError>
    where
        R: Rng + ?Sized,
    {
        check_vector(momentum, self.dimension)?;
        for value in momentum.iter_mut() {
            *value = standard_normal(rng);
        }
        // If G = L L^T is the inverse mass, solve L^T p = z so
        // Cov(p) = G^-1.
        for row in (0..self.dimension).rev() {
            let mut value = momentum[row];
            for (column, momentum_value) in momentum.iter().copied().enumerate().skip(row + 1) {
                value -= self.cholesky[column * self.dimension + row] * momentum_value;
            }
            momentum[row] = value / self.cholesky[row * self.dimension + row];
        }
        Ok(())
    }

    fn velocity(&self, momentum: &[f64], output: &mut [f64]) -> Result<(), McmcError> {
        check_vector(momentum, self.dimension)?;
        check_vector(output, self.dimension)?;
        for (row, output) in output.iter_mut().enumerate() {
            *output = self.inverse_mass[row * self.dimension..(row + 1) * self.dimension]
                .iter()
                .zip(momentum.iter())
                .map(|(matrix, momentum)| matrix * momentum)
                .sum();
        }
        Ok(())
    }

    fn displacement_dot_velocity(
        &self,
        left_position: &[f64],
        right_position: &[f64],
        momentum: &[f64],
    ) -> Result<f64, McmcError> {
        check_vector(left_position, self.dimension)?;
        check_vector(right_position, self.dimension)?;
        check_vector(momentum, self.dimension)?;
        let mut product = 0.0;
        for row in 0..self.dimension {
            let velocity = self.inverse_mass[row * self.dimension..(row + 1) * self.dimension]
                .iter()
                .zip(momentum.iter())
                .map(|(matrix, momentum)| matrix * momentum)
                .sum::<f64>();
            product += (right_position[row] - left_position[row]) * velocity;
        }
        Ok(product)
    }

    fn velocity_dot_momentum(
        &self,
        momentum: &[f64],
        other_momentum: &[f64],
    ) -> Result<f64, McmcError> {
        check_vector(momentum, self.dimension)?;
        check_vector(other_momentum, self.dimension)?;
        let mut product = 0.0;
        for (row, other) in other_momentum.iter().copied().enumerate() {
            let velocity = self.inverse_mass[row * self.dimension..(row + 1) * self.dimension]
                .iter()
                .zip(momentum.iter())
                .map(|(matrix, momentum)| matrix * momentum)
                .sum::<f64>();
            product += velocity * other;
        }
        Ok(product)
    }

    fn velocity_dot_momentum_sum(
        &self,
        momentum: &[f64],
        first: &[f64],
        second: &[f64],
    ) -> Result<f64, McmcError> {
        check_vector(momentum, self.dimension)?;
        check_vector(first, self.dimension)?;
        check_vector(second, self.dimension)?;
        let mut product = 0.0;
        for (row, (first, second)) in first.iter().zip(second.iter()).enumerate() {
            let velocity = self.inverse_mass[row * self.dimension..(row + 1) * self.dimension]
                .iter()
                .zip(momentum.iter())
                .map(|(matrix, momentum)| matrix * momentum)
                .sum::<f64>();
            product += velocity * (first + second);
        }
        Ok(product)
    }

    fn kinetic_energy(&self, momentum: &[f64]) -> Result<f64, McmcError> {
        check_vector(momentum, self.dimension)?;
        let mut quadratic = 0.0;
        for row in 0..self.dimension {
            let velocity = self.inverse_mass[row * self.dimension..(row + 1) * self.dimension]
                .iter()
                .zip(momentum.iter())
                .map(|(matrix, momentum)| matrix * momentum)
                .sum::<f64>();
            quadratic += momentum[row] * velocity;
        }
        if !quadratic.is_finite() || quadratic < 0.0 {
            return Err(McmcError::InvalidConfig(
                "dense metric produced invalid kinetic energy".to_string(),
            ));
        }
        Ok(0.5 * quadratic)
    }

    fn set_diagonal_inverse_mass(&mut self, diagonal: &[f64]) -> Result<(), McmcError> {
        if diagonal.len() != self.dimension {
            return Err(McmcError::DimensionMismatch {
                expected: self.dimension,
                actual: diagonal.len(),
            });
        }
        validate_positive_diagonal(diagonal)?;
        self.inverse_mass.fill(0.0);
        self.cholesky.fill(0.0);
        for (index, value) in diagonal.iter().copied().enumerate() {
            self.inverse_mass[index * self.dimension + index] = value;
            self.cholesky[index * self.dimension + index] = value.sqrt();
        }
        Ok(())
    }

    fn set_dense_inverse_mass(
        &mut self,
        dimension: usize,
        inverse_mass: &[f64],
        jitter: f64,
    ) -> Result<(), McmcError> {
        if dimension != self.dimension {
            return Err(McmcError::DimensionMismatch {
                expected: self.dimension,
                actual: dimension,
            });
        }
        let replacement = Self::from_inverse_mass(dimension, inverse_mass, jitter)?;
        *self = replacement;
        Ok(())
    }

    fn name(&self) -> &'static str {
        "DenseMetric"
    }
}

fn validate_dimension(dimension: usize) -> Result<(), McmcError> {
    if dimension == 0 {
        Err(McmcError::InvalidConfig(
            "metric dimension must be positive".to_string(),
        ))
    } else {
        Ok(())
    }
}

fn check_vector(vector: &[f64], expected: usize) -> Result<(), McmcError> {
    if vector.len() == expected {
        Ok(())
    } else {
        Err(McmcError::DimensionMismatch {
            expected,
            actual: vector.len(),
        })
    }
}

fn validate_positive_diagonal(diagonal: &[f64]) -> Result<(), McmcError> {
    validate_dimension(diagonal.len())?;
    if diagonal
        .iter()
        .all(|value| value.is_finite() && *value > 0.0)
    {
        Ok(())
    } else {
        Err(McmcError::InvalidConfig(
            "inverse-mass diagonal must be finite and positive".to_string(),
        ))
    }
}

fn validate_dense_input(dimension: usize, matrix: &[f64], jitter: f64) -> Result<(), McmcError> {
    validate_dimension(dimension)?;
    let expected = dimension.saturating_mul(dimension);
    if matrix.len() != expected {
        return Err(McmcError::DimensionMismatch {
            expected,
            actual: matrix.len(),
        });
    }
    if matrix.iter().any(|value| !value.is_finite()) {
        return Err(McmcError::InvalidConfig(
            "dense inverse mass must contain only finite values".to_string(),
        ));
    }
    if !jitter.is_finite() || jitter <= 0.0 {
        return Err(McmcError::InvalidConfig(
            "dense metric jitter must be finite and positive".to_string(),
        ));
    }
    Ok(())
}

fn symmetrize(matrix: &[f64], dimension: usize) -> Vec<f64> {
    let mut symmetric = vec![0.0; matrix.len()];
    for row in 0..dimension {
        for column in 0..dimension {
            symmetric[row * dimension + column] =
                0.5 * (matrix[row * dimension + column] + matrix[column * dimension + row]);
        }
    }
    symmetric
}

fn lower_times_transpose(lower: &[f64], dimension: usize) -> Vec<f64> {
    let mut matrix = vec![0.0; dimension.saturating_mul(dimension)];
    for row in 0..dimension {
        for column in 0..dimension {
            matrix[row * dimension + column] = (0..dimension)
                .map(|inner| lower[row * dimension + inner] * lower[column * dimension + inner])
                .sum();
        }
    }
    matrix
}
