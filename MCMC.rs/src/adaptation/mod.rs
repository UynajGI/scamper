mod covariance;
mod dense_covariance;
mod scale;

pub use covariance::DiagonalCovarianceAdaptation;
pub(crate) use dense_covariance::regularized_cholesky;
pub use dense_covariance::DenseCovarianceAdaptation;
pub use scale::RobbinsMonroScale;
