mod covariance;
mod dense_covariance;
mod hmc;
mod scale;

pub use covariance::DiagonalCovarianceAdaptation;
pub(crate) use dense_covariance::regularized_cholesky;
pub use dense_covariance::DenseCovarianceAdaptation;
pub(crate) use hmc::MetricUpdate;
pub use hmc::{DualAveraging, HmcWarmup, MetricAdaptation, WarmupWindowConfig};
pub use scale::RobbinsMonroScale;
