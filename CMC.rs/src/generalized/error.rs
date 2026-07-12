//! Errors returned by generalized-ensemble infrastructure.

use std::fmt::{Display, Formatter};

/// Invalid axis, histogram, checkpoint or refinement configuration.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GeneralizedError {
    detail: String,
}

impl GeneralizedError {
    pub(crate) fn new(detail: impl Into<String>) -> Self {
        Self {
            detail: detail.into(),
        }
    }

    /// Human-readable validation detail.
    #[inline]
    pub fn detail(&self) -> &str {
        &self.detail
    }
}

impl Display for GeneralizedError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        formatter.write_str(&self.detail)
    }
}

impl std::error::Error for GeneralizedError {}
