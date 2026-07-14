//! Errors from dynamic and rejection-free classical Monte Carlo kernels.

use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DynamicsError {
    message: String,
}

impl DynamicsError {
    pub fn new(message: impl Into<String>) -> Self {
        Self {
            message: message.into(),
        }
    }
}

impl fmt::Display for DynamicsError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.message)
    }
}

impl std::error::Error for DynamicsError {}
