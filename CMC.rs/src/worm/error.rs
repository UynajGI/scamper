//! Error type for classical worm configuration and checkpoint validation.

use std::fmt::{Display, Formatter};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WormError {
    message: String,
}

impl WormError {
    pub fn new(message: impl Into<String>) -> Self {
        Self {
            message: message.into(),
        }
    }

    #[inline]
    pub fn message(&self) -> &str {
        &self.message
    }
}

impl Display for WormError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        formatter.write_str(&self.message)
    }
}

impl std::error::Error for WormError {}
