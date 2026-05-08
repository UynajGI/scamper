use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Generic parameter container with string keys.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Params {
    values: HashMap<String, String>,
}

impl Params {
    pub fn new() -> Self {
        Self {
            values: HashMap::new(),
        }
    }

    /// Set a parameter (converts to string).
    pub fn set<T: ToString>(&mut self, key: &str, value: T) {
        self.values.insert(key.to_string(), value.to_string());
    }

    /// Get a parameter (parses from string).
    pub fn get<T: std::str::FromStr>(&self, key: &str) -> Option<T> {
        self.values.get(key).and_then(|v| v.parse::<T>().ok())
    }

    /// Check if a parameter exists.
    pub fn contains(&self, key: &str) -> bool {
        self.values.contains_key(key)
    }

    /// Merge another Params, overwriting duplicates.
    pub fn merge(&mut self, other: &Params) {
        for (k, v) in &other.values {
            self.values.insert(k.clone(), v.clone());
        }
    }
}

impl Default for Params {
    fn default() -> Self {
        Self::new()
    }
}
