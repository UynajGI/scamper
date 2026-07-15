//! JSON serialization for simulation results.

use std::fs;
use std::path::Path;

use crate::{CarloError, Results};

/// Write `results` to `path` as pretty-printed JSON.
pub fn save_json(results: &Results, path: &Path) -> Result<(), CarloError> {
    let json = results.to_json()?;
    fs::write(path, json).map_err(|e| CarloError::InvalidConfig {
        field: "json_output".into(),
        reason: e.to_string(),
    })
}
