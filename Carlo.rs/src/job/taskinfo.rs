//! Task information and progress tracking.

use crate::CarloError;
use std::collections::HashMap;
use std::path::PathBuf;
use std::time::SystemTime;

/// Task parameters with validation.
#[derive(Debug, Clone)]
pub struct TaskInfo {
    name: String,
    params: HashMap<String, String>,
}

impl TaskInfo {
    pub fn new(name: &str, params: HashMap<String, String>) -> Result<Self, CarloError> {
        let required = ["sweeps", "thermalization", "binsize"];
        for key in required {
            if !params.contains_key(key) {
                return Err(CarloError::InvalidConfig {
                    field: key.into(),
                    reason: format!("Task {} missing required parameter {}", name, key),
                });
            }
        }
        Ok(Self {
            name: name.to_string(),
            params,
        })
    }

    pub fn name(&self) -> &str {
        &self.name
    }
    pub fn params(&self) -> &HashMap<String, String> {
        &self.params
    }
    pub fn get<T: std::str::FromStr>(&self, key: &str) -> Option<T> {
        self.params.get(key).and_then(|v| v.parse().ok())
    }
}

pub fn task_name(task_id: u64) -> String {
    format!("task{:04}", task_id)
}

#[derive(Debug, Clone)]
pub struct TaskProgress {
    pub target_sweeps: u64,
    pub sweeps: u64,
    pub num_runs: u64,
    pub thermalization_fraction: f64,
    pub dir: PathBuf,
    /// Last modification time of the task directory.
    pub last_modified: Option<SystemTime>,
    /// Total elapsed time estimate.
    pub elapsed_seconds: f64,
    /// Estimated sweeps per second rate.
    pub sweeps_per_sec: f64,
}

pub fn list_run_files(dir: &PathBuf, pattern: &str) -> Vec<PathBuf> {
    use std::fs;
    let re = regex::Regex::new(pattern).unwrap_or_else(|_| regex::Regex::new(".*").unwrap());
    fs::read_dir(dir)
        .ok()
        .into_iter()
        .flatten()
        .filter_map(|e| e.ok())
        .filter(|e| re.is_match(&e.file_name().to_string_lossy()))
        .map(|e| e.path())
        .collect()
}
