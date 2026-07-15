//! Task information and progress tracking.

use crate::CarloError;
use std::collections::HashMap;
use std::path::PathBuf;
use std::time::SystemTime;

/// One parameter set within a job, with required-key validation.
#[derive(Debug, Clone)]
pub struct TaskInfo {
    name: String,
    params: HashMap<String, String>,
}

impl TaskInfo {
    /// Create a new task.
    ///
    /// Returns [`CarloError::InvalidConfig`] if `sweeps`, `thermalization`,
    /// or `binsize` are absent from `params`.
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

    /// Task name (e.g. `"task0003"`).
    pub fn name(&self) -> &str {
        &self.name
    }
    /// Raw parameter map.
    pub fn params(&self) -> &HashMap<String, String> {
        &self.params
    }
    /// Typed parameter lookup (parses via `FromStr`).
    pub fn get<T: std::str::FromStr>(&self, key: &str) -> Option<T> {
        self.params.get(key).and_then(|v| v.parse().ok())
    }
}

/// Generate a canonical task name from a 1-based task ID.
pub fn task_name(task_id: u64) -> String {
    format!("task{:04}", task_id)
}

/// Progress snapshot for a task.
#[derive(Debug, Clone)]
pub struct TaskProgress {
    /// Configured measurement sweep target.
    pub target_sweeps: u64,
    /// Estimated completed sweeps.
    pub sweeps: u64,
    /// Number of independent runs.
    pub num_runs: u64,
    /// Fraction of sweeps spent in thermalization (0–1).
    pub thermalization_fraction: f64,
    /// Task output directory.
    pub dir: PathBuf,
    /// Last modification time of the task directory.
    pub last_modified: Option<SystemTime>,
    /// Total elapsed time estimate.
    pub elapsed_seconds: f64,
    /// Estimated sweeps per second rate.
    pub sweeps_per_sec: f64,
}

/// List files in `dir` whose names match `pattern` (a regex).
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
