//! Job configuration and time handling.

use crate::{
    job::{TaskInfo, TaskProgress},
    CarloError,
};
use chrono::{DateTime, Utc};
use std::path::PathBuf;
use std::time::Duration;
use std::time::SystemTime;

pub fn parse_duration(s: &str) -> Result<Duration, CarloError> {
    let re = regex::Regex::new(
        r"^((((?P<days>\d+)-)?(?P<hours>\d+):)?(?P<minutes>\d+):)?(?P<seconds>\d+)$",
    )
    .unwrap();
    let caps = re.captures(s).ok_or_else(|| CarloError::InvalidConfig {
        field: "duration".into(),
        reason: format!("{} does not match [[HH:]MM:]SS format", s),
    })?;
    let conv = |name: &str| {
        caps.name(name)
            .map(|m| m.as_str().parse::<u64>().unwrap_or(0))
            .unwrap_or(0)
    };
    Ok(Duration::from_secs(conv("seconds"))
        + Duration::from_secs(conv("minutes") * 60)
        + Duration::from_secs(conv("hours") * 3600)
        + Duration::from_secs(conv("days") * 86400))
}

pub fn run_time_from_slurm(grace_factor: f64, default: Duration) -> Duration {
    if let Some(end_time_str) = std::env::var_os("SLURM_JOB_END_TIME") {
        if let Ok(end_time_unix) = end_time_str.to_string_lossy().parse::<i64>() {
            let now = Utc::now().timestamp();
            let remaining = (end_time_unix - now).max(0) as f64;
            return Duration::from_secs((remaining * grace_factor) as u64);
        }
    }
    default
}

#[derive(Debug, Clone)]
pub struct JobInfo {
    name: String,
    dir: std::path::PathBuf,
    #[allow(dead_code)]
    mc_type: String,
    #[allow(dead_code)]
    rng_type: String,
    pub tasks: Vec<TaskInfo>,
    checkpoint_time: Duration,
    run_time: Duration,
    #[allow(dead_code)]
    ranks_per_run: usize,
}

impl JobInfo {
    pub fn new(
        job_file: &str,
        mc_type: &str,
        rng_type: &str,
        tasks: Vec<TaskInfo>,
        checkpoint_time: Duration,
        run_time: Duration,
        ranks_per_run: usize,
    ) -> Self {
        let expanded = shellexpand::tilde(job_file);
        Self {
            name: std::path::Path::new(expanded.as_ref())
                .file_name()
                .unwrap_or_default()
                .to_string_lossy()
                .to_string(),
            dir: std::path::PathBuf::from(expanded.as_ref()).join(".data"),
            mc_type: mc_type.to_string(),
            rng_type: rng_type.to_string(),
            tasks,
            checkpoint_time,
            run_time,
            ranks_per_run,
        }
    }

    pub fn task_dir(&self, task: &TaskInfo) -> std::path::PathBuf {
        self.dir.join(task.name())
    }
    pub fn name(&self) -> &str {
        &self.name
    }
    pub fn dir(&self) -> &std::path::PathBuf {
        &self.dir
    }
    pub fn is_checkpoint_time(&self, last_checkpoint: DateTime<Utc>) -> bool {
        Utc::now()
            >= last_checkpoint
                + chrono::Duration::from_std(self.checkpoint_time)
                    .expect("checkpoint_time duration out of range")
    }
    pub fn is_end_time(&self, start: DateTime<Utc>) -> bool {
        Utc::now()
            >= start
                + chrono::Duration::from_std(self.run_time).expect("run_time duration out of range")
    }

    /// Returns the filename of the results JSON file.
    pub fn result_filename(&self) -> PathBuf {
        self.dir
            .parent()
            .map(|p| p.join(format!("{}.results.json", self.name)))
            .unwrap_or_else(|| PathBuf::from(format!("{}.results.json", self.name)))
    }

    /// Create the job directory structure.
    pub fn create_directories(&self) -> Result<(), CarloError> {
        std::fs::create_dir_all(&self.dir).map_err(|e| CarloError::IoError {
            path: self.dir.clone(),
            source: e,
        })?;
        for task in &self.tasks {
            let task_dir = self.task_dir(task);
            std::fs::create_dir_all(&task_dir).map_err(|e| CarloError::IoError {
                path: task_dir,
                source: e,
            })?;
        }
        Ok(())
    }

    /// Read progress of all tasks.
    pub fn read_progress(&self) -> Vec<TaskProgress> {
        self.tasks
            .iter()
            .filter_map(|task| {
                let task_dir = self.task_dir(task);
                read_task_progress(&task_dir, task)
            })
            .collect()
    }

    /// Concatenate all task results into a single JSON file.
    pub fn concatenate_results(&self) -> Result<(), CarloError> {
        let results: Vec<serde_json::Value> = self
            .tasks
            .iter()
            .filter_map(|task| {
                let result_file = self.task_dir(task).join("results.json");
                std::fs::read_to_string(&result_file)
                    .ok()
                    .and_then(|s| serde_json::from_str(&s).ok())
            })
            .collect();

        let output =
            serde_json::to_string_pretty(&results).map_err(CarloError::SerializationError)?;

        std::fs::write(self.result_filename(), output).map_err(|e| CarloError::IoError {
            path: self.result_filename(),
            source: e,
        })?;

        Ok(())
    }

    /// Get run time duration.
    pub fn run_time(&self) -> Duration {
        self.run_time
    }

    /// Get checkpoint time interval.
    pub fn checkpoint_time(&self) -> Duration {
        self.checkpoint_time
    }

    /// Build JobInfo from JSON config file.
    pub fn from_config(json: &str) -> Result<Self, CarloError> {
        let config: serde_json::Value =
            serde_json::from_str(json).map_err(CarloError::SerializationError)?;

        let name = config["name"].as_str().unwrap_or("carlo").to_string();
        let mc_type = config["mc_type"].as_str().unwrap_or("unknown").to_string();
        let rng_type = config["rng_type"]
            .as_str()
            .unwrap_or("Xoshiro256PlusPlus")
            .to_string();

        let run_time = config["run_time"]
            .as_str()
            .map(|s| parse_duration(s).unwrap_or(Duration::from_secs(86400)))
            .unwrap_or(Duration::from_secs(86400));

        let checkpoint_time = config["checkpoint_time"]
            .as_str()
            .map(|s| parse_duration(s).unwrap_or(Duration::from_secs(3600)))
            .unwrap_or(Duration::from_secs(3600));

        let ranks_per_run = config["ranks_per_run"].as_u64().unwrap_or(1) as usize;

        // Parse tasks
        let tasks: Vec<TaskInfo> = config["tasks"]
            .as_array()
            .map(|arr| {
                arr.iter()
                    .filter_map(|t| {
                        let name = t["name"].as_str().unwrap_or("task").to_string();
                        let params = t["params"]
                            .as_object()
                            .map(|m| {
                                m.iter()
                                    .map(|(k, v)| {
                                        (k.clone(), v.to_string().trim_matches('"').to_string())
                                    })
                                    .collect()
                            })
                            .unwrap_or_default();
                        TaskInfo::new(&name, params).ok()
                    })
                    .collect()
            })
            .unwrap_or_default();

        Ok(Self {
            name,
            dir: PathBuf::from(".data"),
            mc_type,
            rng_type,
            tasks,
            checkpoint_time,
            run_time,
            ranks_per_run,
        })
    }
}

/// Read progress for a single task.
fn read_task_progress(task_dir: &PathBuf, task: &TaskInfo) -> Option<TaskProgress> {
    if !task_dir.exists() {
        return None;
    }

    let target_sweeps = task.get::<u64>("sweeps").unwrap_or(0);

    // Count run directories
    let run_dirs: Vec<_> = std::fs::read_dir(task_dir)
        .ok()?
        .filter_map(|e| e.ok())
        .filter(|e| e.file_name().to_string_lossy().starts_with("run"))
        .collect();

    let num_runs = run_dirs.len() as u64;

    // Estimate sweeps from run count (simplified)
    let sweeps = if num_runs > 0 {
        target_sweeps / num_runs.max(1)
    } else {
        0
    };

    let thermalization_fraction = if num_runs > 0 { 1.0 } else { 0.0 };

    // Get last modification time and estimate elapsed time and rate
    let (last_modified, elapsed_seconds, sweeps_per_sec) = task_dir
        .metadata()
        .ok()
        .and_then(|m| m.modified().ok())
        .map(|modified| {
            let elapsed = SystemTime::now()
                .duration_since(modified)
                .unwrap_or(Duration::ZERO)
                .as_secs_f64();
            let rate = if elapsed > 0.0 {
                sweeps as f64 / elapsed
            } else {
                0.0
            };
            (Some(modified), elapsed, rate)
        })
        .unwrap_or((None, 0.0, 0.0));

    Some(TaskProgress {
        target_sweeps,
        sweeps,
        num_runs,
        thermalization_fraction,
        dir: task_dir.clone(),
        last_modified,
        elapsed_seconds,
        sweeps_per_sec,
    })
}
