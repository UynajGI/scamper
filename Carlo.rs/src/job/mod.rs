//! Job and task management for Monte Carlo simulations.
//!
//! This module provides tools for defining and managing simulation jobs:
//!
//! - [`JobInfo`]: Overall job configuration (name, directory, tasks, timing)
//! - [`TaskInfo`]: Individual simulation task with parameters
//! - [`TaskMaker`]: Helper for generating parameter sweeps
//!
//! # Job Configuration
//!
//! Jobs are typically defined via JSON configuration files:
//!
//! ```json
//! {
//!     "name": "ising_study",
//!     "mc_type": "Ising",
//!     "run_time": "24:00:00",
//!     "checkpoint_time": "1:00:00",
//!     "tasks": [
//!         {"name": "L100_T2.0", "params": {"L": 100, "beta": 2.0}},
//!         {"name": "L100_T2.5", "params": {"L": 100, "beta": 2.5}}
//!     ]
//! }
//! ```
//!
//! # Time Format
//!
//! Time durations use `[[HH:]MM:]SS` format (e.g., `24:00:00` = 24 hours).
//! For SLURM jobs, [`run_time_from_slurm()`] automatically detects remaining time.

mod jobinfo;
mod taskinfo;
mod taskmaker;

pub use jobinfo::{parse_duration, run_time_from_slurm, JobInfo};
pub use taskinfo::{list_run_files, task_name, TaskInfo, TaskProgress};
pub use taskmaker::TaskMaker;
