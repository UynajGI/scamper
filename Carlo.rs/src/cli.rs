//! Command line interface for Carlo.
//!
//! The CLI provides tools for managing Monte Carlo simulations:
//!
//! # Commands
//!
//! - `carlo-rs run`: Create job infrastructure for an application-owned simulation
//! - `carlo-rs status`: Check simulation progress
//! - `carlo-rs merge`: Combine results from completed runs
//! - `carlo-rs delete`: Remove simulation data
//!
//! # Usage
//!
//! ```bash
//! # Create job infrastructure for an application-owned simulation
//! carlo-rs run --job-dir my_job/
//!
//! # Check progress
//! carlo-rs status --job-dir my_job/
//!
//! # Merge results (requires hdf5 feature)
//! carlo-rs merge --job-dir my_job/
//!
//! # Clean up
//! carlo-rs delete --job-dir my_job/
//! ```
//!
//! # Options
//!
//! - `-s, --single`: Run in single-threaded mode (no parallelization)
//! - `-r, --restart`: Delete existing data and start fresh
//!
//! # Job Directory Structure
//!
//! ```text
//! my_job/
//! ├── my_job.json       # Job configuration
//! └── .data/
//!     ├── task1/
//!     │   ├── run1.meas.h5
//!     │   ├── run2.meas.h5
//!     │   └── results.json
//!     └── task2/
//!         └── ...
//! ```

use crate::{job::JobInfo, progress::print_status_table, CarloError};
use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[cfg(feature = "hdf5")]
use crate::merge::{merge_results, MergeOptions};

#[derive(Parser)]
#[command(name = "carlo-rs")]
#[command(about = "Monte Carlo simulation framework", version)]
struct Cli {
    /// Job directory (contains .data folder)
    #[arg(short, long, default_value = ".")]
    job_dir: PathBuf,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Create job infrastructure for an application-owned simulation
    Run {
        /// Run in single-threaded mode
        #[arg(short, long)]
        single: bool,
        /// Delete existing files and start fresh
        #[arg(short, long)]
        restart: bool,
    },
    /// Check simulation progress
    Status,
    /// Merge results from runs
    Merge,
    /// Delete simulation data
    Delete,
}

pub fn run() -> Result<(), CarloError> {
    let cli = Cli::parse();

    // Build JobInfo from directory
    let job = build_job_info(&cli.job_dir)?;

    match cli.command {
        Commands::Run { single, restart } => cli_run(&job, single, restart),
        Commands::Status => cli_status(&job),
        Commands::Merge => cli_merge(&job),
        Commands::Delete => cli_delete(&job),
    }
}

fn build_job_info(job_dir: &std::path::Path) -> Result<JobInfo, CarloError> {
    let name = job_dir
        .file_name()
        .map(|n| n.to_string_lossy().to_string())
        .unwrap_or_else(|| "carlo".to_string());

    // Try to load job config from JSON
    let config_file = job_dir.with_extension("json");
    if config_file.exists() {
        let content = std::fs::read_to_string(&config_file).map_err(|e| CarloError::IoError {
            path: config_file.clone(),
            source: e,
        })?;
        // Parse and build JobInfo from config
        // For now, create a default job
        JobInfo::from_config(&content)
    } else {
        // Create minimal job info from existing directory
        Ok(JobInfo::new(
            &name,
            "unknown",
            "Xoshiro256PlusPlus",
            vec![],
            std::time::Duration::from_secs(3600),
            std::time::Duration::from_secs(86400),
            1,
        ))
    }
}

fn cli_run(job: &JobInfo, single: bool, restart: bool) -> Result<(), CarloError> {
    if restart {
        cli_delete(job)?;
    }

    // Create directories
    job.create_directories()?;

    let mode = if single {
        "single-threaded"
    } else {
        "MPI parallel"
    };
    tracing::info!("Starting simulation '{}' in {} mode", job.name(), mode);

    // The actual simulation would be run by the user's code
    // This CLI provides infrastructure support
    println!("Job directory: {:?}", job.dir());
    println!("Tasks: {}", job.tasks.len());
    println!("Run time: {:?}", job.run_time());
    println!("Checkpoint interval: {:?}", job.checkpoint_time());

    Ok(())
}

fn cli_status(job: &JobInfo) -> Result<(), CarloError> {
    let tasks = job.read_progress();

    if tasks.is_empty() {
        println!("No progress data found. Has the simulation been started?");
        return Ok(());
    }

    print_status_table(&tasks);

    let all_done = tasks.iter().all(|t| t.sweeps >= t.target_sweeps);
    if all_done {
        std::process::exit(0);
    } else {
        std::process::exit(1);
    }
}

fn cli_merge(job: &JobInfo) -> Result<(), CarloError> {
    println!("Merging results for job '{}'...", job.name());

    #[cfg(feature = "hdf5")]
    {
        let mut success_count = 0;
        let mut error_count = 0;

        for task in &job.tasks {
            let task_dir = job.task_dir(task);
            match merge_results(&task_dir, &MergeOptions::default()) {
                Ok(results) => {
                    // Convert params (String→String) to JSON with typed values
                    let params_json: serde_json::Map<String, serde_json::Value> = task
                        .params()
                        .iter()
                        .map(|(k, v)| {
                            // Try integer, then float, fall back to string
                            let val = v
                                .parse::<i64>()
                                .map(serde_json::Value::from)
                                .or_else(|_| v.parse::<f64>().map(serde_json::Value::from))
                                .unwrap_or_else(|_| serde_json::Value::String(v.clone()));
                            (k.clone(), val)
                        })
                        .collect();

                    // Build Carlo.jl-compatible task result wrapper
                    let task_result = serde_json::json!({
                        "task": task.name(),
                        "parameters": params_json,
                        "results": results,
                    });

                    let result_file = task_dir.join("results.json");
                    let output = serde_json::to_string_pretty(&task_result)
                        .map_err(CarloError::SerializationError)?;
                    std::fs::write(&result_file, output).map_err(|e| CarloError::IoError {
                        path: result_file,
                        source: e,
                    })?;
                    success_count += 1;
                    println!("  ✓ Merged {}", task.name());
                }
                Err(e) => {
                    tracing::warn!("Failed to merge {}: {}", task.name(), e);
                    error_count += 1;
                }
            }
        }

        // Concatenate all results
        job.concatenate_results()?;

        println!(
            "Merged {}/{} tasks successfully",
            success_count,
            job.tasks.len()
        );
        if error_count > 0 {
            println!("{} tasks had errors", error_count);
        }
    }

    #[cfg(not(feature = "hdf5"))]
    {
        println!("Merge requires 'hdf5' feature. Recompile with --features hdf5");
        // Still try to concatenate existing results
        job.concatenate_results()?;
    }

    Ok(())
}

fn cli_delete(job: &JobInfo) -> Result<(), CarloError> {
    let results_file = job.result_filename();
    if results_file.exists() {
        std::fs::remove_file(&results_file).map_err(|e| CarloError::IoError {
            path: results_file,
            source: e,
        })?;
    }

    if job.dir().exists() {
        std::fs::remove_dir_all(job.dir()).map_err(|e| CarloError::IoError {
            path: job.dir().clone(),
            source: e,
        })?;
    }

    println!("Deleted job data for '{}'", job.name());
    Ok(())
}
