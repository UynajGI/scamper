//! MPI backend for distributed Monte Carlo simulations.
//!
//! The implementation has two complementary entry points:
//!
//! - [`MpiBackend`] implements the generic [`Backend`](super::Backend) trait by
//!   assigning logical task IDs to MPI ranks deterministically.
//! - [`run_distributed`] runs parameter tasks through a controller/worker-group
//!   scheduler. Rank 0 is the controller; the remaining ranks are partitioned
//!   into groups of `ranks_per_run` ranks. Ranks inside a group execute the same
//!   [`Run`] and may coordinate through `MonteCarlo::sweep_with_comm` and
//!   `MonteCarlo::measure_with_comm`.
//!
//! MPI is initialized exactly once per entry point and its [`Universe`] is kept
//! alive for the complete lifetime of all communicators. No communicator handle
//! is copied with `unsafe` code.

use rand_core::{Rng, SeedableRng};
use serde::{de::DeserializeOwned, Deserialize, Serialize};
use std::collections::{HashMap, VecDeque};
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

#[cfg(feature = "mpi")]
use mpi::environment::{Threading, Universe};
#[cfg(feature = "mpi")]
use mpi::topology::{Color, SimpleCommunicator};
#[cfg(feature = "mpi")]
use mpi::traits::*;

use crate::{
    CarloError, FromParams, MonteCarlo, Params, Results, RngPhase, RngStreamKey, Run, RunConfig,
    RunId, TaskId,
};

// ============================================================================
// Time limits
// ============================================================================

/// Wall-clock limits used by distributed workers.
#[derive(Debug, Clone)]
pub struct TimeLimits {
    checkpoint_time: Option<Duration>,
    run_time: Option<Duration>,
    last_checkpoint: Instant,
    start_time: Instant,
}

impl TimeLimits {
    /// Create limits with both checkpoint and total runtime enabled.
    pub fn new(checkpoint_time: Duration, run_time: Duration) -> Self {
        Self::optional(Some(checkpoint_time), Some(run_time))
    }

    /// Create limits where either timer may be disabled.
    pub fn optional(checkpoint_time: Option<Duration>, run_time: Option<Duration>) -> Self {
        Self::from_start(checkpoint_time, run_time, Instant::now())
    }

    /// Create limits whose total runtime is measured from a shared job start.
    pub fn from_start(
        checkpoint_time: Option<Duration>,
        run_time: Option<Duration>,
        start_time: Instant,
    ) -> Self {
        Self {
            checkpoint_time,
            run_time,
            last_checkpoint: Instant::now(),
            start_time,
        }
    }

    pub fn should_checkpoint(&self) -> bool {
        self.checkpoint_time
            .is_some_and(|limit| self.last_checkpoint.elapsed() >= limit)
    }

    pub fn should_finish(&self) -> bool {
        self.run_time
            .is_some_and(|limit| self.start_time.elapsed() >= limit)
    }

    pub fn reset_checkpoint(&mut self) {
        self.last_checkpoint = Instant::now();
    }
}

// ============================================================================
// Errors and RNG bounds
// ============================================================================

/// Errors produced by the MPI backend.
#[derive(Debug, thiserror::Error)]
pub enum MpiError {
    #[error("MPI initialization failed or MPI was already initialized by another owner")]
    InitFailed,

    #[error("MPI communication failed: {0}")]
    Communication(String),

    #[error("invalid MPI topology: {0}")]
    InvalidTopology(String),

    #[error("invalid distributed scheduler transition: {0}")]
    InvalidTransition(String),

    #[error("distributed worker failed: {0}")]
    Worker(String),
}

/// RNG requirements for distributed runs.
///
/// With the `hdf5` feature enabled, this additionally requires checkpoint
/// serialization support. Without `hdf5`, every `Rng + SeedableRng + Send`
/// automatically implements this trait.
#[cfg(feature = "hdf5")]
pub trait MpiRng: Rng + SeedableRng + Send + crate::RngCheckpointHdf5 + 'static {}

#[cfg(feature = "hdf5")]
impl<T> MpiRng for T where T: Rng + SeedableRng + Send + crate::RngCheckpointHdf5 + 'static {}

#[cfg(not(feature = "hdf5"))]
pub trait MpiRng: Rng + SeedableRng + Send + 'static {}

#[cfg(not(feature = "hdf5"))]
impl<T> MpiRng for T where T: Rng + SeedableRng + Send + 'static {}

// ============================================================================
// Wire protocol
// ============================================================================

mod tags {
    pub const CONTROLLER_COMMAND: i32 = 0x4351;
    pub const WORKER_REPORT: i32 = 0x4352;
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct Assignment {
    task_id: usize,
    run_id: u64,
    measurement_sweeps: u64,
    #[serde(default)]
    task_target_sweeps: u64,
    #[serde(default)]
    ranks_per_run: i32,
    thermalization_sweeps: u64,
    binsize: usize,
    base_seed: u64,
    params: Params,
}

#[cfg(feature = "hdf5")]
#[derive(Debug, Clone, Serialize, Deserialize)]
struct CheckpointCommit {
    task_id: usize,
    run_id: u64,
    ranks_per_run: i32,
    sweep_count: u64,
    measurement_sweeps_done: u64,
}

#[cfg(feature = "hdf5")]
#[derive(Debug, Clone, Serialize, Deserialize)]
enum CheckpointDecision {
    Fresh,
    Resume(CheckpointCommit),
    Error(String),
}

#[derive(Debug, Serialize, Deserialize)]
enum ControllerCommand {
    Assign(Assignment),
    Stop,
}

#[derive(Debug, Serialize, Deserialize)]
enum WorkerReport {
    Ready,
    Completed {
        task_id: usize,
        run_id: u64,
        measurement_sweeps: u64,
        results: Results,
    },
    Interrupted {
        task_id: usize,
        run_id: u64,
        measurement_sweeps: u64,
    },
    Failed {
        task_id: Option<usize>,
        run_id: Option<u64>,
        message: String,
    },
}

#[cfg(feature = "mpi")]
fn send_json<T: Serialize>(
    comm: &SimpleCommunicator,
    dest: i32,
    tag: i32,
    value: &T,
) -> Result<(), MpiError> {
    let bytes = serde_json::to_vec(value).map_err(|e| MpiError::Communication(e.to_string()))?;
    comm.process_at_rank(dest)
        .send_with_tag(bytes.as_slice(), tag);
    Ok(())
}

#[cfg(feature = "mpi")]
fn recv_json<T: DeserializeOwned>(
    comm: &SimpleCommunicator,
    source: i32,
    tag: i32,
) -> Result<T, MpiError> {
    let (bytes, _) = comm.process_at_rank(source).receive_vec_with_tag::<u8>(tag);
    serde_json::from_slice(&bytes).map_err(|e| MpiError::Communication(e.to_string()))
}

#[cfg(feature = "mpi")]
fn recv_json_any<T: DeserializeOwned>(
    comm: &SimpleCommunicator,
    tag: i32,
) -> Result<(T, i32), MpiError> {
    let (bytes, status) = comm.any_process().receive_vec_with_tag::<u8>(tag);
    let value =
        serde_json::from_slice(&bytes).map_err(|e| MpiError::Communication(e.to_string()))?;
    Ok((value, status.source_rank()))
}

#[cfg(feature = "mpi")]
fn broadcast_json<T>(
    comm: &SimpleCommunicator,
    root_rank: i32,
    value_on_root: Option<&T>,
) -> Result<T, MpiError>
where
    T: Serialize + DeserializeOwned,
{
    let root = comm.process_at_rank(root_rank);
    let mut bytes = if comm.rank() == root_rank {
        serde_json::to_vec(value_on_root.ok_or_else(|| {
            MpiError::Communication("broadcast root did not provide a value".into())
        })?)
        .map_err(|e| MpiError::Communication(e.to_string()))?
    } else {
        Vec::new()
    };

    let mut len = bytes.len() as u64;
    root.broadcast_into(&mut len);
    if comm.rank() != root_rank {
        bytes.resize(len as usize, 0);
    }
    if len > 0 {
        root.broadcast_into(bytes.as_mut_slice());
    }
    serde_json::from_slice(&bytes).map_err(|e| MpiError::Communication(e.to_string()))
}

// ============================================================================
// Tasks and scheduling
// ============================================================================

/// A parameter task and its aggregate target measurement sweeps.
#[derive(Debug, Clone)]
pub struct TaskSpec {
    pub id: usize,
    pub target_sweeps: u64,
    pub thermalization: u64,
    pub params: Params,
}

#[derive(Debug, Clone)]
struct TaskInfo {
    spec: TaskSpec,
    completed: u64,
    reserved: u64,
    next_run_id: u64,
    pending: VecDeque<Assignment>,
}

/// Pure scheduler state used by the rank-0 controller.
///
/// `next()` reserves work immediately, preventing multiple worker groups from
/// being assigned overlapping sweep budgets.
pub struct TaskStream {
    tasks: Vec<TaskInfo>,
    index_by_id: HashMap<usize, usize>,
    cursor: usize,
    chunk_divisor: u64,
}

impl TaskStream {
    pub fn new(specs: Vec<TaskSpec>) -> Self {
        Self::with_parallelism(specs, 1).expect("TaskStream::new received duplicate task IDs")
    }

    pub fn with_parallelism(
        specs: Vec<TaskSpec>,
        parallel_groups: i32,
    ) -> Result<Self, CarloError> {
        let mut index_by_id = HashMap::new();
        let mut tasks = Vec::with_capacity(specs.len());
        for spec in specs {
            let index = tasks.len();
            if index_by_id.insert(spec.id, index).is_some() {
                return Err(CarloError::InvalidConfig {
                    field: "tasks".into(),
                    reason: format!("duplicate MPI task id {}", spec.id),
                });
            }
            tasks.push(TaskInfo {
                spec,
                completed: 0,
                reserved: 0,
                next_run_id: 1,
                pending: VecDeque::new(),
            });
        }
        Ok(Self {
            tasks,
            index_by_id,
            cursor: 0,
            chunk_divisor: parallel_groups.max(1) as u64,
        })
    }

    fn task_mut(&mut self, task_id: usize) -> Option<&mut TaskInfo> {
        let index = *self.index_by_id.get(&task_id)?;
        self.tasks.get_mut(index)
    }

    fn task(&self, task_id: usize) -> Option<&TaskInfo> {
        let index = *self.index_by_id.get(&task_id)?;
        self.tasks.get(index)
    }

    fn assignment_for(&mut self, index: usize, base: &RunConfig) -> Option<Assignment> {
        let task = self.tasks.get_mut(index)?;
        if let Some(pending) = task.pending.pop_front() {
            return Some(pending);
        }

        let unavailable = task.completed.saturating_add(task.reserved);
        if unavailable >= task.spec.target_sweeps {
            return None;
        }
        let remaining = task.spec.target_sweeps - unavailable;
        let nominal = task
            .spec
            .target_sweeps
            .div_ceil(self.chunk_divisor)
            .max(base.binsize as u64)
            .max(1);
        let budget = remaining.min(nominal);
        let run_id = task.next_run_id;
        task.next_run_id = task.next_run_id.saturating_add(1);
        task.reserved = task.reserved.saturating_add(budget);

        Some(Assignment {
            task_id: task.spec.id,
            run_id,
            measurement_sweeps: budget,
            task_target_sweeps: task.spec.target_sweeps,
            ranks_per_run: 0,
            thermalization_sweeps: task.spec.thermalization,
            binsize: base.binsize,
            base_seed: base.base_seed,
            params: task.spec.params.clone(),
        })
    }

    fn next_assignment(&mut self, base: &RunConfig) -> Option<Assignment> {
        if self.tasks.is_empty() {
            return None;
        }
        for offset in 0..self.tasks.len() {
            let index = (self.cursor + offset) % self.tasks.len();
            if let Some(assignment) = self.assignment_for(index, base) {
                self.cursor = (index + 1) % self.tasks.len();
                return Some(assignment);
            }
        }
        None
    }

    /// Compatibility helper returning the next task with unreserved work.
    pub fn next(&mut self, _num_workers: i32) -> Option<&TaskSpec> {
        if self.tasks.is_empty() {
            return None;
        }
        for offset in 0..self.tasks.len() {
            let index = (self.cursor + offset) % self.tasks.len();
            let task = &self.tasks[index];
            if task.completed.saturating_add(task.reserved) < task.spec.target_sweeps
                || !task.pending.is_empty()
            {
                self.cursor = (index + 1) % self.tasks.len();
                return Some(&self.tasks[index].spec);
            }
        }
        None
    }

    pub fn report_progress(&mut self, task_id: usize, sweeps: u64) {
        if let Some(task) = self.task_mut(task_id) {
            task.completed = task.completed.saturating_add(sweeps);
            task.reserved = task.reserved.saturating_sub(sweeps);
        }
    }

    pub fn start_run(&mut self, task_id: usize) {
        if let Some(task) = self.task_mut(task_id) {
            let remaining = task
                .spec
                .target_sweeps
                .saturating_sub(task.completed.saturating_add(task.reserved));
            task.reserved = task.reserved.saturating_add(remaining.min(1));
        }
    }

    pub fn complete_run(&mut self, _task_id: usize) {}

    pub fn sweeps_hint(&self, task_id: usize, num_active: i32) -> u64 {
        self.task(task_id)
            .map(|task| {
                let remaining = task
                    .spec
                    .target_sweeps
                    .saturating_sub(task.completed.saturating_add(task.reserved));
                remaining.div_ceil(num_active.max(1) as u64)
            })
            .unwrap_or(0)
    }

    pub fn is_done(&self, task_id: usize) -> bool {
        self.task(task_id)
            .map(|task| task.completed >= task.spec.target_sweeps)
            .unwrap_or(true)
    }

    fn complete_assignment(&mut self, task_id: usize, budget: u64) -> Result<(), MpiError> {
        let task = self.task_mut(task_id).ok_or_else(|| {
            MpiError::InvalidTransition(format!("completion for unknown task {task_id}"))
        })?;
        task.reserved = task.reserved.saturating_sub(budget);
        task.completed = task.completed.saturating_add(budget);
        Ok(())
    }

    fn interrupt_assignment(&mut self, task_id: usize, budget: u64) -> Result<(), MpiError> {
        let task = self.task_mut(task_id).ok_or_else(|| {
            MpiError::InvalidTransition(format!("interruption for unknown task {task_id}"))
        })?;
        task.reserved = task.reserved.saturating_sub(budget);
        Ok(())
    }

    fn restore_assignment(&mut self, assignment: Assignment) -> Result<(), CarloError> {
        let task = self
            .task_mut(assignment.task_id)
            .ok_or_else(|| CarloError::InvalidConfig {
                field: "job_dir".into(),
                reason: format!(
                    "persisted MPI assignment references unknown task {}",
                    assignment.task_id
                ),
            })?;
        task.next_run_id = task.next_run_id.max(assignment.run_id.saturating_add(1));
        task.reserved = task.reserved.saturating_add(assignment.measurement_sweeps);
        task.pending.push_back(assignment);
        Ok(())
    }

    fn restore_completed(&mut self, task_id: usize, run_id: u64, sweeps: u64) {
        if let Some(task) = self.task_mut(task_id) {
            task.completed = task.completed.saturating_add(sweeps);
            task.next_run_id = task.next_run_id.max(run_id.saturating_add(1));
        }
    }

    fn all_done(&self) -> bool {
        self.tasks
            .iter()
            .all(|task| task.completed >= task.spec.target_sweeps)
    }
}

// ============================================================================
// Results aggregation
// ============================================================================

/// Aggregates independent run results, optionally separated by task ID.
pub struct ResultsAggregator {
    by_task: HashMap<usize, Vec<Results>>,
}

impl ResultsAggregator {
    pub fn new() -> Self {
        Self {
            by_task: HashMap::new(),
        }
    }

    /// Compatibility method: add to task zero.
    pub fn add(&mut self, results: &Results) {
        self.add_for_task(0, results.clone());
    }

    pub fn add_for_task(&mut self, task_id: usize, results: Results) {
        self.by_task.entry(task_id).or_default().push(results);
    }

    pub fn finalize(self) -> Results {
        let all: Vec<Results> = self.by_task.into_values().flatten().collect();
        Results::merge(&all)
    }

    pub fn finalize_ordered(mut self, tasks: &[TaskSpec]) -> Vec<Results> {
        tasks
            .iter()
            .map(|task| {
                let values = self.by_task.remove(&task.id).unwrap_or_default();
                Results::merge(&values)
            })
            .collect()
    }
}

impl Default for ResultsAggregator {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================================
// Type-state compatibility facade
// ============================================================================

pub trait WorkerState: private::Sealed {}
#[derive(Debug, Clone, Copy)]
pub struct Idle;
#[derive(Debug, Clone, Copy)]
pub struct Running {
    pub task_id: usize,
    pub run_id: u64,
}
#[derive(Debug, Clone, Copy)]
pub struct Done;

mod private {
    pub trait Sealed {}
    impl Sealed for super::Idle {}
    impl Sealed for super::Running {}
    impl Sealed for super::Done {}
}

impl WorkerState for Idle {}
impl WorkerState for Running {}
impl WorkerState for Done {}

/// Lightweight worker state marker retained for source compatibility.
/// Communicators are deliberately not stored because rsmpi communicators are
/// neither `Send` nor `Sync` and must not be bitwise copied.
pub struct Worker<S: WorkerState> {
    world_rank: i32,
    state: S,
}

impl Worker<Idle> {
    pub fn new(world_rank: i32) -> Self {
        Self {
            world_rank,
            state: Idle,
        }
    }

    pub fn assign(self, task_id: usize, run_id: u64) -> Worker<Running> {
        Worker {
            world_rank: self.world_rank,
            state: Running { task_id, run_id },
        }
    }

    pub fn finish(self) -> Worker<Done> {
        Worker {
            world_rank: self.world_rank,
            state: Done,
        }
    }
}

impl Worker<Running> {
    pub fn task_id(&self) -> usize {
        self.state.task_id
    }
    pub fn run_id(&self) -> u64 {
        self.state.run_id
    }
    pub fn finish(self) -> Worker<Done> {
        Worker {
            world_rank: self.world_rank,
            state: Done,
        }
    }
}

impl Worker<Done> {
    pub fn reset(self) -> Worker<Idle> {
        Worker {
            world_rank: self.world_rank,
            state: Idle,
        }
    }
}

// ============================================================================
// Public configuration
// ============================================================================

#[derive(Debug, Clone)]
pub struct DistributedConfig {
    pub run_config: RunConfig,
    pub ranks_per_run: i32,
    pub run_time: Option<Duration>,
    pub checkpoint_time: Option<Duration>,
    pub job_dir: PathBuf,
    pub tasks: Vec<TaskSpec>,
}

impl Default for DistributedConfig {
    fn default() -> Self {
        Self {
            run_config: RunConfig::default(),
            ranks_per_run: 1,
            run_time: None,
            checkpoint_time: None,
            job_dir: PathBuf::from("."),
            tasks: Vec::new(),
        }
    }
}

fn validate_config(config: &DistributedConfig, world_size: i32) -> Result<(), CarloError> {
    if world_size < 2 {
        return Err(CarloError::InvalidConfig {
            field: "mpi".into(),
            reason: "distributed execution requires rank 0 plus at least one worker rank".into(),
        });
    }
    if config.ranks_per_run <= 0 {
        return Err(CarloError::InvalidConfig {
            field: "ranks_per_run".into(),
            reason: "must be positive".into(),
        });
    }
    let worker_ranks = world_size - 1;
    if worker_ranks % config.ranks_per_run != 0 {
        return Err(CarloError::InvalidConfig {
            field: "ranks_per_run".into(),
            reason: format!(
                "{worker_ranks} worker ranks are not divisible by {}",
                config.ranks_per_run
            ),
        });
    }
    if config.tasks.is_empty() {
        return Err(CarloError::InvalidConfig {
            field: "tasks".into(),
            reason: "at least one distributed task is required".into(),
        });
    }
    if config.run_config.binsize == 0 {
        return Err(CarloError::InvalidConfig {
            field: "binsize".into(),
            reason: "must be positive".into(),
        });
    }
    if (config.checkpoint_time.is_some() || config.run_config.checkpoint_interval > 0)
        && !cfg!(feature = "hdf5")
    {
        return Err(CarloError::InvalidConfig {
            field: "checkpoint".into(),
            reason: "distributed checkpointing requires the hdf5 feature".into(),
        });
    }
    Ok(())
}

// ============================================================================
// Persistent scheduler files
// ============================================================================

fn run_dir(job_dir: &Path, task_id: usize, run_id: u64) -> PathBuf {
    job_dir.join(format!("task_{task_id:04}/run{run_id:04}"))
}

fn assignment_path(job_dir: &Path, task_id: usize, run_id: u64) -> PathBuf {
    run_dir(job_dir, task_id, run_id).join("mpi-assignment.json")
}

fn result_path(job_dir: &Path, task_id: usize, run_id: u64) -> PathBuf {
    run_dir(job_dir, task_id, run_id).join("result.json")
}

#[cfg(feature = "hdf5")]
fn checkpoint_path(
    job_dir: &Path,
    task_id: usize,
    run_id: u64,
    rank_in_run: i32,
    ranks_per_run: i32,
) -> PathBuf {
    let dir = run_dir(job_dir, task_id, run_id);
    if ranks_per_run == 1 {
        dir.join(format!("run{run_id:04}.dump.h5"))
    } else {
        dir.join(format!("run{run_id:04}.rank{rank_in_run:04}.dump.h5"))
    }
}

#[cfg(feature = "hdf5")]
fn checkpoint_staging_path(path: &Path) -> PathBuf {
    path.with_extension("next.h5")
}

#[cfg(feature = "hdf5")]
fn checkpoint_commit_path(job_dir: &Path, task_id: usize, run_id: u64) -> PathBuf {
    run_dir(job_dir, task_id, run_id).join("mpi-checkpoint.json")
}

fn write_json_atomic<T: Serialize>(path: &Path, value: &T) -> Result<(), CarloError> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).map_err(|source| CarloError::IoError {
            path: parent.to_path_buf(),
            source,
        })?;
    }
    let tmp = path.with_extension("json.tmp");
    let bytes = serde_json::to_vec_pretty(value)?;
    fs::write(&tmp, bytes).map_err(|source| CarloError::IoError {
        path: tmp.clone(),
        source,
    })?;
    fs::rename(&tmp, path).map_err(|source| CarloError::IoError {
        path: path.to_path_buf(),
        source,
    })?;
    Ok(())
}

fn restore_persisted_runs(
    stream: &mut TaskStream,
    aggregator: &mut ResultsAggregator,
    config: &DistributedConfig,
) -> Result<(), CarloError> {
    for task in &config.tasks {
        let task_dir = config.job_dir.join(format!("task_{:04}", task.id));
        let entries = match fs::read_dir(&task_dir) {
            Ok(entries) => entries,
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => continue,
            Err(source) => {
                return Err(CarloError::IoError {
                    path: task_dir,
                    source,
                })
            }
        };

        for entry in entries {
            let entry = entry.map_err(|source| CarloError::IoError {
                path: task_dir.clone(),
                source,
            })?;
            if !entry
                .file_type()
                .map_err(|source| CarloError::IoError {
                    path: entry.path(),
                    source,
                })?
                .is_dir()
            {
                continue;
            }
            let assignment_file = entry.path().join("mpi-assignment.json");
            if !assignment_file.exists() {
                continue;
            }
            let assignment: Assignment =
                serde_json::from_slice(&fs::read(&assignment_file).map_err(|source| {
                    CarloError::IoError {
                        path: assignment_file.clone(),
                        source,
                    }
                })?)?;
            if assignment.task_id != task.id {
                return Err(CarloError::InvalidConfig {
                    field: "job_dir".into(),
                    reason: format!(
                        "persisted assignment in task_{:04} belongs to task {}",
                        task.id, assignment.task_id
                    ),
                });
            }
            if assignment.task_target_sweeps != 0
                && assignment.task_target_sweeps != task.target_sweeps
            {
                return Err(CarloError::InvalidConfig {
                    field: "job_dir".into(),
                    reason: format!(
                        "task {} target changed from {} to {} while restart data exists",
                        task.id, assignment.task_target_sweeps, task.target_sweeps
                    ),
                });
            }
            if assignment.ranks_per_run != 0 && assignment.ranks_per_run != config.ranks_per_run {
                return Err(CarloError::InvalidConfig {
                    field: "ranks_per_run".into(),
                    reason: format!(
                        "task {} was checkpointed with {} ranks per run, not {}",
                        task.id, assignment.ranks_per_run, config.ranks_per_run
                    ),
                });
            }
            if assignment.thermalization_sweeps != task.thermalization
                || assignment.binsize != config.run_config.binsize
                || assignment.base_seed != config.run_config.base_seed
                || assignment.params != task.params
            {
                return Err(CarloError::InvalidConfig {
                    field: "job_dir".into(),
                    reason: format!(
                        "task {} parameters or run configuration changed while restart data exists",
                        task.id
                    ),
                });
            }
            let result_file = entry.path().join("result.json");
            if result_file.exists() {
                let results: Results =
                    serde_json::from_slice(&fs::read(&result_file).map_err(|source| {
                        CarloError::IoError {
                            path: result_file.clone(),
                            source,
                        }
                    })?)?;
                stream.restore_completed(
                    assignment.task_id,
                    assignment.run_id,
                    assignment.measurement_sweeps,
                );
                aggregator.add_for_task(assignment.task_id, results);
            } else {
                stream.restore_assignment(assignment)?;
            }
        }
    }
    Ok(())
}

// ============================================================================
// Generic MPI backend
// ============================================================================

/// Generic task-partitioning MPI backend.
///
/// Every rank must call `spawn_tasks` with the same `n_tasks` and `base_seed`.
/// Logical task `i` is executed by rank `i % world_size`; the final barrier
/// ensures all ranks finish before the method returns.
#[cfg(feature = "mpi")]
struct MpiRuntime {
    universe: Universe,
    // MPI_THREAD_SERIALIZED permits multiple application threads as long as
    // only one of them is inside MPI at a time. Backend methods enforce that
    // rule through this process-local lock.
    call_lock: Mutex<()>,
}

#[cfg(feature = "mpi")]
#[derive(Clone)]
pub struct MpiBackend {
    runtime: Arc<MpiRuntime>,
    rank: i32,
    size: i32,
    ranks_per_run: i32,
}

#[cfg(feature = "mpi")]
impl MpiBackend {
    pub fn new() -> Result<Self, CarloError> {
        Self::with_ranks_per_run(1)
    }

    pub fn with_ranks_per_run(ranks_per_run: i32) -> Result<Self, CarloError> {
        if ranks_per_run <= 0 {
            return Err(CarloError::InvalidConfig {
                field: "ranks_per_run".into(),
                reason: "must be positive".into(),
            });
        }
        let (universe, provided) =
            mpi::initialize_with_threading(Threading::Serialized).ok_or(MpiError::InitFailed)?;
        if provided < Threading::Serialized {
            return Err(MpiError::InvalidTopology(format!(
                "MPI implementation provided {provided:?} thread support; Serialized is required"
            ))
            .into());
        }
        let world = universe.world();
        let rank = world.rank();
        let size = world.size();
        if size > 1 && (size - 1) % ranks_per_run != 0 {
            return Err(CarloError::InvalidConfig {
                field: "ranks_per_run".into(),
                reason: format!(
                    "{} worker ranks are not divisible by {ranks_per_run}",
                    size - 1
                ),
            });
        }
        Ok(Self {
            runtime: Arc::new(MpiRuntime {
                universe,
                call_lock: Mutex::new(()),
            }),
            rank,
            size,
            ranks_per_run,
        })
    }

    pub fn rank(&self) -> i32 {
        self.rank
    }
    pub fn size(&self) -> i32 {
        self.size
    }
    pub fn is_controller(&self) -> bool {
        self.rank == 0
    }
    pub fn is_run_leader(&self) -> bool {
        self.is_controller() || (self.rank > 0 && (self.rank - 1) % self.ranks_per_run == 0)
    }
    pub fn num_workers(&self) -> i32 {
        (self.size - 1).max(0)
    }
    pub fn num_parallel_runs(&self) -> i32 {
        if self.size <= 1 {
            0
        } else {
            self.num_workers() / self.ranks_per_run
        }
    }
    pub fn run_group(&self) -> i32 {
        if self.rank == 0 {
            0
        } else {
            1 + (self.rank - 1) / self.ranks_per_run
        }
    }
    pub fn rank_in_run(&self) -> i32 {
        if self.rank == 0 {
            0
        } else {
            (self.rank - 1) % self.ranks_per_run
        }
    }
    pub fn ranks_per_run(&self) -> i32 {
        self.ranks_per_run
    }
}

#[cfg(feature = "mpi")]
impl Default for MpiBackend {
    fn default() -> Self {
        Self::new().expect("failed to initialize MPI backend")
    }
}

#[cfg(feature = "mpi")]
impl super::Backend for MpiBackend {
    type Rng = rand_xoshiro::Xoshiro256PlusPlus;

    fn spawn_tasks<F>(&self, n_tasks: usize, base_seed: u64, f: F)
    where
        F: Fn(usize, &mut Self::Rng) + Sync,
    {
        let _guard = self
            .runtime
            .call_lock
            .lock()
            .expect("MPI backend call lock was poisoned");
        let world = self.runtime.universe.world();
        for task_id in 0..n_tasks {
            if task_id % self.size as usize != self.rank as usize {
                continue;
            }
            let mut rng: Self::Rng = RngStreamKey::new(base_seed)
                .with_task(task_id as u64)
                .with_replica(self.rank as u64)
                .with_phase(RngPhase::BackendTask)
                .seeded();
            f(task_id, &mut rng);
        }
        world.barrier();
    }

    fn barrier(&self) {
        let _guard = self
            .runtime
            .call_lock
            .lock()
            .expect("MPI backend call lock was poisoned");
        self.runtime.universe.world().barrier();
    }
}

#[cfg(not(feature = "mpi"))]
#[derive(Debug, Clone, Copy)]
pub struct MpiBackend;

#[cfg(not(feature = "mpi"))]
impl MpiBackend {
    pub fn new() -> Result<Self, CarloError> {
        Err(CarloError::InvalidConfig {
            field: "mpi".into(),
            reason: "MPI feature not enabled".into(),
        })
    }

    pub fn with_ranks_per_run(_ranks_per_run: i32) -> Result<Self, CarloError> {
        Self::new()
    }
}

// ============================================================================
// Distributed execution
// ============================================================================

#[cfg(feature = "mpi")]
#[derive(Debug)]
enum GroupOutcome {
    Completed(Results),
    Interrupted,
}

/// Run the distributed controller/worker scheduler.
///
/// Return value semantics:
/// - rank 0 returns one merged [`Results`] value per input task, preserving
///   `config.tasks` order;
/// - worker ranks return an empty vector.
#[cfg(feature = "mpi")]
pub fn run_distributed<MC, R>(config: DistributedConfig) -> Result<Vec<Results>, CarloError>
where
    MC: MonteCarlo<Rng = R> + FromParams<Rng = R>,
    R: MpiRng,
{
    let universe = mpi::initialize().ok_or(MpiError::InitFailed)?;
    let world = universe.world();
    validate_config(&config, world.size())?;
    let job_started = Instant::now();

    let world_rank = world.rank();
    let group_color = if world_rank == 0 {
        Color::undefined()
    } else {
        Color::with_value((world_rank - 1) / config.ranks_per_run)
    };
    let run_comm = world.split_by_color_with_key(group_color, world_rank);

    let is_group_leader = world_rank > 0 && (world_rank - 1) % config.ranks_per_run == 0;
    let leader_color = if world_rank == 0 || is_group_leader {
        Color::with_value(0)
    } else {
        Color::undefined()
    };
    let leader_comm = world.split_by_color_with_key(leader_color, world_rank);

    let result = if world_rank == 0 {
        let leaders = leader_comm.ok_or_else(|| {
            MpiError::InvalidTopology("controller was excluded from leader communicator".into())
        })?;
        run_controller(&leaders, &config)
    } else {
        let group = run_comm.ok_or_else(|| {
            MpiError::InvalidTopology("worker was excluded from run communicator".into())
        })?;
        if is_group_leader {
            let leaders = leader_comm.ok_or_else(|| {
                MpiError::InvalidTopology(
                    "run-group leader was excluded from leader communicator".into(),
                )
            })?;
            run_group_leader::<MC, R>(&leaders, &group, &config, job_started)
        } else {
            run_group_follower::<MC, R>(&group, &config, job_started)
        }
    };

    // Every rank reaches this barrier after its local role has terminated.
    world.barrier();
    result
}

#[cfg(feature = "mpi")]
fn run_controller(
    leader_comm: &SimpleCommunicator,
    config: &DistributedConfig,
) -> Result<Vec<Results>, CarloError> {
    let group_count = leader_comm.size() - 1;
    let initialized = (|| {
        let mut stream = TaskStream::with_parallelism(config.tasks.clone(), group_count)?;
        let mut aggregator = ResultsAggregator::new();
        fs::create_dir_all(&config.job_dir).map_err(|source| CarloError::IoError {
            path: config.job_dir.clone(),
            source,
        })?;
        restore_persisted_runs(&mut stream, &mut aggregator, config)?;
        Ok::<_, CarloError>((stream, aggregator))
    })();

    let (mut stream, mut aggregator) = match initialized {
        Ok(state) => state,
        Err(error) => {
            // Worker groups have already sent Ready and are blocked waiting for
            // a command. Release every group before returning the local error.
            stop_ready_groups(leader_comm, group_count)?;
            return Err(error);
        }
    };

    let mut live_groups = group_count;
    let mut stopping = stream.all_done();
    let mut first_error: Option<String> = None;

    while live_groups > 0 {
        let (report, source) = recv_json_any::<WorkerReport>(leader_comm, tags::WORKER_REPORT)?;

        let report_result: Result<(), CarloError> = match report {
            WorkerReport::Ready => Ok(()),
            WorkerReport::Completed {
                task_id,
                run_id,
                measurement_sweeps,
                results,
            } => (|| {
                stream.complete_assignment(task_id, measurement_sweeps)?;
                write_json_atomic(&result_path(&config.job_dir, task_id, run_id), &results)?;
                #[cfg(feature = "hdf5")]
                cleanup_checkpoint_artifacts(config, task_id, run_id);
                aggregator.add_for_task(task_id, results);
                Ok(())
            })(),
            WorkerReport::Interrupted {
                task_id,
                run_id: _,
                measurement_sweeps,
            } => {
                stopping = true;
                stream
                    .interrupt_assignment(task_id, measurement_sweeps)
                    .map_err(CarloError::from)
            }
            WorkerReport::Failed {
                task_id: _,
                run_id: _,
                message,
            } => {
                stopping = true;
                Err(MpiError::Worker(message).into())
            }
        };

        if let Err(error) = report_result {
            if first_error.is_none() {
                first_error = Some(error.to_string());
            }
            stopping = true;
        }

        let command = if stopping {
            ControllerCommand::Stop
        } else if let Some(mut assignment) = stream.next_assignment(&config.run_config) {
            assignment.ranks_per_run = config.ranks_per_run;
            match write_json_atomic(
                &assignment_path(&config.job_dir, assignment.task_id, assignment.run_id),
                &assignment,
            ) {
                Ok(()) => ControllerCommand::Assign(assignment),
                Err(error) => {
                    if first_error.is_none() {
                        first_error = Some(error.to_string());
                    }
                    stopping = true;
                    ControllerCommand::Stop
                }
            }
        } else {
            ControllerCommand::Stop
        };

        if matches!(command, ControllerCommand::Stop) {
            live_groups -= 1;
        }
        send_json(leader_comm, source, tags::CONTROLLER_COMMAND, &command)?;
    }

    if let Some(message) = first_error {
        return Err(MpiError::Worker(message).into());
    }
    Ok(aggregator.finalize_ordered(&config.tasks))
}

#[cfg(feature = "mpi")]
fn stop_ready_groups(leader_comm: &SimpleCommunicator, group_count: i32) -> Result<(), CarloError> {
    for _ in 0..group_count {
        let (_report, source) = recv_json_any::<WorkerReport>(leader_comm, tags::WORKER_REPORT)?;
        send_json(
            leader_comm,
            source,
            tags::CONTROLLER_COMMAND,
            &ControllerCommand::Stop,
        )?;
    }
    Ok(())
}

#[cfg(feature = "mpi")]
fn run_group_leader<MC, R>(
    leader_comm: &SimpleCommunicator,
    run_comm: &SimpleCommunicator,
    config: &DistributedConfig,
    job_started: Instant,
) -> Result<Vec<Results>, CarloError>
where
    MC: MonteCarlo<Rng = R> + FromParams<Rng = R>,
    R: MpiRng,
{
    send_json(leader_comm, 0, tags::WORKER_REPORT, &WorkerReport::Ready)?;

    loop {
        let command: ControllerCommand = recv_json(leader_comm, 0, tags::CONTROLLER_COMMAND)?;
        let command: ControllerCommand = broadcast_json(run_comm, 0, Some(&command))?;
        match command {
            ControllerCommand::Stop => return Ok(Vec::new()),
            ControllerCommand::Assign(assignment) => {
                let task_id = assignment.task_id;
                let run_id = assignment.run_id;
                let budget = assignment.measurement_sweeps;
                match execute_group_assignment::<MC, R>(run_comm, config, &assignment, job_started)
                {
                    Ok(GroupOutcome::Completed(results)) => {
                        let report = WorkerReport::Completed {
                            task_id,
                            run_id,
                            measurement_sweeps: budget,
                            results,
                        };
                        if let Err(error) = send_json(leader_comm, 0, tags::WORKER_REPORT, &report)
                        {
                            // Results can contain non-finite model output that
                            // JSON refuses to encode. Keep the protocol aligned
                            // by reporting a small serializable failure instead.
                            let fallback = WorkerReport::Failed {
                                task_id: Some(task_id),
                                run_id: Some(run_id),
                                message: format!("failed to serialize worker results: {error}"),
                            };
                            send_json(leader_comm, 0, tags::WORKER_REPORT, &fallback)?;
                        }
                    }
                    Ok(GroupOutcome::Interrupted) => {
                        let report = WorkerReport::Interrupted {
                            task_id,
                            run_id,
                            measurement_sweeps: budget,
                        };
                        send_json(leader_comm, 0, tags::WORKER_REPORT, &report)?;
                    }
                    Err(error) => {
                        let report = WorkerReport::Failed {
                            task_id: Some(task_id),
                            run_id: Some(run_id),
                            message: error.to_string(),
                        };
                        send_json(leader_comm, 0, tags::WORKER_REPORT, &report)?;
                    }
                }
            }
        }
    }
}

#[cfg(feature = "mpi")]
fn run_group_follower<MC, R>(
    run_comm: &SimpleCommunicator,
    config: &DistributedConfig,
    job_started: Instant,
) -> Result<Vec<Results>, CarloError>
where
    MC: MonteCarlo<Rng = R> + FromParams<Rng = R>,
    R: MpiRng,
{
    loop {
        let command: ControllerCommand = broadcast_json(run_comm, 0, None)?;
        match command {
            ControllerCommand::Stop => return Ok(Vec::new()),
            ControllerCommand::Assign(assignment) => {
                // All ranks return the same success/error decision from
                // execute_group_assignment, so the next broadcast stays aligned.
                let _ =
                    execute_group_assignment::<MC, R>(run_comm, config, &assignment, job_started);
            }
        }
    }
}

#[cfg(feature = "mpi")]
fn execute_group_assignment<MC, R>(
    run_comm: &SimpleCommunicator,
    config: &DistributedConfig,
    assignment: &Assignment,
    job_started: Instant,
) -> Result<GroupOutcome, CarloError>
where
    MC: MonteCarlo<Rng = R> + FromParams<Rng = R>,
    R: MpiRng,
{
    let rank_in_run = run_comm.rank();
    let run_config = RunConfig {
        thermalization_sweeps: assignment.thermalization_sweeps,
        measurement_sweeps: assignment.measurement_sweeps,
        binsize: assignment.binsize,
        base_seed: assignment.base_seed,
        progress_interval: config.run_config.progress_interval,
        checkpoint_interval: config.run_config.checkpoint_interval,
    };
    let seed = RngStreamKey::new(assignment.base_seed)
        .with_task(assignment.task_id as u64)
        .with_run(assignment.run_id)
        .with_replica(rank_in_run as u64)
        .with_phase(RngPhase::Initialization)
        .seed();

    #[cfg(feature = "hdf5")]
    let checkpoint = checkpoint_path(
        &config.job_dir,
        assignment.task_id,
        assignment.run_id,
        rank_in_run,
        config.ranks_per_run,
    );
    #[cfg(feature = "hdf5")]
    let checkpoint_staging = checkpoint_staging_path(&checkpoint);
    #[cfg(feature = "hdf5")]
    let checkpoint_commit =
        checkpoint_commit_path(&config.job_dir, assignment.task_id, assignment.run_id);

    #[cfg(feature = "hdf5")]
    let checkpoint_decision: CheckpointDecision = {
        let root_decision = if rank_in_run == 0 {
            Some(read_checkpoint_decision(
                &checkpoint_commit,
                assignment,
                config.ranks_per_run,
            ))
        } else {
            None
        };
        broadcast_json(run_comm, 0, root_decision.as_ref())?
    };

    #[cfg(feature = "hdf5")]
    let run_result: Result<Run<MC, R>, CarloError> = match &checkpoint_decision {
        CheckpointDecision::Error(message) => Err(CarloError::CheckpointCorrupted {
            detail: message.clone(),
        }),
        CheckpointDecision::Fresh => (|| {
            remove_if_exists(&checkpoint)?;
            remove_if_exists(&checkpoint_staging)?;
            Run::new(
                &assignment.params,
                TaskId::new(assignment.task_id),
                RunId::new(assignment.run_id),
                &run_config,
                seed,
            )
        })(),
        CheckpointDecision::Resume(commit) => {
            match Run::<MC, R>::read_checkpoint(&checkpoint, &assignment.params, &run_config, seed)
            {
                Ok(Some(run)) => {
                    if run.sweep_count() != commit.sweep_count
                        || run.sweeps_done() != commit.measurement_sweeps_done
                    {
                        Err(CarloError::CheckpointCorrupted {
                            detail: format!(
                                "rank {rank_in_run} checkpoint counters ({}, {}) do not match committed counters ({}, {})",
                                run.sweep_count(),
                                run.sweeps_done(),
                                commit.sweep_count,
                                commit.measurement_sweeps_done
                            ),
                        })
                    } else {
                        Ok(run)
                    }
                }
                Ok(None) => Err(CarloError::CheckpointCorrupted {
                    detail: format!(
                        "committed checkpoint is missing for rank {rank_in_run}: {}",
                        checkpoint.display()
                    ),
                }),
                Err(error) => Err(error),
            }
        }
    };

    #[cfg(not(feature = "hdf5"))]
    let run_result: Result<Run<MC, R>, CarloError> = Run::new(
        &assignment.params,
        TaskId::new(assignment.task_id),
        RunId::new(assignment.run_id),
        &run_config,
        seed,
    );

    let run_result = run_result.and_then(|run| {
        if run.task_id().as_usize() != assignment.task_id
            || run.run_id().as_u64() != assignment.run_id
            || run.target_sweeps() != assignment.measurement_sweeps
        {
            return Err(CarloError::CheckpointCorrupted {
                detail: format!(
                    "checkpoint identity/target ({}, {}, {}) does not match assignment ({}, {}, {})",
                    run.task_id().as_usize(),
                    run.run_id().as_u64(),
                    run.target_sweeps(),
                    assignment.task_id,
                    assignment.run_id,
                    assignment.measurement_sweeps
                ),
            });
        }
        Ok(run)
    });

    // Construction/checkpoint errors must be turned into a group-wide decision
    // before any rank enters model collectives.
    let local_init_ok = if run_result.is_ok() { 1i32 } else { 0i32 };
    let mut init_ok = vec![0i32; run_comm.size() as usize];
    run_comm.all_gather_into(&local_init_ok, init_ok.as_mut_slice());
    if init_ok.contains(&0) {
        return match run_result {
            Err(error) => Err(error),
            Ok(_) => Err(MpiError::Worker(
                "another rank failed while constructing or restoring the distributed run".into(),
            )
            .into()),
        };
    }
    let mut run = run_result?;

    #[allow(unused_mut)]
    let mut limits = TimeLimits::from_start(config.checkpoint_time, config.run_time, job_started);
    #[allow(unused_mut)]
    let mut last_checkpoint_sweep = run.sweep_count();
    let poll_sweeps = config.run_config.progress_interval.clamp(1, 10_000);

    loop {
        for _ in 0..poll_sweeps {
            if run.is_complete() {
                break;
            }
            if run_comm.size() == 1 {
                run.step();
            } else {
                run.step_with_comm(run_comm);
            }
        }

        let mut control = [0u64; 2]; // action: 0=continue, 1=complete, 2=interrupt
        if rank_in_run == 0 {
            control[0] = if run.is_complete() {
                1
            } else if limits.should_finish() {
                2
            } else {
                0
            };
            let sweep_checkpoint_due = config.run_config.checkpoint_interval > 0
                && run.sweep_count().saturating_sub(last_checkpoint_sweep)
                    >= config.run_config.checkpoint_interval;
            control[1] = if limits.should_checkpoint() || sweep_checkpoint_due || control[0] == 2 {
                1
            } else {
                0
            };
        }
        run_comm
            .process_at_rank(0)
            .broadcast_into(control.as_mut_slice());

        if control[1] != 0 && !run.is_complete() {
            #[cfg(feature = "hdf5")]
            {
                write_group_checkpoint(
                    run_comm,
                    &mut run,
                    &checkpoint,
                    &checkpoint_staging,
                    &checkpoint_commit,
                    assignment,
                    config.ranks_per_run,
                )?;
                limits.reset_checkpoint();
                last_checkpoint_sweep = run.sweep_count();
            }
        }

        match control[0] {
            0 => continue,
            1 => {
                let results = run.finalize(assignment.base_seed);
                return if rank_in_run == 0 {
                    Ok(GroupOutcome::Completed(results))
                } else {
                    // Followers discard rank-local results. A genuinely
                    // distributed model should reduce global observables in
                    // measure_with_comm before they are recorded on rank 0.
                    Ok(GroupOutcome::Completed(Results::new()))
                };
            }
            2 => return Ok(GroupOutcome::Interrupted),
            other => {
                return Err(MpiError::InvalidTransition(format!(
                    "unknown run-group control action {other}"
                ))
                .into())
            }
        }
    }
}

#[cfg(feature = "hdf5")]
fn read_checkpoint_decision(
    commit_path: &Path,
    assignment: &Assignment,
    ranks_per_run: i32,
) -> CheckpointDecision {
    if !commit_path.exists() {
        return CheckpointDecision::Fresh;
    }
    let result = (|| {
        let bytes = fs::read(commit_path).map_err(|error| error.to_string())?;
        let commit: CheckpointCommit =
            serde_json::from_slice(&bytes).map_err(|error| error.to_string())?;
        if commit.task_id != assignment.task_id
            || commit.run_id != assignment.run_id
            || commit.ranks_per_run != ranks_per_run
        {
            return Err(format!(
                "checkpoint commit {:?} does not match task {}, run {}, ranks_per_run {}",
                commit, assignment.task_id, assignment.run_id, ranks_per_run
            ));
        }
        Ok(commit)
    })();
    match result {
        Ok(commit) => CheckpointDecision::Resume(commit),
        Err(message) => CheckpointDecision::Error(message),
    }
}

#[cfg(feature = "hdf5")]
fn remove_if_exists(path: &Path) -> Result<(), CarloError> {
    match fs::remove_file(path) {
        Ok(()) => Ok(()),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(source) => Err(CarloError::IoError {
            path: path.to_path_buf(),
            source,
        }),
    }
}

#[cfg(feature = "hdf5")]
fn cleanup_checkpoint_artifacts(config: &DistributedConfig, task_id: usize, run_id: u64) {
    let _ = remove_if_exists(&checkpoint_commit_path(&config.job_dir, task_id, run_id));
    for rank_in_run in 0..config.ranks_per_run {
        let path = checkpoint_path(
            &config.job_dir,
            task_id,
            run_id,
            rank_in_run,
            config.ranks_per_run,
        );
        let _ = remove_if_exists(&path);
        let _ = remove_if_exists(&checkpoint_staging_path(&path));
    }
}

#[cfg(feature = "hdf5")]
fn write_group_checkpoint<MC, R>(
    run_comm: &SimpleCommunicator,
    run: &mut Run<MC, R>,
    checkpoint: &Path,
    staging: &Path,
    commit_path: &Path,
    assignment: &Assignment,
    ranks_per_run: i32,
) -> Result<(), CarloError>
where
    MC: MonteCarlo<Rng = R> + FromParams<Rng = R>,
    R: MpiRng,
{
    let rank_in_run = run_comm.rank();

    // Invalidate the previous generation before any final file is replaced.
    // A crash from this point until the new marker is written causes a clean
    // restart from the assignment, never a mixed-rank resume.
    let mut marker_removed = 1i32;
    let mut marker_error = None;
    if rank_in_run == 0 {
        if let Err(error) = remove_if_exists(commit_path) {
            marker_removed = 0;
            marker_error = Some(error);
        }
    }
    run_comm
        .process_at_rank(0)
        .broadcast_into(&mut marker_removed);
    if marker_removed == 0 {
        return match marker_error {
            Some(error) => Err(error),
            None => Err(CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: "run-group leader could not invalidate the previous checkpoint".into(),
            }),
        };
    }

    let local_write: Result<(), CarloError> = (|| {
        if let Some(parent) = staging.parent() {
            fs::create_dir_all(parent).map_err(|source| CarloError::IoError {
                path: parent.to_path_buf(),
                source,
            })?;
        }
        remove_if_exists(staging)?;
        run.write_checkpoint(staging)
    })();
    let local_write_ok = if local_write.is_ok() { 1i32 } else { 0i32 };
    let mut write_ok = vec![0i32; run_comm.size() as usize];
    run_comm.all_gather_into(&local_write_ok, write_ok.as_mut_slice());
    if write_ok.contains(&0) {
        let _ = remove_if_exists(staging);
        return match local_write {
            Err(error) => Err(error),
            Ok(()) => Err(CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: "another MPI rank failed to stage its checkpoint".into(),
            }),
        };
    }

    let local_publish: Result<(), CarloError> = (|| {
        remove_if_exists(checkpoint)?;
        fs::rename(staging, checkpoint).map_err(|source| CarloError::IoError {
            path: checkpoint.to_path_buf(),
            source,
        })
    })();
    let local_publish_ok = if local_publish.is_ok() { 1i32 } else { 0i32 };
    let mut publish_ok = vec![0i32; run_comm.size() as usize];
    run_comm.all_gather_into(&local_publish_ok, publish_ok.as_mut_slice());
    if publish_ok.contains(&0) {
        return match local_publish {
            Err(error) => Err(error),
            Ok(()) => Err(CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: "another MPI rank failed to publish its checkpoint".into(),
            }),
        };
    }

    let mut commit_written = 1i32;
    let mut commit_error = None;
    if rank_in_run == 0 {
        let commit = CheckpointCommit {
            task_id: assignment.task_id,
            run_id: assignment.run_id,
            ranks_per_run,
            sweep_count: run.sweep_count(),
            measurement_sweeps_done: run.sweeps_done(),
        };
        if let Err(error) = write_json_atomic(commit_path, &commit) {
            commit_written = 0;
            commit_error = Some(error);
        }
    }
    run_comm
        .process_at_rank(0)
        .broadcast_into(&mut commit_written);
    if commit_written == 0 {
        return match commit_error {
            Some(error) => Err(error),
            None => Err(CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: "run-group leader could not commit the checkpoint generation".into(),
            }),
        };
    }

    Ok(())
}

// ============================================================================
// Backward-compatible configuration and entry point
// ============================================================================

#[derive(Debug, Clone)]
pub struct MpiRunConfig {
    pub run_config: RunConfig,
    pub ranks_per_run: i32,
    pub run_time: Option<Duration>,
    pub checkpoint_time: Option<Duration>,
    pub job_dir: PathBuf,
    pub tasks: Vec<Params>,
}

impl Default for MpiRunConfig {
    fn default() -> Self {
        Self {
            run_config: RunConfig::default(),
            ranks_per_run: 1,
            run_time: None,
            checkpoint_time: None,
            job_dir: PathBuf::from("."),
            tasks: Vec::new(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct SchedulerTask {
    pub target_sweeps: u64,
    pub sweeps: u64,
    pub thermalization: u64,
    pub dir: PathBuf,
    pub scheduled_runs: u64,
    pub max_scheduled_runs: u64,
}

impl SchedulerTask {
    pub fn new(target_sweeps: u64, thermalization: u64, dir: PathBuf) -> Self {
        Self {
            target_sweeps,
            sweeps: 0,
            thermalization,
            dir,
            scheduled_runs: 0,
            max_scheduled_runs: u64::MAX,
        }
    }

    pub fn is_done(&self) -> bool {
        self.sweeps >= self.target_sweeps
    }
}

#[cfg(feature = "mpi")]
pub fn run_distributed_compat<
    MC: FromParams<Rng = rand_xoshiro::Xoshiro256PlusPlus>
        + MonteCarlo<Rng = rand_xoshiro::Xoshiro256PlusPlus>,
>(
    config: MpiRunConfig,
) -> Result<Vec<Results>, CarloError> {
    let tasks = config
        .tasks
        .iter()
        .enumerate()
        .map(|(id, params)| TaskSpec {
            id,
            target_sweeps: config.run_config.measurement_sweeps,
            thermalization: config.run_config.thermalization_sweeps,
            params: params.clone(),
        })
        .collect();

    run_distributed::<MC, rand_xoshiro::Xoshiro256PlusPlus>(DistributedConfig {
        run_config: config.run_config,
        ranks_per_run: config.ranks_per_run,
        run_time: config.run_time,
        checkpoint_time: config.checkpoint_time,
        job_dir: config.job_dir,
        tasks,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn task(id: usize, target_sweeps: u64) -> TaskSpec {
        TaskSpec {
            id,
            target_sweeps,
            thermalization: 10,
            params: Params::new(),
        }
    }

    #[test]
    fn task_stream_reservations_never_overschedule() {
        let base = RunConfig {
            measurement_sweeps: 100,
            binsize: 4,
            ..RunConfig::default()
        };
        let mut stream = TaskStream::with_parallelism(vec![task(7, 10)], 3).unwrap();
        let mut assignments = Vec::new();
        while let Some(assignment) = stream.next_assignment(&base) {
            assignments.push(assignment);
        }
        assert_eq!(
            assignments
                .iter()
                .map(|assignment| assignment.measurement_sweeps)
                .sum::<u64>(),
            10
        );
        assert!(assignments
            .iter()
            .all(|assignment| assignment.measurement_sweeps > 0));
    }

    #[test]
    fn interrupted_assignment_releases_its_reservation() {
        let base = RunConfig {
            binsize: 2,
            ..RunConfig::default()
        };
        let mut stream = TaskStream::with_parallelism(vec![task(3, 8)], 2).unwrap();
        let first = stream.next_assignment(&base).unwrap();
        stream
            .interrupt_assignment(first.task_id, first.measurement_sweeps)
            .unwrap();
        let replacement = stream.next_assignment(&base).unwrap();
        assert_eq!(replacement.measurement_sweeps, first.measurement_sweeps);
    }

    #[test]
    fn duplicate_task_ids_are_rejected() {
        let error = TaskStream::with_parallelism(vec![task(1, 1), task(1, 2)], 1)
            .err()
            .expect("duplicate ids must fail");
        assert!(error.to_string().contains("duplicate MPI task id"));
    }
}
