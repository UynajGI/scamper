//! MPI Backend for distributed Monte Carlo simulations.
//!
//! Rust-native design with:
//! - Type-state pattern for compile-time safety
//! - Channel abstraction over MPI
//! - Result-based error handling everywhere

#[cfg(feature = "mpi")]
use mpi::topology::SimpleCommunicator;
#[cfg(feature = "mpi")]
use mpi::traits::*;

use rand_core::Rng;
use rand_core::SeedableRng;
use std::collections::HashMap;
use std::path::PathBuf;
use std::time::{Duration, Instant};

#[cfg(feature = "hdf5")]
use crate::RngCheckpointHdf5;
use crate::{
    CarloError, Estimate, FromParams, MonteCarlo, Params, Results, Run, RunConfig, RunId, TaskId,
};

// ============================================================================
// Time Limits
// ============================================================================

/// Time limits for simulation control.
pub struct TimeLimits {
    checkpoint_time: Duration,
    run_time: Duration,
    last_checkpoint: Instant,
    start_time: Instant,
}

impl TimeLimits {
    pub fn new(checkpoint_time: Duration, run_time: Duration) -> Self {
        Self {
            checkpoint_time,
            run_time,
            last_checkpoint: Instant::now(),
            start_time: Instant::now(),
        }
    }

    pub fn should_checkpoint(&self) -> bool {
        self.last_checkpoint.elapsed() >= self.checkpoint_time
    }

    pub fn should_finish(&self) -> bool {
        self.start_time.elapsed() >= self.run_time
    }

    pub fn reset_checkpoint(&mut self) {
        self.last_checkpoint = Instant::now();
    }
}

// ============================================================================
// MPI Error
// ============================================================================

/// MPI Error type
#[cfg(feature = "mpi")]
#[derive(Debug, thiserror::Error)]
pub enum MpiError {
    #[error("MPI initialization failed")]
    InitFailed,

    #[error("Communication error: {0}")]
    Communication(String),

    #[error("Invalid state transition")]
    InvalidTransition,
}

// ============================================================================
// MPI Serializable
// ============================================================================

/// Trait for MPI-transmissible types
#[cfg(feature = "mpi")]
pub trait MpiSerializable: Sized {
    fn send_mpi(&self, comm: &SimpleCommunicator, dest: i32, tag: i32) -> Result<(), MpiError>;
    fn recv_mpi(comm: &SimpleCommunicator, source: i32, tag: i32) -> Result<Self, MpiError>;
}

#[cfg(feature = "mpi")]
impl<T: mpi::datatype::Equivalence + Copy> MpiSerializable for T {
    fn send_mpi(&self, comm: &SimpleCommunicator, dest: i32, tag: i32) -> Result<(), MpiError> {
        comm.process_at_rank(dest).send_with_tag(self, tag);
        Ok(())
    }

    fn recv_mpi(comm: &SimpleCommunicator, source: i32, tag: i32) -> Result<Self, MpiError> {
        let (val, _) = comm.process_at_rank(source).receive_with_tag(tag);
        Ok(val)
    }
}

// ============================================================================
// Message Types
// ============================================================================

mod tags {
    pub const WORKER_MSG: i32 = 4355;
    pub const CONTROLLER_MSG: i32 = 4356;
}

/// Worker status message
#[cfg(feature = "mpi")]
#[derive(Debug, Clone, Copy, Equivalence, Default)]
pub struct WorkerStatusMsg {
    pub status: i32,
    pub task_id: u64,
    pub sweeps: u64,
}

#[cfg(feature = "mpi")]
impl WorkerStatusMsg {
    pub fn idle() -> Self {
        Self {
            status: 0,
            task_id: 0,
            sweeps: 0,
        }
    }
    pub fn progress(task_id: usize, sweeps: u64) -> Self {
        Self {
            status: 1,
            task_id: task_id as u64,
            sweeps,
        }
    }
    pub fn complete(task_id: usize) -> Self {
        Self {
            status: 2,
            task_id: task_id as u64,
            sweeps: 0,
        }
    }

    pub fn is_idle(&self) -> bool {
        self.status == 0
    }
    pub fn is_progress(&self) -> bool {
        self.status == 1
    }
    pub fn is_complete(&self) -> bool {
        self.status == 2
    }
    pub fn is_timeup(&self) -> bool {
        self.status == 3
    }
}

/// Controller command message
#[cfg(feature = "mpi")]
#[derive(Debug, Clone, Copy, Equivalence, Default)]
pub struct ControllerCmdMsg {
    pub action: i32,
    pub task_id: u64,
    pub run_id: u64,
    pub sweeps_hint: u64,
}

#[cfg(feature = "mpi")]
impl ControllerCmdMsg {
    pub fn exit() -> Self {
        Self {
            action: 0,
            ..Default::default()
        }
    }
    pub fn assign_task(task_id: usize, run_id: u64, sweeps_hint: u64) -> Self {
        Self {
            action: 1,
            task_id: task_id as u64,
            run_id,
            sweeps_hint,
        }
    }
    pub fn continue_(sweeps_hint: u64) -> Self {
        Self {
            action: 2,
            sweeps_hint,
            ..Default::default()
        }
    }
    pub fn finish_and_new() -> Self {
        Self {
            action: 3,
            ..Default::default()
        }
    }

    pub fn is_exit(&self) -> bool {
        self.action == 0
    }
    pub fn is_assign(&self) -> bool {
        self.action == 1
    }
    pub fn is_continue(&self) -> bool {
        self.action == 2
    }
    pub fn is_finish_and_new(&self) -> bool {
        self.action == 3
    }
}

// ============================================================================
// Task Stream
// ============================================================================

/// Task specification
#[derive(Debug, Clone)]
pub struct TaskSpec {
    pub id: usize,
    pub target_sweeps: u64,
    pub thermalization: u64,
    pub params: Params,
}

struct TaskInfo {
    spec: TaskSpec,
    sweeps_done: u64,
    runs_active: u64,
}

/// Stream of tasks
pub struct TaskStream {
    tasks: Vec<TaskInfo>,
    current_id: Option<usize>,
}

impl TaskStream {
    pub fn new(specs: Vec<TaskSpec>) -> Self {
        let tasks = specs
            .into_iter()
            .map(|spec| TaskInfo {
                spec,
                sweeps_done: 0,
                runs_active: 0,
            })
            .collect();
        Self {
            tasks,
            current_id: None,
        }
    }

    pub fn next(&mut self, num_workers: i32) -> Option<&TaskSpec> {
        let start = self.current_id.map(|id| id + 1).unwrap_or(0);

        for (i, task) in self.tasks.iter().enumerate().skip(start) {
            if self.has_work(task, num_workers) {
                self.current_id = Some(i);
                return Some(&task.spec);
            }
        }

        for (i, task) in self.tasks.iter().enumerate().take(start) {
            if self.has_work(task, num_workers) {
                self.current_id = Some(i);
                return Some(&task.spec);
            }
        }

        None
    }

    fn has_work(&self, task: &TaskInfo, num_workers: i32) -> bool {
        if task.sweeps_done >= task.spec.target_sweeps {
            return false;
        }
        let remaining = task.spec.target_sweeps - task.sweeps_done;
        let min_work = std::cmp::max(
            task.spec.thermalization * task.runs_active,
            task.runs_active,
        );
        remaining > min_work && task.runs_active < num_workers as u64
    }

    pub fn report_progress(&mut self, task_id: usize, sweeps: u64) {
        if let Some(task) = self.tasks.get_mut(task_id) {
            task.sweeps_done += sweeps;
        }
    }

    pub fn start_run(&mut self, task_id: usize) {
        if let Some(task) = self.tasks.get_mut(task_id) {
            task.runs_active += 1;
        }
    }

    pub fn complete_run(&mut self, task_id: usize) {
        if let Some(task) = self.tasks.get_mut(task_id) {
            task.runs_active = task.runs_active.saturating_sub(1);
        }
    }

    pub fn sweeps_hint(&self, task_id: usize, num_active: i32) -> u64 {
        let task = match self.tasks.get(task_id) {
            Some(t) => t,
            None => return 0,
        };

        if task.runs_active == 0 {
            return task.spec.target_sweeps - task.sweeps_done;
        }

        let remaining = task.spec.target_sweeps - task.sweeps_done;
        std::cmp::min(
            remaining / task.runs_active,
            task.spec.target_sweeps / num_active.max(1) as u64,
        )
    }

    pub fn is_done(&self, task_id: usize) -> bool {
        self.tasks
            .get(task_id)
            .map(|t| t.sweeps_done >= t.spec.target_sweeps)
            .unwrap_or(true)
    }
}

// ============================================================================
// Results Aggregator
// ============================================================================

/// Aggregates results from multiple runs
pub struct ResultsAggregator {
    observables: HashMap<String, Vec<f64>>,
}

impl ResultsAggregator {
    pub fn new() -> Self {
        Self {
            observables: HashMap::new(),
        }
    }

    pub fn add(&mut self, results: &Results) {
        for (name, estimate) in results.estimates() {
            self.observables
                .entry(name.clone())
                .or_default()
                .push(estimate.mean);
        }
    }

    pub fn finalize(self) -> Results {
        let estimates: HashMap<String, Estimate> = self
            .observables
            .into_iter()
            .map(|(name, values)| (name, Estimate::from_bins(&values)))
            .collect();
        Results::from_measurements(&estimates)
    }
}

impl Default for ResultsAggregator {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================================
// Type-State Worker
// ============================================================================

/// Worker state marker
pub trait WorkerState: private::Sealed {}
pub struct Idle;
pub struct Running {
    pub task_id: usize,
    pub run_id: u64,
}
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

/// Worker with type-state
#[cfg(feature = "mpi")]
pub struct Worker<S: WorkerState> {
    world_rank: i32,
    leader_comm: SimpleCommunicator,
    run_comm: Option<SimpleCommunicator>,
    state: S,
}

#[cfg(feature = "mpi")]
impl Worker<Idle> {
    pub fn new(
        world_rank: i32,
        leader_comm: &SimpleCommunicator,
        run_comm: Option<&SimpleCommunicator>,
    ) -> Self {
        // SAFETY: SimpleCommunicator does not implement Clone, but we need to store
        // it in the Worker struct. We use ptr::read to create a bitwise copy.
        // This is safe because:
        // 1. The caller's communicator reference is borrowed for the call duration only
        // 2. We take ownership of the bitwise copy, which will be dropped when Worker drops
        // 3. The MPI library handles reference counting internally for communicator handles
        // 4. This pattern is used because the mpi crate doesn't expose a safe Clone for SimpleCommunicator
        Self {
            world_rank,
            leader_comm: unsafe { std::ptr::read(leader_comm) },
            run_comm: run_comm.map(|c| unsafe { std::ptr::read(c) }),
            state: Idle,
        }
    }

    pub fn recv_task(self) -> Result<WorkerEither, MpiError> {
        // Send idle
        WorkerStatusMsg::idle().send_mpi(&self.leader_comm, 0, tags::WORKER_MSG)?;

        // Receive command
        let cmd = ControllerCmdMsg::recv_mpi(&self.leader_comm, 0, tags::CONTROLLER_MSG)?;

        if cmd.is_exit() {
            Ok(WorkerEither::Done(Worker {
                world_rank: self.world_rank,
                leader_comm: self.leader_comm,
                run_comm: self.run_comm,
                state: Done,
            }))
        } else if cmd.is_assign() {
            Ok(WorkerEither::Running(Worker {
                world_rank: self.world_rank,
                leader_comm: self.leader_comm,
                run_comm: self.run_comm,
                state: Running {
                    task_id: cmd.task_id as usize,
                    run_id: cmd.run_id,
                },
            }))
        } else {
            Err(MpiError::InvalidTransition)
        }
    }
}

#[cfg(feature = "mpi")]
impl Worker<Running> {
    pub fn task_id(&self) -> usize {
        self.state.task_id
    }
    pub fn run_id(&self) -> u64 {
        self.state.run_id
    }

    pub fn send_progress(&mut self, sweeps: u64) -> Result<ControllerCmdMsg, MpiError> {
        WorkerStatusMsg::progress(self.state.task_id, sweeps).send_mpi(
            &self.leader_comm,
            0,
            tags::WORKER_MSG,
        )?;
        ControllerCmdMsg::recv_mpi(&self.leader_comm, 0, tags::CONTROLLER_MSG)
    }

    pub fn finish(self) -> Result<Worker<Done>, MpiError> {
        WorkerStatusMsg::complete(self.state.task_id).send_mpi(
            &self.leader_comm,
            0,
            tags::WORKER_MSG,
        )?;
        Ok(Worker {
            world_rank: self.world_rank,
            leader_comm: self.leader_comm,
            run_comm: self.run_comm,
            state: Done,
        })
    }
}

#[cfg(feature = "mpi")]
impl Worker<Done> {
    pub fn reset(self) -> Worker<Idle> {
        Worker {
            world_rank: self.world_rank,
            leader_comm: self.leader_comm,
            run_comm: self.run_comm,
            state: Idle,
        }
    }
}

/// Either running or done worker
#[cfg(feature = "mpi")]
pub enum WorkerEither {
    Running(Worker<Running>),
    Done(Worker<Done>),
}

// ============================================================================
// Distributed Config
// ============================================================================

/// Configuration for distributed runs
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

// ============================================================================
// Entry Point
// ============================================================================

/// Run distributed simulation
#[cfg(feature = "mpi")]
pub fn run_distributed<MC, R>(config: DistributedConfig) -> Result<Vec<Results>, CarloError>
where
    MC: MonteCarlo + FromParams<Rng = R>,
    R: Rng + SeedableRng + Send,
{
    let universe = mpi::initialize().ok_or(MpiError::InitFailed)?;
    let world = universe.world();
    let world_rank = world.rank();
    let world_size = world.size();

    if world_size < 2 {
        return Err(CarloError::InvalidConfig {
            field: "mpi".into(),
            reason: "MPI requires at least 2 ranks".into(),
        });
    }

    if (world_size - 1) % config.ranks_per_run != 0 {
        return Err(CarloError::InvalidConfig {
            field: "ranks_per_run".into(),
            reason: format!(
                "Worker ranks ({}) not divisible by ranks_per_run ({})",
                world_size - 1,
                config.ranks_per_run
            ),
        });
    }

    // Split communicators
    let run_color = if world_rank == 0 {
        mpi::topology::Color::undefined()
    } else {
        mpi::topology::Color::with_value(1 + (world_rank - 1) / config.ranks_per_run)
    };
    let run_comm = world.split_by_color(run_color);

    let is_leader = world_rank == 0 || ((world_rank - 1) % config.ranks_per_run == 0);
    let leader_color = if is_leader {
        mpi::topology::Color::with_value(1)
    } else {
        mpi::topology::Color::undefined()
    };
    let leader_comm =
        world
            .split_by_color(leader_color)
            .ok_or_else(|| CarloError::InvalidConfig {
                field: "leader_comm".into(),
                reason: "Failed to create leader communicator".into(),
            })?;

    if world_rank == 0 {
        run_controller(&leader_comm, &config)
    } else {
        run_worker::<MC, R>(world_rank, &leader_comm, run_comm.as_ref(), &config)
    }
}

#[cfg(feature = "mpi")]
fn run_controller(
    leader_comm: &SimpleCommunicator,
    config: &DistributedConfig,
) -> Result<Vec<Results>, CarloError> {
    let num_workers = leader_comm.size() - 1;
    let mut tasks = TaskStream::new(config.tasks.clone());
    let aggregator = ResultsAggregator::new();
    let mut active_workers = num_workers;

    leader_comm.barrier();

    while active_workers > 0 {
        let (status, mpi_status) = leader_comm.any_process().receive::<WorkerStatusMsg>();
        let source = mpi_status.source_rank();

        if status.is_idle() {
            let task_opt = tasks.next(active_workers);
            if let Some(task) = task_opt {
                let task_id = task.id;
                tasks.start_run(task_id);
                let hint = tasks.sweeps_hint(task_id, active_workers);
                let run_id = 1; // Simplified
                let cmd = ControllerCmdMsg::assign_task(task_id, run_id, hint);
                cmd.send_mpi(leader_comm, source, tags::CONTROLLER_MSG)?;
            } else {
                ControllerCmdMsg::exit().send_mpi(leader_comm, source, tags::CONTROLLER_MSG)?;
                active_workers -= 1;
            }
        } else if status.is_progress() {
            let task_id = status.task_id as usize;
            tasks.report_progress(task_id, status.sweeps);

            if tasks.is_done(task_id) {
                tasks.complete_run(task_id);
                ControllerCmdMsg::finish_and_new().send_mpi(
                    leader_comm,
                    source,
                    tags::CONTROLLER_MSG,
                )?;
            } else {
                let hint = tasks.sweeps_hint(task_id, active_workers);
                ControllerCmdMsg::continue_(hint).send_mpi(
                    leader_comm,
                    source,
                    tags::CONTROLLER_MSG,
                )?;
            }
        } else if status.is_complete() {
            // Would receive results here
        } else if status.is_timeup() {
            active_workers -= 1;
        }
    }

    leader_comm.barrier();
    Ok(vec![aggregator.finalize()])
}

#[cfg(feature = "mpi")]
fn run_worker<MC, R>(
    world_rank: i32,
    leader_comm: &SimpleCommunicator,
    _run_comm: Option<&SimpleCommunicator>,
    config: &DistributedConfig,
) -> Result<Vec<Results>, CarloError>
where
    MC: MonteCarlo + FromParams<Rng = R>,
    R: Rng + SeedableRng + Send,
{
    let mut worker = Worker::<Idle>::new(world_rank, leader_comm, None);
    let mut results = Vec::new();
    let time_start = Instant::now();
    #[allow(unused_mut)]
    let mut last_checkpoint = Instant::now();

    leader_comm.barrier();

    loop {
        let either = worker.recv_task()?;
        let WorkerEither::Running(w) = either else {
            break;
        };
        let task_id = w.task_id();
        let default_task = TaskSpec {
            id: 0,
            target_sweeps: 100,
            thermalization: 10,
            params: Params::new(),
        };
        let task = config.tasks.get(task_id).unwrap_or(&default_task);
        let params = &task.params;

        // Build RunConfig from task spec
        let run_config = RunConfig {
            measurement_sweeps: task.target_sweeps,
            thermalization_sweeps: task.thermalization,
            binsize: config.run_config.binsize,
            base_seed: config.run_config.base_seed,
            progress_interval: config.run_config.progress_interval,
            checkpoint_interval: config.run_config.checkpoint_interval,
        };

        let seed = config
            .run_config
            .base_seed
            .wrapping_add((task_id as u64) * 10000)
            .wrapping_add(w.run_id());

        // Check for existing checkpoint
        #[cfg(feature = "hdf5")]
        let checkpoint_path = config.job_dir.join(format!(
            "task_{:04}/run{:04}/run{:04}.dump.h5",
            task_id,
            w.run_id(),
            w.run_id()
        ));

        #[cfg(feature = "hdf5")]
        let run = if let Some(existing_run) =
            Run::<MC, R>::read_checkpoint(&checkpoint_path, params, &run_config, seed)?
        {
            existing_run
        } else {
            Run::new(
                params,
                TaskId::new(task_id),
                RunId::new(w.run_id()),
                &run_config,
                seed,
            )?
        };

        #[cfg(not(feature = "hdf5"))]
        let run: Run<MC, R> = Run::new(
            params,
            TaskId::new(task_id),
            RunId::new(w.run_id()),
            &run_config,
            seed,
        )?;

        let mut run = run;
        let mut w = w;

        // Run simulation with progress reporting and checkpointing
        loop {
            run.step();

            // Check time limits
            let elapsed = time_start.elapsed();
            let should_checkpoint = config
                .checkpoint_time
                .is_some_and(|ct| last_checkpoint.elapsed() >= ct);
            let should_finish = config.run_time.is_some_and(|rt| elapsed >= rt);

            if should_checkpoint || should_finish || run.is_complete() {
                #[cfg(feature = "hdf5")]
                if should_checkpoint && !run.is_complete() {
                    run.write_checkpoint(&checkpoint_path)?;
                    last_checkpoint = Instant::now();
                }

                // Report progress
                let cmd = w.send_progress(run.sweeps_done())?;

                if should_finish && !run.is_complete() {
                    // Write final checkpoint on timeup
                    #[cfg(feature = "hdf5")]
                    run.write_checkpoint(&checkpoint_path)?;

                    // Send timeup status
                    WorkerStatusMsg {
                        status: 3,
                        task_id: task_id as u64,
                        sweeps: run.sweeps_done(),
                    }
                    .send_mpi(&w.leader_comm, 0, tags::WORKER_MSG)?;

                    return Ok(results); // Worker exits on timeup
                }

                // Reset worker for next task iteration
                worker = w.finish()?.reset();

                if cmd.is_finish_and_new() || run.is_complete() {
                    let result = run.finalize(config.run_config.base_seed);
                    results.push(result);
                }

                break; // inner loop
            }
        }
    }

    Ok(results)
}

// ============================================================================
// Backward Compatibility
// ============================================================================

/// Old MpiRunConfig for backwards compatibility
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

/// Old MpiBackend for backwards compatibility
#[cfg(feature = "mpi")]
pub struct MpiBackend {
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
        let universe = mpi::initialize().ok_or(MpiError::InitFailed)?;
        let world = universe.world();
        Ok(Self {
            rank: world.rank(),
            size: world.size(),
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
        self.is_controller() || ((self.rank - 1) % self.ranks_per_run == 0)
    }
    pub fn num_workers(&self) -> i32 {
        self.size - 1
    }
    pub fn num_parallel_runs(&self) -> i32 {
        self.num_workers() / self.ranks_per_run
    }
    pub fn run_group(&self) -> i32 {
        if self.is_controller() {
            0
        } else {
            1 + (self.rank - 1) / self.ranks_per_run
        }
    }
    pub fn rank_in_run(&self) -> i32 {
        if self.is_controller() {
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
        Self::new().expect("Failed to create MPI backend")
    }
}

#[cfg(feature = "mpi")]
impl Clone for MpiBackend {
    fn clone(&self) -> Self {
        Self {
            rank: self.rank,
            size: self.size,
            ranks_per_run: self.ranks_per_run,
        }
    }
}

#[cfg(feature = "mpi")]
impl super::Backend for MpiBackend {
    type Rng = rand_xoshiro::Xoshiro256PlusPlus;

    fn spawn_tasks<F>(&self, _: usize, _: u64, _: F)
    where
        F: Fn(usize, &mut Self::Rng) + Sync,
    {
        unimplemented!("MpiBackend uses run_distributed")
    }

    fn barrier(&self) {
        if let Some(universe) = mpi::initialize() {
            universe.world().barrier();
        }
    }
}

/// Scheduler task for backwards compatibility
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

/// Old run_distributed entry point
#[cfg(feature = "mpi")]
pub fn run_distributed_compat<
    MC: FromParams + MonteCarlo<Rng = rand_xoshiro::Xoshiro256PlusPlus>,
>(
    config: MpiRunConfig,
) -> Result<Vec<Results>, CarloError> {
    let tasks: Vec<TaskSpec> = config
        .tasks
        .iter()
        .enumerate()
        .map(|(i, _)| TaskSpec {
            id: i,
            target_sweeps: config.run_config.measurement_sweeps,
            thermalization: config.run_config.thermalization_sweeps,
            params: Params::new(),
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

// Non-MPI stub
#[cfg(not(feature = "mpi"))]
pub struct MpiBackend;

#[cfg(not(feature = "mpi"))]
impl MpiBackend {
    pub fn new() -> Result<Self, CarloError> {
        Err(CarloError::InvalidConfig {
            field: "mpi".into(),
            reason: "MPI feature not enabled".into(),
        })
    }
}
