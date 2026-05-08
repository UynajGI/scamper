---
name: Carlo.jl to Carlo.rs Complete Migration
description: Complete migration of Carlo.jl functionality to Carlo.rs with Rust-native design advantages
type: project
---

# Carlo.jl → Carlo.rs Complete Migration Design Spec

**Date**: 2026-04-08
**Status**: Approved for implementation

## 1. Project Overview

### Goal
Complete migration of all Carlo.jl functionality to Carlo.rs, leveraging Rust's advantages while maintaining 100% functional equivalence.

### Scope
- All 25 Julia source modules → Rust equivalents
- Preserve Rust-native advantages already implemented
- Ensure numerical consistency between implementations

### Rust Advantages to Preserve (No Changes)
| Julia Module | Rust Replacement | Reason |
|--------------|------------------|--------|
| `random_wrap.jl` | `rand` crate + `SeedableRng` trait | Rust's RNG trait system is superior |
| `log.jl` | `log` + `tracing` crates | Better structured logging ecosystem |
| `tinyargparse.jl` | `clap` derive macros | Type-safe CLI with compile-time validation |

### Modules to Implement

**Layer 1 - Core Analysis (no dependencies):**
- `merge.jl` → `src/merge.rs`
- `evaluable.jl` → `src/evaluable.rs`

**Layer 2 - Run Lifecycle (depends on Layer 1):**
- `run.jl` → `src/run.rs`
- `jobtools/taskinfo.jl` → `src/job/taskinfo.rs`
- `jobtools/jobinfo.jl` → `src/job/jobinfo.rs`
- `jobtools/taskmaker.jl` → `src/job/taskmaker.rs`

**Layer 3 - Scheduler Completion (depends on Layer 2):**
- MPI backend completion (communicator splitting, time limits, checkpoint)
- `cli.jl` → `src/cli.rs`
- `parallel_tempering.jl` → `src/parallel_tempering.rs`

---

## 2. Architecture

### File Structure After Completion

```
Carlo.rs/src/
├── lib.rs              # Public exports
├── error.rs            # ✓ CarloError enum
├── estimate.rs         # ✓ Estimate struct
├── measurements.rs     # ✓ Accumulator, Measurements
├── context.rs          # ✓ Context - needs checkpoint support
├── monte_carlo.rs      # ✓ MonteCarlo trait
├── params.rs           # ✓ Params struct
├── results.rs          # ✓ Results - needs merge support
├── lattice.rs          # ✓ LatticeParams
├── scheduler.rs        # ✓ Scheduler, RunConfig
│
├── merge.rs            # NEW - ResultObservable, rebin analysis
├── evaluable.rs        # NEW - Evaluable, Evaluator, jackknife
├── run.rs              # NEW - Run struct, step!, checkpoint lifecycle
│
├── job/
│   ├── mod.rs          # NEW - JobTools module
│   ├── taskinfo.rs     # NEW - TaskInfo, TaskProgress
│   ├── jobinfo.rs      # NEW - JobInfo, duration parsing
│   └── taskmaker.rs    # NEW - TaskMaker builder
│
├── backend/
│   ├── mod.rs          # ✓ Backend trait, exports
│   ├── rayon.rs        # ✓ RayonBackend
│   └── mpi.rs          # ✓ Worker type-state - needs completion
│
├── output/
│   ├── mod.rs          # ✓ save_hdf5, save_json
│   ├── hdf5.rs         # ✓ - needs checkpoint format
│   └── json.rs         # ✓
│
├── cli.rs              # NEW - CLI commands (run/status/merge/delete)
└── parallel_tempering.rs  # NEW - ParallelTemperingMC
```

---

## 3. Layer 1: Core Analysis

### 3.1 `merge.rs` - Result Merging

**Purpose**: Merge measurements from multiple runs, perform rebinning analysis, compute autocorrelation times.

**Key Types**:
```rust
/// Observable metadata from HDF5 files
pub struct ObservableType<T, N> {
    internal_bin_length: u64,
    shape: [usize; N],
    total_sample_count: u64,
}

/// Merged observable with statistics
pub struct ResultObservable<T, N> {
    internal_bin_length: u64,
    rebin_length: u64,
    mean: Array<T, N>,
    error: Array<T, N>,
    covariance: Option<Array<T, N+N>>,
    autocorrelation_time: Array<f64, N>,
    rebin_means: Array<T, N+1>,  // shape + bin_count
}
```

**Key Functions**:
```rust
/// Calculate optimal rebin count
pub fn calc_rebin_count(sample_count: u64, min_bin_count: u64 = 10) -> u64;

/// Merge results from multiple HDF5 measurement files
pub fn merge_results(
    filenames: &[PathBuf],
    rebin_length: Option<u64>,
    sample_skip: u64,
    estimate_covariance: bool,
) -> Result<HashMap<String, ResultObservable>, CarloError>;

/// Iterate over observables in HDF5 files
pub fn iterate_measfile_observables<F, T>(
    filenames: &[PathBuf],
    f: F,
) -> Result<HashMap<String, T>, CarloError>;
```

**Autocorrelation Time Computation**:
- Regular: `τ = 0.5 * ((σ_rebin / σ_no_rebin)^2 - 1)`
- Decorrelated (with covariance): eigenvalue decomposition of Σ

**Dependencies**: `hdf5`, `ndarray`

---

### 3.2 `evaluable.rs` - Jackknife Analysis

**Purpose**: Compute derived quantities with correct error propagation via jackknife resampling.

**Key Types**:
```rust
/// Derived observable with jackknife error
pub struct Evaluable<T, R, N, C> {
    internal_bin_length: u64,
    rebin_length: u64,
    rebin_count: u64,
    mean: Array<T, N>,
    error: Array<R, N>,
    covariance: Option<C>,
}

/// Evaluator for defining derived observables
pub struct Evaluator {
    observables: HashMap<String, ResultObservable>,
    evaluables: HashMap<String, Evaluable>,
    estimate_covariance: bool,
}
```

**Key Functions**:
```rust
/// Jackknife resampling for error propagation
pub fn jackknife<F, T, N>(
    func: F,
    sample_sets: &[Array<T, N+1>],  // N dims + sample_count
    estimate_covariance: bool,
) -> Result<(Array<T, N>, Array<f64, N>, Option<Array<T, N+N>>), CarloError>;

/// Define an evaluable from observables
pub fn evaluate<E>(
    evaluator: &mut Evaluator,
    name: &str,
    ingredients: &[&str],
    evaluation: E,
) -> Result<(), CarloError>;
```

**Jackknife Algorithm**:
1. Compute complete evaluation from all samples
2. Remove one sample at a time, recompute
3. Bias-corrected mean: `n * complete - (n-1) * jacked_mean`
4. Error: `sqrt((n-1)/n * Σ(jacked_eval - jacked_mean)^2)`

**Dependencies**: `merge.rs`

---

## 4. Layer 2: Run Lifecycle

### 4.1 `run.rs` - Simulation Run

**Purpose**: Single run lifecycle - initialization, stepping, checkpointing.

**Key Types**:
```rust
/// A single Monte Carlo run
pub struct Run<MC: MonteCarlo, R: Rng + SeedableRng> {
    context: Context<R>,
    implementation: MC,
}

/// Run state for checkpoint
pub struct RunCheckpoint {
    sweeps: u64,
    thermalization_sweeps: u64,
    rng_state: Vec<u8>,
}
```

**Key Functions**:
```rust
/// Create a new run
pub fn new_run<MC, R>(
    params: &Params,
    seed_variation: u64,
) -> Result<Run<MC, R>, CarloError>;

/// Perform one MC step, return thermalized sweep count
pub fn step<MC, R>(
    run: &mut Run<MC, R>,
) -> Result<u64, CarloError>;

/// Write checkpoint to HDF5
pub fn write_checkpoint<MC, R>(
    run: &Run<MC, R>,
    path: &Path,
) -> Result<(), CarloError>;

/// Read checkpoint from HDF5
pub fn read_checkpoint<MC, R>(
    path: &Path,
    params: &Params,
) -> Result<Option<Run<MC, R>>, CarloError>;
```

**Checkpoint Format (HDF5)**:
```
/run.dump.h5
├── context/
│   ├── 0001/   # rank 0
│   │   ├── sweeps
│   │   ├── thermalization_sweeps
│   │   └── rng_state
│   ├── 0002/   # rank 1 (parallel runs)
│   └── ...
├── simulation/
│   └── [MC-specific data]
└── version/
    ├── carlo_version
    └── mc_version
```

**Dependencies**: `merge.rs`, `evaluable.rs`, HDF5

---

### 4.2 `job/taskinfo.rs` - Task Information

**Purpose**: Define a single parameter set for MC calculation.

**Key Types**:
```rust
/// Task parameters with validation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TaskInfo {
    name: String,
    params: HashMap<String, Value>,
}

/// Task progress tracking
#[derive(Debug, Clone)]
pub struct TaskProgress {
    target_sweeps: u64,
    sweeps: u64,
    num_runs: u64,
    thermalization_fraction: f64,
    dir: PathBuf,
}
```

**Required Parameters**:
- `sweeps`: Minimum measurement sweeps
- `thermalization`: Thermalization sweeps
- `binsize`: Internal bin size for observables

**Optional Parameters**:
- `seed`: Fixed seed for debugging
- `rebin_length`: Override automatic rebinning
- `rebin_sample_skip`: Skip first N bins
- `estimate_covariance`: Estimate covariance matrices
- `max_runs_per_task`: Limit runs

**Key Functions**:
```rust
pub fn task_name(task_id: u64) -> String;  // "task%04d"
pub fn list_run_files(dir: &Path, pattern: &str) -> Vec<PathBuf>;
pub fn read_dump_progress(dir: &Path) -> Result<Vec<(u64, u64)>, CarloError>;
```

---

### 4.3 `job/jobinfo.rs` - Job Configuration

**Purpose**: Complete job definition with tasks and timing.

**Key Types**:
```rust
/// Job configuration
#[derive(Debug, Clone)]
pub struct JobInfo {
    name: String,
    dir: PathBuf,
    mc_type: String,  // Type name for dynamic dispatch
    rng_type: String,
    tasks: Vec<TaskInfo>,
    checkpoint_time: Duration,
    run_time: Duration,
    ranks_per_run: usize,  // 0 = all
}
```

**Key Functions**:
```rust
/// Parse duration "[[[days-]hours:]minutes:]seconds"
pub fn parse_duration(s: &str) -> Result<Duration, CarloError>;

/// Get run time from SLURM environment
pub fn run_time_from_slurm(grace_factor: f64, default: Duration) -> Duration;

/// Task directory path
pub fn task_dir(job: &JobInfo, task: &TaskInfo) -> PathBuf;

/// Read all task progress
pub fn read_progress(job: &JobInfo) -> Result<Vec<TaskProgress>, CarloError>;

/// Check timing conditions
pub fn is_checkpoint_time(job: &JobInfo, last_checkpoint: DateTime) -> bool;
pub fn is_end_time(job: &JobInfo, start_time: DateTime) -> bool;
```

---

### 4.4 `job/taskmaker.rs` - Task Builder

**Purpose**: Fluent builder for generating multiple tasks.

**Key Types**:
```rust
/// Builder for task list
pub struct TaskMaker {
    tasks: Vec<TaskInfo>,
    current_params: HashMap<String, Value>,
}
```

**Key Methods**:
```rust
impl TaskMaker {
    pub fn new() -> Self;

    /// Set a parameter
    pub fn set(&mut self, key: &str, value: Value) -> &mut Self;

    /// Create task with current params + overrides
    pub fn task(&mut self, overrides: HashMap<String, Value>) -> &mut Self;

    /// Finalize and return tasks
    pub fn make_tasks(self) -> Vec<TaskInfo>;

    /// Current task name (task%04d format)
    pub fn current_task_name(&self) -> String;
}
```

**Example Usage**:
```rust
let mut tm = TaskMaker::new();
tm.set("sweeps", 10000)
   .set("thermalization", 2000)
   .set("binsize", 500);

tm.task(Params::new().insert("T", 0.04));

tm.set("sweeps", 5000);
for T in (0.1..10.0).step_by(2) {
    tm.task(Params::new().insert("T", T));
}

let tasks = tm.make_tasks();
```

---

## 5. Layer 3: Scheduler Completion

### 5.1 MPI Backend Completion

**Current Status**: Type-state Worker pattern implemented, needs:

**Missing Features**:

1. **Communicator Splitting** (for parallel runs):
```rust
/// Split MPI world into run communicators
pub fn split_communicators(
    world: &SimpleCommunicator,
    ranks_per_run: usize,
) -> Result<(SimpleCommunicator, Option<SimpleCommunicator>), MpiError>;
```

2. **Time Limits**:
```rust
pub struct TimeLimits {
    checkpoint_time: Duration,
    run_time: Duration,
    last_checkpoint: Instant,
    start_time: Instant,
}
```

3. **Result Collection**:
```rust
/// Gather results from all workers
pub fn gather_results(
    worker: &Worker<Done>,
) -> Result<HashMap<String, ResultObservable>, MpiError>;
```

4. **Checkpoint Integration**:
- Worker<Running>::checkpoint() → Worker<Running>
- Worker<Idle>::resume_from_checkpoint() → WorkerEither

---

### 5.2 `cli.rs` - Command Line Interface

**Purpose**: CLI entry point with subcommands.

**Subcommands**:
- `run` (-s single, -r restart)
- `status`
- `merge`
- `delete`

**Key Types**:
```rust
#[derive(Parser)]
#[command(name = "carlo")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Start simulation
    Run {
        #[arg(short, long)]
        single: bool,
        #[arg(short, long)]
        restart: bool,
    },
    /// Check progress
    Status,
    /// Merge results
    Merge,
    /// Delete job data
    Delete,
}
```

**Entry Point**:
```rust
pub fn start(job: JobInfo, args: &[String]) -> Result<(), CarloError>;
```

---

### 5.3 `parallel_tempering.rs` - Parallel Tempering MC

**Purpose**: Run multiple chains with parameter exchange.

**Key Types**:
```rust
/// Parallel tempering wrapper
pub struct ParallelTemperingMC<T, MC: MonteCarlo> {
    parameter_name: String,
    parameter_values: Vec<T>,
    tempering_interval: u64,
    chain_idx: usize,
    child_mc: MC,
}

/// Parallel measurements queue
pub struct ParallelMeasurements {
    queue: Vec<(String, ArrayD<Value>)>,
}
```

**Key Traits**:
```rust
/// Required trait for PT-compatible MC
pub trait ParallelTemperingCompatible: MonteCarlo {
    fn log_weight_ratio(&self, param: &str, new_value: Value) -> f64;
    fn change_parameter(&mut self, param: &str, new_value: Value);
}
```

**MPI Protocol**:
- Tag `PT_WEIGHT_MSG = 4573792`: Exchange log weight ratios
- Tag `PT_SWITCH_MSG = 4573793`: Accept/reject switch

---

## 6. Dependencies

### Cargo.toml Additions
```toml
[dependencies]
hdf5 = "0.8"          # HDF5 for checkpoint/measurements
ndarray = "0.15"      # Array operations
chrono = "0.4"        # Time handling
clap = { version = "4", features = ["derive"] }  # CLI
serde_json = "1.0"    # Already present
mpi = { version = "0.7", optional = true }  # Already present

[features]
mpi = ["dep:mpi"]
```

---

## 7. Implementation Order

### Phase 1: Core Analysis (Layer 1)
1. `merge.rs` - ObservableType, ResultObservable, rebinning
2. `evaluable.rs` - Evaluable, Evaluator, jackknife
3. Tests for numerical correctness

### Phase 2: Run Lifecycle (Layer 2)
4. `job/taskinfo.rs` - TaskInfo, TaskProgress
5. `job/jobinfo.rs` - JobInfo, duration parsing
6. `job/taskmaker.rs` - TaskMaker builder
7. `run.rs` - Run struct, checkpoint I/O
8. Tests for checkpoint persistence

### Phase 3: Scheduler Completion (Layer 3)
9. MPI communicator splitting
10. Time limits in scheduler
11. Result collection
12. `cli.rs` - CLI with clap
13. Tests for distributed runs

### Phase 4: Advanced Features
14. `parallel_tempering.rs` - PT algorithm
15. Integration tests
16. Numerical consistency tests with Carlo.jl

---

## 8. Testing Strategy

### Unit Tests
- Each module has dedicated test file
- Numerical correctness tests with known inputs
- Edge case handling (empty files, missing observables)

### Integration Tests
- Full run lifecycle with checkpoint recovery
- MPI distributed runs (requires `mpi` feature)
- CLI command tests

### Numerical Consistency Tests
- Generate same random seeds in Julia and Rust
- Compare observable means/errors within tolerance
- Compare autocorrelation times

---

## 9. Success Criteria

1. **Functional Equivalence**: All Carlo.jl operations work identically
2. **Numerical Consistency**: Results match within floating-point tolerance
3. **MPI Correctness**: Distributed runs produce same results as single runs
4. **Checkpoint Recovery**: Can resume from Julia checkpoints (HDF5 format compatible)
5. **Clippy Clean**: No warnings with `-D warnings`
6. **Test Coverage**: All modules have unit tests

---

## 10. Why This Design

**Why:** User requested complete functional migration from Carlo.jl to Carlo.rs after three incomplete attempts. The migration needs to cover all 25 Julia modules while preserving Rust advantages already implemented.

**How to apply:**
- Implement in layers (1→2→3→4) to manage complexity
- Each layer tests independently before next layer
- Numerical consistency tests verify correctness at each phase
- Final integration tests confirm full system works