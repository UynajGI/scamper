# Rust-Native MPI Backend Design

## Overview

Complete rewrite of MPI backend leveraging Rust's type system, ownership model, and zero-cost abstractions. Not a translation of Carlo.jl, but a design that would be impossible in Julia.

## Design Principles

1. **Type-state pattern** - Invalid states are unrepresentable at compile time
2. **Channel abstraction** - MPI communication wrapped in mpsc-style API
3. **Generic constraints** - Traits encode capabilities, monomorphization removes runtime cost
4. **Result everywhere** - No panics, no unimplemented!()
5. **Ownership-based safety** - One owner per communicator, transferred explicitly

## Architecture

### 1. Type-State State Machine

```
┌─────────┐     recv_task()      ┌──────────┐
│  Idle   │ ──────────────────► │  Running │
└─────────┘                     └──────────┘
     ▲                               │
     │                               │ finish()
     │                               ▼
     │                          ┌──────────┐
     └───────────────────────── │   Done   │
           reset()              └──────────┘
```

**Rust encoding:**
```rust
pub struct Worker<S: MpiState> {
    comm: MpiComm,
    rank: i32,
    _state: PhantomData<S>,
}

pub trait MpiState {}
pub struct Idle;
pub struct Running { task_id: usize, run_id: u64 }
pub struct Done;

// Compile-time enforced transitions
impl Worker<Idle> {
    pub fn recv_task(self) -> Result<Worker<Running>, MpiError>;
}

impl Worker<Running> {
    pub fn send_progress(&mut self, sweeps: u64) -> Result<(), MpiError>;
    pub fn finish(self, results: Results) -> Result<Worker<Done>, MpiError>;
}

impl Worker<Done> {
    pub fn reset(self) -> Worker<Idle>;
}
```

**Benefit:** Cannot accidentally send progress while idle, cannot finish twice.

### 2. Channel Abstraction over MPI

```rust
pub struct MpiReceiver<T> { comm: MpiComm, tag: i32, _marker: PhantomData<T> }
pub struct MpiSender<T> { comm: MpiComm, tag: i32, _marker: PhantomData<T> }

impl<T: MpiSerializable> MpiSender<T> {
    pub fn send(&self, value: &T) -> Result<(), MpiError>;
}

impl<T: MpiSerializable> MpiReceiver<T> {
    pub fn recv(&self) -> Result<T, MpiError>;
    pub fn try_recv(&self) -> Result<Option<T>, MpiError>;
}

// Factory
pub fn mpi_channel<T: MpiSerializable>(tag: i32) -> (MpiSender<T>, MpiReceiver<T>);
```

**Benefit:** Testable without MPI (mock channels), familiar API.

### 3. Task as Stream

```rust
pub struct TaskStream {
    tasks: VecDeque<TaskSpec>,
    current: Option<usize>,
}

impl TaskStream {
    pub fn next(&mut self) -> Option<TaskSpec>;
    pub fn report_progress(&mut self, task_id: usize, sweeps: u64);
    pub fn is_complete(&self) -> bool;
}

pub struct TaskSpec {
    pub id: usize,
    pub target_sweeps: u64,
    pub thermalization: u64,
    pub params: Params,
}
```

### 4. Results Aggregation

```rust
pub struct ResultsAggregator {
    estimates: HashMap<String, BinsAccumulator>,
}

impl ResultsAggregator {
    pub fn add(&mut self, results: &Results);
    pub fn finalize(self) -> Results;
}

struct BinsAccumulator {
    bins: Vec<f64>,
    bin_size: usize,
}
```

### 5. Run Lifecycle

```rust
pub struct Simulation<MC: MonteCarlo> {
    context: Context<MC::Rng>,
    mc: MC,
    spec: TaskSpec,
    sweeps_done: u64,
}

impl<MC: MonteCarlo> Simulation<MC> {
    pub fn new(params: &Params, spec: TaskSpec, seed: u64) -> Result<Self, CarloError>
    where MC: FromParams;
    
    pub fn step(&mut self) -> StepResult;
    pub fn is_thermalized(&self) -> bool;
    pub fn finalize(self) -> Results;
}

pub enum StepResult {
    Thermalizing { remaining: u64 },
    Measuring { thermalized_sweeps: u64 },
    Complete,
}
```

### 6. Controller-Worker Protocol (Typed)

```rust
// Controller -> Worker messages
#[derive(MpiSerializable)]
pub enum ControllerMsg {
    AssignTask { task: TaskSpec, run_id: u64, sweeps_hint: u64 },
    Continue { sweeps_hint: u64 },
    Exit,
}

// Worker -> Controller messages  
#[derive(MpiSerializable)]
pub enum WorkerMsg {
    Idle,
    Progress { task_id: usize, sweeps: u64 },
    Complete { task_id: usize, results: Results },
    Timeup,
}
```

### 7. Main Entry Point

```rust
pub fn run_distributed<MC: MonteCarlo + FromParams>(
    config: DistributedConfig,
) -> Result<Vec<Results>, CarloError> {
    let world = MpiWorld::init()?;
    
    match world.rank() {
        0 => run_controller(world, config),
        _ => run_worker::<MC>(world, config),
    }
}

fn run_controller(world: MpiWorld, config: DistributedConfig) -> Result<Vec<Results>, CarloError>;
fn run_worker<MC: MonteCarlo + FromParams>(world: MpiWorld, config: DistributedConfig) -> Result<Vec<Results>, CarloError>;
```

## Error Handling

```rust
#[derive(Debug, thiserror::Error)]
pub enum MpiError {
    #[error("MPI initialization failed")]
    InitFailed,
    
    #[error("Communication error: {0}")]
    Communication(String),
    
    #[error("Serialization error: {0}")]
    Serialization(#[from] bincode::Error),
    
    #[error("Invalid state transition")]
    InvalidTransition,
}
```

## Testing Strategy

1. **Unit tests** - Mock channels, no MPI required
2. **Integration tests** - `mpirun -np 4 cargo test`
3. **Property tests** - Proptest for state machine invariants

## Implementation Order

1. `MpiSerializable` trait + derive macro
2. Channel abstraction (`MpiSender`, `MpiReceiver`)
3. Type-state `Worker` and `Controller`
4. `TaskStream` and `ResultsAggregator`
5. `Simulation` struct
6. Message types
7. `run_distributed` entry point
8. Tests

## Dependencies

- `mpi` (existing)
- `bincode` (for serialization) or custom MPI derive
- `crossbeam` (for channel-like API inspiration, optional)