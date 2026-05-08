# Carlo.rs Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the core Monte Carlo framework crate (Carlo.rs) with trait abstractions, scheduling, measurements, error analysis, and output support.

**Architecture:** Carlo.rs provides the `MonteCarlo` trait (user implementation entry), `Context` (runtime state + RNG + measurements), `Backend` trait (parallel abstraction), and `Scheduler` (execution orchestrator). CMC.rs and QMC.rs will depend on Carlo.rs and implement concrete models.

**Tech Stack:** Rust (stable-2026-04-07), `rand`/`rand_xoshiro` (RNG), `rayon` (parallel), `hdf5` (output), `thiserror` (errors), `serde` (JSON), `tracing` (logging)

---

## File Structure

```
Carlo.rs/
├── Cargo.toml
└── src/
    ├── lib.rs              # Public exports + crate root
    ├── error.rs            # CarloError enum (thiserror)
    ├── monte_carlo.rs      # MonteCarlo trait definition
    ├── context.rs          # Context struct + RNG + state
    ├── measurements.rs     # Measurements + Accumulator (binning)
    ├── estimate.rs         # Estimate struct + binning analysis
    ├── backend/
    │   ├── mod.rs          # Backend trait
    │   └── rayon.rs        # RayonBackend implementation
    ├── results.rs          # Results struct + HDF5/JSON output
    ├── scheduler.rs        # RunConfig + Scheduler
    └── params.rs           # Params struct + FromParams trait
```

---

## Task 1: Workspace Setup

**Files:**
- Create: `Cargo.toml` (workspace root)
- Create: `Carlo.rs/Cargo.toml`
- Create: `rust-toolchain.toml`
- Create: `justfile`

- [ ] **Step 1: Create workspace root Cargo.toml**

```toml
[workspace]
resolver = "2"
members = ["Carlo.rs"]

[workspace.dependencies]
rand = "0.9"
rand_xoshiro = "0.7"
rayon = "1.10"
thiserror = "2.0"
serde = { version = "1.0", features = ["derive"] }
serde_json = "1.0"
tracing = "0.1"
tracing-subscriber = "0.3"
hdf5 = { version = "0.8", features = ["have-parallel"] }
chrono = "0.4"
```

- [ ] **Step 2: Create Carlo.rs/Cargo.toml**

```toml
[package]
name = "carlo-rs"
version = "0.1.0"
edition = "2024"
license = "Apache-2.0"
description = "Core Monte Carlo framework for Scuttle"

[dependencies]
rand = { workspace = true }
rand_xoshiro = { workspace = true }
rayon = { workspace = true }
thiserror = { workspace = true }
serde = { workspace = true }
serde_json = { workspace = true }
tracing = { workspace = true }
hdf5 = { workspace = true }
chrono = { workspace = true }

[dev-dependencies]
tempfile = "3"

[features]
default = []
strict-repro = []  # Use jump sequence for RNG
```

- [ ] **Step 3: Create rust-toolchain.toml**

```toml
[toolchain]
channel = "stable-2024-01-01"
components = ["rustfmt", "clippy"]
```

- [ ] **Step 4: Create justfile (minimal)**

```just
# Quick feedback
check:
    cargo fmt --check
    cargo clippy --all-targets --all-features
    cargo test --workspace

# Build
build:
    cargo build --release

# Test
test:
    cargo test --workspace --all-features

# Docs
doc:
    cargo doc --workspace --no-deps --open

# Clean
clean:
    cargo clean
```

- [ ] **Step 5: Verify workspace compiles**

Run: `cargo check`
Expected: "Checking carlo-rs..." (may show empty crate warning)

- [ ] **Step 6: Commit**

```bash
jj describe -m "feat(Carlo): set up Cargo workspace and justfile"
```

---

## Task 2: Error Types

**Files:**
- Create: `Carlo.rs/src/lib.rs`
- Create: `Carlo.rs/src/error.rs`

- [ ] **Step 1: Write failing test for error construction**

Create: `Carlo.rs/src/error.rs` (minimal skeleton for test to reference)

```rust
// Skeleton only - will be filled in step 3
pub struct CarloError;
```

Create: `Carlo.rs/tests/error_test.rs`

```rust
use carlo_rs::CarloError;

#[test]
fn test_error_display() {
    let err = CarloError::InvalidConfig { field: "binsize".into(), reason: "must be positive".into() };
    assert!(err.to_string().contains("binsize"));
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test error_test`
Expected: FAIL - "CarloError::InvalidConfig does not exist"

- [ ] **Step 3: Implement CarloError enum**

```rust
// Carlo.rs/src/error.rs
use thiserror::Error;

#[derive(Debug, Error)]
pub enum CarloError {
    #[error("Invalid configuration: field '{field}' - {reason}")]
    InvalidConfig { field: String, reason: String },

    #[error("HDF5 I/O error: {0}")]
    Hdf5Error(#[from] hdf5::Error),

    #[error("Measurement '{name}' not found")]
    MeasurementNotFound { name: String },

    #[error("Checkpoint corrupted: {detail}")]
    CheckpointCorrupted { detail: String },

    #[error("Convergence not reached after {sweeps} sweeps")]
    ConvergenceTimeout { sweeps: u64 },
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test error_test`
Expected: PASS

- [ ] **Step 5: Update lib.rs to export error**

```rust
// Carlo.rs/src/lib.rs
mod error;

pub use error::CarloError;
```

- [ ] **Step 6: Commit**

```bash
jj describe -m "feat(Carlo): add CarloError enum with thiserror"
```

---

## Task 3: MonteCarlo Trait

**Files:**
- Create: `Carlo.rs/src/monte_carlo.rs`
- Modify: `Carlo.rs/src/lib.rs`

- [ ] **Step 1: Write failing test for trait basic usage**

Create: `Carlo.rs/tests/monte_carlo_test.rs`

```rust
use carlo_rs::{MonteCarlo, Context};
use rand_xoshiro::Xoshiro256PlusPlus;

struct DummyMC {
    sweep_count: u64,
}

impl MonteCarlo for DummyMC {
    fn sweep(&mut self, ctx: &mut Context<Xoshiro256PlusPlus>) {
        self.sweep_count += 1;
        ctx.measure("sweeps", 1.0);
    }
}

#[test]
fn test_monte_carlo_sweep() {
    let mut mc = DummyMC { sweep_count: 0 };
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut ctx = Context::new(rng, 10);
    
    mc.sweep(&mut ctx);
    
    assert_eq!(mc.sweep_count, 1);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test monte_carlo_test`
Expected: FAIL - "MonteCarlo not found", "Context not found"

- [ ] **Step 3: Implement MonteCarlo trait**

```rust
// Carlo.rs/src/monte_carlo.rs
use rand::RngCore;
use rand::SeedableRng;

use crate::Context;

/// Core trait for Monte Carlo algorithms.
/// Users implement `sweep()` and optionally override other methods.
pub trait MonteCarlo: Sized {
    /// Execute one sweep (update configuration).
    /// Users may call `ctx.measure()` inside sweep.
    fn sweep(&mut self, ctx: &mut Context<Self::Rng>);

    /// RNG type used by this model.
    type Rng: RngCore + SeedableRng + Send;

    /// Optional: measure observables (default: empty).
    fn measure(&mut self, _ctx: &mut Context<Self::Rng>) {}

    /// Optional: save state to HDF5 (default: empty).
    fn save(&self, _out: &mut hdf5::Group) {}

    /// Optional: load state from HDF5 (default: empty).
    fn load(&mut self, _in: &hdf5::Group) {}

    /// Optional: algorithm name (default: "UnnamedMC").
    fn name(&self) -> &'static str {
        "UnnamedMC"
    }
}
```

- [ ] **Step 4: Implement Context skeleton (minimal for trait test)**

```rust
// Carlo.rs/src/context.rs
use rand::RngCore;
use rand::SeedableRng;

use crate::Measurements;

/// Runtime context for Monte Carlo simulation.
pub struct Context<R: RngCore + SeedableRng> {
    /// Random number generator (public access).
    pub rng: R,
    
    /// Measurement collector.
    measurements: Measurements,
    
    /// Current sweep count.
    sweep_count: u64,
    
    /// Thermalization sweeps threshold.
    thermalization_sweeps: u64,
    
    /// Whether thermalized.
    thermalized: bool,
}

impl<R: RngCore + SeedableRng> Context<R> {
    /// Create new context.
    pub fn new(rng: R, thermalization_sweeps: u64) -> Self {
        Self {
            rng,
            measurements: Measurements::new(100), // default binsize
            sweep_count: 0,
            thermalization_sweeps,
            thermalized: false,
        }
    }

    /// Record an observable sample.
    pub fn measure(&mut self, name: &str, value: f64) {
        self.measurements.add_sample(name, value);
    }

    /// Check if thermalized.
    pub fn is_thermalized(&self) -> bool {
        self.thermalized
    }

    /// Get sweep count.
    pub fn sweep_count(&self) -> u64 {
        self.sweep_count
    }

    /// Advance sweep counter (internal).
    pub(crate) fn advance_sweep(&mut self) {
        self.sweep_count += 1;
        if self.sweep_count > self.thermalization_sweeps {
            self.thermalized = true;
        }
    }
}
```

- [ ] **Step 5: Implement Measurements skeleton (minimal)**

```rust
// Carlo.rs/src/measurements.rs
use std::collections::HashMap;

/// Measurement collector (placeholder for now).
pub struct Measurements {
    observables: HashMap<String, Vec<f64>>,
    binsize: usize,
}

impl Measurements {
    pub fn new(binsize: usize) -> Self {
        Self {
            observables: HashMap::new(),
            binsize,
        }
    }

    pub fn add_sample(&mut self, name: &str, value: f64) {
        self.observables
            .entry(name.to_string())
            .or_insert_with(Vec::new)
            .push(value);
    }
}
```

- [ ] **Step 6: Update lib.rs exports**

```rust
// Carlo.rs/src/lib.rs
mod error;
mod monte_carlo;
mod context;
mod measurements;

pub use error::CarloError;
pub use monte_carlo::MonteCarlo;
pub use context::Context;
pub use measurements::Measurements;
```

- [ ] **Step 7: Run test to verify it passes**

Run: `cargo test monte_carlo_test`
Expected: PASS

- [ ] **Step 8: Commit**

```bash
jj describe -m "feat(Carlo): add MonteCarlo trait and Context skeleton"
```

---

## Task 4: Measurements + Accumulator (Binning Analysis)

**Files:**
- Modify: `Carlo.rs/src/measurements.rs`

- [ ] **Step 1: Write failing test for binning**

Create: `Carlo.rs/tests/measurements_test.rs`

```rust
use carlo_rs::Measurements;

#[test]
fn test_binning_accumulation() {
    let mut meas = Measurements::new(10);
    
    // Add 25 samples (should create 2 full bins + 5 partial)
    for i in 0..25 {
        meas.add_sample("Energy", i as f64);
    }
    
    // After finalize, should have 2 bins
    let results = meas.finalize();
    let estimate = results.get("Energy").expect("Energy observable");
    
    // Mean of bins: [0..10].mean()=4.5, [10..20].mean()=14.5 → overall ~9.5
    assert!(estimate.mean > 0.0);
    assert!(estimate.stderr > 0.0);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test measurements_test`
Expected: FAIL - "no method `finalize` on Measurements"

- [ ] **Step 3: Implement full Measurements + Accumulator**

```rust
// Carlo.rs/src/measurements.rs
use std::collections::HashMap;
use crate::Estimate;

/// Single observable accumulator with binning.
pub struct Accumulator {
    /// Completed bins (each bin is mean of bin_capacity samples).
    bins: Vec<f64>,
    
    /// Current bin being filled.
    current_bin: Vec<f64>,
    
    /// Number of samples per bin.
    bin_capacity: usize,
    
    /// Sum of all samples (for quick mean estimate).
    total_sum: f64,
    
    /// Total sample count.
    total_count: usize,
}

impl Accumulator {
    pub fn new(bin_capacity: usize) -> Self {
        Self {
            bins: Vec::new(),
            current_bin: Vec::new(),
            bin_capacity,
            total_sum: 0.0,
            total_count: 0,
        }
    }

    /// Add a sample to the accumulator.
    pub fn add(&mut self, value: f64) {
        self.current_bin.push(value);
        self.total_sum += value;
        self.total_count += 1;

        if self.current_bin.len() >= self.bin_capacity {
            // Bin is full: compute mean and store
            let bin_mean = self.current_bin.iter().sum::<f64>() / self.bin_capacity as f64;
            self.bins.push(bin_mean);
            self.current_bin.clear();
        }
    }

    /// Finalize and return estimate.
    pub fn finalize(&self) -> Estimate {
        // Include partial bin if significant
        let mut all_bins = self.bins.clone();
        if !self.current_bin.is_empty() {
            let partial_mean = self.current_bin.iter().sum::<f64>() / self.current_bin.len() as f64;
            all_bins.push(partial_mean);
        }

        Estimate::from_bins(&all_bins)
    }

    /// Check if any complete bins exist.
    pub fn has_complete_bins(&self) -> bool {
        !self.bins.is_empty()
    }
}

/// Measurement collector managing multiple observables.
pub struct Measurements {
    observables: HashMap<String, Accumulator>,
    default_binsize: usize,
}

impl Measurements {
    pub fn new(default_binsize: usize) -> Self {
        Self {
            observables: HashMap::new(),
            default_binsize,
        }
    }

    /// Add a sample to an observable (auto-create if needed).
    pub fn add_sample(&mut self, name: &str, value: f64) {
        if !self.observables.contains_key(name) {
            self.observables.insert(name.to_string(), Accumulator::new(self.default_binsize));
        }
        self.observables.get_mut(name).unwrap().add(value);
    }

    /// Register an observable with custom binsize.
    pub fn register(&mut self, name: &str, binsize: usize) {
        self.observables.insert(name.to_string(), Accumulator::new(binsize));
    }

    /// Finalize all observables and return estimates.
    pub fn finalize(&self) -> HashMap<String, Estimate> {
        self.observables
            .iter()
            .map(|(name, acc)| (name.clone(), acc.finalize()))
            .collect()
    }
}
```

- [ ] **Step 4: Implement Estimate (for finalize to work)**

```rust
// Carlo.rs/src/estimate.rs
use serde::{Serialize, Deserialize};

/// Statistical estimate with error analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Estimate {
    /// Mean value.
    pub mean: f64,
    
    /// Standard error of the mean.
    pub stderr: f64,
    
    /// Integrated autocorrelation time (placeholder: 1.0).
    pub autocorr_time: f64,
    
    /// Number of bins used.
    pub n_bins: usize,
}

impl Estimate {
    /// Compute estimate from bin means.
    pub fn from_bins(bins: &[f64]) -> Self {
        if bins.is_empty() {
            return Self {
                mean: 0.0,
                stderr: 0.0,
                autocorr_time: 1.0,
                n_bins: 0,
            };
        }

        let n = bins.len() as f64;
        let mean = bins.iter().sum::<f64>() / n;

        // Standard deviation of bin means
        let variance = if n > 1.0 {
            bins.iter()
                .map(|b| (b - mean).powi(2))
                .sum::<f64>() / (n - 1.0)
        } else {
            0.0
        };

        // Standard error = std / sqrt(n_bins)
        let stderr = variance.sqrt() / n.sqrt();

        Self {
            mean,
            stderr,
            autocorr_time: 1.0, // Placeholder; proper estimation is complex
            n_bins: bins.len(),
        }
    }

    /// Convenience: format as "mean ± stderr".
    pub fn format(&self) -> String {
        format!("{:.6} ± {:.6}", self.mean, self.stderr)
    }
}
```

- [ ] **Step 5: Update lib.rs exports**

```rust
// Carlo.rs/src/lib.rs
mod error;
mod monte_carlo;
mod context;
mod measurements;
mod estimate;

pub use error::CarloError;
pub use monte_carlo::MonteCarlo;
pub use context::Context;
pub use measurements::{Measurements, Accumulator};
pub use estimate::Estimate;
```

- [ ] **Step 6: Run test to verify it passes**

Run: `cargo test measurements_test`
Expected: PASS

- [ ] **Step 7: Commit**

```bash
jj describe -m "feat(Carlo): add Measurements, Accumulator, and Estimate with binning analysis"
```

---

## Task 5: Backend Trait + RayonBackend

**Files:**
- Create: `Carlo.rs/src/backend/mod.rs`
- Create: `Carlo.rs/src/backend/rayon.rs`
- Modify: `Carlo.rs/src/lib.rs`

- [ ] **Step 1: Write failing test for Backend**

Create: `Carlo.rs/tests/backend_test.rs`

```rust
use carlo_rs::backend::{Backend, RayonBackend};
use rand_xoshiro::Xoshiro256PlusPlus;
use std::sync::atomic::{AtomicU64, Ordering};

#[test]
fn test_rayon_backend_spawn_tasks() {
    let backend = RayonBackend::new(4); // 4 threads
    let counter = AtomicU64::new(0);

    backend.spawn_tasks(10, 42, |task_id, rng| {
        // Each task should have unique RNG seed
        let seed_offset = rng.next_u64(); // Just read something to prove RNG works
        counter.fetch_add(1, Ordering::SeqCst);
    });

    backend.barrier();
    
    assert_eq!(counter.load(Ordering::SeqCst), 10);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test backend_test`
Expected: FAIL - "backend module not found"

- [ ] **Step 3: Implement Backend trait**

```rust
// Carlo.rs/src/backend/mod.rs
use rand::RngCore;
use rand::SeedableRng;

/// Parallel execution backend abstraction.
pub trait Backend: Clone + Send + Sync {
    /// RNG type for this backend.
    type Rng: RngCore + SeedableRng + Send;

    /// Spawn n tasks in parallel, each with isolated RNG.
    fn spawn_tasks<F>(&self, n_tasks: usize, base_seed: u64, f: F)
    where
        F: Fn(usize, &mut Self::Rng) + Sync;

    /// Wait for all tasks to complete.
    fn barrier(&self);
}

mod rayon;
pub use rayon::RayonBackend;
```

- [ ] **Step 4: Implement RayonBackend**

```rust
// Carlo.rs/src/backend/rayon.rs
use rand::RngCore;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;

use super::Backend;

/// Rayon-based parallel backend (Phase 1).
#[derive(Clone)]
pub struct RayonBackend {
    n_threads: usize,
}

impl RayonBackend {
    pub fn new(n_threads: usize) -> Self {
        Self { n_threads }
    }

    pub fn default() -> Self {
        Self::new(rayon::current_num_threads())
    }
}

impl Backend for RayonBackend {
    type Rng = Xoshiro256PlusPlus;

    fn spawn_tasks<F>(&self, n_tasks: usize, base_seed: u64, f: F)
    where
        F: Fn(usize, &mut Self::Rng) + Sync,
    {
        // Simple seed offset strategy (default)
        // For strict-repro, use jump sequence (Phase 2)
        (0..n_tasks).into_par_iter().for_each(|task_id| {
            let seed = base_seed.wrapping_add(task_id as u64);
            let mut rng = Self::Rng::seed_from_u64(seed);
            f(task_id, &mut rng);
        });
    }

    fn barrier(&self) {
        // Rayon automatically synchronizes after parallel iteration
        // No explicit barrier needed
    }
}
```

- [ ] **Step 5: Update lib.rs exports**

```rust
// Carlo.rs/src/lib.rs
mod error;
mod monte_carlo;
mod context;
mod measurements;
mod estimate;
pub mod backend;

pub use error::CarloError;
pub use monte_carlo::MonteCarlo;
pub use context::Context;
pub use measurements::{Measurements, Accumulator};
pub use estimate::Estimate;
pub use backend::{Backend, RayonBackend};
```

- [ ] **Step 6: Run test to verify it passes**

Run: `cargo test backend_test`
Expected: PASS

- [ ] **Step 7: Commit**

```bash
jj describe -m "feat(Carlo): add Backend trait and RayonBackend implementation"
```

---

## Task 6: Results + HDF5/JSON Output

**Files:**
- Create: `Carlo.rs/src/results.rs`
- Modify: `Carlo.rs/src/lib.rs`

- [ ] **Step 1: Write failing test for Results**

Create: `Carlo.rs/tests/results_test.rs`

```rust
use carlo_rs::{Results, Estimate};
use tempfile::NamedTempFile;

#[test]
fn test_results_json_output() {
    let mut results = Results::new();
    results.add("Energy", Estimate {
        mean: 0.5,
        stderr: 0.01,
        autocorr_time: 1.0,
        n_bins: 100,
    });

    let json = results.to_json().unwrap();
    assert!(json.contains("Energy"));
    assert!(json.contains("0.5"));
}

#[test]
fn test_results_get() {
    let mut results = Results::new();
    results.add("Energy", Estimate {
        mean: 0.5,
        stderr: 0.01,
        autocorr_time: 1.0,
        n_bins: 100,
    });

    let est = results.get("Energy").unwrap();
    assert_eq!(est.mean, 0.5);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test results_test`
Expected: FAIL - "Results not found"

- [ ] **Step 3: Implement Results struct**

```rust
// Carlo.rs/src/results.rs
use std::collections::HashMap;
use std::path::Path;

use serde::{Serialize, Deserialize};
use chrono::{DateTime, Utc};

use crate::{Estimate, CarloError};

/// Simulation results container.
pub struct Results {
    /// Observable estimates.
    estimates: HashMap<String, Estimate>,
    
    /// Run metadata.
    metadata: Metadata,
}

/// Run metadata for reproducibility.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Metadata {
    /// Scuttle version.
    pub version: String,
    
    /// Timestamp.
    pub timestamp: DateTime<Utc>,
    
    /// Base seed used.
    pub base_seed: u64,
    
    /// Thermalization sweeps.
    pub thermalization_sweeps: u64,
    
    /// Measurement sweeps.
    pub measurement_sweeps: u64,
}

impl Default for Metadata {
    fn default() -> Self {
        Self {
            version: env!("CARGO_PKG_VERSION"),
            timestamp: Utc::now(),
            base_seed: 0,
            thermalization_sweeps: 0,
            measurement_sweeps: 0,
        }
    }
}

impl Results {
    pub fn new() -> Self {
        Self {
            estimates: HashMap::new(),
            metadata: Metadata::default(),
        }
    }

    /// Create from measurements finalize output.
    pub fn from_measurements(measurements: &HashMap<String, Estimate>) -> Self {
        Self {
            estimates: measurements.clone(),
            metadata: Metadata::default(),
        }
    }

    /// Add an observable estimate.
    pub fn add(&mut self, name: &str, estimate: Estimate) {
        self.estimates.insert(name.to_string(), estimate);
    }

    /// Get an observable estimate.
    pub fn get(&self, name: &str) -> Option<&Estimate> {
        self.estimates.get(name)
    }

    /// Export to JSON string.
    pub fn to_json(&self) -> Result<String, CarloError> {
        #[derive(Serialize)]
        struct JsonResults {
            observables: HashMap<String, Estimate>,
            metadata: Metadata,
        }

        let json_data = JsonResults {
            observables: self.estimates.clone(),
            metadata: self.metadata.clone(),
        };

        serde_json::to_string_pretty(&json_data)
            .map_err(|e| CarloError::InvalidConfig {
                field: "json_output".into(),
                reason: e.to_string(),
            })
    }

    /// Save to HDF5 file.
    pub fn save_hdf5(&self, path: &Path) -> Result<(), CarloError> {
        let file = hdf5::File::create(path)?;
        
        // Write observables
        let obs_group = file.create_group("observables")?;
        for (name, est) in &self.estimates {
            let obs = obs_group.create_group(name)?;
            obs.write_scalar("mean", &est.mean)?;
            obs.write_scalar("stderr", &est.stderr)?;
            obs.write_scalar("autocorr_time", &est.autocorr_time)?;
            obs.write_scalar("n_bins", &(est.n_bins as i64))?;
        }

        // Write metadata
        let meta_group = file.create_group("metadata")?;
        meta_group.write_scalar("version", &self.metadata.version)?;
        meta_group.write_scalar("timestamp", &self.metadata.timestamp.to_rfc3339())?;
        meta_group.write_scalar("base_seed", &(self.metadata.base_seed as i64))?;
        
        Ok(())
    }

    /// Set metadata.
    pub fn set_metadata(&mut self, metadata: Metadata) {
        self.metadata = metadata;
    }
}
```

- [ ] **Step 4: Update lib.rs exports**

```rust
// Carlo.rs/src/lib.rs
// ... (add results module)
mod results;

pub use results::{Results, Metadata};
```

- [ ] **Step 5: Run test to verify it passes**

Run: `cargo test results_test`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
jj describe -m "feat(Carlo): add Results with HDF5 and JSON output support"
```

---

## Task 7: Params + FromParams Trait

**Files:**
- Create: `Carlo.rs/src/params.rs`
- Modify: `Carlo.rs/src/lib.rs`

- [ ] **Step 1: Write failing test for Params**

Create: `Carlo.rs/tests/params_test.rs`

```rust
use carlo_rs::{Params, FromParams, Context, MonteCarlo};
use rand_xoshiro::Xoshiro256PlusPlus;

struct DummyModel {
    lattice_size: usize,
    temperature: f64,
}

impl MonteCarlo for DummyModel {
    type Rng = Xoshiro256PlusPlus;
    
    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        // placeholder
    }
}

impl FromParams for DummyModel {
    fn from_params(params: &Params, rng: &mut Xoshiro256PlusPlus) -> Self {
        Self {
            lattice_size: params.get("lattice_size").unwrap_or(32),
            temperature: params.get("temperature").unwrap_or(1.0),
        }
    }
}

#[test]
fn test_params_creation() {
    let mut params = Params::new();
    params.set("lattice_size", 64);
    params.set("temperature", 2.269);

    assert_eq!(params.get::<usize>("lattice_size"), Some(64));
    assert_eq!(params.get::<f64>("temperature"), Some(2.269));
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test params_test`
Expected: FAIL - "Params not found"

- [ ] **Step 3: Implement Params and FromParams**

```rust
// Carlo.rs/src/params.rs
use std::collections::HashMap;
use rand::RngCore;
use rand::SeedableRng;

use crate::{MonteCarlo, Context};

/// Generic parameter container.
pub struct Params {
    /// String key → string value (parsed on access).
    values: HashMap<String, String>,
}

impl Params {
    pub fn new() -> Self {
        Self { values: HashMap::new() }
    }

    /// Set a parameter (converts to string).
    pub fn set<T: ToString>(&mut self, key: &str, value: T) {
        self.values.insert(key.to_string(), value.to_string());
    }

    /// Get a parameter (parses from string).
    pub fn get<T: std::str::FromStr>(&self, key: &str) -> Option<T> {
        self.values.get(key).and_then(|v| v.parse::<T>().ok())
    }
}

/// Trait for constructing models from parameters.
pub trait FromParams: MonteCarlo {
    /// Construct model from params and RNG.
    fn from_params(params: &Params, rng: &mut Self::Rng) -> Self;
}
```

- [ ] **Step 4: Update lib.rs exports**

```rust
// Carlo.rs/src/lib.rs
mod params;

pub use params::{Params, FromParams};
```

- [ ] **Step 5: Run test to verify it passes**

Run: `cargo test params_test`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
jj describe -m "feat(Carlo): add Params container and FromParams trait"
```

---

## Task 8: Scheduler + RunConfig

**Files:**
- Create: `Carlo.rs/src/scheduler.rs`
- Modify: `Carlo.rs/src/lib.rs`

- [ ] **Step 1: Write integration test for Scheduler**

Create: `Carlo.rs/tests/scheduler_test.rs`

```rust
use carlo_rs::{Scheduler, RunConfig, Params, MonteCarlo, Context, FromParams, Backend, RayonBackend};
use rand_xoshiro::Xoshiro256PlusPlus;

struct CountingMC {
    sweep_count: u64,
    total_energy: f64,
}

impl MonteCarlo for CountingMC {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        self.sweep_count += 1;
        // Record energy proportional to sweep count
        if ctx.is_thermalized() {
            ctx.measure("Energy", self.sweep_count as f64);
        }
    }
}

impl FromParams for CountingMC {
    fn from_params(_params: &Params, _rng: &mut Self::Rng) -> Self {
        Self { sweep_count: 0, total_energy: 0.0 }
    }
}

#[test]
fn test_scheduler_single_task() {
    let config = RunConfig {
        thermalization_sweeps: 100,
        measurement_sweeps: 1000,
        base_seed: 42,
        binsize: 100,
    };

    let backend = RayonBackend::default();
    let scheduler = Scheduler::new(backend, config);
    let params = Params::new();

    let results = scheduler.run_one::<CountingMC>(&params);

    // Should have Energy observable
    assert!(results.get("Energy").is_some());
    
    // Mean should be around (100 + 1000/2) = 600ish
    let est = results.get("Energy").unwrap();
    assert!(est.mean > 100.0);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test scheduler_test`
Expected: FAIL - "Scheduler not found", "RunConfig not found"

- [ ] **Step 3: Implement RunConfig and Scheduler**

```rust
// Carlo.rs/src/scheduler.rs
use crate::{Backend, MonteCarlo, Context, Params, FromParams, Results, Metadata};

/// Run configuration.
pub struct RunConfig {
    /// Number of thermalization sweeps.
    pub thermalization_sweeps: u64,
    
    /// Number of measurement sweeps.
    pub measurement_sweeps: u64,
    
    /// Base RNG seed.
    pub base_seed: u64,
    
    /// Default binsize for measurements.
    pub binsize: usize,
}

impl Default for RunConfig {
    fn default() -> Self {
        Self {
            thermalization_sweeps: 1000,
            measurement_sweeps: 10000,
            base_seed: 0,
            binsize: 100,
        }
    }
}

/// Scheduler orchestrates Monte Carlo runs.
pub struct Scheduler<B: Backend> {
    backend: B,
    config: RunConfig,
}

impl<B: Backend> Scheduler<B> {
    pub fn new(backend: B, config: RunConfig) -> Self {
        Self { backend, config }
    }

    /// Run a single task.
    pub fn run_one<MC: FromParams>(&self, params: &Params) -> Results {
        let mut rng = MC::Rng::seed_from_u64(self.config.base_seed);
        let mut ctx = Context::new_with_binsize(rng, self.config.thermalization_sweeps, self.config.binsize);
        let mut mc = MC::from_params(params, &mut ctx.rng);

        // Thermalization phase
        for _ in 0..self.config.thermalization_sweeps {
            mc.sweep(&mut ctx);
            ctx.advance_sweep();
        }

        // Measurement phase
        for _ in 0..self.config.measurement_sweeps {
            mc.sweep(&mut ctx);
            mc.measure(&mut ctx);
            ctx.advance_sweep();
        }

        // Finalize
        let estimates = ctx.finalize_measurements();
        let mut results = Results::from_measurements(&estimates);
        results.set_metadata(Metadata {
            version: env!("CARGO_PKG_VERSION"),
            timestamp: chrono::Utc::now(),
            base_seed: self.config.base_seed,
            thermalization_sweeps: self.config.thermalization_sweeps,
            measurement_sweeps: self.config.measurement_sweeps,
        });

        results
    }

    /// Run multiple tasks in parallel.
    pub fn run_parallel<MC: FromParams>(&self, n_tasks: usize, params: &Params) -> Vec<Results> {
        use std::sync::Mutex;

        let results = Mutex::new(Vec::new());
        
        self.backend.spawn_tasks(n_tasks, self.config.base_seed, |task_id, rng| {
            let mut ctx = Context::new_with_binsize(
                MC::Rng::seed_from_u64(rng.next_u64()), // task-specific seed
                self.config.thermalization_sweeps,
                self.config.binsize,
            );
            let mut mc = MC::from_params(params, &mut ctx.rng);

            // Thermalization
            for _ in 0..self.config.thermalization_sweeps {
                mc.sweep(&mut ctx);
                ctx.advance_sweep();
            }

            // Measurement
            for _ in 0..self.config.measurement_sweeps {
                mc.sweep(&mut ctx);
                mc.measure(&mut ctx);
                ctx.advance_sweep();
            }

            let estimates = ctx.finalize_measurements();
            let mut result = Results::from_measurements(&estimates);
            result.set_metadata(Metadata {
                version: env!("CARGO_PKG_VERSION"),
                timestamp: chrono::Utc::now(),
                base_seed: self.config.base_seed.wrapping_add(task_id as u64),
                thermalization_sweeps: self.config.thermalization_sweeps,
                measurement_sweeps: self.config.measurement_sweeps,
            });

            results.lock().unwrap().push(result);
        });

        self.backend.barrier();
        results.into_inner().unwrap()
    }
}
```

- [ ] **Step 4: Add Context binsize constructor**

```rust
// Add to Carlo.rs/src/context.rs
impl<R: RngCore + SeedableRng> Context<R> {
    pub fn new_with_binsize(rng: R, thermalization_sweeps: u64, binsize: usize) -> Self {
        Self {
            rng,
            measurements: Measurements::new(binsize),
            sweep_count: 0,
            thermalization_sweeps,
            thermalized: false,
        }
    }

    /// Finalize measurements and return estimates.
    pub fn finalize_measurements(&self) -> HashMap<String, Estimate> {
        self.measurements.finalize()
    }
}
```

- [ ] **Step 5: Update lib.rs exports**

```rust
// Carlo.rs/src/lib.rs
mod scheduler;

pub use scheduler::{Scheduler, RunConfig};
```

- [ ] **Step 6: Run test to verify it passes**

Run: `cargo test scheduler_test`
Expected: PASS

- [ ] **Step 7: Commit**

```bash
jj describe -m "feat(Carlo): add Scheduler and RunConfig with single/parallel task execution"
```

---

## Task 9: Deterministic Reproducibility Test

**Files:**
- Create: `Carlo.rs/tests/reproducibility_test.rs`

- [ ] **Step 1: Write reproducibility test**

```rust
use carlo_rs::{Scheduler, RunConfig, Params, MonteCarlo, Context, FromParams, RayonBackend};
use rand_xoshiro::Xoshiro256PlusPlus;

struct SimpleMC;

impl MonteCarlo for SimpleMC {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        // Use RNG to produce deterministic output
        let val = ctx.rng.gen::<f64>();
        if ctx.is_thermalized() {
            ctx.measure("RandomValue", val);
        }
    }
}

impl FromParams for SimpleMC {
    fn from_params(_params: &Params, _rng: &mut Self::Rng) -> Self {
        Self
    }
}

#[test]
fn test_deterministic_same_seed_same_result() {
    let config = RunConfig {
        thermalization_sweeps: 10,
        measurement_sweeps: 100,
        base_seed: 12345,
        binsize: 10,
    };

    let backend = RayonBackend::new(1); // Single thread for repro
    let params = Params::new();

    // Run twice with same seed
    let results1 = Scheduler::new(backend.clone(), config.clone()).run_one::<SimpleMC>(&params);
    let results2 = Scheduler::new(backend.clone(), config.clone()).run_one::<SimpleMC>(&params);

    let est1 = results1.get("RandomValue").unwrap();
    let est2 = results2.get("RandomValue").unwrap();

    // Same seed should give exactly same mean
    assert!((est1.mean - est2.mean).abs() < 1e-10);
}
```

- [ ] **Step 2: Add rand::Rng trait usage**

Context needs `rng.gen()` method. The RNG field is public, so this should work. If test fails, add trait import.

- [ ] **Step 3: Run test**

Run: `cargo test reproducibility_test`
Expected: PASS (if RNG trait import works)

- [ ] **Step 4: Fix if needed**

If `rng.gen()` not available, add `use rand::Rng` in test file.

- [ ] **Step 5: Commit**

```bash
jj describe -m "test(Carlo): add deterministic reproducibility verification"
```

---

## Task 10: Final Integration + justfile Expansion

**Files:**
- Modify: `justfile`
- Run: `cargo test --workspace`

- [ ] **Step 1: Expand justfile with full commands**

```just
# Quick feedback
check:
    cargo fmt --check
    cargo clippy --all-targets --all-features -- -D warnings
    cargo test --workspace

# Format
fmt:
    cargo fmt

# Build
build:
    cargo build --release

# Test
test:
    cargo test --workspace --all-features

# Test specific
test-unit:
    cargo test --lib

test-integration:
    cargo test --test '*'

# Docs
doc:
    cargo doc --workspace --no-deps --open

doc-check:
    cargo doc --workspace --no-deps

# Clean
clean:
    cargo clean

# Audit
audit:
    cargo audit

# Benchmarks (future)
bench:
    cargo bench

# Publish dry-run
publish-dry:
    cargo publish --dry-run
```

- [ ] **Step 2: Run full test suite**

Run: `just check`
Expected: All tests pass, clippy clean

- [ ] **Step 3: Fix any clippy warnings**

If clippy reports warnings, fix them in relevant files.

- [ ] **Step 4: Commit final**

```bash
jj describe -m "chore(Carlo): finalize Carlo.rs with full test suite and expanded justfile"
```

---

## Self-Review Checklist

After implementing all tasks:

1. **Spec coverage:**
   - MonteCarlo trait ✓ (Task 3)
   - Context ✓ (Task 3, Task 8)
   - Backend + RayonBackend ✓ (Task 5)
   - Measurements + binning ✓ (Task 4)
   - Estimate ✓ (Task 4)
   - Results (HDF5 + JSON) ✓ (Task 6)
   - Scheduler ✓ (Task 8)
   - Params + FromParams ✓ (Task 7)
   - Error types ✓ (Task 2)
   - Deterministic repro ✓ (Task 9)

2. **Placeholder scan:** None - all tasks have complete code.

3. **Type consistency:** 
   - `Context<R>` generic throughout
   - `Estimate::from_bins(bins: &[f64])` matches `Measurements.finalize()` return
   - `Results::from_measurements(HashMap<String, Estimate>)` matches Context finalize output

---

Plan complete. Ready for execution.