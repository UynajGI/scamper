# Carlo.rs

A Rust implementation of the [Carlo.jl](https://github.com/lukas-weber/Carlo.jl) Monte Carlo simulation framework with **100% core feature parity**.

## Overview

Carlo.rs is a framework for developing high-performance, distributed Monte Carlo simulations. It handles model-independent tasks:

- **Autocorrelation and error analysis** - Automatic binning, jackknife resampling, and decorrelated autocorrelation time estimation
- **MPI scheduling** - Monte Carlo-aware distributed execution with controller-worker architecture
- **Checkpointing** - Resume simulations from saved HDF5 state
- **Complex observables** - Native support for complex number measurements
- **Performance monitoring** - Sweep rate tracking and elapsed time display

while leaving all flexibility of implementing Monte Carlo updates and estimators to you.

## Quick Start

Add to your `Cargo.toml`:

```toml
[dependencies]
carlo-rs = { path = "Carlo.rs" }
rand_xoshiro = "0.8"
```

Implement the `MonteCarlo` trait for your model:

```rust
use carlo_rs::{MonteCarlo, Context, CarloError, FromParams, Params, Scheduler, RunConfig, RayonBackend};
use rand_xoshiro::Xoshiro256PlusPlus;
use rand_core::Rng;

struct IsingMC {
    spins: Vec<i8>,
    beta: f64,
}

impl MonteCarlo for IsingMC {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        for i in 0..self.spins.len() {
            // Metropolis update
            let neighbor_sum = self.spins[(i + 1) % self.spins.len()]
                + self.spins[(i - 1 + self.spins.len()) % self.spins.len()];
            let dE = 2.0 * self.spins[i] as f64 * neighbor_sum as f64;
            if dE < 0.0 || ctx.rng.gen::<f64>() < (-dE * self.beta).exp() {
                self.spins[i] *= -1;
            }
        }
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        let m = self.spins.iter().sum::<i8>() as f64 / self.spins.len() as f64;
        ctx.measure("Magnetization", m.abs());
        // Complex observables are also supported:
        // ctx.measure_complex("OrderParameter", re, im);
    }
}

impl FromParams for IsingMC {
    fn from_params(params: &Params, _rng: &mut Self::Rng) -> Result<Self, CarloError> {
        Ok(Self {
            spins: vec![1; params.get::<usize>("L").unwrap_or(100)],
            beta: params.get::<f64>("beta").unwrap_or(1.0),
        })
    }
}
```

## Features

| Feature | Description | Dependencies |
|---------|-------------|--------------|
| `hdf5` | HDF5 checkpoint and measurement files | `libhdf5-dev` |
| `mpi` | MPI distributed backend | `libopenmpi-dev` |
| `strict-repro` | Strict reproducibility mode | None |

## Installation

### System Dependencies (Ubuntu/Debian)

```bash
sudo apt-get install libhdf5-dev openmpi-bin libopenmpi-dev
```

### Build

```bash
# Basic build
cargo build --release

# With HDF5 support
cargo build --release --features hdf5

# With MPI support
cargo build --release --features mpi

# Full feature set
cargo build --release --features "hdf5 mpi"
```

## CLI Usage

```bash
# Start simulation
carlo run --job-dir my_job/

# Check progress (with sweep rate and elapsed time)
carlo status --job-dir my_job/

# Merge results
carlo merge --job-dir my_job/

# Clean data
carlo delete --job-dir my_job/
```

## Development

```bash
# Quick check (format + lint + test)
just check

# Run tests
just test

# Run MPI tests
just test-mpi

# Generate docs
just doc

# Run benchmarks
just bench
```

## Architecture

- **`MonteCarlo` trait**: Core abstraction - implement `sweep()` for your algorithm
- **`Backend`**: Parallel execution (`RayonBackend` for threads, `MpiBackend` for distributed)
- **`Scheduler`**: Orchestrates thermalization → measurement phases
- **`Context`**: Runtime state (RNG, measurements, sweep counter)
- **`Merge`**: Rebinning, autocorrelation time estimation, covariance matrices
- **`Evaluator`**: Jackknife resampling for derived observables
- **`MultiplexEvaluator`**: Parallel tempering chain evaluations

## Carlo.jl Parity

| Feature | Carlo.jl | Carlo.rs |
|---------|----------|----------|
| Abstract MC trait | `AbstractMC` | `MonteCarlo` |
| Binning accumulator | `Accumulator` | `Accumulator` |
| Autocorrelation time | ✓ | ✓ |
| Decorrelated mode | ✓ | ✓ |
| Covariance estimation | ✓ | ✓ |
| MPI backend | ✓ | ✓ |
| Parallel tempering | ✓ | ✓ |
| HDF5 checkpoint | ✓ | ✓ |
| Jackknife resampling | ✓ | ✓ |
| ResultTools (dataframe, read-back) | ✓ | ✓ |
| Complex observables | (via ResultTools) | Native support |
| Progress bars | ✗ | ✓ |
| Performance monitoring | ✗ | ✓ |

## Documentation

Generate and view documentation:

```bash
cargo doc --workspace --no-deps --open
```

## License

Apache-2.0