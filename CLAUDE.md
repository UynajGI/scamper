# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Scuttle is a Monte Carlo simulation framework for (quantum) Monte Carlo algorithms. It consists of:
- **Carlo.rs** (Rust): Core framework implementation being developed here
- **Carlo.jl** (Julia): Reference implementation from [lukas-weber/Carlo.jl](https://github.com/lukas-weber/Carlo.jl)
- **StochasticSeriesExpansion.jl** (Julia): SSE QMC algorithm built on Carlo.jl

The framework handles model-independent tasks: autocorrelation/error analysis, MPI scheduling, checkpointing, while leaving MC update/estimator implementation to users.

## Build Commands

```bash
# Quick check (format + lint + test)
just check

# Build release
just build

# Build with MPI support (requires MPI installed)
just build-mpi

# Run all tests
just test

# Run unit tests only
just test-unit

# Run integration tests only
just test-integration

# Run MPI tests (requires mpirun)
just test-mpi

# Run benchmarks
just bench

# Install system dependencies (Ubuntu/Debian)
just install-deps
```

## Architecture

### Core Trait: MonteCarlo

The `MonteCarlo` trait ([Carlo.rs/src/monte_carlo.rs](Carlo.rs/src/monte_carlo.rs)) is the central abstraction:

```rust
pub trait MonteCarlo: Sized {
    type Rng: Rng + SeedableRng + Send;
    fn sweep(&mut self, ctx: &mut Context<Self::Rng>);
    fn measure(&mut self, _ctx: &mut Context<Self::Rng>) {}
}
```

Users implement `sweep()` for configuration updates and optionally `measure()` for observables. The `FromParams` trait constructs models from parameter dictionaries.

### Execution Pipeline

1. **Backend** ([Carlo.rs/src/backend/mod.rs](Carlo.rs/src/backend/mod.rs)): Parallel execution abstraction
   - `RayonBackend`: Thread-parallel (default)
   - `MpiBackend`: Distributed MPI (requires `mpi` feature)

2. **Scheduler** ([Carlo.rs/src/scheduler.rs](Carlo.rs/src/scheduler.rs)): Orchestrates runs
   - Thermalization phase → Measurement phase
   - `run_one()` for single task, `run_parallel()` for multiple

3. **Context** ([Carlo.rs/src/context.rs](Carlo.rs/src/context.rs)): Runtime state
   - Holds RNG, measurements accumulator, sweep counter
   - `ctx.measure(name, value)` records observables
   - `ctx.is_thermalized()` checks thermalization status

### Results and Analysis

- **Measurements** ([Carlo.rs/src/measurements.rs](Carlo.rs/src/measurements.rs)): Binned accumulation during simulation
- **Merge** ([Carlo.rs/src/merge.rs](Carlo.rs/src/merge.rs)): Rebinning and autocorrelation time estimation after simulation
- **Evaluable** ([Carlo.rs/src/evaluable.rs](Carlo.rs/src/evaluable.rs)): Jackknife analysis for derived observables

### CLI Commands

The CLI ([Carlo.rs/src/cli.rs](Carlo.rs/src/cli.rs)) provides:
- `carlo run`: Start simulation
- `carlo status`: Check progress
- `carlo merge`: Combine results
- `carlo delete`: Clean data

## Features

- `hdf5`: HDF5 checkpoint/measurement files (requires `libhdf5-dev`)
- `mpi`: MPI distributed backend (requires `libopenmpi-dev`)
- `strict-repro`: Use jump sequence for RNG (strict reproducibility mode)

## Key Patterns

When implementing a new Monte Carlo model:

1. Define struct holding configuration state
2. Implement `MonteCarlo` trait with `sweep()` method
3. Implement `FromParams` trait for construction from params
4. Optionally implement `MonteCarloCheckpoint` for HDF5 checkpointing (requires `hdf5` feature)
5. Call `ctx.measure("ObservableName", value)` in `sweep()` or `measure()` methods

<!-- code-review-graph MCP tools -->
## MCP Tools: code-review-graph

**IMPORTANT: This project has a knowledge graph. ALWAYS use the
code-review-graph MCP tools BEFORE using Grep/Glob/Read to explore
the codebase.** The graph is faster, cheaper (fewer tokens), and gives
you structural context (callers, dependents, test coverage) that file
scanning cannot.

### When to use graph tools FIRST

- **Exploring code**: `semantic_search_nodes` or `query_graph` instead of Grep
- **Understanding impact**: `get_impact_radius` instead of manually tracing imports
- **Code review**: `detect_changes` + `get_review_context` instead of reading entire files
- **Finding relationships**: `query_graph` with callers_of/callees_of/imports_of/tests_for
- **Architecture questions**: `get_architecture_overview` + `list_communities`

Fall back to Grep/Glob/Read **only** when the graph doesn't cover what you need.

### Key Tools

| Tool | Use when |
|------|----------|
| `detect_changes` | Reviewing code changes — gives risk-scored analysis |
| `get_review_context` | Need source snippets for review — token-efficient |
| `get_impact_radius` | Understanding blast radius of a change |
| `get_affected_flows` | Finding which execution paths are impacted |
| `query_graph` | Tracing callers, callees, imports, tests, dependencies |
| `semantic_search_nodes` | Finding functions/classes by name or keyword |
| `get_architecture_overview` | Understanding high-level codebase structure |
| `refactor_tool` | Planning renames, finding dead code |

### Workflow

1. The graph auto-updates on file changes (via hooks).
2. Use `detect_changes` for code review.
3. Use `get_affected_flows` to understand impact.
4. Use `query_graph` pattern="tests_for" to check coverage.
