# CLAUDE.md

Scuttle — Monte Carlo framework. Carlo.rs (Rust, primary) ports Carlo.jl (Julia, reference).

## Commands

```bash
just check          # fmt + clippy + test
just test           # cargo test --workspace
just build          # cargo build --release
just bench          # cargo bench
just install-deps   # apt-get libhdf5-dev openmpi-bin libopenmpi-dev
```

Pre-commit hook (`.githooks/pre-commit`) runs `fmt --check` → `clippy -- -D warnings` → `test` on staged `.rs` files. Enable: `git config core.hooksPath .githooks`.

## Architecture

| Module | File | Purpose |
|--------|------|---------|
| `MonteCarlo` trait | `monte_carlo.rs` | Core: `sweep(ctx)`, `measure(ctx)`, `Rng` type |
| `FromParams` trait | `monte_carlo.rs` | Construct model from `Params` dict |
| `Context` | `context.rs` | RNG, measurements, sweep counter, `ctx.measure(name, val)` |
| `Run` | `run.rs` | Single run lifecycle, `step()`, checkpoint/restart |
| `Scheduler` | `scheduler.rs` | Thermalization → measurement loop, `run_one` / `run_parallel` |
| `Backend` | `backend/` | `RayonBackend` (threads), `MpiBackend` (MPI) |
| `Measurements` | `measurements.rs` | Binned `Accumulator`, complex observables |
| `Merge` | `merge.rs` | Rebinning, autocorr time, `merge_results`, `merge_task_results` |
| `Evaluable` | `evaluable.rs` | Jackknife resampling, `Evaluator`, `MultiplexEvaluator` |
| `ResultTools` | `output/resulttools.rs` | `dataframe()`, `measurement_from_obs()` — read-back `results.json` |
| `ParallelTempering` | `parallel_tempering.rs` | PT MC with chain scheduling |
| `CLI` | `cli.rs` | `carlo run/status/merge/delete` |
| `Job` | `job/` | `JobInfo`, `TaskInfo`, `TaskMaker`, progress tracking |

## Features

- `hdf5` — checkpoint/measurement files (needs `libhdf5-dev`)
- `mpi` — distributed backend (needs `libopenmpi-dev`)

## Key Patterns

1. Define struct with config state
2. `impl MonteCarlo` — `type Rng`, `fn sweep()`, optional `fn measure()`
3. `impl FromParams` — construct from `Params`
4. Call `ctx.measure("Name", value)` in `sweep()` or `measure()`
5. Optional: `impl MonteCarloCheckpoint` for HDF5 save/load

## MCP Tools

Use `code-review-graph` MCP **before** Grep/Glob/Read: `semantic_search_nodes`, `query_graph` (callers_of/callees_of/tests_for), `get_impact_radius`, `detect_changes`. The graph auto-updates on file changes.
