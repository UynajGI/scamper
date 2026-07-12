# CLAUDE.md

Scuttle — Monte Carlo framework with three Rust crates: Carlo.rs (core), QMC.rs (quantum), CMC.rs (classical).


Hooks via **lefthook** (`lefthook.yml` + `.lefthook/` scripts). pre-commit runs `fmt --check` → `cargo check` → `clippy` → `typos` on staged `.rs` (affected crates only); commit-msg enforces Conventional Commits; post-commit runs `codegraph sync` (non-blocking); pre-push runs `cargo deny` only (test moved to CI for fast pushes — run `just test` manually before pushing). Install: `lefthook install` or `just hooks`. Skip: `LEFTHOOK=0 git commit`. Lint policy (incl. `unsafe_code = "deny"`) is codified in `[workspace.lints]` (`Cargo.toml`) — no `RUSTFLAGS` needed; the sole exception is `Carlo.rs/src/backend/mpi.rs` (`#![allow(unsafe_code)]` for MPI FFI). CI (`.github/workflows/ci.yml`) runs fmt + clippy + test + deny (`--all-features`) as parallel jobs.

## Workspace

| Crate | Role | Description |
|-------|------|-------------|
| Carlo.rs | Core framework | `MonteCarlo` trait, `Scheduler`, `Context`, `Measurements`, `Merge`, `Backend` |
| QMC.rs | Quantum MC | General continuous-time lattice QMC (`LatticeSpinQmc` implements `MonteCarlo` + `FromParams`), spin-boson wormhole QMC
| CMC.rs | Classical MC | Layered: `Lattice` → `System` → `Model` → `Algorithm` → `ClassicalMC` wrapper |

## Carlo.rs Architecture

MonteCarlo trait → Scheduler.run_one() → Results flow:

| Module              | File                    | Purpose                                                            |
| ------------------- | ----------------------- | ------------------------------------------------------------------ |
| `MonteCarlo` trait  | `monte_carlo.rs`        | Core: `sweep(ctx)`, `measure(ctx)`, `Rng` type, lifecycle hooks    |
| `FromParams` trait  | `monte_carlo.rs`        | Construct model from `Params` dict                                 |
| `Context`           | `context.rs`            | RNG, measurements, sweep counter, `RunPhase`, checkpoint state     |
| `Run`               | `run.rs`                | Single run lifecycle, `step()`, checkpoint/restart                 |
| `Scheduler`         | `scheduler.rs`          | Thermalization → measurement loop, `run_one`/`run_parallel`/`run_controlled` |
| `Backend`           | `backend/`              | `RayonBackend` (threads), `MpiBackend` (MPI)                       |
| `Measurements`      | `measurements.rs`       | Binned `Accumulator`, complex observables                          |
| `Merge`             | `merge.rs`              | Rebinning, autocorr time, `merge_results`, `merge_task_results`    |
| `Evaluable`         | `evaluable.rs`          | Jackknife resampling, `Evaluator`, `MultiplexEvaluator`            |
| `ResultTools`       | `output/resulttools.rs` | `dataframe()`, `measurement_from_obs()` — read-back `results.json` |
| `ParallelTempering` | `parallel_tempering.rs` | PT MC with chain scheduling                                        |
| `RunPhase`          | `phase.rs`              | Explicit lifecycle: Initialization→Thermalization→Measurement→Finished |
| `AdaptiveRunControl`| `run_control.rs`        | Algorithm-driven phase transitions, `RunDecision`                  |
| `CLI`               | `cli.rs`                | `carlo run/status/merge/delete`                                    |
| `Job`               | `job/`                  | `JobInfo`, `TaskInfo`, `TaskMaker`, progress tracking              |

## CMC.rs Architecture

Orthogonal traits instead of a monolithic trait — each concern is a separate trait.
Directories group related modules; public API is re-exported flat from `lib.rs`.

| Directory | Files | Purpose |
|-----------|-------|---------|
| `core/` | `move.rs`, `cache.rs`, `trial.rs`, `ensemble.rs`, `acceptance.rs`, `visit.rs` | Move types, incremental patches, `TrialEvaluator`, `MetropolisHastingsAcceptance`, visit schedules |
| `lattice/` | `graph.rs`, `state.rs`, `interaction.rs`, `models.rs`, `proposal.rs` | `CsrLattice` + builders, `System`, `Hamiltonian` + capability traits, built-in models, `ProposalStrategy` |
| `algorithms/` | `metropolis.rs`, `wolff.rs`, `swendsen_wang.rs`, `heat_bath.rs`, `microcanonical.rs`, `hybrid.rs`, `common.rs` | `Algorithm<H>` trait + 6 kernels, `SimulationPhase`, `checked_probability` |
| `observables/` | `energy.rs`, `magnetization.rs`, `correlation.rs`, `common.rs` | `Observable<H>`, `DefaultObservableSet`, `TotalEnergy`, `Magnetization`, `compute_correlation_1d` |
| Top-level | `classical_mc.rs`, `multi_spin.rs`, `postprocess.rs` | `ClassicalMC` Carlo.rs adapter, `MultiSpinIsing`, derived observables |

Key patterns:
- `ClassicalMC<IsingModel, MetropolisCore>` → `Scheduler.run_one()` → `Results`
- Models are stateless — temperature (β) lives in `System`, not in model structs
- Algorithm trait: `sweep(&mut self, system, model, rng)` — directly mutates system.energy
- CSR lattice: cache-friendly neighbor iteration via `lattice.neighbors(site)` returning `&[usize]`
- `SmallVec<[f64; 3]>` for spin proposals (stack-allocated, no heap for spin_dim ≤ 3)
- Heat-bath (Glauber): discrete (Ising/Potts → `HeatBathable`) + continuous (XY/Heisenberg → `ContinuousHeatBathable`, vMF/Best-Fisher sampling)
- Over-relaxation: `MicrocanonicalCore` reflects spins across local field (ΔE=0), no acceptance
- Derived observables measured post-run from E²/M²/M⁴ moments stored in `Results`
- `lattice_type` param: `"chain"`, `"square"`, `"triangular"`, `"honeycomb"`, `"kagome"`
- Users can ignore `ClassicalMC` and compose manually for custom behavior

## QMC.rs Architecture

Boundary: Carlo.rs owns runtime (scheduling, RNG, measurements, I/O); QMC.rs owns physics (representations, kernels, estimators).

Two production backends:

### Lattice QMC (primary)

| Module | File | Purpose |
|--------|------|---------|
| `graph` | `graph.rs` | `CsrGraph` — typed/weighted CSR adjacency + unique edge table; builders (chain, square, hypercubic, adjacency, CSR) |
| `local_space` | `local_space.rs` | `LocalHilbertSpace` trait, `SpinSpace` — site-resolved arbitrary `S` with exact integer `2S` algebra |
| `lattice::model` | `lattice/model.rs` | `SpinModelBuilder` → `SpinLatticeModel` — sparse positive `OperatorTerm` catalog, Marshall `Z2` gauge, Heisenberg/XY/XXZ/XYZ/tfim/generic |
| `lattice::vertex` | `lattice/vertex.rs` | `VertexKind` (positive local matrix elements), `Vertex` (sampled insertion), `Event` (worldline endpoint) |
| `lattice::configuration` | `lattice/configuration.rs` | `LatticeConfiguration` — product state + unsorted vertex vector; `WorldlineIndex` — time-ordered leg links |
| `lattice::scattering` | `lattice/scattering.rs` | `ScatteringTable` — exact local-detailed-balance directed-loop routing (LowBounce + Metropolis policies) |
| `lattice::updates` | `lattice/updates.rs` | `ContinuousLatticeEngine<M>` — diagonal add/remove + directed-loop blocks; journal-based rollback |
| `lattice::observables` | `lattice/observables.rs` | Magnetization, staggered magnetization, susceptibility, energy, vertex orders, edge SzSz correlation |
| `lattice::mc` | `lattice/mc.rs` | `LatticeSpinQmc` — Carlo.rs adapter (`MonteCarlo` + `FromParams`); warmup schedule adaptation |

### Spin-Boson QMC

| Module | File | Purpose |
|--------|------|---------|
| `spin_boson` | `spin_boson/` | Continuous-time retarded-interaction wormhole QMC (spin-boson impurity, bath samplers, scattering table, diagonal/loop updates, observables) |

### Shared

| Module | File | Purpose |
|--------|------|---------|
| `algorithm` | `algorithm.rs` | `QmcKernel<C,R>` trait, `UpdateSchedule` — reusable kernel contract |

Key patterns:
- `LatticeSpinQmc` wraps `ContinuousLatticeEngine<SpinLatticeModel>` and implements `MonteCarlo` + `FromParams` — drop-in for `Scheduler.run_one()`
- Models are compiled from physical couplings into positive `K=C-H` operator catalogs; the engine has no model-name branches
- `CsrGraph` is data, not an algorithm — `from_csr` mirrors CMC's layout without crate coupling
- Marshall gauge: BFS solves bipartite `Z2` phases; rejects frustrated/non-stoquastic models
- `rand` 0.10: use `RngExt` for `random()`/`random_range()`, `Rng` for trait bounds (`R: Rng + ?Sized`)
- `UpdateSchedule` (diagonal proposals, directed loops, max loop steps) fixed during measurement; adaptation allowed during thermalization
- `from_sparse()` on `OperatorTerm` is the extension point for future bosonic/fermionic catalogs

## Features

- `hdf5` — checkpoint/measurement files (needs `libhdf5-dev`)
- `mpi` — distributed backend (needs `libopenmpi-dev`)

## Key Patterns

1. Define struct with config state
2. `impl MonteCarlo` — `type Rng`, `fn sweep()`, optional `fn measure()`
3. `impl FromParams` — construct from `Params`
4. Call `ctx.measure("Name", value)` in `sweep()` or `measure()`
5. Optional: `impl MonteCarloCheckpoint` for HDF5 save/load


## Behavioral Guidelines

**Tradeoff:** These guidelines bias toward caution over speed. For trivial tasks, use judgment.

## 1. Think Before Coding

**Don't assume. Don't hide confusion. Surface tradeoffs.**

Before implementing:
- State your assumptions explicitly. If uncertain, ask.
- If multiple interpretations exist, present them - don't pick silently.
- If a simpler approach exists, say so. Push back when warranted.
- If something is unclear, stop. Name what's confusing. Ask.

## 2. Simplicity First

**Minimum code that solves the problem. Nothing speculative.**

- No features beyond what was asked.
- No abstractions for single-use code.
- No "flexibility" or "configurability" that wasn't requested.
- No error handling for impossible scenarios.
- If you write 200 lines and it could be 50, rewrite it.

Ask yourself: "Would a senior engineer say this is overcomplicated?" If yes, simplify.

## 3. Surgical Changes

**Touch only what you must. Clean up only your own mess.**

When editing existing code:
- Don't "improve" adjacent code, comments, or formatting.
- Don't refactor things that aren't broken.
- Match existing style, even if you'd do it differently.
- If you notice unrelated dead code, mention it - don't delete it.

When your changes create orphans:
- Remove imports/variables/functions that YOUR changes made unused.
- Don't remove pre-existing dead code unless asked.

The test: Every changed line should trace directly to the user's request.

## 4. Goal-Driven Execution

**Define success criteria. Loop until verified.**

Transform tasks into verifiable goals:
- "Add validation" → "Write tests for invalid inputs, then make them pass"
- "Fix the bug" → "Write a test that reproduces it, then make it pass"
- "Refactor X" → "Ensure tests pass before and after"

For multi-step tasks, state a brief plan:
```
1. [Step] → verify: [check]
2. [Step] → verify: [check]
3. [Step] → verify: [check]
```

Strong success criteria let you loop independently. Weak criteria ("make it work") require constant clarification.

---

**These guidelines are working if:** fewer unnecessary changes in diffs, fewer rewrites due to overcomplication, and clarifying questions come before implementation rather than after mistakes.

## CodeGraph

This project has a CodeGraph index (`.codegraph/`). ALWAYS use
`codegraph explore` (CLI) or `codegraph_explore` (MCP) BEFORE Grep/Glob/Read
to explore the codebase. One call returns verbatim source + call paths +
blast-radius summary — faster and cheaper than reading files yourself.

Syncs automatically on post-commit via lefthook.
