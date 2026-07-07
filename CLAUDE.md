# CLAUDE.md

Scuttle — Monte Carlo framework with three Rust crates: Carlo.rs (core), QMC.rs (quantum), CMC.rs (classical).


Hooks via **lefthook** (`lefthook.yml` + `.lefthook/` scripts). pre-commit runs `fmt --check` → `cargo check` → `clippy` → `typos` on staged `.rs` (affected crates only); commit-msg enforces Conventional Commits; pre-push runs `cargo deny` only (test moved to CI for fast pushes — run `just test` manually before pushing). Install: `lefthook install` or `just hooks`. Skip: `LEFTHOOK=0 git commit`. Lint policy (incl. `unsafe_code = "deny"`) is codified in `[workspace.lints]` (`Cargo.toml`) — no `RUSTFLAGS` needed; the sole exception is `Carlo.rs/src/backend/mpi.rs` (`#![allow(unsafe_code)]` for MPI FFI). CI (`.github/workflows/ci.yml`) runs fmt + clippy + test + deny (`--all-features`) as parallel jobs.

## Workspace

| Crate | Role | Description |
|-------|------|-------------|
| Carlo.rs | Core framework | `MonteCarlo` trait, `Scheduler`, `Context`, `Measurements`, `Merge`, `Backend` |
| QMC.rs | Quantum MC | Worldline objects (continuous/discrete) — pure toolbox, no `MonteCarlo` impl |
| CMC.rs | Classical MC | Layered: `Lattice` → `System` → `Model` → `Algorithm` → `ClassicalMC` wrapper |

## Carlo.rs Architecture

MonteCarlo trait → Scheduler.run_one() → Results flow:

| Module              | File                    | Purpose                                                            |
| ------------------- | ----------------------- | ------------------------------------------------------------------ |
| `MonteCarlo` trait  | `monte_carlo.rs`        | Core: `sweep(ctx)`, `measure(ctx)`, `Rng` type                     |
| `FromParams` trait  | `monte_carlo.rs`        | Construct model from `Params` dict                                 |
| `Context`           | `context.rs`            | RNG, measurements, sweep counter, `ctx.measure(name, val)`         |
| `Run`               | `run.rs`                | Single run lifecycle, `step()`, checkpoint/restart                 |
| `Scheduler`         | `scheduler.rs`          | Thermalization → measurement loop, `run_one` / `run_parallel`      |
| `Backend`           | `backend/`              | `RayonBackend` (threads), `MpiBackend` (MPI)                       |
| `Measurements`      | `measurements.rs`       | Binned `Accumulator`, complex observables                          |
| `Merge`             | `merge.rs`              | Rebinning, autocorr time, `merge_results`, `merge_task_results`    |
| `Evaluable`         | `evaluable.rs`          | Jackknife resampling, `Evaluator`, `MultiplexEvaluator`            |
| `ResultTools`       | `output/resulttools.rs` | `dataframe()`, `measurement_from_obs()` — read-back `results.json` |
| `ParallelTempering` | `parallel_tempering.rs` | PT MC with chain scheduling                                        |
| `CLI`               | `cli.rs`                | `carlo run/status/merge/delete`                                    |
| `Job`               | `job/`                  | `JobInfo`, `TaskInfo`, `TaskMaker`, progress tracking              |

## CMC.rs Architecture

Orthogonal traits instead of a monolithic trait — each concern is a separate trait:

| Layer | File | Purpose |
|-------|------|---------|
| Lattice | `lattice.rs` | `CsrLattice` — flat CSR arrays (offsets + neighbors), `BondType`, builders (chain, square, hypercubic, triangular, honeycomb, kagome) |
| System | `system.rs` | `System { lattice, spins, energy, beta }` — pub fields, β moved here from model structs |
| Traits | `hamiltonian.rs` | `Hamiltonian`, `ClusterModel`, `Proposable`, `Measurable`, `HeatBathable`, `ContinuousHeatBathable` — orthogonal model traits |
| Models | `models.rs` | `IsingModel`, `PottsModel`, `XYModel`, `HeisenbergModel` — implement traits above |
| Algorithm | `algorithm.rs` | `Algorithm<H>` trait — `MetropolisCore`, `WolffCore`, `SWCore`, `HeatBathCore`, `MicrocanonicalCore`, `ContinuousHeatBathCore` |
| Proposal | `proposal.rs` | `ProposalStrategy<H>` — Standard, OPSS (adaptive over-relaxation) |
| Wrapper | `classical_mc.rs` | `ClassicalMC<H, A>` — `MonteCarlo`, `FromParams`, `ParallelTemperingCompatible`, JSON checkpoint |
| Multi-spin | `multi_spin.rs` | `MultiSpinIsing` — bit-parallel Ising with 64 replicas, impl `MonteCarlo` + `FromParams` + PT |
| Observables | `observables.rs` | Pluggable `Observable<H>` + `DefaultObservableSet` (Energy, Magnetization) |
| Postprocess | `postprocess.rs` | Derived observables: `susceptibility()`, `specific_heat()`, `binder_cumulant()`, `compute_correlation_1d()` |

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

| Tool                        | Use when                                               |
| --------------------------- | ------------------------------------------------------ |
| `detect_changes`            | Reviewing code changes — gives risk-scored analysis    |
| `get_review_context`        | Need source snippets for review — token-efficient      |
| `get_impact_radius`         | Understanding blast radius of a change                 |
| `get_affected_flows`        | Finding which execution paths are impacted             |
| `query_graph`               | Tracing callers, callees, imports, tests, dependencies |
| `semantic_search_nodes`     | Finding functions/classes by name or keyword           |
| `get_architecture_overview` | Understanding high-level codebase structure            |
| `refactor_tool`             | Planning renames, finding dead code                    |

### Workflow

1. The graph auto-updates on file changes (via hooks).
2. Use `detect_changes` for code review.
3. Use `get_affected_flows` to understand impact.
4. Use `query_graph` pattern="tests_for" to check coverage.
