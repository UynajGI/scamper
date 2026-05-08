# CLAUDE.md

Scuttle — Monte Carlo framework with three Rust crates: Carlo.rs (core), QMC.rs (quantum), CMC.rs (classical). Carlo.jl and StochasticSeriesExpansion.jl are reference libraries, not part of this implementation.


Pre-commit hook (`.githooks/pre-commit`) runs `fmt --check` → `clippy -- -D warnings` → `test` on staged `.rs` files. Enable: `git config core.hooksPath .githooks`.

## Architecture

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
