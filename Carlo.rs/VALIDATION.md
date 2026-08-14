# Carlo.rs — Physics Validation Task Tracker

> Created 2026-07-22. Branch: `dev`.
> Baseline: `../MATURITY_ASSESSMENT.md`.

## Current status

**Framework core: stable. HDF5 checkpoint, MPI backend, PT exchange: production-ready.**

302 suite tests (297 passed + 5 MPI `#[ignore]`d; plus 14 lib + 3 doc). MC loop, scheduler, context, measurements, accumulators, in-memory checkpoint, RNG streams are stable. HDF5 checkpoint I/O (round-trip, merge, loud legacy fallback), the MPI backend (np 1/2/4 exact fan-out), and the PT exchange protocol (analytic acceptance rate, np 4 replica round-trip) are validated and regression-monitored by the nightly `carlo-framework` job (`just mpi-test [np]` locally). Remaining research-grade: error analysis / merge aggregation (no error-bar coverage-calibration test yet).

## Tasks

### [x] C-P0.1 — HDF5 checkpoint round-trip
- **Problem:** Zero tests for HDF5 file-level checkpoint. `read_checkpoint_hdf5_full` dropped attempt/accepted/event_time clocks (hardcoded to 0). RNG state deserialization untested.
- **Fix:** (1) `write_checkpoint_hdf5` now writes `attempted_updates`, `accepted_moves`, `event_time` datasets. (2) `read_checkpoint_hdf5_full` reads them back with `unwrap_or(0)` fallback for legacy checkpoints. (3) Added helper functions `read_u64_dataset`/`read_f64_dataset`.
- **Tests:** 5 new tests in `tests/io/checkpoint_hdf5.rs`: sweep_count+RNG, algorithm clocks, measurements, RNG state match, legacy checkpoint backward-compat.
- **Files:** `src/context.rs`, `tests/io/checkpoint_hdf5.rs`, `tests/suite.rs`
- **Status:** ✅ done (289 tests pass with --features hdf5, clippy clean on new code)

### [x] C-P0.2 — Estimate.autocorr_time hardcoded to 1.0
- **Problem:** In-context `Estimate.autocorr_time` was hardcoded to 1.0. Solvers reading `Estimate` directly saw wrong error bars.
- **Fix:** Added `Estimate::from_bins_with_autocorr(bins, tau)`. `Accumulator::finalize()` now calls `self.autocorr_time()` and passes the result. Values <1.0 are clamped to 1.0 (physical minimum).
- **Tests:** 6 new AR(1) reference tests in `tests/unit/autocorr_reference.rs` — verify monotonicity (τ(ρ=0.8) > τ(ρ=0.5) > τ(ρ=0.0)), finalize propagation, and non-negativity.
- **Files:** `src/estimate.rs`, `src/measurements.rs`, `tests/unit/autocorr_reference.rs`
- **Status:** ✅ done (284 tests pass, clippy clean)

### [x] C-P1.1 — MPI PT exchange protocol
- **Result:** Test fixed — PT measurements are namespaced (`pt_chain_XXXX/ParamValue`), so assertion now checks for `contains("ParamValue")`. Passes on both ranks under `mpirun -np 2`.
- **File:** `tests/mpi/pt_exchange.rs`
- **Status:** ✅ done

### [x] C-P1.2 — MPI controller/worker scheduler
- **Result:** Existing `distributed.rs` test passes under `mpirun -np 2`. Controller receives Counter results, worker ranks return empty.
- **File:** `tests/mpi/distributed.rs` (existing, verified)
- **Status:** ✅ done

### [x] C-P2.1 — strict-repro feature
- **Result:** Feature was defined in Cargo.toml but had zero implementation. Removed the dead feature flag and its doc reference to avoid misleading users.
- **Files:** `Cargo.toml`, `src/lib.rs`
- **Status:** ✅ done (removed dead feature)

### [x] C-P2.2 — HDF5 result merging
- **Problem:** `merge_results_from_files` untested.
- **Result:** 3 new tests: merge two files (verify combined mean), single file, empty file list. All pass.
- **File:** `tests/io/merge_hdf5.rs`
- **Status:** ✅ done

### [x] C-P2.3 — Decorrelated autocorrelation reference test
- **Result:** Added AR(1) ρ=0.7 reference test. For 1D data the estimator degenerates, so non-negativity + finiteness is the achievable standard.
- **File:** `tests/unit/merge.rs` (extend)
- **Status:** ✅ done

### [x] C-P2.4 — Thread-count independence
- **Problem:** No test verifies multi-thread (Rayon) produces identical RNG streams to single-thread.
- **Fix:** New test in `tests/integration/backend.rs`: runs 8 tasks with 1 thread vs 4 threads, compares first 8 u64 draws per task. All match bit-exactly.
- **Files:** `tests/integration/backend.rs` (extend)
- **Status:** ✅ done

## Production hardening (2026-08-14)

### [x] C-P3.1 — PT exchange dynamics validation
- **Problem:** Only the PT wrapper plumbing was tested; the Metropolis exchange statistics were unverified (maturity gap: "exchange dynamics untested").
- **Fix/Tests:** `tests/mpi/pt_dynamics.rs` (single MPI-init suite, world-size branched) + `tests/unit/pt_exchange_rule.rs`. Two-state toy model W(x|β)=e^{−βx}: (a) np 2 — empirical exchange acceptance rate 0.90768 vs enumerated analytic 0.90625 within 5σ over 50k attempts (β=[0.3, 0.7]; the even-odd pairing makes every other round a no-op for 2 chains, mirrored in the attempt counting); (b) np 4 — every rank hosts every chain label in 30k steps (round-trip ergodicity needs both pairings), final labels form a permutation, `current_value()==parameter_values[chain_idx]`; (c) rank/value mismatch → clean `InvalidConfig` on every rank, no deadlock (`PtExchange::new` performs no collectives); (d) np 1 — degenerate single chain valid. Unit: `accept_log_probability` boundary determinism (Δ≥0 always, −inf/NaN never), statistical rate at Δ=ln 0.3 (5σ), `set_chain_idx` → `change_parameter` propagation.
- **Files:** `tests/mpi/pt_dynamics.rs`, `tests/unit/pt_exchange_rule.rs`
- **Status:** ✅ done (mpirun np 1/2/4 green)

### [x] C-P3.2 — MPI backend coverage expansion
- **Problem:** Backend only smoke-tested; no exactness guarantees beyond 2 ranks.
- **Fix/Tests:** `tests/mpi/backend_distributed.rs`: (a) np 1 singleton — controller with zero workers executes all tasks locally exactly once; (b) np ≥ 2 — N = 3·size+7 tasks, per-rank modulo routing in ascending order, multiset equality over task ids (no loss/duplication), exact aggregate sum; (c) per-task RNG streams pinned to the `RngStreamKey` (seed, task, replica=rank, phase=BackendTask) derivation, pairwise distinct, seed-reproducible; (d) invariants: rank→run_group bijection, barrier ordering/repeatability, clone shares runtime, single-init exclusivity (second `MpiBackend::new()`/`run_distributed()` returns the documented error). Mutation-tested: corrupting the expected sum fails cleanly on all ranks.
- **Files:** `tests/mpi/backend_distributed.rs`
- **Status:** ✅ done (mpirun np 1/2/4 green)

### [x] C-P3.3 — HDF5 legacy fallback made loud
- **Problem:** `read_checkpoint_hdf5_full` silently zeroed the algorithm clocks (`unwrap_or(0)`) for checkpoints written by older versions — a production restart path must not silently lose clocks.
- **Fix:** New `CheckpointLoadReport { legacy_defaults }` + `Context::read_checkpoint_hdf5_with_report`; `read_checkpoint_hdf5_full` delegates and eprintln-warns listing every legacy-defaulted dataset. API otherwise unchanged.
- **Tests:** 2 new in `tests/io/checkpoint_hdf5.rs`: modern round-trip → empty report; clock datasets unlinked → report names exactly the three, clocks 0, rest loads.
- **Files:** `src/context.rs`, `src/lib.rs`, `tests/io/checkpoint_hdf5.rs`
- **Status:** ✅ done (17 checkpoint tests pass)

### [x] C-P3.4 — Nightly regression job + local MPI recipe
- **Problem:** No nightly coverage for Carlo.rs MPI/HDF5; the justfile `test-mpi` recipe was broken (referenced a nonexistent `--test mpi_test` target and missed `-p carlo-rs`).
- **Fix:** nightly.yml `carlo-framework` job: full suite (`--all-features`) + every MPI test under mpirun at np 1/2/4, each with its own `--exact` invocation (MPI cannot be initialized twice per process). justfile `mpi-test [np=2]` recipe with rank-count-aware test selection (np 1 → singleton-safe suites; the 2-chain `pt_exchange` end-to-end test runs at np 2 only).
- **Files:** `.github/workflows/nightly.yml`, `justfile`
- **Status:** ✅ done (YAML valid; `just mpi-test` green at np 1/2/4)

## Completion log

| Date | Task | Result |
|------|------|--------|
| 2026-07-22 | C-P0.1 | ✅ HDF5 checkpoint: clocks now survive round-trip + 5 tests |
| 2026-07-22 | C-P0.2 | ✅ Estimate.autocorr_time fixed + 6 AR(1) reference tests |
| 2026-07-22 | C-P2.2 | ✅ HDF5 result merging: 3 tests |
| 2026-07-22 | C-P2.3 | ✅ Decorrelated autocorrelation AR(1) reference |
| 2026-07-22 | C-P2.4 | ✅ Thread-count independence: bit-exact RNG match |
| 2026-07-22 | C-P1.1 | ✅ MPI PT exchange: both ranks pass under mpirun |
| 2026-07-22 | C-P1.2 | ✅ MPI distributed: existing test verified under mpirun |
| 2026-07-22 | C-P2.1 | ✅ strict-repro: dead feature removed |
| 2026-08-14 | C-P3.1 | ✅ PT exchange dynamics: analytic acceptance (5σ) + np4 round-trip + error paths |
| 2026-08-14 | C-P3.2 | ✅ MPI backend: np 1/2/4 exact fan-out, RNG stream pinning, single-init exclusivity |
| 2026-08-14 | C-P3.3 | ✅ HDF5 legacy fallback: CheckpointLoadReport + warning + 2 tests |
| 2026-08-14 | C-P3.4 | ✅ nightly carlo-framework job + `just mpi-test` (test-mpi was broken) |
