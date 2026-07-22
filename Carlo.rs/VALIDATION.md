# Carlo.rs — Physics Validation Task Tracker

> Created 2026-07-22. Branch: `dev`.
> Baseline: `../MATURITY_ASSESSMENT.md`.

## Current status

**Framework core: stable. HDF5/MPI: experimental.**

294 tests. MC loop, scheduler, context, measurements, accumulators, in-memory checkpoint, RNG streams are solid. HDF5 checkpoint I/O, MPI backend, and PT exchange protocol are the critical gaps.

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

### [ ] C-P1.1 — MPI PT exchange protocol
- **Problem:** MPI parallel-tempering exchange untested (2 ignored smoke tests only). Controller/worker protocol, chain permutation, measurement synchronization all unverified.
- **Plan:** `#[cfg(feature = "mpi")] #[ignore]` test: 4-replica ladder, verify exchange acceptance ratio, chain ordering after swaps, measurement gather correctness.
- **File:** `Carlo.rs/tests/mpi/pt_exchange.rs` (new)
- **Status:** not started

### [ ] C-P1.2 — MPI controller/worker scheduler
- **Problem:** `MpiBackend` controller/worker task partitioning, checkpoint two-phase commit (`*.next.h5` staging, `mpi-checkpoint.json` commit marker), restart validation — all untested.
- **Plan:** `#[cfg(feature = "mpi")] #[ignore]` test: 8 tasks across 2 ranks. Verify task partitioning (`task_id % size == rank`), results aggregation, checkpoint commit.
- **File:** `Carlo.rs/tests/mpi/distributed.rs` (extend existing)
- **Status:** not started

### [ ] C-P2.1 — strict-repro feature test
- **Problem:** The `strict-repro` feature (jump-sequence RNG for exact reproducibility across task counts) has zero tests.
- **Plan:** Run same model with `strict-repro` on, different task counts. Verify identical RNG streams.
- **File:** `Carlo.rs/tests/integration/reproducibility.rs` (extend)
- **Status:** not started

### [ ] C-P2.2 — HDF5 result merging
- **Problem:** `merge_results` / `merge_results_from_files` for HDF5 measurement files is untested. This is the primary production analysis path.
- **Plan:** Create 2 HDF5 result files, merge, verify rebinned estimates match manual calculation.
- **File:** `Carlo.rs/tests/io/merge_hdf5.rs` (new)
- **Feature gate:** `#[cfg(feature = "hdf5")]`
- **Status:** not started

### [ ] C-P2.3 — Decorrelated autocorrelation time reference test
- **Problem:** `compute_decorrelated_autocorr_time` validated only for non-negativity/finiteness. No closed-form reference.
- **Plan:** Feed AR(1) multivariate data with known covariance structure. Assert per-component τ matches analytic value.
- **File:** `Carlo.rs/tests/unit/merge.rs` (extend)
- **Status:** not started

### [ ] C-P2.4 — Thread-count independence
- **Problem:** No test verifies multi-thread (Rayon) produces statistically equivalent results to single-thread.
- **Plan:** Run 4-task job with 1 thread vs 4 threads. Assert mean estimates agree within 2σ.
- **File:** `Carlo.rs/tests/integration/backend.rs` (extend)
- **Status:** not started

## Completion log

| Date | Task | Result |
|------|------|--------|
| 2026-07-22 | C-P0.1 | ✅ HDF5 checkpoint: clocks now survive round-trip + 5 tests |
| 2026-07-22 | C-P0.2 | ✅ Estimate.autocorr_time fixed + 6 AR(1) reference tests |
