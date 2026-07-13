# Stage 4.1 cross-cutting hardening report

This revision reviews and hardens Stage 4 only. Stage 5 worm functionality is intentionally absent.

## Stage 4 review fixes

- Registered `generalized_bench` with `harness = false`, so the existing Criterion target now builds and runs correctly.
- Added exact bin-count and accepted-macrostate audit checks to generalized-ensemble kernels.
- Preserved the last accepted macrostate bin in Wang-Landau snapshots with backward-compatible optional restore.
- Enabled serde_json `float_roundtrip`. Real testing exposed that JSON checkpoint traces otherwise changed by roughly one ULP after restore, violating the exact-future-trajectory contract.
- Removed a stale generalized audit helper and retained the canonical-kernel gate that prevents Wang-Landau or multicanonical kernels from using canonical beta exchange.

No structural error was found in the Stage 4 Wang-Landau lifecycle, DOS updates, frozen production path, canonical reweighting, or checkpoint versioning.

## Cross-cutting additions

- Unified `cache-audit` feature and common lattice, particle and macrostate validation policy.
- Criterion coverage for all requested throughput paths plus integrated autocorrelation time and ESS/s pilot reports.
- Domain-separated `RngStreamKey` covering task, run, chain, replica, logical thread, phase and substream identities.
- Schedule-independent Rayon task streams: physical worker assignment is deliberately not included in a task seed.
- Shared `accept_log_probability` helper and migration of Metropolis-Hastings decisions in Carlo.rs, CMC.rs, MCMC.rs and QMC.rs to log-domain comparisons.
- MCMC and Carlo checkpoint paths continue from serialized concrete RNG state.

## Validation environment

- rustc 1.90.0 (1159e78c4 2025-09-14)
- cargo 1.90.0 (840b83a10 2025-07-30)

## Successful validation

```bash
cargo fmt --all --check
cargo check --workspace --all-targets
cargo clippy --workspace --all-targets -- -D warnings
cargo test -p carlo-rs --lib --tests
cargo test -p cmc-rs --lib --tests
cargo test -p mcmc-rs --lib --tests
cargo test -p qmc-rs --lib --tests
cargo test -p carlo-rs --doc
cargo test -p cmc-rs --doc
cargo test -p mcmc-rs --doc
cargo test -p qmc-rs --doc
cargo check -p cmc-rs --all-targets --features cache-audit
cargo test -p cmc-rs --lib --tests --features cache-audit
cargo bench -p cmc-rs --no-run
cargo bench -p cmc-rs --bench performance_bench -- --test
```

Library and integration tests: **343 passed, 0 failed**.

- Carlo.rs: 134
- CMC.rs: 140
- MCMC.rs: 29
- QMC.rs: 40

Doctests: **3 passed, 0 failed, 7 ignored**.

The `cache-audit` build repeated all 140 CMC library/integration tests successfully. All four CMC benchmark executables were generated: library benches, `generalized_bench`, `particle_bench`, and `performance_bench`.

## Benchmark smoke-run statistics

These fixed-seed pilot values only validate the reporting path and are not portable performance claims:

```text
STAT_EFF ising_metropolis tau_int=11.750457 ess=43.573 ess_per_second=768.869
STAT_EFF ising_wolff tau_int=2.503790 ess=204.490 ess_per_second=4413.621
STAT_EFF ising_swendsen_wang tau_int=2.569426 ess=99.633 ess_per_second=2420.828
STAT_EFF lj_translation tau_int=10.261391 ess=6.237 ess_per_second=50.490
```

The smoke run also executed successfully for:

- Ising Metropolis attempted updates;
- Wolff cluster sites;
- Swendsen-Wang physical edges;
- Lennard-Jones trial translations;
- cell-list neighbor queries;
- batch delta-energy;
- Wang-Landau JSON checkpoint serialization.

## Optional native features

`cargo check --workspace --all-targets --all-features` was attempted. It stopped in the external `hdf5-sys 0.8.1` build script because the host advertises HDF5 `1.14.5`, which that build script rejects as an invalid version. Default features and the new `cache-audit` feature compile cleanly; no Rust source error was reported before the native HDF5 build failure. MPI/HDF5 dependency modernization is outside this Stage 4.1 scope.

See `CMC.rs/CACHE_AUDIT.md`, `CMC.rs/BENCHMARKS.md`, and `Carlo.rs/RNG_STREAMS.md` for operating contracts.
