# MCMC.rs v0.2 Validation Report

## Review performed

- Confirmed the uploaded workspace contains the v0.1.1 target module, atomic
  component/slice fixes and state-invariant regression tests.
- Reviewed every MCMC source and test module after the target-typed
  `TransitionKernel<T>` refactor.
- Checked balanced Rust delimiters, comments and strings across all 190
  workspace Rust source files.
- Checked that every declared Rust module resolves, every struct has unique
  field names, the workspace includes MCMC.rs, and both the package manifest
  and root lockfile identify `mcmc-rs` as 0.2.0.
- Checked changed Rust files for trailing whitespace and lines longer than 100
  columns.
- Numerically checked dense covariance/Cholesky reconstruction, correlated
  proposal geometry, simplex normalization and finite-difference Jacobians,
  and the generic replica-exchange acceptance ratio independently.
- Audited v0.1 checkpoint compatibility for every newly serialized field.
- Added tests for Gibbs error atomicity, composition report aggregation,
  legacy/empty report normalization, dense covariance, correlated proposal
  generation, transform round trips, transformed-density correction,
  replica-exchange diagnostics and deterministic replay.
- Rechecked component-wise Metropolis and slice sampling input-state validation
  and consistent one-transition iteration semantics.
- Rechecked accepted and rejected replica exchanges so state/log-density caches
  remain synchronized and iteration counts advance deterministically.

## Important validation limitation

The execution container used for this revision has no `cargo`, `rustc`,
`rustfmt`, `clippy`, local Cargo registry or network access for installing crate
dependencies. Consequently this report does **not** claim an actual Rust build
or test pass. Run the following in a Rust-enabled checkout before merging:

```bash
cargo fmt --all --check
cargo clippy -p mcmc-rs --all-targets
cargo test -p mcmc-rs
cargo clippy -p mcmc-rs --all-targets --features hdf5
cargo test -p mcmc-rs --features hdf5
```

## Workspace note

MCMC.rs's HDF5 code has been updated for hdf5 0.8. Carlo.rs still contains its
pre-existing `create_dataset_simple` calls behind Carlo's own `hdf5` feature.
Therefore a full workspace `--all-features` build may still require a separate
Carlo.rs HDF5 migration; the default workspace and MCMC-only HDF5 feature do
not enable Carlo's HDF5 feature.
