# MCMC.rs v0.5 Validation Report

## v0.4 review result

The uploaded v0.4 workspace had no blocker requiring a repair release. NUTS
candidate selection, accepted-state atomicity, warmup freeze, checkpoint
continuation, E-BFMI and backward-compatible trace loading all passed review and
were suitable foundations for v0.5. The detailed review is in
`REVIEW_V04.md`.

The review identified three production-readiness improvements rather than
correctness blockers:

- retain summed trajectory momentum and apply merged-tree plus both
  cross-subtree generalized U-turn checks;
- find a reasonable initial HMC/NUTS step size before dual averaging, with an
  option to repeat after metric-window updates;
- expose a finite-difference gradient checker for differentiable targets.

## v0.5 implementation checks

### Generalized multinomial NUTS termination

- Every tree node retains the sum of all momenta in that subtree.
- A merged trajectory checks the full-tree criterion and both cross-subtree
  joins.
- Unit, diagonal and dense metrics provide allocation-free summed-momentum
  products in the NUTS hot path.
- The existing corrected multinomial candidate rule remains: a subtree that
  terminated internally contributes diagnostics but does not add candidate mass
  to the outer multinomial pool.
- Custom `Metric` implementations remain source-compatible through default
  trait methods.

### Dynamic initial step-size search

- `StepSizeSearch` performs a bounded doubling/halving search using one-step
  Hamiltonian proposals.
- Static HMC and NUTS can opt into the same search implementation.
- Search is restricted to warmup and rejects invalid bounds or iteration
  limits.
- By default it can repeat after a metric-window update, then restarts dual
  averaging from the new scale.
- Search completion state is serialized so checkpoint continuation is exact.
- New fields use serde defaults; v0.4 HMC/NUTS JSON loads with the search
  disabled, preserving old behavior.

### Gradient validation

- `check_gradient` compares the analytic gradient with central finite
  differences.
- It uses coordinate-scaled perturbations and independent absolute/relative
  tolerances.
- The report contains analytic, numerical, absolute and relative error for each
  component and identifies mismatches without panicking.
- Invalid configurations, non-finite evaluations and support-crossing probes
  return structured errors.

## Actual Rust validation

The reconstructed implementation workspace was validated with Rust 1.90.0:

```bash
cargo fmt --all --check
cargo check -p mcmc-rs --all-targets
cargo clippy -p mcmc-rs --all-targets -- -D warnings
cargo test -p mcmc-rs
cargo check --workspace --all-targets
```

Results:

- MCMC.rs: **61 passed, 0 failed**;
- strict Clippy: passed with warnings denied;
- default-feature workspace check for Carlo.rs, CMC.rs, QMC.rs and MCMC.rs:
  passed;
- `MCMC.rs/Cargo.toml` and the root `Cargo.lock` both contain version `0.5.0`.

The v0.5-specific suite covers:

- generalized summed-momentum and cross-subtree U-turn criteria;
- metric summed products against explicit velocity calculations;
- HMC and NUTS recovery from a deliberately poor initial scale;
- re-search after metric updates;
- exact checkpoint continuation with step-size search state;
- loading v0.4 HMC/NUTS JSON without the new fields;
- two-dimensional correlated-Gaussian NUTS mean, variance and covariance;
- correct-gradient acceptance, wrong-component reporting and invalid-input
  rejection.

All v0.1-v0.4 regression tests also passed, including state atomicity,
checkpoint round trips, replica exchange, transforms, HMC, NUTS, E-BFMI and
multi-chain diagnostics.

## Optional HDF5 feature

The command

```bash
cargo check -p mcmc-rs --all-targets --features hdf5
```

again stopped before MCMC.rs source compilation. `hdf5-sys 0.8.1` found the
system HDF5 installation but rejected its 1.14.5 header with
`Invalid H5_VERSION: "1.14.5"`. This is an external dependency/system-library
compatibility limitation; default-feature validation is unaffected.

## Final archive requirements

The final workspace archive must:

- contain the new step-size search, gradient checker, v0.5 tests, examples and
  review/validation documents;
- declare `mcmc-rs` version 0.5.0 in both `Cargo.toml` and `Cargo.lock`;
- exclude `target/`, `.git/`, Cargo registries, toolchains and authentication
  configuration;
- contain no internal registry credentials;
- pass ZIP CRC/integrity verification after creation.
