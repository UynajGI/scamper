# MCMC.rs v0.3 Validation Report

## v0.2 review result

The uploaded workspace contained the v0.2 kernel composition, Gibbs update,
dense covariance, transform and replica-exchange implementation. No
architecture-level blocker was found. The HMC-readiness gaps documented in
`REVIEW_V02.md` were small enough to fix while proceeding directly to v0.3.

## v0.3 implementation checks

- Added static HMC with private trajectory workspace and atomic state commit.
- Added unit, diagonal and dense inverse-mass metrics.
- Added a reusable leapfrog integrator.
- Added dual averaging and expanding fast/slow/fast warmup windows with
  diagonal or dense metric adaptation.
- Added accepted-state gradient synchronization and cache reuse after
  rejection.
- Added invalid-trajectory and absolute energy-error divergence handling.
- Added analytic gradient pullback and log-Jacobian gradients for every v0.2
  built-in transform.
- Made `TransformedTarget` implement `DifferentiableLogDensity` when both the
  constrained target and bijector are differentiable.
- Added HMC scalar measurements to the Carlo.rs adapter.
- Excluded reconstructible HMC trajectory workspaces from serde checkpoints so
  non-finite divergent paths cannot poison JSON serialization.
- Enabled serde_json `float_roundtrip` for exact f64 checkpoint continuation.
- Added v0.3 tests for HMC moments, divergence atomicity, gradient-cache reuse,
  warmup checkpoint trajectory equivalence, post-divergence serialization,
  metrics, leapfrog reversibility and energy scaling, dual averaging, windowed
  metric adaptation, incomplete warmup rejection and transformed gradients.
- Updated package and lockfile `mcmc-rs` version to `0.3.0`.

## Actual Rust validation

Validation used Rust 1.90.0. The following commands completed successfully:

```bash
cargo fmt --all --check
cargo check -p mcmc-rs --all-targets
cargo clippy -p mcmc-rs --all-targets
cargo test -p mcmc-rs
```

The default-feature test run completed with **41 passed, 0 failed**. This
includes all retained v0.1/v0.2 tests, the new v0.3 tests and exact checkpoint
continuation tests.

## Structural and numerical validation

- Parsed all 211 Rust source files in the workspace with the tree-sitter Rust
  grammar: no syntax-error nodes.
- Resolved every `mod name;` declaration in MCMC.rs to a source file/module.
- Checked Rust struct declarations for duplicate field names: none found.
- Checked the workspace still includes `MCMC.rs` and Carlo.rs still exposes
  `Run::from_parts` / `Run::finalize_with_mc`.
- Dense Cholesky reconstruction maximum error: approximately `4.44e-16`.
- Dense metric velocity/kinetic-energy reference check passed.
- Empirical dense-momentum covariance maximum error: approximately `1.96e-3`.
- Leapfrog reversibility error: approximately `5.55e-17`.
- Simplex transformed-gradient finite-difference maximum error: approximately
  `2.58e-10`.
- Independent fixed-HMC standard-normal run produced mean approximately
  `0.00895` and variance approximately `1.0033`.

## Optional HDF5 feature

The command below did not reach MCMC.rs source compilation:

```bash
cargo check -p mcmc-rs --all-targets --features hdf5
```

`hdf5-sys 0.8.1` rejected the installed HDF5 `1.14.5` header because that
crate version's build-time parser recognizes the 1.8/1.10/1.12 version lines.
Therefore the default feature set is fully validated, while the optional HDF5
feature remains unverified in this environment due to dependency/system-library
compatibility.
