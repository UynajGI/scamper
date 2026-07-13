# MCMC.rs v0.4 Validation Report

## v0.3 review result

The uploaded v0.3 implementation had no blocker requiring a separate repair
release. Static HMC, metric geometry, windowed warmup, differentiable
transforms, state atomicity and checkpoint boundaries were suitable foundations
for NUTS. The detailed review is in `REVIEW_V03.md`.

## v0.4 implementation checks

- Added multinomial `Nuts<M>` for unit, diagonal and dense metrics.
- Added signed leapfrog integration for bidirectional trajectories.
- Added a metric-aware, allocation-free built-in U-turn hot path.
- Added maximum tree-depth and energy-error stopping.
- Kept all trajectory construction in private workspace and committed state
  atomically.
- Reused dual averaging and windowed diagonal/dense metric adaptation.
- Added transition acceptance statistic, Hamiltonian energy and depth-limit
  reporting.
- Extended `MemoryTrace`, Carlo measurements and multi-chain diagnostics.
- Added per-chain E-BFMI.
- Added backward-compatible serde defaults for new report, trace and diagnostic
  fields.
- Verified NUTS warmup checkpoint continuation step by step.
- Verified that divergent workspaces remain JSON-serializable because they are
  reconstructible and skipped.

During statistical testing, an early implementation incorrectly merged a
candidate from a subtree that had already terminated internally on a U-turn.
That produced an inflated standard-normal variance near 1.24. The subtree merge
rule was corrected so terminated subtree diagnostics remain counted but its
candidate mass is excluded from the outer multinomial pool. A 30,000-draw
independent numerical check then produced variance approximately 1.002.

## Actual Rust validation

The implementation workspace was validated with Rust 1.90.0:

```bash
cargo fmt --all --check
cargo check -p mcmc-rs --all-targets
cargo clippy -p mcmc-rs --all-targets -- -D warnings
cargo test -p mcmc-rs
cargo check --workspace --all-targets
```

Results:

- MCMC.rs: **51 passed, 0 failed**;
- strict Clippy: passed;
- default-feature workspace check for Carlo.rs, CMC.rs, QMC.rs and MCMC.rs:
  passed.

The suite covers all v0.1–v0.3 regressions plus:

- standard-normal NUTS moments;
- depth-limit reporting;
- divergence atomicity and post-divergence serialization;
- warmup checkpoint trajectory equivalence;
- metric displacement–velocity algebra;
- E-BFMI edge cases;
- dynamic-HMC trace columns and old JSON trace extension.

## Optional HDF5 feature

```bash
cargo check -p mcmc-rs --all-targets --features hdf5
```

did not reach MCMC.rs source compilation in the validation environment.
`hdf5-sys 0.8.1` rejected the installed HDF5 1.14.5 header with
`Invalid H5_VERSION: "1.14.5"`. Default-feature validation is unaffected.

## Final archive checks

The final workspace archive is required to:

- contain `MCMC.rs/src/kernel/nuts.rs`, v0.4 tests and documentation;
- declare `mcmc-rs` version 0.4.0 in both `Cargo.toml` and `Cargo.lock`;
- exclude build caches, registries and authentication configuration;
- pass ZIP CRC/integrity verification after creation.
