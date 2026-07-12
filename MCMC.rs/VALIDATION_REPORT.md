# MCMC.rs v0.1 Validation Report

## Implemented checks

- Static delimiter/string/comment balance over every Rust source file.
- Workspace membership and local `Cargo.lock` package entry checked.
- JSON checkpoint fields audited to avoid serializing NaN/Infinity sentinels.
- Trace layout invariants checked by `MemoryTrace::validate`.
- Independent numerical cross-check of the implemented rank-normalized R-hat
  and ESS formulas on synthetic data:
  - four IID normal chains: R-hat approximately 1.0001, bulk ESS 8000/8000;
  - one chain shifted by 1.5 standard deviations: R-hat approximately 1.215.
- Tests included for Gaussian moments, adaptation freeze, slice sampling,
  multi-chain reproducibility, diagnostics, JSON checkpoint continuation, and
  Carlo.rs trace recovery.

## Environment limitation

The artifact-generation container does not contain `cargo`, `rustc`, or
`rustfmt`, and outbound DNS is unavailable, so the Rust toolchain could not be
installed. Consequently, these commands could not be executed here:

```bash
cargo fmt --all --check
cargo clippy -p mcmc-rs --all-targets --all-features
cargo test -p mcmc-rs --all-features
```

Run the commands above in a Rust environment before merging. This report does
not claim successful compilation.
