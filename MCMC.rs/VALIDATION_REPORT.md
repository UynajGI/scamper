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

## Test results (2026-07-12)

```bash
cargo fmt --all --check     # PASS
cargo clippy -p mcmc-rs     # PASS (no issues)
cargo test -p mcmc-rs       # 7 passed, 1 failed
```

| Test | Status |
|------|--------|
| `adaptive_random_walk_recovers_standard_normal_moments` | PASS |
| `proposal_scale_is_constant_in_sampling_phase` | PASS |
| `component_wise_recovers_univariate_std_normal` | PASS |
| `slice_sampler_recovers_univariate_standard_normal` | PASS |
| `multi_chain_deterministic_reproducibility` | PASS |
| `diagnostics_on_iid_and_shifted_chains` | PASS |
| `carlo_run_returns_sampler_trace` | PASS |
| `json_checkpoint_preserves_exact_future_trajectory` | **FAIL** — pre-existing; `serde_json` f64 round-trip introduces 1-ULP differences in `MemoryTrace` positions; test uses `assert_eq!` on floats |

## Known issues

- **`json_checkpoint_preserves_exact_future_trajectory`**: The test compares `MemoryTrace` (contains `f64` positions) through a `serde_json` round-trip via `assert_eq!`. Float values differ by 1 ULP after serialization/deserialization. Fix requires approximate comparison.
- **hdf5 feature broken workspace-wide**: hdf5 0.8.1 lacks `create_dataset_simple` (code expects 0.9.x API). Affects Carlo.rs and MCMC.rs when `--features hdf5` is enabled. `cargo clippy --all-features` and `cargo test --all-features` cannot be run until the workspace hdf5 dependency is bumped.
