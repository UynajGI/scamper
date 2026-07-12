# MCMC.rs v0.1.1 Validation Report

## v0.1.1 changes (2026-07-12)

- Flattened `target/` directory into single `MCMC.rs/src/target.rs`
- `ComponentWiseMetropolis`: `#[serde(default)] proposed_position` workspace, atomic swap on sweep completion, single iteration per transition, O(d) copy + O(1) swap, `log_acceptance` NaN guard
- `SliceSampler`: `#[serde(default)] working_position` workspace, atomic swap, single iteration per transition, `max_shrink_steps == 0` guard, bracket interval validity check
- `EuclideanState::validate()`: finite position/density, gradient cache dimension and content checks; called by `MemoryTrace::record()` and `ChainCheckpoint::validate_format()`
- `MemoryTrace::validate()` strengthened: expected draw count, discrete column ranges, `chain_id` consistency
- 6 new state-invariant tests: atomic error recovery (component + slice), unified iteration counting, invalid slice limits, old checkpoint backward compat

## Implemented checks

- Static delimiter/string/comment balance over every Rust source file.
- Workspace membership and local `Cargo.lock` package entry checked.
- JSON checkpoint fields audited to avoid serializing NaN/Infinity sentinels.
- Trace layout invariants checked by `MemoryTrace::validate`.
- State invariants checked by `EuclideanState::validate()` on every `MemoryTrace::record()` and checkpoint load.
- Independent numerical cross-check of the implemented rank-normalized R-hat
  and ESS formulas on synthetic data:
  - four IID normal chains: R-hat approximately 1.0001, bulk ESS 8000/8000;
  - one chain shifted by 1.5 standard deviations: R-hat approximately 1.215.

## Test results (2026-07-12)

```bash
cargo fmt --all --check     # PASS
cargo clippy -p mcmc-rs     # PASS (no issues)
cargo test -p mcmc-rs       # 13 passed, 1 failed (pre-existing)
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
| `fallible_component_transition_leaves_accepted_state_unchanged` | PASS (new) |
| `fallible_slice_transition_leaves_accepted_state_unchanged` | PASS (new) |
| `every_kernel_advances_iteration_once_per_transition` | PASS (new) |
| `invalid_slice_limits_fail_without_mutating_state` | PASS (new) |
| `legacy_component_checkpoint_without_workspace_remains_usable` | PASS (new) |
| `legacy_slice_checkpoint_without_workspace_remains_usable` | PASS (new) |
| `json_checkpoint_preserves_exact_future_trajectory` | **FAIL** — pre-existing; `serde_json` f64 round-trip introduces 1-ULP differences in `MemoryTrace` positions; test uses `assert_eq!` on floats |

## Known issues

- **`json_checkpoint_preserves_exact_future_trajectory`**: The test compares `MemoryTrace` (contains `f64` positions) through a `serde_json` round-trip via `assert_eq!`. Float values differ by 1 ULP after serialization/deserialization. Fix requires approximate comparison.
- **hdf5 feature broken workspace-wide**: hdf5 0.8.1 lacks `create_dataset_simple` (code expects 0.9.x API). Affects Carlo.rs and MCMC.rs when `--features hdf5` is enabled. `cargo clippy --all-features` and `cargo test --all-features` cannot be run until the workspace hdf5 dependency is bumped.
