# Validation report

Date: 2026-07-12

## Phase 0 completion status

The sampling-core-v2 stabilization (Phase 0) is complete. All local verification passes:

```bash
cargo fmt --all --check     # PASS
cargo check --workspace      # PASS
cargo clippy --workspace     # PASS (no warnings)
cargo test --workspace       # 245 passed, 7 ignored
```

Note: `--all-features` cannot run locally because `mpi` requires `libopenmpi-dev` and `hdf5` has a crate API mismatch — both are pre-existing issues tracked separately, not introduced by sampling-core-v2.

## Bug fixes applied

- Snapshot format tag validated on load (rejects non-`cmc-rs-snapshot-v2`)
- `BondType::as_label()` / `from_label()` stable labels replace `Debug`-based serialization
- `CarloError::CheckpointCorrupted` used for snapshot errors; `InvalidConfig` retained for parameter errors
- Unused `energy` field removed from save_snapshot (always recomputed on load)
- Duplicate test removed from `Carlo.rs/tests/checkpoint_test.rs`

## New test coverage

### Statistical correctness (7 tests, `statistical_correctness_test.rs`)
- Exact Ising N=2,3,4 energy vs Boltzmann enumeration (within 3σ)
- Potts q=3, N=4 exact energy mean (within 3σ)
- Algorithm consistency: Metropolis, Wolff, SW at 8×8 (pairwise 3σ)
- PT energy consistency under `change_parameter` (physical energy invariant)
- Fixed seed reproducibility (bitwise-identical results)

### Detailed balance (6 tests, `detailed_balance_test.rs`)
- Asymmetric Hastings proposal (custom `ProposalStrategy`, 80% bias)
- Batch move (all-spin flip with p=0.3 on N=3)
- Parallel edges (two bonds between same sites)
- Self-loop + normal bond
- Heat bath conditional sampler
- Wolff cluster rejection-free transitions

### Checkpoint persistence (6 tests, `checkpoint_test.rs`)
- Split-run state identity (200→save→200 vs 400 continuous)
- 1000-sweep split (400→save→restore→600 vs 1000 continuous)
- Format tag validation (rejects "cmc-rs-snapshot-v1")
- Edge kind corruption detection
- Topology mismatch detection
- Energy recomputed on load

## API stability

No API changes were made during Phase 0. The `sampling-core-v2` git tag marks the stabilized interface point.

## Contract document

`docs/SAMPLING_CORE_CONTRACT.md` documents all interface invariants: trial evaluation protocol, proposal ratio convention, phase lifecycle, cache invariants, ensemble independence, snapshot format, and error types.
