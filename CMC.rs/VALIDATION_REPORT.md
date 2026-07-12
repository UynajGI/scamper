# Validation report

Date: 2026-07-12

## Phase 2 implementation status

Phase 2 adds the first continuous-system backend without changing the lattice public API:

- const-generic periodic `OrthorhombicCell<D>` geometry and minimum-image displacement;
- AoS `ParticleConfiguration<D>` with `u16` species labels;
- species-aware `PairPotential` and Lennard-Jones 12-6 with truncated, shifted-potential and shifted-force cutoffs;
- packed `CellList<D>` buckets, precomputed periodic neighbor-cell stencils and O(1) accepted membership patches;
- transactional `ParticleTranslation<D>` evaluation using reusable `ParticleEnergyPatch` scratch;
- warmup-only adaptive `TranslateParticle<D>` proposal and frozen production kernel;
- scheduler-ready `LennardJonesNvt<D>` adapter with energy, energy-per-particle and density measurements;
- Criterion benchmarks for cell-list translation attempts and the O(N²) reference energy.

New particle tests cover analytic pair energy, periodic minimum image, boundary-crossing commits, rejection immutability, cell-list candidate completeness, repeated cache patches, cutoff limits, deterministic trajectories, adaptation freeze, scheduler integration, and a two-particle canonical energy distribution checked against deterministic midpoint quadrature.

Validation completed for this handoff:

```text
cargo fmt --all --check                                      PASS
cargo check --workspace --all-targets                        PASS
cargo test -p cmc-rs                                         PASS: 95 passed, 1 ignored
cargo clippy -p cmc-rs --lib --no-deps -- -D warnings        PASS
cargo clippy (particle integration tests and benchmark)       PASS
```

The Criterion benchmark target is type-checked and linted. Run measurements locally with:

```bash
cargo bench -p cmc-rs --bench particle_bench
```

A full `cargo clippy --workspace --all-targets -- -D warnings` under Rust 1.88 also reports pre-existing `uninlined_format_args` warnings in Phase 0/1 Carlo.rs and CMC.rs test code; no Stage-2 target emits a Clippy warning.

## Phase 1 completion status

Phase 1 (CMC.rs module reorganization) is complete. The 15 flat source files have been
reorganized into `core/`, `lattice/`, `algorithms/`, and `observables/` subdirectories
(30 files total), with a new `AcceptanceRule<D>` trait extracted from the trial layer.
Public API is preserved via flat re-exports from `lib.rs`.

```bash
cargo fmt --all --check     # PASS
cargo check --workspace      # PASS
cargo clippy --workspace     # PASS (no warnings)
cargo test --workspace       # 245 passed, 7 ignored
```

### Phase 1 changes

- **Module tree**: `core/` (6 files), `lattice/` (5 files), `algorithms/` (8 files), `observables/` (4 files), plus 3 top-level modules
- **`AcceptanceRule<D>` trait**: extracted to `core/acceptance.rs` with `MetropolisHastingsAcceptance`; `metropolis_hastings_step` now accepts `&impl AcceptanceRule`
- **Algorithm split**: 776-line `algorithm.rs` → `algorithms/{common,metropolis,wolff,swendsen_wang,heat_bath,microcanonical,hybrid}.rs`
- **Cache split**: `EnergyPatch`/`BatchEnergyWorkspace` → `core/cache.rs`; move types stay in `core/move.rs`
- **Observables split**: `observables.rs` → `observables/{energy,magnetization,correlation}.rs`; `compute_correlation_1d` moved from `postprocess.rs`

## Phase 0 completion status

The sampling-core-v2 stabilization (Phase 0) was completed before Phase 1. All local verification passes:

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
