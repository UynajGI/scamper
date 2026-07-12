# Validation report

Date: 2026-07-12

## Scope validated

This package contains the Carlo.rs lifecycle/run-control extension and the CMC.rs sampling-core refactor. It deliberately does not add particle models, Wang-Landau, density-of-states sampling, or worm sectors.

## Static checks completed in this environment

- Parsed all 103 Rust source/test files in the workspace with the Tree-sitter Rust grammar: 0 files with syntax errors or missing syntax nodes.
- Parsed all 4 workspace `Cargo.toml` files with Python `tomllib`.
- Counted 134 Carlo.rs and 53 CMC.rs `#[test]` functions after the refactor.
- Confirmed the only `system.recompute_energy(model)` call left in CMC update kernels is the optional Metropolis audit interval; Wolff and Swendsen-Wang no longer use unconditional full-energy repair.
- Confirmed removed legacy cluster APIs (`fk_bond_probability`, `random_cluster_spin`, and full `membership.fill`) have no residual uses.
- Ran an independent affected-edge reference check over 5,000 randomized weighted multigraph cases, including parallel edges, self-loops, onsite fields, empty batches, and all-site batches: incremental and full energy differences agreed in every case.
- Parsed the final archive with `unzip -t` after packaging.

## Added regression coverage

Carlo.rs tests now cover:

- exact fixed warmup boundary semantics;
- explicit `RunPhase` transitions and zero-warmup lifecycle;
- lifecycle start/end hooks;
- legacy checkpoint phase inference;
- preservation of an explicit adaptation phase across checkpoint state;
- algorithm-driven `AdaptiveRunControl` transition and stop decisions;
- explicit thermalization remaining authoritative beyond a fixed counter threshold.

CMC.rs tests now cover:

- transactional site trials (no mutation before commit);
- generic Metropolis-Hastings rejection and Hastings correction;
- canonical beta application exactly once;
- affected-edge batch energy versus full recomputation;
- direct multi-body Hamiltonian scratch-backed batch fallback;
- parallel-edge and self-loop batch correctness;
- cached energy consistency for Metropolis, Wolff, SW, heat-bath-compatible paths, and microcanonical reflection;
- SW independent per-cluster assignments;
- adaptive proposal freezing outside thermalization;
- visit-order workspace restoration after random traversal.

## Toolchain limitation

The execution container does not provide `rustc`, `cargo`, `rustfmt`, or `clippy`, and outbound access needed to install the Rust toolchain is unavailable. Therefore the following commands could not be executed here and must be the first local/CI verification step:

```bash
cargo fmt --all --check
cargo check --workspace --all-targets --all-features
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo test --workspace --all-features
```

No claim of successful Rust compilation or runtime test execution is made in this report.
