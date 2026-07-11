# Spin-boson wormhole implementation status

## Completed in this delivery

- QMC-side kernel abstraction and fixed measured-sweep schedule;
- continuous-time retarded four-leg configuration;
- circular worldline indexing and invariant validation;
- normalized single-mode, power-law, and tabulated bath proposals;
- diagonal insertion/removal with proposal cancellation;
- generic exact local-detailed-balance low-bounce scattering table;
- closed directed-loop update with wormhole endpoint traversal;
- Jaynes-Cummings, XXZ, XYZ, and rotated original spin-boson/Rabi catalogs;
- Carlo.rs `MonteCarlo` / `FromParams` adapter;
- warmup-only schedule adaptation;
- longitudinal observables and update diagnostics;
- unit and integration test source;
- runnable Carlo.rs example and architecture documentation.

## Deliberate boundary of the impurity engine

The delivered engine covers positive, sign-free retarded four-leg impurity
vertices. The following are not silently approximated:

- spatial lattice indices and exchange vertices;
- variational/projector algorithms;
- complex or sign-indefinite matrix spectral functions;
- the special single-time in-plane magnetic-field vertex construction;
- off-diagonal improved-loop estimators;
- bath-observable reconstruction estimators;
- a proof-optimal analytic/linear-programmed scattering policy (the delivered
  residual-flow policy is exact and low-bounce, but not claimed LP-optimal for
  every custom graph).

These are extension modules, not missing pieces of the four-leg impurity
sampler. The delivered residual-flow table is exact and low-bounce; the
reference Metropolis policy is retained for debugging and cross-checks.

## Validation status in this build environment

Static source audit performed:

- all newly added Rust source files have balanced delimiters and module paths;
- model catalogs were mirrored and checked for probability-row normalization
  and local detailed balance;
- mirrored mixed diagonal/loop chains for JC, XXZ, XYZ, and rotated Rabi stayed
  periodic and closed over thousands of sweeps;
- the source contains unit/integration tests for the same invariants.

A Rust toolchain was not present in the execution container, and package
installation was unavailable, so `cargo fmt`, `cargo clippy`, and `cargo test`
could not be executed here. Run `QMC.rs/scripts/verify.sh` in a Rust-enabled
environment before merging.
