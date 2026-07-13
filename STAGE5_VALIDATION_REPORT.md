# Stage 5 validation report

## Scope

This delivery reviews the Stage 4.1 cross-cut hardening and adds Stage 5 only:
a persistent classical worm extended-configuration-space framework plus a
ferromagnetic, zero-field Ising high-temperature graph implementation. Stage 6
kinetic/event-chain algorithms are not included.

## Stage 4.1 review

The review found no structural correctness problem in the generalized-ensemble,
cache-audit, deterministic RNG-stream, checkpoint, or log-domain acceptance
changes. The current tree retains the intended boundaries:

- generalized and worm kernels do not acquire canonical beta-tempering APIs;
- cache audits validate before any recomputation or repair;
- concrete RNG state remains part of scheduler/checkpoint state;
- Metropolis-Hastings decisions use log probabilities;
- Wang-Landau DOS/histogram state remains checkpoint validated.

No Stage 4.1 algorithm was rolled back. The only changes adjacent to Stage 4
are documentation/API integration for the new Stage 5 module.

## Stage 5 design

The reusable layer consists of:

- `WormSector` and `WormState<Configuration, Defect>`;
- a model-owned `WormModel` trait for local proposals, trial patches, commits,
  cache validation, and endpoint bins;
- `WormKernel<Model>` with persistent physical/worm sectors;
- exact open, close, local-step branch, and local-proposal Hastings factors;
- log-domain acceptance through `carlo_rs::accept_log_probability`;
- transition diagnostics and optional endpoint-pair histograms;
- versioned JSON snapshots with runtime consistency validation.

The first backend is the ferromagnetic zero-field Ising high-temperature graph
representation on an arbitrary loop-free `CsrLattice`. Weighted and parallel
edges are supported. Negative effective couplings and self-loops are rejected
because this first backend does not implement signed weights or self-loop defect
semantics.

The chain may remain in the worm sector across sweep and checkpoint boundaries.
This is deliberate: the extended distribution is stationary and endpoint
observables have an explicit sampling measure.

## Correctness validation

The Stage 5 tests cover:

- generic sector invariants;
- transactional and reversible edge toggles;
- reciprocal open/close acceptance ratios;
- local branch and irregular-degree Hastings corrections;
- weighted parallel-edge graph cache auditing;
- zero-coupling finite behavior;
- exact high-temperature graph enumeration on a four-site ring;
- physical-sector occupied-edge and energy estimators versus exact enumeration;
- endpoint-count correlation versus the analytic two-point function;
- exact future trajectory after JSON snapshot restoration;
- rejection of inconsistent checkpoint counters;
- Carlo.rs scheduler construction and measurements.

The benchmark smoke run reported:

```text
STAT_EFF classical_worm tau_int=1.647775 ess=2485.776 ess_per_second=8294.032
```

This value is a deterministic smoke-run diagnostic on the benchmark setup, not
a universal performance claim.

## Commands executed

The following commands were run with Rust/Cargo 1.90.0:

```text
cargo fmt --all --check
cargo check --workspace --all-targets
cargo clippy --workspace --all-targets -- -D warnings
cargo test -p cmc-rs --lib --tests
cargo test -p cmc-rs --lib --tests --features cache-audit
cargo test -p carlo-rs --lib --tests
cargo test -p mcmc-rs --lib --tests
cargo test -p qmc-rs --lib --tests
cargo test --workspace --doc
cargo bench -p cmc-rs --no-run
cargo bench -p cmc-rs --bench worm_bench -- --test
```

Results:

- formatting: passed;
- workspace all-target compilation: passed;
- workspace all-target Clippy with warnings denied: passed;
- library and integration tests: 367 passed, 0 failed;
- CMC.rs tests: 152 passed, 0 failed;
- CMC.rs with `cache-audit`: 152 passed, 0 failed;
- doctests: 3 passed, 7 ignored, 0 failed;
- all five CMC.rs benchmark executables compiled;
- all four classical-worm benchmark smoke groups passed.

A single monolithic `cargo test --workspace --lib --tests` invocation exceeded
the command time window during test-binary generation. The same workspace tests
were then run per crate, using the identical final source tree, and all passed.

## Scope limits

This release does not claim a universal worm algorithm. Integer-current, dimer,
loop-gas, nonzero-field, frustrated/sign-problem, and quantum directed-loop
representations require their own model-specific extended spaces. QMC.rs remains
separate, and Stage 6 was not started.
