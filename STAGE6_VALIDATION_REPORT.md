# Stage 5 review and Stage 6 validation report

## Decision

The Stage 5 persistent classical worm implementation did not require a 5.1
correctness repair.  Review and real tests confirmed:

- physical/worm sector invariants;
- transactional local edge toggles;
- reciprocal open/close and head-step Hastings terms;
- exact energy/log-graph-weight cache audits;
- exact four-site graph expansion observables;
- endpoint-correlation agreement with the analytic Ising ring;
- exact future trajectory after JSON snapshot restoration;
- rejection of self loops, negative effective couplings and inconsistent
  checkpoint transition counters.

The review therefore proceeded to Stage 6.

## Stage 6 scope

### Carlo.rs clocks

- `SimulationClock::{Sweeps, Attempts, AcceptedMoves, EventTime}`;
- explicit attempt, executed-move and event-time counters in `Context`;
- serde-defaulted clock fields in `ContextCheckpoint`;
- typed clock access and JSON checkpoint round trip;
- legacy JSON checkpoint compatibility.

The existing optional HDF5 backend is not claimed as Stage 6 validated; see the
known workspace issue below.

### Conserved dynamics

`KawasakiCore` performs nearest-edge exchange of unlike Ising spins through the
existing two-site transactional batch-energy path and log-domain canonical
acceptance.  The scheduler adapter records attempted and accepted exchanges
without assigning them a physical time.

### Rejection-free continuous time

- generic finite-catalog `RejectionFreeModel`;
- direct `GillespieKernel` event selection and exponential waiting times;
- exact fixed-event-time observation windows;
- zero-rate absorbing-state handling;
- stable Glauber and continuous-time Metropolis rate laws;
- `KineticIsingModel` reference backend.

### BKL / n-fold way

`BklIsingKernel` stores per-site rates in a Fenwick tree.  It provides
`O(log N)` event selection and local rate updates on arbitrary weighted CSR
Ising graphs.  Versioned snapshots include the exact rate vector and Fenwick
accumulation structure so checkpoint restoration preserves the future event
trajectory bit-for-bit for a fixed RNG state.

### Event-chain Monte Carlo

`HardSphereEventChain<D>` implements straight lifted chains for identical hard
spheres in periodic orthorhombic cells.  The first backend uses exact `O(N)`
collision search as a correctness baseline.  It validates overlap, geometry,
wrapping, lifting and versioned snapshots.  No cell-list acceleration claim is
made in this stage.

### Cross-cutting work

- Stage 6 kernels participate in the existing `cache-audit` policy;
- all stochastic acceptance/rate calculations avoid unstable positive
  exponentials;
- Criterion coverage includes Kawasaki attempts, direct Gillespie events,
  Fenwick BKL events, event-chain lifted distance and BKL JSON serialization;
- the benchmark prints BKL integrated autocorrelation time, ESS and ESS/s.

## Real validation

The following commands were run with Rust 1.90.0:

```text
cargo fmt --all --check
cargo check --workspace --all-targets
cargo clippy --workspace --all-targets -- -D warnings
cargo test -p cmc-rs --lib --tests
cargo test -p cmc-rs --features cache-audit --lib --tests
cargo test -p carlo-rs --lib --tests
cargo test -p mcmc-rs --lib --tests
cargo test -p qmc-rs --lib --tests
cargo test --workspace --doc
cargo bench -p cmc-rs --bench dynamics_bench --no-run
cargo bench -p cmc-rs --bench dynamics_bench -- \
  --warm-up-time 0.01 --measurement-time 0.02 --sample-size 10
```

Results:

```text
fmt                                                   passed
workspace all-target check                            passed
workspace all-target clippy -D warnings               passed
CMC default library/integration tests                 163 passed
CMC cache-audit library/integration tests             163 passed
Carlo library/integration tests                       136 passed
MCMC library/integration tests                         52 passed
QMC library/integration tests                          54 passed
Total default library/integration tests               405 passed
doctests                                                3 passed, 7 ignored
dynamics benchmark optimized build                    passed
dynamics benchmark smoke groups                         5 passed
```

The short smoke benchmark produced the following diagnostic values on this
container.  They validate execution and units, not stable comparative
performance because the measurement window was intentionally tiny:

```text
BKL statistical pilot: tau_int=91.917287, ESS=22.281, ESS/s=1038.549
Kawasaki attempted exchanges: approximately 14.8-19.4 million/s
Direct Gillespie events: approximately 7.0-62.1 thousand/s
Fenwick BKL events: approximately 1.15-1.44 million/s
Hard-sphere lifted distance: approximately 2.63-2.82 million units/s
BKL JSON serialization: approximately 0.78-0.89 ms
```

## Stage 6 physics tests

The dedicated test target verifies:

1. rate-proportional direct Gillespie event frequencies;
2. the exponential waiting-time mean;
3. exact fixed event-time observation boundaries without overshoot;
4. Kawasaki magnetization conservation and energy-cache consistency;
5. BKL rate/Fenwick cache consistency;
6. exact BKL future trajectory after snapshot restoration;
7. fixed-time BKL equilibrium energy against exact enumeration of a four-site
   periodic Ising chain;
8. analytic hard-sphere collision location and lifting transfer;
9. periodic wrapping and exact event-chain snapshot continuation;
10. distinct scheduler clocks for kinetic and geometric algorithms.

Three internal dynamics tests additionally validate rate-law detailed balance,
Fenwick selection/update equivalence to linear weights, and constructor energy
cache repair.

## Known pre-existing optional-feature issue

`cargo check -p carlo-rs --features hdf5` cannot currently validate any HDF5
path.  The unmodified dependency `hdf5-sys 0.8.1` first rejects the installed
HDF5 1.14.5 header.  Relaxing that dependency parser locally for diagnosis
exposes older, repository-wide Carlo.rs HDF5 API mismatches such as the removed
`hdf5::H5` import and `create_dataset_simple` calls.  These failures predate
Stage 6 and occur throughout the existing HDF5 implementation, not in the
default Stage 6 code path.  The package therefore does not claim
`--all-features` success.

## Files added or materially changed

- `Carlo.rs/src/clock.rs`
- `Carlo.rs/src/context.rs`
- `Carlo.rs/src/lib.rs`
- `Carlo.rs/tests/context_test.rs`
- `CMC.rs/src/dynamics/*`
- `CMC.rs/src/algorithms/mod.rs`
- `CMC.rs/src/lib.rs`
- `CMC.rs/tests/dynamics_stage6_test.rs`
- `CMC.rs/benches/dynamics_bench.rs`
- `CMC.rs/CLASSICAL_DYNAMICS.md`
- relevant README, benchmark and migration documentation

No Stage 7 work was started.
