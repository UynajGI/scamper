# CMC.rs physical validation tests

This directory is the physics-correctness gate for CMC.rs. It is included by
`../physics_validation.rs` because Cargo integration tests require a top-level
harness.

## Placement policy

- **`tests/physics/`**: physical correctness. Default tests are deterministic:
  finite-state enumeration, exact detailed balance, analytical identities,
  cache transactions and geometric invariants.
- **Ignored tests**: long stochastic convergence checks against an exact target.
  They are secondary evidence and use explicit error budgets.
- **`examples/`**: user-facing demonstrations only. An example must never be
  treated as a correctness gate.
- **`benches/`**: throughput and statistical efficiency only. A faster wrong
  algorithm must still fail the tests.

## Commands

```bash
cargo test -p cmc-rs --test physics_validation
cargo test -p cmc-rs --features cache-audit --test physics_validation
cargo test -p cmc-rs --test physics_validation -- --ignored --test-threads=1
```

## Strictness rules

1. Exact finite state spaces are enumerated whenever feasible.
2. Detailed balance is checked on transition probabilities, not inferred from
   a visually plausible time series.
3. Floating tolerances are tied to roundoff scale for deterministic formulas.
4. Stochastic tests compare to analytical/exact targets and are ignored by
   default to avoid flaky CI.
5. Event-driven dynamics are sampled in event time, not event count.
6. Continuous-energy DOS reweighting is not called exact when bin centers are
   only a discretization approximation.
