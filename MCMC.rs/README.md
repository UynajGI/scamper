# MCMC.rs v0.5

`mcmc-rs` is Scuttle's statistical-inference MCMC layer. Carlo.rs owns the
generic run lifecycle; MCMC.rs owns target densities, kernels, warmup,
transforms, traces, replica exchange and convergence diagnostics.

## New in v0.5

- generalized multinomial NUTS termination using summed trajectory momentum;
- full merged-tree and both cross-subtree U-turn checks;
- allocation-free summed-momentum products for unit, diagonal and dense metrics;
- bounded reasonable-step-size search shared by static HMC and NUTS;
- optional re-search after metric-window updates and exact checkpoint recovery;
- public central finite-difference gradient validation with component reports;
- backward-compatible v0.4 HMC/NUTS JSON defaults.

All v0.1–v0.4 kernels, adaptations, transforms, multi-chain diagnostics,
replica exchange, JSON trace/checkpoint support and optional HDF5 export remain.

## Automatic initial step size

```rust
use mcmc_rs::{Nuts, StepSizeSearch};

let warmup = 1_000;
let kernel = Nuts::diagonal(vec![1.0, 1.0], 1.0, 10)?
    .with_diagonal_adaptation(warmup, 0.8, 1.0e-3)?
    .with_step_size_search(StepSizeSearch::default())?;
# Ok::<(), mcmc_rs::McmcError>(())
```

The search is warmup-only, doubles or halves a one-step scale within configured
bounds, and by default repeats after metric updates before restarting dual
averaging.

## Gradient validation

```rust
use mcmc_rs::{check_gradient, GradientCheckConfig};

let report = check_gradient(
    &mut target,
    &[0.3, -1.2],
    GradientCheckConfig::default(),
)?;
assert!(report.passed, "{report:?}");
# Ok::<(), mcmc_rs::McmcError>(())
```

See `examples/nuts_gaussian.rs`, `examples/gradient_check.rs`, `REVIEW_V04.md`
and `VALIDATION_REPORT.md`.
