# MCMC.rs v0.4

`mcmc-rs` is Scuttle's statistical-inference MCMC layer. It remains independent
of CMC.rs and QMC.rs:

```text
Carlo.rs
   ↑
   ├── CMC.rs
   ├── QMC.rs
   └── MCMC.rs
```

Carlo.rs owns the generic run lifecycle, RNG contexts, scheduling and scalar
online measurements. MCMC.rs owns target densities, transition kernels,
warmup adaptation, constrained-parameter transforms, posterior traces,
replica exchange and multi-chain diagnostics.

## v0.4 capabilities

v0.4 retains all v0.1–v0.3 functionality and adds dynamic Hamiltonian Monte
Carlo:

- `Nuts<M>` with unit, diagonal or dense Euclidean metrics;
- binary trajectory doubling with a generalized metric-aware U-turn test;
- multinomial candidate selection in the log domain;
- configurable maximum tree depth and energy-error divergence threshold;
- atomic state commit: invalid or divergent trajectories never corrupt the
  accepted position, log density or gradient cache;
- dual-averaged step-size adaptation and fast/slow/fast windowed metric
  adaptation shared with static HMC;
- per-transition Hamiltonian energy, tree depth, depth-limit hit, leapfrog,
  acceptance-statistic and divergence reporting;
- trace storage for the new dynamic-HMC diagnostics;
- per-chain E-BFMI and aggregate maximum-tree-depth-hit diagnostics;
- exact JSON checkpoint continuation, including warmup state, metric state and
  RNG, while reconstructible trajectory workspaces remain un-serialized.

The package also keeps:

- random-walk, component-wise, slice and Gibbs kernels;
- static HMC;
- static kernel composition;
- diagonal and dense covariance adaptation;
- parameter transforms and differentiable transformed targets;
- parallel independent chains and replica exchange;
- contiguous in-memory traces, JSON export and optional HDF5 export;
- rank-normalized split R-hat, bulk/tail ESS and MCSE.

## NUTS example

```rust
use mcmc_rs::{
    run_multichain, DifferentiableLogDensity, LogDensity, McmcConfig, Nuts,
};

#[derive(Clone, Copy)]
struct Gaussian;

impl LogDensity<[f64]> for Gaussian {
    fn log_density(&mut self, x: &[f64]) -> f64 {
        -0.5 * x.iter().map(|value| value * value).sum::<f64>()
    }
}

impl DifferentiableLogDensity for Gaussian {
    fn log_density_and_gradient(
        &mut self,
        x: &[f64],
        gradient: &mut [f64],
    ) -> f64 {
        for (gradient, value) in gradient.iter_mut().zip(x.iter().copied()) {
            *gradient = -value;
        }
        self.log_density(x)
    }
}

fn main() -> Result<(), mcmc_rs::McmcError> {
    let warmup = 1_000;
    let output = run_multichain(
        |_| Gaussian,
        |_| {
            Nuts::diagonal(vec![1.0, 1.0], 0.2, 10)
                .and_then(|kernel| {
                    kernel.with_diagonal_adaptation(
                        warmup,
                        0.8,
                        1.0e-3,
                    )
                })
                .expect("valid NUTS configuration")
        },
        vec![
            vec![-2.0, -1.0],
            vec![2.0, 1.0],
            vec![-1.0, 2.0],
            vec![1.0, -2.0],
        ],
        McmcConfig {
            chains: 4,
            warmup,
            samples: 2_000,
            parameter_names: vec!["x".into(), "y".into()],
            ..McmcConfig::default()
        },
    )?;

    println!("E-BFMI: {:?}", output.diagnostics.chain_ebfmi);
    println!(
        "maximum-tree-depth hits: {}",
        output.diagnostics.max_tree_depth_hits
    );
    Ok(())
}
```

A complete correlated-Gaussian example is in
`examples/nuts_gaussian.rs`.

## NUTS diagnostics

`TransitionReport` records:

```text
acceptance_statistic
energy
energy_error
leapfrog_steps
tree_depth
max_tree_depth_reached
divergent
target_evaluations
gradient_evaluations
```

`MemoryTrace` stores retained energy, tree depth and depth-limit flags.
`diagnose()` reports per-chain E-BFMI and aggregate depth-limit hits in addition
to R-hat, ESS, MCSE, divergence count and mean acceptance.

E-BFMI is returned as `None` when the retained energy series is incomplete,
non-finite, too short or constant.

## Checkpoint contract

Persistent state includes:

- accepted position, log density and valid gradient cache;
- RNG state;
- metric and adaptation state;
- current step size and warmup windows;
- trace and transition report.

Private leapfrog/tree workspaces are skipped during serialization and rebuilt
from dimension metadata. Consequently, a divergent temporary phase point cannot
poison JSON serialization.

## Feature flags

Default builds are pure Rust apart from normal crate dependencies.

The optional `hdf5` feature enables post-run HDF5 trace export. In the validation
environment, `hdf5-sys 0.8.1` did not accept the installed HDF5 1.14.5 header,
so that feature remains dependent on a compatible system HDF5 installation.

## Validation

The v0.4 implementation was exercised with:

```bash
cargo fmt --all --check
cargo check -p mcmc-rs --all-targets
cargo clippy -p mcmc-rs --all-targets -- -D warnings
cargo test -p mcmc-rs
cargo check --workspace --all-targets
```

The MCMC.rs suite completed with 51 passing tests. See
`VALIDATION_REPORT.md` for the review and validation record.
