# MCMC.rs v0.1

`mcmc-rs` is Scuttle's statistical-inference MCMC layer. It is intentionally
independent of CMC.rs and QMC.rs:

```text
Carlo.rs
   ↑
   ├── CMC.rs
   ├── QMC.rs
   └── MCMC.rs
```

Carlo.rs owns generic execution lifecycle, RNG contexts and scalar online
measurements. MCMC.rs owns log densities, Markov kernels, warmup adaptation,
posterior traces and multi-chain diagnostics.

## v0.1 scope

- `LogDensity<[f64]>` and closure adapter
- synchronized Euclidean chain state
- isotropic/diagonal random-walk Metropolis
- component-wise Metropolis
- coordinate-wise slice sampling
- warmup-only Robbins-Monro scale adaptation
- online diagonal covariance adaptation
- contiguous row-major `MemoryTrace` with thinning and JSON/HDF5 export
- rank-normalized split R-hat, folded R-hat, bulk/tail ESS and MCSE
- deterministic Rayon multi-chain runner
- generic JSON checkpoints containing state, kernel, RNG and trace
- `McmcSampler` adapter for Carlo.rs

HMC and NUTS are deliberately deferred. Their future implementation can reuse
`DifferentiableLogDensity`, `EuclideanCache` and `TransitionReport` without
changing the v0.1 target/trace contracts.

## Example

```rust
use mcmc_rs::{run_multichain, FnLogDensity, McmcConfig, RandomWalkMetropolis};

let config = McmcConfig {
    chains: 4,
    warmup: 2_000,
    samples: 5_000,
    parameter_names: vec!["x".into(), "y".into()],
    ..McmcConfig::default()
};
let initial = vec![vec![-2.0, 1.0], vec![2.0, -1.0], vec![0.5, 2.0], vec![-0.5, -2.0]];
let kernel = RandomWalkMetropolis::isotropic(2, 0.5)?
    .with_scale_adaptation(0.234)?
    .with_diagonal_covariance_adaptation(1e-3)?;
let output = run_multichain(
    |_| FnLogDensity::new(|x: &[f64]| -0.5 * x.iter().map(|v| v * v).sum::<f64>()),
    |_| kernel.clone(),
    initial,
    config,
)?;
println!("R-hat = {}", output.diagnostics.parameters[0].rhat);
# Ok::<(), mcmc_rs::McmcError>(())
```

## Trace policy

Posterior draws are not stored through `Context::measure_array()`. A parameter
vector is not one scalar observable, and flattening it would make R-hat, ESS,
quantiles and per-parameter summaries impossible. Carlo measurements are used
only for online scalars such as acceptance and target-evaluation counts.

## Reproducibility

Each chain owns a target evaluator, kernel, RNG and trace. Chain seeds are
stable functions of `(base_seed, chain_id)`, so Rayon scheduling order does not
change trajectories.

## Optional features

Enable `hdf5` for post-run trace import/export:

```bash
cargo test -p mcmc-rs --features hdf5
```
