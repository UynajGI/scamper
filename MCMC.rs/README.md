# MCMC.rs v0.2

`mcmc-rs` is Scuttle's statistical-inference MCMC layer. It remains independent
of CMC.rs and QMC.rs:

```text
Carlo.rs
   ↑
   ├── CMC.rs
   ├── QMC.rs
   └── MCMC.rs
```

Carlo.rs owns generic execution lifecycle, RNG contexts and scalar online
measurements. MCMC.rs owns target densities, statistical transition kernels,
warmup adaptation, constrained-parameter transforms, posterior traces,
replica exchange and multi-chain diagnostics.

## v0.2 capabilities

v0.2 retains the v0.1.1 random-walk, component-wise and slice kernels and adds:

- target-parameterized `TransitionKernel<T>` for model-specific kernels;
- static `Then`, `Repeat` and `Mixture` kernel composition;
- atomic target-specific Gibbs/block updates through `GibbsUpdate<T>`;
- dense Gaussian proposals backed by row-major lower Cholesky factors;
- online dense covariance adaptation with regularized Cholesky finalization;
- positive, interval, ordered, simplex, identity and product transforms;
- `TransformedTarget<T, B>` with automatic log-Jacobian correction;
- local Rayon replica exchange with generic cross-target exchange ratios;
- fixed-slot trace and per-edge exchange diagnostics;
- HDF5 trace writing updated to the hdf5 0.8 dataset-builder API.

HMC and NUTS remain deliberately deferred to v0.3/v0.4. Their future
implementation can reuse `DifferentiableLogDensity`, `EuclideanCache`, dense
geometry, transforms and `TransitionReport`.

## Independent chains

```rust
use mcmc_rs::{run_multichain, FnLogDensity, McmcConfig, RandomWalkMetropolis};

let config = McmcConfig {
    chains: 4,
    warmup: 2_000,
    samples: 5_000,
    parameter_names: vec!["x".into(), "y".into()],
    ..McmcConfig::default()
};
let initial = vec![
    vec![-2.0, 1.0],
    vec![2.0, -1.0],
    vec![0.5, 2.0],
    vec![-0.5, -2.0],
];
let kernel = RandomWalkMetropolis::isotropic(2, 0.5)?
    .with_scale_adaptation(0.234)?
    .with_dense_covariance_adaptation(1e-4)?;
let output = run_multichain(
    |_| FnLogDensity::new(|x: &[f64]| -0.5 * x.iter().map(|v| v * v).sum::<f64>()),
    |_| kernel.clone(),
    initial,
    config,
)?;
println!("R-hat = {}", output.diagnostics.parameters[0].rhat);
# Ok::<(), mcmc_rs::McmcError>(())
```

## Target-specific Gibbs update

`TransitionKernel<T>` is parameterized by the concrete target type. A Gibbs
updater can therefore use model-specific conditional parameters without
runtime type erasure:

```rust
use mcmc_rs::proposal::standard_normal;
use mcmc_rs::{
    EuclideanState, GibbsKernel, GibbsUpdate, GibbsUpdateResult, LogDensity,
    McmcError, SamplingPhase,
};
use rand::Rng;

struct Model { conditional_mean: f64 }
impl LogDensity<[f64]> for Model {
    fn log_density(&mut self, x: &[f64]) -> f64 {
        -0.5 * (x[0] - self.conditional_mean).powi(2)
    }
}

struct ExactConditionalNormal;
impl GibbsUpdate<Model> for ExactConditionalNormal {
    fn update<R: Rng + ?Sized>(
        &mut self,
        target: &mut Model,
        _current: &EuclideanState,
        proposal: &mut [f64],
        rng: &mut R,
        _phase: SamplingPhase,
    ) -> Result<GibbsUpdateResult, McmcError> {
        proposal[0] = target.conditional_mean + standard_normal(rng);
        Ok(GibbsUpdateResult::requiring_target_evaluation())
    }
}

let kernel = GibbsKernel::new(ExactConditionalNormal);
# let _ = kernel;
```

The updater writes only to a private proposal workspace. If it returns an error
or produces an invalid state, the accepted chain state remains unchanged.

## Constrained targets

```rust
use mcmc_rs::{LogDensity, Positive, TransformedTarget};

struct Exponential;
impl LogDensity<[f64]> for Exponential {
    fn log_density(&mut self, x: &[f64]) -> f64 {
        if x[0] > 0.0 { -x[0] } else { f64::NEG_INFINITY }
    }
}

let target = TransformedTarget::new(Exponential, Positive)?;
# let _ = target;
# Ok::<(), mcmc_rs::McmcError>(())
```

Kernels see the unconstrained coordinate `z`; the wrapped target evaluates
`x = exp(z)` and adds `log |dx/dz|` automatically.

## Replica exchange

`run_parallel_tempering` keeps targets and kernels attached to fixed ladder
slots and swaps accepted states. Each exchange ratio uses cross-evaluation:

```text
log α = log π_i(x_j) + log π_j(x_i) - log π_i(x_i) - log π_j(x_j)
```

The ladder value is only metadata passed to the factories. This supports full
posterior tempering, likelihood-only tempering and arbitrary model-parameter
ladders without hard-coding a beta convention into the runtime.

Traces are attached to fixed ladder slots. A production draw is recorded after
its local transition; when that step completes an exchange interval, the
neighbor exchange is attempted immediately afterward. The post-exchange state
therefore appears through the next local transition, while `final_position`
always reflects the last accepted local or exchange transition.

## Trace policy

Posterior draws are not stored through `Context::measure_array()`. A parameter
vector is not one scalar observable, and flattening it would make R-hat, ESS,
quantiles and per-parameter summaries impossible. Carlo measurements remain
for online scalars such as acceptance and target-evaluation counts.
