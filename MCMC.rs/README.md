# MCMC.rs v0.3

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

## v0.3 capabilities

v0.3 retains every v0.2 kernel and adds a production-oriented static HMC slice:

- `StaticHmc<M>` with fixed trajectory length and atomic accept/reject commit;
- `UnitMetric`, `DiagonalMetric` and `DenseMetric` inverse-mass geometries;
- allocation-conscious `LeapfrogIntegrator` and reusable `PhasePoint` storage;
- accepted-state gradient caching through `EuclideanCache`;
- finite-trajectory and configurable energy-error divergence detection;
- Nesterov dual averaging for step size;
- fast/slow/fast windowed warmup with diagonal or dense metric adaptation;
- differentiable positive, interval, ordered, simplex and product transforms;
- `TransformedTarget<T, B>: DifferentiableLogDensity` with analytic pullback and
  log-Jacobian gradients;
- HMC scalar measurements in the Carlo.rs adapter: energy error, divergence,
  gradient evaluations, leapfrog steps and frozen step size.

NUTS and dynamic trajectory building remain deliberately deferred to v0.4.

## Static HMC

```rust
use mcmc_rs::{
    run_multichain, DifferentiableLogDensity, LogDensity, McmcConfig, StaticHmc,
};

#[derive(Clone, Copy)]
struct Gaussian;

impl LogDensity<[f64]> for Gaussian {
    fn log_density(&mut self, x: &[f64]) -> f64 {
        -0.5 * x.iter().map(|value| value * value).sum::<f64>()
    }
}

impl DifferentiableLogDensity for Gaussian {
    fn log_density_and_gradient(&mut self, x: &[f64], gradient: &mut [f64]) -> f64 {
        for (gradient, value) in gradient.iter_mut().zip(x.iter().copied()) {
            *gradient = -value;
        }
        self.log_density(x)
    }
}

let warmup = 1_000;
let output = run_multichain(
    |_| Gaussian,
    |_| {
        StaticHmc::diagonal(vec![1.0, 1.0], 0.15, 8)
            .and_then(|kernel| kernel.with_diagonal_adaptation(warmup, 0.8, 1e-3))
            .expect("valid HMC configuration")
    },
    vec![vec![-2.0, 1.0], vec![2.0, -1.0], vec![0.5, 2.0], vec![-0.5, -2.0]],
    McmcConfig {
        chains: 4,
        warmup,
        samples: 2_000,
        parameter_names: vec!["x".into(), "y".into()],
        ..McmcConfig::default()
    },
)?;
println!("R-hat = {}", output.diagnostics.parameters[0].rhat);
# Ok::<(), mcmc_rs::McmcError>(())
```

When HMC adaptation is enabled, the `total_warmup` supplied to the kernel must
match the number of warmup transitions run by `ChainRunner`, `McmcConfig`, or
the Carlo.rs scheduler. Entering sampling early returns an error instead of
silently freezing an incomplete metric.

## Metric convention

Metrics store the inverse mass matrix `G = M^-1`:

```text
p ~ Normal(0, M)
K(p) = 0.5 pᵀ G p
q̇ = G p
```

Windowed adaptation estimates the position covariance and installs it directly
as `G`. This makes the linearized frequencies approximately isotropic for a
Gaussian target. Dense metrics cache a lower Cholesky factor and keep momentum
sampling, velocity and kinetic energy synchronized to the same regularized
matrix.

## Constrained HMC

```rust
use mcmc_rs::{
    DifferentiableLogDensity, LogDensity, Positive, TransformedTarget,
};

struct Exponential;
impl LogDensity<[f64]> for Exponential {
    fn log_density(&mut self, x: &[f64]) -> f64 {
        if x[0] > 0.0 { -x[0] } else { f64::NEG_INFINITY }
    }
}
impl DifferentiableLogDensity for Exponential {
    fn log_density_and_gradient(&mut self, x: &[f64], gradient: &mut [f64]) -> f64 {
        gradient[0] = -1.0;
        self.log_density(x)
    }
}

let target = TransformedTarget::new(Exponential, Positive)?;
# let _ = target;
# Ok::<(), mcmc_rs::McmcError>(())
```

The HMC state is the unconstrained coordinate `z`; `TransformedTarget` computes
`x = transform(z)`, adds `log |dx/dz|`, pulls the constrained gradient back to
`z`, and adds the analytic log-Jacobian gradient.

## Existing v0.2 functionality

The following remain available and compatible:

- random-walk, component-wise, slice and Gibbs kernels;
- `Then`, `Repeat` and `Mixture` static composition;
- diagonal and dense random-walk covariance adaptation;
- `MemoryTrace`, thinning, JSON/HDF5 trace persistence;
- rank-normalized split R-hat, bulk/tail ESS and MCSE;
- deterministic Rayon multi-chain execution;
- fixed-slot generic replica exchange;
- generic serde chain checkpoints including persistent kernel/adaptation and RNG state;
- exact f64 JSON round trips, with transient HMC trajectory workspaces rebuilt after restore.

## Trace policy

Posterior vectors remain in `MemoryTrace`, not `Context::measure_array()`.
`MemoryTrace` records divergence and energy error for HMC. Carlo.rs measurements
provide online scalar summaries including gradient evaluations and leapfrog
work. Static HMC has a fixed leapfrog count and a frozen production step size,
so those values do not require additional per-draw vector columns.
