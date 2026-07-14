# MCMC.rs Architecture — v0.5

## Core invariants

1. `ChainState.position`, cached `log_density`, and a valid gradient cache always
   describe the same accepted position after a completed transition.
2. Accepted positions, log densities and gradients are finite. `NEG_INFINITY`
   is allowed only as a target response for a proposed point outside support.
3. Every kernel mutates private proposal/trajectory workspace first and commits
   accepted state atomically.
4. Adaptation is legal only in `SamplingPhase::Warmup`; production sampling
   freezes step size and geometry.
5. Metrics store inverse mass `G = M^-1`. Momentum sampling, kinetic energy,
   velocity and U-turn checks all use that same geometry.
6. Leapfrog and NUTS tree construction never mutate the accepted chain state.
7. Reconstructible trajectory workspaces are omitted from checkpoints.
8. One completed kernel call increments chain iteration exactly once, regardless
   of the number of leapfrog or component substeps.
9. Raw posterior draws and Hamiltonian diagnostics live in MCMC.rs traces;
   Carlo.rs measurements remain scalar online summaries.
10. Independent chains and temperature slots own independent target workspaces,
    kernels, adaptation state, RNG streams and traces.

## Target interfaces

```rust,ignore
pub trait LogDensity<S: ?Sized>: Send {
    fn log_density(&mut self, state: &S) -> f64;
}

pub trait DifferentiableLogDensity: LogDensity<[f64]> {
    fn log_density_and_gradient(
        &mut self,
        position: &[f64],
        gradient: &mut [f64],
    ) -> f64;
}
```

The combined differentiable interface avoids a duplicate model evaluation per
leapfrog step and permits target-owned reusable workspaces.

## Metric abstraction

`Metric` supplies:

- momentum sampling;
- kinetic energy;
- velocity `G p`;
- displacement–velocity products;
- velocity products against one momentum or a sum of momenta;
- optional diagonal or dense inverse-mass updates.

All new product primitives have correct default implementations for external
metric implementations. Unit, diagonal and dense built-ins override them
without allocating temporary velocity vectors in the NUTS tree hot path.

## Static HMC

`StaticHmc<M>` samples momentum, executes a fixed number of leapfrog steps,
checks the final Hamiltonian error and performs a Metropolis correction.
Candidate position, log density and gradient are committed together only after
a successful acceptance decision.

## NUTS

`Nuts<M>` uses iterative outer doubling and recursive binary subtree building.

For each transition:

1. synchronize or evaluate the accepted-state gradient;
2. sample momentum and compute the initial Hamiltonian;
3. repeatedly choose a direction and build a subtree of depth `d`;
4. combine valid subtree candidates with log-domain multinomial weighting;
5. stop on a generalized U-turn, divergence or configured depth limit;
6. atomically commit the selected candidate or mark a rejection;
7. feed the mean trajectory acceptance statistic to warmup adaptation.

A subtree that terminates internally because of a U-turn or divergence
contributes evaluation and diagnostic counts, but its candidate mass is not
merged into the outer trajectory. This preserves the valid trajectory set used
for multinomial selection.

Each subtree retains its summed momentum `rho`. The generalized termination
test checks the full merged trajectory and both cross-subtree joins. In compact
form, every checked segment uses its endpoint displacement together with the
metric velocity of that segment's summed momentum. This catches a turn that can
be hidden by checking only the outermost endpoint momenta.

Tree depth is capped at 20 internally to bound the largest possible trajectory.
Users normally configure a lower limit such as 8–12.

## Warmup

Static HMC and NUTS share `HmcWarmup`:

- an optional bounded doubling/halving search finds a reasonable one-step
  initial scale before dual averaging;
- the search can repeat after metric-window updates;
- dual averaging tunes step size toward the requested acceptance statistic;
- slow windows accumulate diagonal or dense covariance estimates;
- geometry updates occur only at slow-window boundaries;
- dual averaging restarts after a geometry update;
- the terminal buffer retunes the final step size;
- entering sampling before configured warmup completes is an error.

## Gradient validation

`check_gradient` compares `DifferentiableLogDensity` output with central finite
differences at a supplied point. It reports component-wise analytic, numerical,
absolute and relative error and treats non-finite/support-crossing probes as
structured errors. It is a model-development tool and is not called from the
transition hot path.

## Trace and diagnostics

`MemoryTrace` uses contiguous row-major position storage and parallel scalar
columns. v0.4 added backward-compatible serde-default columns for:

- Hamiltonian energy;
- tree depth;
- maximum-tree-depth reached.

When an old v0.3 JSON trace receives a new draw, missing columns are backfilled
with `None` or `0` before append.

Cross-chain diagnostics include:

- rank-normalized and folded split R-hat;
- bulk and tail ESS;
- MCSE;
- posterior moments and quantiles;
- divergence count and mean acceptance;
- per-chain E-BFMI;
- aggregate maximum-tree-depth hits.

## Checkpoints

A checkpoint serializes accepted state, valid gradient cache, persistent kernel
configuration, metric/adaptation state, RNG, trace and target fingerprint.
NUTS `PhasePoint` and integrator workspaces are serde-skipped because they are
fully reconstructible and may contain non-finite values after a rejected
numerical trajectory.

JSON uses serde_json float round-tripping so resumed trajectories preserve
future random decisions and floating-point state exactly.

## Extension boundary

v0.5 deliberately does not add:

- automatic differentiation framework bindings;
- Riemannian metrics;
- generalized NUTS for discrete parameters;
- dynamic trait-object kernel graphs;
- SMC, particle MCMC or reversible-jump methods.

Those can be layered over the stabilized target, metric, kernel, trace and
checkpoint contracts without changing accepted-state invariants.
