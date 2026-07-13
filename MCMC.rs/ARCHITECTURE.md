# MCMC.rs Architecture — v0.3

## Core invariants

1. `ChainState.position`, cached `log_density`, and valid gradient cache describe
   the same accepted position after every completed transition.
2. `NEG_INFINITY` denotes outside support. Accepted states must have finite log
   density and finite position; HMC accepted gradients must be finite.
3. Every built-in kernel commits a proposal only after all target, geometry and
   numerical checks succeed. A divergent HMC trajectory is rejected atomically.
4. Adaptation is legal only in `SamplingPhase::Warmup`. Entering sampling with
   incomplete configured HMC warmup is an error.
5. HMC metrics store inverse mass `G = M^-1`; momentum sampling, kinetic energy
   and velocity use one synchronized geometry.
6. Leapfrog integration mutates private `PhasePoint` workspace, never the
   accepted chain state.
7. Windowed metric adaptation observes post-transition accepted positions,
   updates geometry only at slow-window boundaries, resets dual averaging after
   a geometry update, and retunes step size in the terminal buffer.
8. Transform kernels operate in unconstrained coordinates. Differentiable
   transforms provide analytic target-gradient pullback and log-Jacobian
   gradients.
9. Replica-exchange targets, kernels, RNGs and traces remain attached to fixed
   ladder slots. Only accepted states move between slots.
10. Checkpoints serialize RNG, persistent kernel/adaptation state, accepted
    state, gradient cache and trace. Reconstructible HMC trajectory workspaces
    are skipped and rebuilt; targets are reconstructed by the caller and
    fingerprint-checked.

## Hamiltonian target contract

`DifferentiableLogDensity` combines value and gradient evaluation:

```rust,ignore
pub trait DifferentiableLogDensity: LogDensity<[f64]> {
    fn log_density_and_gradient(
        &mut self,
        position: &[f64],
        gradient: &mut [f64],
    ) -> f64;
}
```

The combined interface avoids duplicate model evaluation at every leapfrog
step. `StaticHmc` lazily evaluates the accepted-state gradient once and stores it
in `EuclideanCache`. Rejected trajectories reuse that cache; accepted proposals
commit position, log density and proposed gradient together.

## Metric boundary

`Metric` exposes four operations:

```text
dimension
sample p ~ Normal(0, M)
velocity = G p
kinetic energy = 0.5 pᵀ G p
```

where `G = M^-1`. Built-ins:

- `UnitMetric`: identity geometry;
- `DiagonalMetric`: positive inverse-mass diagonal;
- `DenseMetric`: symmetric positive-definite row-major inverse mass plus cached
  lower Cholesky factor.

For dense `G = L Lᵀ`, momentum sampling solves `Lᵀ p = z`, `z ~ Normal(0, I)`,
which gives `Cov(p) = G^-1`. Dense construction symmetrizes the supplied matrix,
uses increasing diagonal jitter when necessary, and stores the exact matrix
reconstructed from the accepted factor.

## Leapfrog and HMC transaction

One `StaticHmc` transition performs:

```text
validate accepted state
→ obtain/copy accepted gradient
→ sample private momentum
→ integrate private PhasePoint for L leapfrog steps
→ compute ΔH = H_proposed - H_current
→ reject on invalid trajectory or |ΔH| above threshold
→ otherwise Metropolis accept with min(1, exp(-ΔH))
→ atomically commit position/log-density/gradient or mark rejection
→ update warmup controller from the resulting accepted state
```

`TransitionReport` records target/gradient evaluations, completed leapfrog
steps, energy error, divergence, acceptance and the step size used by that
trajectory.

## Warmup architecture

`HmcWarmup` combines:

- `DualAveraging` for step size;
- `WarmupWindowConfig` with initial buffer, expanding slow windows and terminal
  buffer;
- diagonal or dense Welford covariance accumulation;
- serialized iteration/window position for checkpoint continuation.

At a slow-window boundary:

1. covariance is regularized;
2. the matching metric geometry is installed;
3. the accumulator is reset;
4. dual averaging restarts around the current step size.

The terminal buffer contains no metric updates and only retunes the step size.
The final production step size is the dual-averaged value.

## Differentiable transform boundary

`DifferentiableBijector` extends `Bijector` with:

```text
pullback: (d log π / dx) → (d log π / dz)
log_jacobian_gradient: d/dz log |det(dx/dz)|
```

Implemented analytically for:

- identity;
- positive exponential;
- bounded logistic interval;
- ordered vectors;
- simplex stick breaking;
- static binary products.

`TransformedTarget<T, B>` owns constrained value/gradient workspaces and
implements `DifferentiableLogDensity` when both wrapped components support the
required derivative contract.

## Composition and lifecycle

`TransitionKernel<T>` remains the common boundary. `StaticHmc` composes with
`Then`, `Repeat` and `Mixture` on differentiable targets. Phase hooks receive the
target and state; HMC uses them to reject incomplete warmup and freeze the final
step size before production.

The Carlo.rs adapter records scalar HMC diagnostics without changing Carlo's
measurement model. Raw posterior positions and per-draw divergence/energy error
remain owned by `MemoryTrace`.

## Deferred to v0.4

v0.4 should add NUTS on top of the v0.3 metric, gradient cache and leapfrog
contracts:

- bidirectional tree expansion;
- slice or multinomial trajectory sampling;
- generalized U-turn criterion for all built-in metrics;
- maximum tree depth and saturation diagnostics;
- E-BFMI and trajectory-level energy storage;
- optional reasonable-step-size search.

The v0.4 work should not change trace ownership, independent-chain diagnostics
or the `TransitionKernel<T>` target boundary.
