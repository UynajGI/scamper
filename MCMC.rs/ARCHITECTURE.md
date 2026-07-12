# MCMC.rs Architecture — v0.1.1

## Invariants

1. `ChainState.position`, cached `log_density`, and cache validity describe the
   same accepted state at every completed transition.
2. Component-wise and slice kernels use private `#[serde(default)]` workspace
   buffers. A transition either fully succeeds (atomic `swap_position`) or leaves
   every byte of `ChainState` unchanged.
3. Every built-in kernel advances `ChainState::iteration` exactly once per
   `transition()` call.
4. `NEG_INFINITY` is a valid support result. `NaN` and positive infinity are
   errors.
5. Hastings convention is `ln q(old|new) - ln q(new|old)`. Built-in Gaussian
   random walks are symmetric and contribute zero.
6. Adaptation is legal only in `SamplingPhase::Warmup`. Entering sampling
   freezes all adaptive parameters.
7. Production traces are row-major contiguous buffers; one trace never mixes
   chain IDs.
8. Multi-chain execution creates independent mutable target and kernel values.
9. Checkpoints include RNG state and adaptive kernel state. The target is
   reconstructed by the caller and checked with a stable fingerprint.
10. `EuclideanState::validate()` is called by `MemoryTrace::record()` and
    `ChainCheckpoint::validate_format()` as a consistency gate.

## Module ownership

- `target.rs`: density evaluation contracts (single file, was `target/` directory in v0.1)
- `state`: accepted state, gradient cache, and `EuclideanState::validate()`
- `kernel`: transition algorithms with private workspace and fixed-layout reports
- `adaptation`: reusable warmup estimators
- `trace`: posterior draw storage and views
- `diagnostics`: rank-normalized cross-chain diagnostics
- `multichain`: deterministic independent-chain parallelism
- `checkpoint`: generic serde checkpoint envelope
- `carlo_adapter`: Carlo.rs lifecycle integration

## Carlo.rs integration

`RunPhase::Thermalization` maps to `SamplingPhase::Warmup` and
`RunPhase::Measurement` maps to `SamplingPhase::Sampling`. The small additions
`Run::from_parts` and `Run::finalize_with_mc` permit closure-based targets and
return the sampler so its trace can be recovered after finalization.

## Deferred interfaces

HMC/NUTS should add metric and integrator modules, not modify `LogDensity`,
`MemoryTrace`, or multi-chain diagnostics. Transforms should wrap targets and
include Jacobian terms before kernels see the unconstrained density.
