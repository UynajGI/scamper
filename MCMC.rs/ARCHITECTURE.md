# MCMC.rs Architecture — v0.1

## Invariants

1. `ChainState.position`, cached `log_density`, and cache validity describe the
   same accepted state at every completed transition.
2. `NEG_INFINITY` is a valid support result. `NaN` and positive infinity are
   errors.
3. Hastings convention is `ln q(old|new) - ln q(new|old)`. Built-in Gaussian
   random walks are symmetric and contribute zero.
4. Adaptation is legal only in `SamplingPhase::Warmup`. Entering sampling
   freezes all adaptive parameters.
5. Production traces are row-major contiguous buffers; one trace never mixes
   chain IDs.
6. Multi-chain execution creates independent mutable target and kernel values.
7. Checkpoints include RNG state and adaptive kernel state. The target is
   reconstructed by the caller and checked with a stable fingerprint.

## Module ownership

- `target`: density evaluation contracts
- `state`: accepted state and future gradient cache
- `kernel`: transition algorithms and fixed-layout reports
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
