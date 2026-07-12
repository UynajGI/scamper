# MCMC.rs Architecture — v0.2

## Core invariants

1. `ChainState.position`, cached `log_density`, and cache validity describe the
   same accepted state at every completed elementary transition.
2. `NEG_INFINITY` is a valid support result. `NaN` and positive infinity are
   errors.
3. Adaptation is legal only in `SamplingPhase::Warmup`. Entering sampling
   freezes every adaptive parameter.
4. Gibbs/block updaters write into private workspaces. The accepted state is
   committed only after update and density validation succeed.
5. Dense Gaussian geometry stores a row-major lower Cholesky factor. Entries
   above the diagonal are zero and diagonal entries are strictly positive.
6. Transform kernels operate in unconstrained Euclidean coordinates;
   `TransformedTarget` owns constrained workspaces and Jacobian correction.
7. Replica-exchange targets, kernels, RNGs and traces remain attached to fixed
   ladder slots. Only accepted states move between slots.
8. Production traces are row-major contiguous buffers; one trace never mixes
   slot/chain IDs.
9. Checkpoints include RNG state and adaptive kernel state. The target is
   reconstructed by the caller and checked with a stable fingerprint.

## Kernel type boundary

v0.2 changes the central trait from a target-generic method to a target-typed
trait:

```rust,ignore
pub trait TransitionKernel<T>
where
    T: LogDensity<[f64]> + ?Sized,
{
    fn transition<R: rand::Rng + ?Sized>(
        &mut self,
        target: &mut T,
        state: &mut EuclideanState,
        rng: &mut R,
        phase: SamplingPhase,
    ) -> Result<TransitionReport, McmcError>;

    fn name(&self, target: &T) -> &'static str;
}
```

Built-in kernels still implement the trait for every `LogDensity`. The target
parameter permits exact Gibbs and blocked kernels to implement the same trait
for one concrete statistical model and access its conditional-distribution
state.

Phase hooks also receive `&mut T`. This removes target-type inference ambiguity
for universal built-in kernels and lets model-specific kernels prepare or freeze
target-coupled workspaces at phase boundaries.

## Composition semantics

- `Then<A, B>` executes two elementary kernels sequentially.
- `Repeat<K>` executes a positive fixed number of elementary transitions.
- `Mixture<A, B>` chooses exactly one kernel per transition.
- `TransitionReport.subtransitions` records the elementary count.
- Sequential composition is intentionally not transactional across child
  kernels: if the second child fails, the already completed first child remains
  a valid accepted transition.

## Dense adaptation

`DenseCovarianceAdaptation` uses a matrix Welford recurrence. At warmup end it:

1. symmetrizes the sample covariance;
2. adds positive diagonal regularization;
3. computes a lower Cholesky factor;
4. retries with increasing jitter if numerical roundoff prevents factorization;
5. freezes and installs the factor into `RandomWalkMetropolis`.

Diagonal and dense covariance adaptation are mutually exclusive on one kernel.
The global Robbins-Monro multiplier remains independent of proposal geometry.

## Transform layer

`Bijector::forward` returns `log |det(dx/dz)|`; `inverse` returns its inverse
counterpart. Implemented transforms are:

- `Identity`;
- scalar `Positive` and `Interval`;
- vector `Ordered`;
- stick-breaking `Simplex`;
- static binary `Product<A, B>`.

Dynamic trait-object transform graphs are deferred. Static products preserve
inlining and straightforward ownership.

## Replica exchange

Local transitions run with Rayon. At each exchange interval, neighboring fixed
slots use alternating even/odd pairing. For slots `i` and `j`, the runtime
cross-evaluates both states under both targets and applies the exact generic
exchange ratio. No assumption is made about how the ladder value modifies the
target.

Production draws are recorded after local transitions and before any exchange
triggered at that interval boundary. Consequently traces remain fixed-slot
observations, and `final_position` may include a final exchange not represented
as a retained draw. Cross-temperature traces are not fed into ordinary
multi-chain R-hat/ESS, because different ladder slots generally target
different distributions.

## Module ownership

- `target`: density evaluation contracts
- `state`: accepted state and future gradient cache
- `kernel`: built-in, composed and target-specific transitions
- `adaptation`: scale, diagonal and dense covariance warmup estimators
- `proposal`: Gaussian proposal geometry
- `transform`: constrained/unconstrained bijectors and target wrapper
- `trace`: posterior draw storage and views
- `diagnostics`: rank-normalized independent-chain diagnostics
- `multichain`: deterministic independent-chain parallelism
- `tempering`: local fixed-slot replica exchange
- `checkpoint`: generic serde checkpoint envelope
- `carlo_adapter`: Carlo.rs lifecycle integration

## Deferred interfaces

v0.3 should add static HMC, metrics, leapfrog integration, dual averaging and
windowed warmup. v0.4 should add NUTS and HMC-specific diagnostics. These should
build on the v0.2 target, dense geometry and transform boundaries rather than
modify trace or independent-chain diagnostics.
