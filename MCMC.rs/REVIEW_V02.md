# MCMC.rs v0.2 Review

## Decision

v0.2 did not contain an architecture-level blocker, so development proceeded to
v0.3. The v0.2 target-typed kernel boundary, state transaction model, contiguous
trace ownership, dense random-walk geometry, transform value maps and generic
replica-exchange runtime were retained.

## Areas reviewed

- crate/workspace membership and package metadata;
- every v0.2 MCMC source and test module;
- state/log-density atomicity in RW, component, slice and Gibbs kernels;
- `TransitionReport` composition and fixed-layout diagnostics;
- dense Welford covariance and regularized Cholesky construction;
- transform round trips and simplex log-Jacobian formula;
- fixed-slot replica-exchange cross-target acceptance ratio;
- serde defaults for v0.1/v0.1.1 checkpoint compatibility;
- trace dimensions, thinning, chain-ID isolation and HDF5 schema;
- multi-chain deterministic seed derivation and diagnostic inputs.

## Findings

### No blocker

No defect comparable to the v0.1 missing target module or v0.1 component/slice
state corruption was found. The v0.2 public boundaries can support HMC without
reworking `TransitionKernel<T>`, `EuclideanState`, trace ownership or Carlo.rs.

### Small issues/gaps fixed in v0.3

1. `Bijector` documentation described a differentiable mapping, but v0.2 did
   not expose gradient pullback or log-Jacobian derivatives. Consequently a
   `TransformedTarget` could not implement `DifferentiableLogDensity` and could
   not be used directly with HMC. v0.3 adds `DifferentiableBijector` and analytic
   derivatives for all v0.2 built-in transforms.
2. `EuclideanCache` was reserved but lacked an atomic API for synchronizing a
   combined value/gradient evaluation with the accepted state. v0.3 adds
   controlled synchronization and Hamiltonian proposal commit methods.
3. The Carlo adapter recorded target evaluations and divergence but omitted
   gradient evaluations, leapfrog work and energy error already represented by
   `TransitionReport`. v0.3 records these scalar measurements.
4. The v0.2 documentation referred to the gradient cache as future
   infrastructure. v0.3 updates the contract now that the cache is active.

## Existing limitations not treated as v0.2 blockers

- The final v0.3 validation compiled all retained v0.2 code and tests; the
  default-feature MCMC.rs suite completed with 41 passed and 0 failed.
- The optional HDF5 build is blocked in this environment before MCMC.rs source
  compilation because `hdf5-sys 0.8.1` does not recognize the installed HDF5
  1.14.5 header format.
- NUTS, dynamic trajectory building and E-BFMI are intentionally v0.4 scope.
