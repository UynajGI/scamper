# v0.4 Review

The uploaded v0.4 workspace had no blocker requiring a repair release. State
atomicity, NUTS candidate selection, warmup freeze, diagnostics, trace
compatibility and checkpoint continuation were suitable foundations for v0.5.

Small production-readiness gaps addressed in v0.5:

- v0.4 used the valid endpoint-displacement U-turn criterion, but did not retain
  summed trajectory momentum or check both cross-subtree joins;
- static HMC and NUTS required a manually chosen initial step size before dual
  averaging;
- differentiable targets had no public finite-difference gradient validator.

v0.5 resolves these without changing the accepted-state or
`TransitionKernel<T>` model. New serialized fields use serde defaults, so v0.4
HMC and NUTS JSON remains readable.
