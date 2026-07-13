# v0.3 Review

## Decision

The uploaded v0.3 workspace was suitable for direct development of v0.4.
No architecture-level defect was found in the static-HMC, metric, warmup,
transform, checkpoint or trace boundaries.

The following were treated as small v0.4-readiness gaps rather than blockers:

- `Metric` lacked a direct displacement–velocity inner-product primitive, so a
  naïve NUTS U-turn implementation would allocate a temporary velocity vector
  at every internal tree node.
- `TransitionReport` did not expose the trajectory acceptance statistic,
  initial Hamiltonian energy or a maximum-tree-depth flag.
- `MemoryTrace` did not retain energy or tree-depth diagnostics.
- multi-chain diagnostics did not compute E-BFMI.
- the leapfrog integrator accepted only positive step sizes, while bidirectional
  NUTS needs signed integration steps.

## Review coverage

The review examined:

- accepted state and gradient-cache synchronization;
- static-HMC Metropolis correction and divergence rejection;
- leapfrog reversibility and metric conventions;
- dual averaging and windowed adaptation freeze semantics;
- transformed-target gradient pullback;
- JSON checkpoint continuation;
- Carlo.rs measurement integration;
- trace backward compatibility;
- public API extension points for a dynamic trajectory kernel.

## Resolution in v0.4

v0.4 adds a source-compatible default metric primitive, signed non-zero
leapfrog step sizes, NUTS-specific report fields with serde defaults, compatible
trace columns, E-BFMI, and the `Nuts<M>` kernel. Existing v0.3 JSON report and
trace payloads remain readable.
