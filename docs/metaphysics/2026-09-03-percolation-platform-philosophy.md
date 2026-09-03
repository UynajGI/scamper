# CMC.rs Percolation Platform Philosophy

> This document defines the ontology, success criteria, methodology, and
> boundaries of the CMC.rs percolation domain. It is the decision standard for
> future architecture, implementation, review, and validation work.

## I. Ontology: A Layered Percolation Platform

CMC.rs percolation is a layered scientific-computing platform, not one model,
one algorithm, or one universal trait.

Its shared core consists of stable concepts that recur across percolation
families:

1. **Topology or neighborhood access** describes what may connect. For
   discrete systems this is a read-only graph view; for continuous systems it
   is a spatial neighborhood defined by particle geometry and search
   structures.
2. **Configuration** records the current occupied, removed, active, or
   connected entities. Its storage is independent of the rule that produced
   it.
3. **Law or process** defines how configurations are generated or evolved:
   independent Bernoulli sampling, heterogeneous probabilities, ordered
   attack/removal, threshold cascades, growth, directed spreading, or
   correlated-cluster sampling.
4. **Connectivity semantics** define what counts as connected: undirected
   components, directed reachability, strong/weak components, geometric
   overlap, winding classes, or model-specific correlated clusters.
5. **Analysis and observables** turn a configuration into scientifically
   meaningful quantities: component moments, giant-component fraction,
   susceptibility excluding the giant component, boundary crossing, wrapping,
   degree-resolved statistics, robustness curves, or process-time records.
6. **Runtime adaptation** connects a complete experiment to Carlo.rs through
   `MonteCarlo`, `FromParams` where appropriate, `Context`, and `Results`.
   Carlo.rs owns scheduling, RNG contexts, accumulation, metadata, and output;
   CMC.rs owns the physical semantics.

These concepts are related by composition, not inheritance. Static undirected
Bernoulli percolation may use union-find; directed percolation may use graph
reachability; bootstrap and k-core percolation may use decremental queues;
continuum percolation may use particle cell lists. They belong to one platform
because they share scientific language and selected infrastructure, not
because they execute through one algorithm.

### Domain ownership

- `percolation/` owns discrete-graph static, network, directed, and evolving
  percolation laws and their analyses.
- `lattice/` owns concrete lattice topology and the existing `CsrLattice`.
- `particle/` owns continuous geometry, periodic cells, coordinates, and
  cell-list neighborhood search; continuum percolation reuses these rather
  than duplicating them.
- `worm/` and `lattice/` own correlated FK/Ising representations and update
  rules; percolation exposes compatible cluster analysis without reimplementing
  the physical sampler.
- Carlo.rs remains the runtime boundary. No percolation design may require
  Carlo.rs to understand percolation-specific states or algorithms.

## II. Teleology: Scientific Correctness First

The first success criterion is scientific correctness. Breadth, throughput,
and API convenience matter only after the physical semantics are explicit and
validated.

A percolation family is production-ready only when all applicable gates pass:

1. **Semantic specification:** the occupied entities, transition/generation
   law, connectivity relation, boundary conditions, and measured ensemble are
   stated without relying on implementation names.
2. **Exact or analytic anchors:** finite-size enumeration, closed-form limits,
   duality where rigorous, reduction identities, conservation laws, or known
   thresholds are tested where available.
3. **Independent reference:** optimized paths are checked against a second
   implementation that does not share the same algorithmic failure mode.
4. **Statistical validation:** stochastic outputs pass multi-seed z-score or
   confidence-interval gates, scaled by `SCUTTLE_ZSCORE_SEEDS` for nightly
   monitoring.
5. **Invalid-domain rejection:** unsupported topology, parameter combinations,
   or undefined observables return explicit errors rather than plausible but
   meaningless values.
6. **Reproducibility:** fixed seeds and domain-separated streams produce stable
   results independent of scheduling order where the Carlo.rs contract
   promises determinism.
7. **Performance evidence:** asymptotic complexity, peak memory, allocation
   behavior, and representative CPU throughput are measured for the intended
   scale. Optimizations must retain a reference/audit path.

### CPU success criteria

The current platform is CPU-first. A high-performance CPU implementation:

- accepts borrowed topology views and does not require copying an external CSR
  merely to run an analysis;
- uses static dispatch in hot loops unless measurement proves otherwise;
- reuses explicit workspaces and scratch buffers across samples;
- supports compact occupancy storage when scale justifies it, without exposing
  storage layout as physical semantics;
- separates one-sample analysis from parameter scans and Carlo.rs scheduling;
- reports the actual algorithm and topology semantics used in result metadata;
- establishes benchmarks before introducing parallelism or specialized fast
  paths.

GPU support is not a current success criterion and must not shape the CPU API
prematurely.

### Compatibility

Compatibility means compatibility with the established Carlo.rs framework:
`MonteCarlo` lifecycle, RNG through `Context`, measurement accumulation, and
`Results`. It does not mean preserving the provisional `PercolationMC`,
`PercolationMode`, `OccupancyState`, parameter names, or observable schema on
the unmerged feature branch. Those may be replaced to prevent long-term debt.

## III. Methodology: Stable Layers, Specialized Algorithms

### 1. Abstract the stable read-only topology boundary

Discrete algorithms target a minimal, zero-cost, read-only graph-view trait
whose semantics include vertex count, physical-edge access, and incidence or
out-neighbor traversal. Algorithms do not own the graph.

`CsrLattice` implements this view. Future directed CSR, borrowed external CSR,
memory-mapped graph data, and generated topology views may implement their
appropriate contracts without conversion into `CsrLattice`.

Directedness is a semantic distinction, not a boolean flag hidden in an
undirected graph abstraction. Shared traits may expose common facts, while
separate capabilities express undirected incidence and directed in/out
adjacency.

### 2. Separate configuration, law, analysis, and workspace

- A configuration contains data, not sampling policy.
- A law or process generates or evolves configurations, not observables.
- An analyzer computes observables from a topology and configuration.
- A workspace owns reusable temporary memory and makes allocation behavior
  explicit.
- A Carlo.rs adapter composes these layers into a scheduled experiment; it
  does not define their core APIs.

This separation is mandatory for heterogeneous probabilities, network attacks,
directed reachability, cascading processes, compact occupancy formats, and
alternative cluster algorithms.

### 3. Prefer capability-specific contracts over a universal model trait

Use small contracts for stable capabilities, such as:

- read-only undirected/directed graph access;
- vertex or edge activity queries;
- one-sample generators;
- evolving processes with explicit event/step semantics;
- component, reachability, boundary, and wrapping analyses;
- geometry-based neighbor predicates.

Do not force independent sampling, irreversible dynamics, directed spreading,
and correlated FK clusters into one `PercolationModel` trait. Shared code must
follow shared semantics, not shared naming.

### 4. Keep data ownership and storage choices replaceable

Public physics APIs must not require `Vec<bool>`, owned `Vec`-backed CSR, or a
specific union-find layout. Dense booleans, bitsets, sparse active sets, and
borrowed data are storage strategies selected by scale and workload.

Concrete fast paths may exist behind generic contracts, but every specialized
path requires an equivalent reference path and parity tests.

### 5. Treat observables as first-class scientific definitions

Raw cluster moments and derived physical observables are distinct:

- expose `sum(s_i^k)` with explicit inclusion/exclusion policy;
- compute susceptibility with its stated normalization and giant-component
  exclusion per sample, not from insufficient aggregate means;
- distinguish largest-cluster size from giant-component fraction;
- distinguish boundary crossing from periodic wrapping;
- support multiple simultaneous boundary/wrapping queries without rerunning
  connectivity unnecessarily;
- define degree-resolved and attack-curve observables at the network layer,
  where their conditioning is known.

Observable names, units, normalization, and finite-size conventions belong in
validation documentation and result metadata.

### 6. Build validation before optimization

For each family:

1. specify the finite-state or analytic reference;
2. implement the simplest correct reference algorithm;
3. implement the production path with reusable workspaces;
4. prove parity on exhaustive small cases and adversarial topologies;
5. add multi-seed statistical gates;
6. benchmark and optimize only measured bottlenecks;
7. retain the reference/audit path.

Reduction identities between families are especially valuable: mixed
site-bond must reduce to pure site/bond limits; network random failure must
reduce to heterogeneous independent site occupation; continuum algorithms
must agree with brute-force overlap checks at small N.

### 7. Let neighboring domains keep their expertise

Continuum percolation uses `particle::SimulationCell`, particle coordinates,
and cell lists. FK/Ising correlated clusters use lattice/worm state and update
rules. The percolation platform may define shared analyzers or adapters, but it
must not copy these representations merely to make the module tree look
uniform.

## IV. Boundaries: What This Platform Is Not

1. **Not one universal percolation trait.** Similar output names do not imply
   identical state, evolution, or connectivity semantics.
2. **Not an extension of `CsrLattice` into every graph type.** `CsrLattice`
   remains an owned undirected multigraph; directed and external graphs retain
   their own representations behind view contracts.
3. **Not union-find as the domain model.** Union-find is a static-undirected
   component algorithm, not the definition of directed, dynamic, wrapping,
   continuum, or correlated percolation.
4. **Not a GPU-driven design.** The current target is a correct, high-performance
   CPU platform compatible with Carlo.rs. GPU concerns are explicitly deferred
   to Carlo.rs backend architecture work.
5. **Not a duplicate particle or worm framework.** Continuous geometry stays
   in `particle/`; FK sampling stays in `worm/`/`lattice/`.
6. **Not a promise to preserve unmerged provisional APIs.** The current feature
   branch may be rewritten before merge if its types bind physical semantics to
   temporary storage or adapter choices.
7. **Not silent extrapolation.** A model is not advertised as supported merely
   because a generic graph routine returns a number. Directed, temporal,
   weighted, correlated, and finite-size semantics require explicit contracts
   and validation.
8. **Not immediate coverage of every named family.** Families enter in
   dependency order after their shared foundation is sound. Documentation must
   keep an honest supported/validated/not-implemented matrix.

## V. Decision Criteria

When reviewing a percolation design, ask in order:

1. What physical ensemble or process is this, exactly?
2. Is its connectivity/reachability definition explicit and distinct from its
   storage or algorithm?
3. Does it reuse the correct existing domain boundary (`lattice`, `particle`,
   `worm`, or Carlo.rs) without copying responsibilities?
4. Is the shared abstraction based on a stable semantic fact, or only on two
   current implementations looking similar?
5. Can topology and configuration be borrowed without conversion or ownership
   transfer?
6. Can hot-loop memory be reused through an explicit workspace?
7. Is static dispatch available on the CPU production path?
8. Are normalization, giant-cluster policy, boundaries, and wrapping semantics
   visible in observable definitions?
9. What exact, independent, statistical, and reduction tests would falsify an
   incorrect implementation?
10. Does unsupported input fail loudly before simulation?
11. Does the Carlo.rs adapter compose an existing physical core, rather than
    forcing the core to mirror scheduler or `Params` concerns?
12. Is the new API still appropriate for a million-edge borrowed CSR, a
    directed graph, and the existing lattice case without making all three pay
    for capabilities they do not use?
13. Is complexity justified by a current family and measured need? If not,
    defer it.

A proposal that fails the first nine questions is not ready for implementation.
A proposal that passes correctness but fails performance questions remains a
reference implementation, not the production path.

## Implementation Roadmap

Execution is staged in the
[CMC.rs Percolation Platform Implementation Roadmap](../plans/2026-09-03-percolation-platform-roadmap.md).
The roadmap is subordinate to this philosophy: implementation phases may be
reordered as evidence changes, but the semantic and validation criteria above
remain the decision standard.

---

*This philosophy was clarified through Socratic dialogue and established on
2026-09-03.*
