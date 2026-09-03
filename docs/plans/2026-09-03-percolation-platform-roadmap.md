# CMC.rs Percolation Platform Implementation Roadmap

> This roadmap implements the principles in the
> [CMC.rs Percolation Platform Philosophy](../metaphysics/2026-09-03-percolation-platform-philosophy.md).
> The philosophy is normative; this roadmap may change as benchmarks and
> validation evidence accumulate.

**Status:** Approved implementation roadmap  
**Established:** 2026-09-03  
**Target:** CPU-first, Carlo.rs-compatible layered percolation platform

## 1. Current Branch and PR #4

PR #4 is a scientifically useful reference implementation, not the final
production API.

Retain and migrate its validation assets:

- exact enumeration and hand-derived finite-size formulas;
- the independent flood-fill reference;
- one-dimensional and 2x2 closed forms;
- site-bond reduction identities;
- multi-seed z-score validation;
- cross-family topology tests, random-graph tests, and three-dimensional
  critical-region checks.

Do not treat these provisional interfaces as compatibility commitments:

- `PercolationMode`;
- `OccupancyState` backed by public `Vec<bool>` fields;
- `cluster_stats(&CsrLattice, ...)`;
- the current `PercolationMC` constructor, parameter schema, or observable
  names.

Foundation work starts from `dev` in reviewable PRs. Validation migrates as the
new layers become available. After the production static adapter lands, PR #4
is closed as superseded rather than rebased into the new architecture.

The philosophy and roadmap should be reviewed independently from provisional
implementation code.

## 2. Dependency Graph

```text
F0 semantics and benchmark baseline
├── F1 read-only topology capabilities
│   ├── F2 configuration/activity + Bernoulli laws
│   │   └── F3 undirected component analyzer + reusable workspace
│   │       └── F4 observables + StaticPercolationMC
│   │           ├── N1 heterogeneous probabilities
│   │           │   └── N2 network failure + degree observables
│   │           │       └── N3 attacks + robustness curves
│   │           └── Z1 Newman-Ziff parameter scans
│   ├── D1 directed topology
│   │   ├── D2 reachability / weak components / SCC
│   │   └── D3 layered space-time directed percolation
│   ├── E1 queue/counter process workspace
│   │   ├── E2 bootstrap percolation
│   │   └── E3 k-core percolation
│   ├── G1 incremental component process
│   │   ├── G2 Achlioptas processes
│   │   └── G3 invasion percolation
│   ├── W1 periodic displacement/embedding capability
│   │   └── W2 wrapping/winding analyzer
│   └── R1 external bond-configuration analysis
│       ├── R2 q=1 random-cluster reduction
│       └── R3 q=2 Edwards-Sokal / SW exposure
├── C0 Carlo.rs vector estimates + extensible provenance
│   └── production curve output for N3/D3/E2/E3/G2/G3
└── P1 particle candidate-pair adapter
    └── P2 continuum/Boolean percolation
```

The serial foundation is `F1 -> F2 -> F3 -> F4`. Static site, bond, mixed, and
arbitrary undirected graph percolation become production-ready only after F4.
Directed, cascade, growth, continuum, random-cluster, and wrapping work may
proceed as parallel tracks once their shared prerequisites are stable. Any
process that emits a curve or profile as one scientific realization must wait
for C0.

## 3. F0: Freeze Semantics and Measure the Baseline

**PR type:** documentation and benchmark baseline  
**Blocks:** no production code, but establishes review criteria for every
subsequent PR

### Deliverables

1. Add a support matrix with four states:
   - **Supported:** the API and semantics exist;
   - **Validated:** all family-specific production gates pass;
   - **Experimental:** usable but not production-ready;
   - **Not implemented:** no supported claim.
2. Define the scientific quantities used throughout the platform:
   - component size;
   - largest component and giant-component fraction;
   - raw cluster moments;
   - finite-cluster susceptibility and its denominator;
   - boundary crossing;
   - periodic wrapping/winding;
   - degree-conditioned network observables.
3. Record the current PR #4 results as a behavioral baseline.
4. Add a Criterion benchmark entry for the current implementation without
   optimizing it.
5. Record topology, `V`, `E`, probabilities, bytes per vertex/edge,
   samples/second, and steady-state allocations/sample.

### Exit Gate

- Every advertised model maps to an explicit ensemble or process.
- A generic routine returning a value is not evidence that a physical model is
  supported.
- Baseline commands and machine/compiler metadata are reproducible.

## 4. F1: Read-Only Topology Capabilities

**PR 1**  
**Files:**

- `CMC.rs/src/topology/id.rs`
- `CMC.rs/src/topology/view.rs`
- `CMC.rs/src/topology/borrowed_csr.rs`
- `CMC.rs/src/topology/error.rs`

### API

Use fixed dense IDs:

```rust
#[repr(transparent)]
pub struct VertexId(usize);

#[repr(transparent)]
pub struct EdgeId(usize);
```

IDs represent validated `[0, n)` indices. Do not use associated ID types:
union-find, bit masks, CSR arrays, and workspaces all require dense indexing,
and associated IDs would spread conversion constraints through every law and
analyzer.

Define capability-specific read-only traits:

- `GraphView::vertex_count()`;
- `UndirectedGraphView::edge_count()`;
- `UndirectedGraphView::edge_endpoints(edge)`;
- `UndirectedGraphView::incidences(vertex)`;
- separate `DirectedOutGraphView` and `DirectedInGraphView` capabilities.

Do not use a `directed: bool` flag. Directedness changes connectivity semantics,
not only storage.

The CPU hot path uses static dispatch and return-position `impl Iterator` in
traits where suitable. Object safety is not a goal.

`CsrLattice` receives trait implementations only. Existing fields and inherent
methods remain unchanged so current CMC.rs users do not pay a migration cost.

`BorrowedUndirectedCsr<'a>` borrows:

- offsets;
- neighbors;
- physical edge IDs per incidence;
- physical edge endpoints.

Construction validates the representation once. Analysis performs no copy.
An offsets/neighbors-only symmetric CSR cannot safely synthesize stable
physical edge IDs in the presence of parallel edges or self-loops; such input
must either provide the missing arrays or implement a weaker adjacency-only
capability.

### Correctness Gate

- Owned/borrowed parity.
- Empty graph, isolated vertices, disconnected components.
- Parallel edges and self-loops.
- Invalid offsets, endpoints, edge IDs, incidence multiplicity, and asymmetric
  undirected incidence.
- No `unsafe` code.
- Existing CMC.rs tests pass without call-site changes.

### Performance Gate

Benchmark direct field access against generic trait access. A reproducible
regression above 5% blocks review; this is a local review gate, not a noisy CI
hard threshold.

## 5. F2: Configuration, Activity, and Bernoulli Laws

**PR 2**  
**Files:**

- `CMC.rs/src/percolation/activity.rs`
- `CMC.rs/src/percolation/configuration.rs`
- `CMC.rs/src/percolation/probability.rs`
- `CMC.rs/src/percolation/law/bernoulli.rs`

### API

Configuration stores activity only. It does not store mode, probability,
sampling law, topology, analyzer, or observable definitions.

Use small capabilities:

- `VertexActivity`;
- `EdgeActivity`.

`StaticConfiguration` uses an internal replaceable `ActivityMask`. The first
production mask is a dense byte buffer (`Vec<u8>` or a transparent equivalent)
for predictable indexing. Public physics APIs do not expose the layout.
`All(n)` and `None(n)` representations may avoid filling buffers at exact
probability limits.

Probability fields support uniform and heterogeneous data without copying:

```rust
ProbabilityField<'a> = Uniform(f64) | Borrowed(&'a [f64]) | Owned(Vec<f64>)
```

Provide independent laws:

- `SiteBernoulli`;
- `BondBernoulli`;
- `MixedBernoulli`.

Do not add more variants to `PercolationMode`. Laws fill an existing
configuration and neither analyze connectivity nor allocate in the steady
state.

### Correctness Gate

- Exact `p = 0` and `p = 1` behavior.
- NaN, infinity, out-of-range probability, and heterogeneous length errors.
- Fixed-seed bitwise reproducibility.
- Mixed activity at `p_edge = 1` reduces to site activity; at `p_site = 1` it
  reduces to bond activity.

### Performance Gate

- Steady-state allocations/sample: zero.
- Record vertices/s, edges/s, bytes/entity, and `p = 0.01`, `0.5`, `0.99`.
- Compare byte masks with a private bitset prototype before deciding whether a
  second production storage backend is justified.

## 6. F3: Undirected Component Analyzer and Reusable Workspace

**PR 3**  
**Files:**

- `CMC.rs/src/percolation/analyzer/components.rs`
- `CMC.rs/src/percolation/workspace/components.rs`
- `CMC.rs/src/percolation/query/boundary.rs`
- `CMC.rs/src/percolation/summary.rs`
- a test-only or audit-only flood-fill reference

### API

The production analyzer is generic over a borrowed
`G: UndirectedGraphView` and activity capabilities. It owns neither topology nor
configuration.

`ComponentWorkspace` owns reusable temporary storage:

- parent and component size;
- root/generation stamps;
- boundary flags;
- optional component labels;
- query scratch.

After `reserve/prepare`, repeated samples do not allocate.

Scan the physical edge table once in edge order. Do not scan duplicate CSR
incidences and deduplicate them.

`ComponentSummary` contains at least:

- active vertex count;
- active edge count;
- component count;
- largest component identity and size;
- raw second moment.

Labels, higher moments, and boundary queries are explicitly requested through
an analysis plan so callers do not pay for unused outputs. One component pass
may answer several crossing queries.

### Scientific Gate

Migrate the PR #4 validation assets:

- independent flood-fill parity;
- exhaustive small configurations;
- chain, square, cubic, triangular, honeycomb, kagome, and random graphs;
- site, bond, and mixed laws;
- self-loops, parallel edges, disconnected graphs, and overlapping boundary
  sets.

### Performance Gate

- Steady-state allocations/sample: zero.
- Benchmarks: chain, square, 3D cubic, sparse Erdos-Renyi, and power-law
  networks.
- Record edges/s, samples/s, workspace capacity, and peak resident memory.

## 7. F4: Scientific Observables and Carlo.rs Adapter

**PR 4**  
**Files:**

- `CMC.rs/src/percolation/observable.rs`
- `CMC.rs/src/percolation/carlo.rs`

### Observable Definitions

`ObservablePlan` explicitly requests moments, labels, boundary queries,
normalization, and largest-component exclusion policy.

Keep these distinct:

- `LargestSize`;
- `GiantFraction`;
- `RawSecondMoment`;
- `FiniteClusterSusceptibility`.

Finite-cluster susceptibility is computed per sample after removing the chosen
largest component. It cannot be reconstructed from separately aggregated mean
largest size and mean second moment.

Tie behavior for several equal largest components is explicit. The initial
default removes one component, selecting the lowest root ID deterministically.
Alternative conventions require separate names and validation.

### Carlo.rs Integration

`StaticPercolationMC` is an owned convenience adapter composed from:

- `CsrLattice`;
- a concrete static-law enum used only at the adapter boundary;
- static configuration;
- component workspace;
- observable plan.

External borrowed CSR and user-defined laws call the physical core directly or
compose an adapter with `Run::from_parts()`. Do not force borrowed data or RNG-
generic closures through `FromParams`.

Replace the provisional PR #4 public API in one migration. Do not keep a
compatibility wrapper for an API that has never merged.

### Exit Gate

- All PR #4 exact, reduction, scheduler, and z-score tests pass after migration.
- README and VALIDATION define every observable name, unit, denominator,
  exclusion convention, and boundary semantics.
- `thermalization_sweeps = 0` is explicit in adapter docs and examples.
- After merge, ordinary site, bond, mixed, and arbitrary static undirected graph
  percolation are production-ready.
- PR #4 is closed as superseded.

## 8. C0: Carlo.rs Vector Estimates and Extensible Provenance

**Independent parallel PR; required before scientific curve output**

### Current Blocker

`Accumulator::add_array()` stores vector-valued bins, but `finalize()` averages
all elements of each bin into one scalar. `Results` stores only scalar
`Estimate` values. A robustness curve, survival profile, or cascade trajectory
would therefore lose its coordinate structure.

CMC.rs must not work around this with dynamically generated scalar observable
names.

### Deliverables

- A generic `ArrayEstimate { shape, elements: Vec<Estimate> }`, or an equivalent
  element-preserving result representation.
- Shape-preserving JSON, readback, merge, and rebin operations.
- Existing scalar APIs remain valid.
- Shape mismatch becomes an explicit error, not only a `debug_assert`.
- Generic extensible provenance attributes; no percolation-specific fields in
  Carlo.rs.

### Exit Gate

- Every element has independent mean, stderr, and applicable autocorrelation.
- JSON round-trip and old result fixtures pass.
- ResultTools/readback and merge preserve shapes.
- Scalar behavior remains byte/semantic compatible where promised.

## 9. N Track: Heterogeneous and Network Percolation

### N1: Heterogeneous Bernoulli Laws — PR 5

Use borrowed or owned `ProbabilityField` values per vertex and per edge. The
analyzer remains unchanged.

Gates:

- exact non-identically-distributed small-graph enumeration;
- uniform field reduction, bitwise where the RNG draw order is identical and
  statistical otherwise;
- mixed reduction identities;
- length and invalid-probability rejection;
- uniform versus heterogeneous sampling benchmark.

### N2: Random Network Failure and Degree Observables — PR 6

Add a `NetworkProfile` that caches original-graph degree. Random failure is a
domain facade over heterogeneous site occupation, not a duplicate sampler.

Observables:

- giant-component fraction;
- finite-cluster susceptibility;
- active fraction by original degree;
- giant membership by original degree;
- component moments conditioned on original degree.

Do not silently substitute residual degree for original degree.

Gates:

- random failure reduces to the site law;
- path, cycle, star, and complete graph hand calculations;
- graph relabeling invariance;
- stable array shapes and explicit empty-degree-bin semantics.

### N3: Targeted Attacks and Robustness Curves — PR 7, depends on C0

The first production attack orders vertices by original degree. Tie policy is
explicit:

- stable vertex ID; or
- seeded random permutation.

One complete removal trajectory is one realization. Its curve is one vector
measurement, not a collection of independent scalar samples.

Gates:

- each step matches naive component recomputation;
- active count is strictly monotone;
- ties reproduce for fixed seeds;
- curve coordinates and provenance record removal fraction and attack policy.

Adaptive residual-degree attack is a separate later PR because it is a
different process.

### Network Input — small PR parallel with N2

- edge-list input builds an owned CSR;
- complete external CSR arrays create a borrowed zero-copy view;
- parsing reports line/endpoint/duplicate policy errors explicitly.

Do not choose a custom binary or memory-mapped format before million-edge
memory and ingestion benchmarks justify one.

## 10. Z Track: Newman-Ziff Parameter Scans

**Independent algorithm PR after F3/F4**

Randomly permute vertices or edges, activate them incrementally, record the
occupation-count curve, then apply binomial reweighting.

It is not a fast mode of Bernoulli sampling and initially does not support:

- heterogeneous probabilities;
- deletion processes;
- a two-dimensional mixed `(p_site, p_bond)` surface.

Gates:

- every occupation-count step matches naive recomputation;
- reweighted estimates agree with direct Bernoulli estimates within confidence
  intervals;
- benchmarks demonstrate an actual parameter-scan advantage before merge.

## 11. D Track: Directed Percolation

### D1: Directed Topology

Add owned and borrowed directed CSR views with separate outgoing and incoming
capabilities. Validate arc multiplicity, endpoints, and optional transpose
arrays.

### D2: Static Directed Analysis

Keep result types distinct:

- source-to-target reachability;
- weak components;
- strongly connected components.

Validate Tarjan against Kosaraju, weak components against an undirected
flood-fill, directed chains/trees against closed forms, and SCC partitions
against graph transpose.

### D3: Layered Space-Time Directed Percolation — depends on C0

Represent causal layers explicitly. One sweep generates one complete
space-time realization and emits survival/arrival profiles. Do not model each
time layer as an independent measurement.

Validate layered DAG dynamic programming against explicit path enumeration.
Use two frontier buffers or bitsets so memory scales with layer width, not
space-time volume when history is not requested.

## 12. E Track: Bootstrap and k-Core Percolation

### E1: Queue/Counter Process Workspace

Share allocation and queue/counter primitives only. Do not create a common
physical transition type.

### E2: Bootstrap Percolation

Specify initial seeds, homogeneous or heterogeneous activation thresholds, and
whether recorded time means synchronous rounds or queue closure order.

Gates:

- exhaustive initial states on small graphs versus repeated full scans;
- `r = 0` activates all vertices;
- `r > max_degree` preserves only initial seeds;
- paths, cycles, stars, and complete graphs;
- update-order independence of the final closure where mathematically valid.

### E3: k-Core Percolation

Track residual degrees, deletion order, and optional peeling rounds.

Gates:

- exhaustive small-graph initial subgraphs versus repeated full scans;
- `k = 0` removes nothing;
- `k > max_degree` removes everything;
- paths, cycles, stars, and complete graphs;
- deterministic final core independent of queue ordering.

Production curve output from E2/E3 depends on C0.

## 13. G Track: Achlioptas and Invasion Processes

### G1: Incremental Component Process

Provide incremental union operations with online component count, largest size,
and second moment. Record process coordinates without forcing a specific law.

### G2: Achlioptas Processes

Physical parameters include:

- candidate count;
- sampling with or without replacement;
- edge universe;
- product or sum rule;
- tie policy.

A generic `score(size_a, size_b)` callback is insufficient to define the
process.

Gates:

- deterministic replay of a fixed candidate trace;
- every step against naive BFS component recomputation;
- candidate count one reduces to ordinary random edge addition;
- the terminal state contains every permitted edge.

### G3: Invasion Percolation

Freeze site or bond thresholds for each realization, grow from explicit seeds
through a minimum-threshold frontier, and make tie and stop policies explicit.

Gates:

- fixed-weight hand-worked graphs;
- binary-heap implementation versus scanning the complete frontier each step;
- tree uniqueness;
- non-trapping bond invasion agrees with the applicable Prim/MST prefix
  relation when thresholds are unique.

Achlioptas and invasion share infrastructure, not a fabricated common law
trait. Their curve output depends on C0.

## 14. P Track: Continuum / Boolean Percolation

### P1: Particle Candidate-Pair Adapter

Borrow coordinates and `SimulationCell`; reuse `CellList` candidate search and
minimum-image geometry. Do not materialize a persistent CSR for every particle
snapshot.

### P2: Boolean Percolation

Initial scope:

- orthorhombic cells;
- disks/spheres;
- monodisperse and finite-polydisperse radii;
- binomial and Poisson point processes.

For an existing equilibrium particle simulation, connectivity is a measurement
of that chain's accepted coordinates, not a second geometry sampler.

Gates:

- candidate-pair analysis versus all-pairs `O(N^2)` overlap for small `N`;
- translation, particle-permutation, and periodic representative invariance;
- zero-radius and fully-overlapping limits;
- Poisson particle-count mean and variance;
- known continuum thresholds used only as statistical brackets, never labeled
  exact.

Arbitrary cell geometry, arbitrary grains, and multilevel cell lists remain
experimental until measured workloads justify them.

## 15. R Track: FK / Random-Cluster Percolation

### R1: External Bond-Configuration Analysis

Analyze a caller-supplied bond configuration using the shared component
analyzer. Do not couple analysis to the sampler that produced the bonds.

### R2: q = 1

Validate the random-cluster model's q=1 limit against independent bond
Bernoulli percolation configuration-by-configuration and statistically.

### R3: q = 2 Edwards-Sokal / Swendsen-Wang

Expose the actual FK bond draw used by a dedicated Edwards-Sokal or
Swendsen-Wang adapter, then analyze that same configuration.

Gates:

- exact bond-subset enumeration on small graphs;
- partition-function agreement with exact zero-field ferromagnetic Ising;
- `P_FK(i connected to j) = <s_i s_j>`;
- independent consistency with spin observables and worm correlation results.

The current worm's high-temperature even subgraph has `tanh(beta J)` weights
and is not an FK configuration. It must never be relabeled as one.

Generic real `q` samplers (Sweeny or Chayes-Machta) require a separate plan
after q=1 and q=2 are production-ready.

## 16. W Track: Crossing Versus Wrapping/Winding

### W1: Periodic Embedding Capability

Each directed incidence must expose an explicit integer cell displacement or
cocycle. Do not infer it from `BondType` or vertex IDs; that fails for small
periodic cells, parallel edges, and non-Bravais lattices.

### W2: Winding Analyzer

Use a union-find with integer potentials. A cycle whose accumulated potential
is nonzero supplies a winding generator. Results report axes/generators or
winding rank, not only one Boolean.

Gates:

- non-wrapping loop;
- single-axis and multi-axis wrapping;
- independent finite-cover lift-BFS reference;
- invariance under origin translation, representative choice, and edge
  direction reversal;
- wrapping query on an open/non-embedded graph returns an explicit error.

Continuum and FK wrapping may reuse this analysis contract only after their
own geometry and physical configurations are defined.

## 17. Uniform Definition of Done for Every PR

A family or capability PR is complete only when all applicable items pass:

1. Physical ensemble/process, connectivity, and boundary semantics are explicit.
2. Exact, analytic, or rigorous limiting anchors exist.
3. An independent reference implementation checks the production path.
4. Multi-seed statistical gates exist for stochastic outputs.
5. Invalid-domain tests reject unsupported input before simulation.
6. Fixed-seed reproducibility is demonstrated.
7. `cargo fmt --all --check`, workspace clippy with `-D warnings`, relevant
   full tests, and `cargo deny` pass.
8. Benchmarks report throughput, peak capacity, bytes/entity, and steady-state
   allocations; optimization data appears in the same PR.
9. README and VALIDATION update Supported/Validated/Experimental/Not
   implemented claims.
10. Crate-root re-exports expose only stable user entry points.

A correct but unbenchmarked optimized path remains experimental. A fast path
without independent parity evidence does not merge.

## 18. Milestones

### M1: Foundation

F0 through F3. Borrowed topology, activity/law separation, and a reusable
allocation-free undirected analyzer are stable.

### M2: Static Production

F4, N1, W2, and Z1. Ordinary and heterogeneous site/bond/mixed percolation,
scientific observables, crossing/wrapping, and efficient parameter scans are
production-ready.

### M3: Network Production

C0, N2/N3, and network input. Real network failure and targeted-attack curves
have honest vector results and degree-conditioned observables.

### M4: Directed and Cascading Processes

D and E tracks. Directed reachability/space-time processes, bootstrap, and
k-core pass their family-specific gates.

### M5: Growth and Neighboring Domains

G, P, and R tracks. Achlioptas/invasion, continuum Boolean, and q=1/q=2
random-cluster analyses are production-ready within documented boundaries.

### M6: Coverage Audit

Audit the full domain matrix, cross-module benchmarks, nightly statistical
monitoring, API documentation, and provenance. Any family that has not passed
its definition of done remains Experimental or Not implemented; it is not
advertised as production-ready.

## 19. Stage Status

| Stage | Status | Evidence / next gate |
|-------|--------|----------------------|
| Philosophy and roadmap | **Established** | Normative philosophy and this staged roadmap are committed |
| F0 semantics and benchmark baseline | **Complete** | Support matrix, scientific definitions, isolated throughput benchmark, allocation/storage probe, independent review and F0 verification are recorded in `CMC.rs/PERCOLATION.md` |
| F1 read-only topology capabilities | **Next** | Review exact trait signatures against `CsrLattice`, a borrowed CSR fixture, and a directed CSR sketch before implementation |
| F2-F4 static production foundation | **Not started** | Execute only after the preceding foundation stage passes review and knowledge synchronization |
| C0 Carlo.rs vector estimates | **Not started (parallel track)** | May proceed independently; process-curve adapters remain blocked until it lands |
| N/Z/D/E/G/P/R/W tracks | **Not started** | Their dependency gates are defined above |

PR #4 remains open as a validation source until F4 migrates all applicable
gates. It is not a stable API baseline and must not be merged as the final
platform architecture.
