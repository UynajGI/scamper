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

Foundation work is delivered as successive independently reviewed commits
stacked on the `feat/cmc-percolation` branch (PR #4). The provisional API that
PR #4 originally introduced is deleted in-branch before merge — superseded in
place under the never-merged-API rule — so the merged diff contains only the
final architecture; PR #4 itself is updated and serves as the delivery vehicle
rather than being closed.

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
- PR #4 carries the final architecture; its provisional API is deleted in-branch
  before merge.

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

### N1: Heterogeneous Bernoulli Laws — delivered in-branch on PR #4

Use borrowed or owned `ProbabilityField` values per vertex and per edge. The
analyzer remains unchanged.

Gates:

- exact non-identically-distributed small-graph enumeration;
- uniform field reduction, bitwise where the RNG draw order is identical and
  statistical otherwise;
- mixed reduction identities;
- length and invalid-probability rejection;
- uniform versus heterogeneous sampling benchmark.

### N2: Random Network Failure and Degree Observables — delivered in-branch, awaiting independent review

`NetworkProfile` caches original-graph degree once per graph. Random failure
(`RandomFailureMC`) is a domain facade that delegates sampling to the F2
uniform `SiteBernoulli` law — not a duplicate sampler — and records six C0
element-preserving degree arrays plus the F4 scalar plan. Literature anchors
(Callaway et al. PRL 85, 5468; Cohen et al. PRL 85, 4626) are archived in
`docs/ref/` with their comparisons recorded in `CMC.rs/PERCOLATION.md`.

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

### N3: Targeted Attacks and Robustness Curves — delivered in-branch, awaiting independent review

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

### D1: Directed Topology — delivered in-branch, awaiting independent review

Add owned and borrowed directed CSR views with separate outgoing and incoming
capabilities. Validate arc multiplicity, endpoints, and optional transpose
arrays.

### D2: Static Directed Analysis — delivered in-branch, awaiting independent review

Keep result types distinct:

- source-to-target reachability;
- weak components;
- strongly connected components.

Validate Tarjan against Kosaraju, weak components against an undirected
flood-fill, directed chains/trees against closed forms, and SCC partitions
against graph transpose.

### D3: Layered Space-Time Directed Percolation — delivered in-branch, awaiting independent review

Represent causal layers explicitly. One sweep generates one complete
space-time realization and emits survival/arrival profiles. Do not model each
time layer as an independent measurement.

Validate layered DAG dynamic programming against explicit path enumeration.
Use two frontier buffers or bitsets so memory scales with layer width, not
space-time volume when history is not requested.

## 12. E Track: Bootstrap and k-Core Percolation

### E1: Queue/Counter Process Workspace — delivered, reviewed

Share allocation and queue/counter primitives only. Do not create a common
physical transition type.

### E2: Bootstrap Percolation — delivered in-branch, awaiting independent review

Specify initial seeds, homogeneous or heterogeneous activation thresholds, and
whether recorded time means synchronous rounds or queue closure order.

Gates:

- exhaustive initial states on small graphs versus repeated full scans;
- `r = 0` activates all vertices;
- `r > max_degree` preserves only initial seeds;
- paths, cycles, stars, and complete graphs;
- update-order independence of the final closure where mathematically valid.

### E3: k-Core Percolation — delivered in-branch, awaiting independent review

Track residual degrees, deletion order, and optional peeling rounds.

Gates:

- exhaustive small-graph initial subgraphs versus repeated full scans;
- `k = 0` removes nothing;
- `k > max_degree` removes everything;
- paths, cycles, stars, and complete graphs;
- deterministic final core independent of queue ordering.

Production curve output from E2/E3 depends on C0.

## 13. G Track: Achlioptas and Invasion Processes

### G1: Incremental Component Process — delivered in-branch, awaiting independent review

Provide incremental union operations with online component count, largest size,
and second moment. Record process coordinates without forcing a specific law.

### G2: Achlioptas Processes — delivered in-branch, awaiting independent review

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

### G3: Invasion Percolation — delivered in-branch, awaiting independent review

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

### P1: Particle Candidate-Pair Adapter — delivered in-branch, awaiting independent review

Borrow coordinates and `SimulationCell`; reuse `CellList` candidate search and
minimum-image geometry. Do not materialize a persistent CSR for every particle
snapshot. Delivered as `percolation/continuum.rs`: a borrowed flat-coordinate
`CandidatePairGraph` adapter plus zero-copy and prepared-list entry points,
emitting boundary-inclusive candidate pairs in a documented cell-major
canonical order into an ephemeral `PairGraph` (transient per-analysis
`UndirectedGraphView`, never a persistent per-snapshot CSR) that feeds the
F3 analyzer. Gates and measured numbers in `CMC.rs/PERCOLATION.md` (P1).

### P2: Boolean Percolation — delivered in-branch, awaiting independent review

Delivered as `percolation/boolean.rs`: `Radii` (monodisperse or
finite-polydisperse, typed validation), a `PointProcess` sampler (binomial
fixed count; Poisson with a documented chunked-Knuth draw contract), the
equilibrium **measurement** path over caller-supplied coordinates
(`overlap_pairs` / `overlap_pairs_from_configuration` / frozen
`BooleanPercolationMC::measuring` — no resampling), and the
`BooleanPercolationMC` Carlo.rs adapter recording the F4 scalar plan per
realization. Overlap pairs compose P1 at the maximum contact cutoff
`2 r_max` with an in-place per-pair filter `d^2 <= (r_i + r_j)^2`
(monodisperse radii pass through unfiltered — bit-identical criterion).
All gates below green; the 2D disk threshold gate brackets
`eta_c = 1.12808737(6)` (Mertens & Moore 2012, archived) as a statistical
bracket, never exact. Evidence in `CMC.rs/PERCOLATION.md` (P2).

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

### R1: External Bond-Configuration Analysis — delivered in-branch, awaiting independent review

Analyze a caller-supplied bond configuration using the shared component
analyzer. Do not couple analysis to the sampler that produced the bonds.

### R2: q = 1 — delivered in-branch, awaiting independent review

Validate the random-cluster model's q=1 limit against independent bond
Bernoulli percolation configuration-by-configuration and statistically.
Delivered as `percolation/random_cluster.rs`: the `RandomClusterMC`
adapter with an explicit `ClusterQ` pinned to `q = 1` (any other q is a
typed rejection citing this roadmap's separate-plan requirement for
generic real-q samplers), sampling by delegation to the F2 `BondBernoulli`
law — at q = 1 the FK measure `p^{|A|} (1-p)^{E-|A|} q^{k(A)}/Z` collapses
exactly to the product Bernoulli measure because `1^{k(A)} = 1` and
`Z(p,1) = 1` — and analysis through the R1 external-bond path, exposing
each realization as a `BondSample`. Gates in `CMC.rs/PERCOLATION.md` (R2).

### R3: q = 2 Edwards-Sokal / Swendsen-Wang — delivered in-branch, awaiting independent review

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
10. Crate-root re-exports expose only independently reviewed public entry points
    with explicit maturity; crate-private sketches and unreviewed types remain
    unexported.

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
| F1 read-only topology capabilities | **Complete** | Dense IDs, generic undirected views, `CsrLattice` integration, validated zero-copy borrowed CSR, independent review, correctness tests, allocation probe, and repeated <5% trait-overhead benchmark gate complete |
| F2 configuration/activity + Bernoulli laws | **Complete** | Private reusable activity storage, `StaticConfiguration`, validated uniform/heterogeneous probability fields, independent site/bond/mixed laws, atomic validation, six-domain aggregate statistical gates over 1-4096 seeds, 14-ID Criterion data, and 25 zero-allocation probe records complete |
| F3 undirected analyzer + reusable workspace | **Complete** | Mode-free generic component analysis, physical-edge scan, canonical identity, multi-query crossing, independent flood-fill parity, typed errors, workspace reuse/stamp-wrap tests, 15 zero-allocation records, measured workspace/RSS, and 48-ID equivalent-work Criterion evidence passed independent review |
| F4 scientific observables + production adapter | **Complete** | `ObservablePlan` observables with per-sample canonical-largest exclusion and undefined-sample indicator, owned `StaticPercolationMC` adapter, provisional API deleted in-branch with every validation asset migrated, zero-allocation adapter probe, 63-ID Criterion evidence, and independent review passed; the static family is Validated (production on the `dev` merge). Definitions and data in `CMC.rs/PERCOLATION.md` |
| C0 Carlo.rs vector estimates | **Complete** | All §8 gates green on `feat/cmc-percolation`: `ArrayEstimate { shape, elements }` with per-element mean/stderr/variance-ratio autocorrelation; `finalize_array`/`finalize_arrays`/`finalize_split` with legacy scalar collapse preserved; `array_observables` JSON with `Results::from_json` readback accepting legacy fixtures; `Results::merge` element-wise with the scalar weighting and shape-mismatch skip (warn, never flatten); shape mismatches are explicit `CarloError::ShapeMismatch` at accumulator/measurements/context level in every profile (silent zip-truncation eliminated; all in-workspace callers updated); generic `Metadata::attributes` provenance (BTreeMap<String, Value>, stored/merged/round-tripped, never interpreted); 18 new tests incl. scalar-accumulator parity, end-to-end Run curve, and HDF5 checkpoint round-trip; `benches/array_finalize.rs` Criterion data. Evidence in `Carlo.rs/VALIDATION.md`; Experimental until review passes. Process-curve adapters N3/D3/E2/E3/G2/G3 are unblocked; Independent review passed 2026-10-22 ([report](../../reviews/2026-10-22-stage-reviews.md)): the reviewer hand-derived the per-element statistics (bit-exact vs scalar accumulators), parsed a hand-written legacy JSON, and proved the old release build silently zip-truncated; one major follow-up applied — merged stderr now uses the pooled-bin variance (the pre-existing formula underestimated merged error bars by sqrt(n_runs)) and per-element tau now matches the scalar completed-bins convention, both with regression tests |
| N1 heterogeneous Bernoulli laws | **Complete** | All §9 N1 gates pass: owned heterogeneous `StaticLaw` variants on `StaticPercolationMC` with typed adapter-level length validation (analyzer unchanged); exact non-iid small-graph enumeration vs scheduler (\|z\| < 4, reference independently reproduced), bitwise uniform and mixed-pure reductions, endpoint determinism, length/invalid-probability rejection, 16-seed scheduler z-scores for site/bond/mixed (enumeration + Poisson-binomial references, verified at 64 seeds), 16-ID uniform-vs-heterogeneous Criterion comparison, and 15 zero-allocation heterogeneous adapter probe records. Evidence in `CMC.rs/PERCOLATION.md` and `CMC.rs/VALIDATION.md`; `FromParams` stays uniform-only by design |
| N2 random failure + degree observables | **Complete** | All §9 N2 gates green: `NetworkProfile` original-degree cache; `RandomFailureMC` facade delegating to the F2 site law (bitwise reduction to `StaticLaw::site` verified); F4 scalar plan plus six C0 element-preserving degree arrays with exact-zero empty bins, occupancy markers, and a `V^3 <= 2^53` exactness gate; independent-BFS exhaustive-enumeration z-scores on path/cycle/star/complete at keep 0.3/0.5/0.7 (16 seeds, \|z\| < 4, \|z̄\| < 2); statistical relabeling invariance; 0-allocation probe on five topologies; Criterion `percolation_network` data; literature anchors — Cohen et al. `q_c = 1/(kappa-1)` bracket (empirical 0.3996 vs Poisson ideal 0.4) and Callaway et al. Eq. (5) giant fraction (\|ΔS\| <= 0.0005 at q = 0.5/0.65/0.8), PDFs archived in `docs/ref/`. Evidence in `CMC.rs/PERCOLATION.md`; Experimental until review passes. N3 (attacks, now unblocked by C0) builds on this profile/facade substrate; Independent review passed 2026-10-22 ([report](../../reviews/2026-10-22-stage-reviews.md)): the reviewer independently re-derived kappa/q_c and the Callaway fixed points, hand-enumerated the star-graph degree observables to exact fractions (including the tie-sensitive 47/128), and recomputed the V^3 gate boundary; follow-ups applied (bracket transcription, 208063/208064 boundary test) |
| N3 targeted attacks + robustness curves | **Complete** | All §9 N3 gates green: `DegreeAttackMC` with explicit tie policies (stable-ID deterministic; seeded shuffle bitwise-reproducible, Fisher-Yates over tie groups); one trajectory recorded as one C0 vector realization (`RobustnessCurve`, defined-indicator susceptibility, constant coordinate array); Full/Every/At grids; per-step independent-BFS exact parity and strict survivor monotonicity on star/path/ER for both policies; coordinates and policy provenance through `Results` attributes with JSON round-trip; literature anchors — Callaway et al. Eq. (14) degree-cutoff curves on ER and power-law configuration graphs (exact, empirical `p_k`) and the Cohen et al. breakdown bracket (`f_c(power-law) < 0.5 f_c(ER)`; survives `0.3 f_c`, collapsed at `1.5 f_c`), PDFs archived in `docs/ref/`; 0-allocation trajectory probe on five topologies; Criterion data. Adaptive residual-degree attack remains a separate future PR by design. Evidence in `CMC.rs/PERCOLATION.md`; Experimental until review passes |

| Network edge-list input | **Complete** | `csr_from_edge_list` builds owned CSR from text edge lists with an explicit `EdgeListSpec` (vertex count, duplicate Reject/AllowParallel, self-loop Reject/Allow); every rejection typed with the 1-based line (parse/endpoint/duplicate/self-loop). Zero-copy path is the reviewed F1 `BorrowedUndirectedCsr`. Gates: parsed-vs-builder structural parity, per-variant line-number rejections, policy semantics (parallel edges distinct; self-loop = two incidences), parsed graphs through the N2 facade and F1 borrowed view with identical analysis; ingestion benchmark 5.4 M edges/s at 16k and 1.23 M edges/s at 1M (text stays sufficient; no binary format introduced). Evidence in `CMC.rs/PERCOLATION.md`. Independent review passed 2026-10-22 ([report](../../reviews/2026-10-22-network-input-full-report.md)): throwaway-crate parser attack and empirical CSR-semantics checks; follow-ups applied (typed `InvalidVertexCount` replacing a zero-vertex panic, duplicate detection extended to repeated self-loops, doc alignment) |

| D1 directed topology | **Complete** | Public `DirectedGraphView` + separate out/in capabilities with dense `ArcId`s; owned `DirectedCsr` (public fields, `validate`, parallel arcs and self-loops first-class) and zero-copy `BorrowedDirectedOutCsr`/`BorrowedDirectedCsr` (transpose required for in-adjacency). Gates: owned/borrowed parity both directions; empty/isolated/parallel/self-loop cases; 12 typed rejections covering offset shape, lengths, ranges, per-arc sense on both sides, and exact arc multiplicity (missing/duplicated); `validate` catches mutated fields; F1-style trait-overhead bench — adjacency scans within the F1 band (generic 14-16% faster than direct), arc-table ratio +17% attributed to direct-fold auto-vectorization (generic per-element below the recorded F1 generic baseline). Evidence in `CMC.rs/PERCOLATION.md`. Independent review passed 2026-10-22 ([report](../../reviews/2026-10-22-d1-full-report.md)): builder hand-traced on a parallel-arc + self-loop multigraph, all validator attacks predicted then confirmed empirically in debug and release, trait separation compile-checked, SIMD explanation verified from struct-layout evidence; follow-ups applied (checked arithmetic for a magic-value declared count that could overflow a validation error path, stale F1 sentence, per-element figure corrected to the recorded 0.72 ns) |

| D2 static directed analysis | **Complete** | Three distinct result types: `analyze_weak_components` (union-find over the arc table), `analyze_strong_components` (iterative Tarjan, canonical lowest-ID labels, `acyclic` distinguishes self-loop cycles), `reachable` (stamped early-exit BFS). Gates: exhaustive enumeration of ALL digraphs on 3 (2^9 = 512) and 4 (2^16 = 65,536) vertices, self-loops included against transitive-closure ground truth for SCC + weak + reachability; independently written Kosaraju parity (exhaustive + 40 random multigraphs); weak vs undirected flood-fill; SCC transpose invariance; closed forms for chain/tree/star/ring/DAG/self-loop; 0-allocation probe on five topologies (weak+strong+reachability per round); Criterion data (grid and random digraph at V=16384). Tarjan 1972 pre-dates arXiv (no PDF archivable); the definitional exhaustive anchor is strictly stronger than published-number reproduction. Evidence in `CMC.rs/PERCOLATION.md`. Independent review passed 2026-10-22 ([report](../../reviews/2026-10-22-d2-full-report.md)): frame-by-frame Tarjan traces, a third independent SCC implementation on 300 fresh multigraphs through reused workspaces, >2^32-query stamp rollover, and 1M-vertex depth; follow-ups applied (exhaustive counts corrected to 512/65,536 with self-loops, probe-target clippy, u32 bound guard, typed queue growth) |

| D3 layered space-time directed percolation | **Complete** | `LayeredDirectedPercolationMC`: one sweep = one complete T-layer realization through a double frontier bitset (memory O(L/8) words, not space-time volume), recorded as ONE C0 vector sample (`DpLayerActive`, `DpSurvival`, `DpExtinctionTime`) — layers are never independent measurements. Directed bond DP with a documented bitwise draw contract; single-seed and full-occupation initial conditions. Gates: frontier DP equals explicit-DAG path enumeration via the D2 reachability engine (oracle replays the draw contract); one-realization-one-vector; Hinrichsen (Adv. Phys. 49, 815 (2000), archived PDF) anchors — p_c = 0.6447001 with P(t) ~ t^-0.159464 for single-seed survival AND full-seed density decay (two independent anchors, same exponent); sub/super-critical brackets with an honest transient-plateau analysis; 0-allocation realization probe on five topologies; Criterion data. "Arrival profiles" are realized as the per-layer activity/density profile plus extinction time (no separate first-arrival-per-site observable — a defensible reading, recorded here explicitly). Evidence in `CMC.rs/PERCOLATION.md`. Independent review passed 2026-10-22 ([report](../../reviews/2026-10-22-d3-full-report.md)): from-scratch xoshiro256++ model reproduced production traces bit-exactly (draw-per-incidence count, L=70 trailing-bit mask, p=0/1, layers=1, determinism); Hinrichsen values verified in the archived text and delta_hat converged to 0.1604/0.1590 at 8000 realizations; follow-ups applied (real error bars in the statistical gates, throughput figure corrected to the bench's own ~401 M site-layers/s, ceil(L/64) wording, xi_parallel ~330) |
| E3 k-core percolation | **Complete** | `kcore_decompose` + `KCoreWorkspace` + `KCoreMC` on the E1 primitives: residual-degree FIFO peeling with exact synchronous-round accounting (`KCoreSurvivorCount` ONE C0 vector sample per realization + `KCoreFinalSize`/`KCoreRounds`/`KCoreDeletedTotal`), uniform or per-vertex thresholds, arbitrary initial activity, deterministic deletion order/waves exposed post-run, mathematical order-independence argued (union-closed constraint family) and verified by 90 shuffled relabeling orders. Gates: exhaustive 2^V initial subsets on path/cycle/star/complete vs an independent repeated-full-scan reference (1456 runs incl. heterogeneous thresholds); k=0 removes nothing / k>max_degree removes everything in one round; hand-worked heterogeneous cascades; 1-core-of-connected-graph and empty-tree-k2 closed forms; bitwise site-law reduction of the activity draw; Janson-Luczak (RSA 30, 50-62 (2007), arXiv:math/0508453, archived PDF) 3-core emergence bracket with an own fixed-point solver — lambda_3 = 3.350919 reproduced, G(n=5e4) below c=lambda_3-0.15 extinguished (0.00000 ± 0.00000), above c=lambda_3+0.15 measured 0.45309 ± 0.00134 vs theory 0.45167; 0-allocation realization probe on five topologies; Criterion data. Evidence in `CMC.rs/PERCOLATION.md` (E3). Independent review passed 2026-10-23 ([report](../../reviews/2026-10-23-e-track-reviews.md)): FIFO↔synchronous-rounds equivalence re-proved and verified on 3200 random multigraph runs; λ₃ and ψ₃ independently recomputed (measured 1.06σ from theory); the probe now asserts zero allocations; an in-repo multigraph gate added |
| G1 incremental component process | **Complete** | `IncrementalUnionWorkspace` + `IncrementalSummary` + `ProcessTrajectory` (`percolation/process/incremental.rs`): the shared substrate for G2/G3 — generic incremental union over any `UndirectedGraphView` with online observables (component count, largest size with canonical lowest-index identity per the F3 convention, raw second moment, active vertex/edge counts) returned by every `activate_vertex`/`activate_edge` and recorded as `(step, summary)` trajectory entries whose budget `V + E` is exact; the step coordinate is caller-defined and never interpreted ("record process coordinates without forcing a specific law"). F3 semantics (edge active iff activated with both endpoints active, self-loop once, parallel edges separately, order independence); typed refusals leave state untouched; reset/re-prepare reuse with zero steady-state allocation; `component_size` as the G2 scoring primitive. Gates: random-interleaving prefix parity vs fresh F3 `analyze` on chain/square/random/multigraph (4 seeds); exhaustive activation orderings (5040+5040+720) vs an independent flood-fill reference at EVERY prefix; canonical identity vs F3 labels (labels and `component_size` cross-checked); random vertex replay reproduces `scan_site` prefix summaries exactly (test-only Z1 demonstration); reuse/reset/grow-shrink-grow prepare cycles; asserted 0-allocation probe on the five shared topologies (V = 4096, 128 realizations each); Criterion `incremental_process` unions/s data (25.4-27.8 M ops/s). Evidence in `CMC.rs/PERCOLATION.md` (G1). Independent review passed 2026-10-23 ([report](../../reviews/2026-10-23-gpr-entry-reviews.md)): hand-traced through a multigraph incl. edge-before-endpoints and self-loop-inside-merge; a 3,384-op drift audit found zero aggregate drift; the u128 width and trajectory budget were proven; and complete Achlioptas + invasion drivers were written against the public substrate unchanged — the law-agnostic contract verified by construction |
| G2 Achlioptas processes | **Complete** | `AchlioptasProcessMC` + `AchlioptasRule`/`CandidateSampling`/`EdgeUniverse`/`TiePolicy`/`AchlioptasGrid` (`percolation/process/achlioptas.rs`, integrator-verified 2026-10-23): rule-driven competitive edge addition on the G1 substrate — all physical parameters explicit closed enums (roadmap: "a generic `score(size_a, size_b)` callback is insufficient"), candidates scored through `component_size`, minimum-score winner, runs to complete occupation of the universe. Null steps under with-replacement are recorded, never filtered. Bitwise RNG draw contract documented (draw counts, ranges, Fisher-Yates tail arrangement, winner compaction) and replayable via `step_with_candidates` (no RNG); records ONE C0 vector sample `AchlioptasLargestSize`/`AchlioptasOccupiedEdges` over the occupation-count axis + `AchlioptasGiantReached`/`AchlioptasStepsToGiant` (explicit `2*S_max > V` criterion)/`AchlioptasSteps`/`AchlioptasNullSteps` scalars. Gates: contract-driver replay reproduces sweeps bit-for-bit + pure-trace re-replay; EVERY step vs independent flood-fill on chain/square/random/multigraph x both rules x both sampling modes; m=1 == test-side uniform permutation bitwise and == `scan_bond` statistically (|z| <= 1.1), m=1 with-replacement == coupon-collector statistics (z = -0.87); terminal = every permitted edge exactly (full + subset universes); hand-worked product-vs-sum split; asserted 0-allocation probe (five topologies, 128 realizations); Criterion `achlioptas_process` ~9-12 M steps/s; Riordan-Warnke (Science 335, 1189 (2012), arXiv:1102.5306, archived) continuity anchor: finite-cluster susceptibility peak grows 363.4 +/- 16.3 (N=2^12) -> 1898.1 +/- 116.0 (N=2^15), x5.22 at 13 sigma, and max-window S_max/N jumps shrink under window refinement (0.36/0.56 ratios). Evidence in `CMC.rs/PERCOLATION.md` (G2) |
| G3 invasion percolation | **Complete** | `invasion_run` + `InvasionWorkspace` + `InvasionPercolationMC` (`percolation/process/invasion.rs`, independently reviewed 2026-10-23 (PASS)): frozen per-realization threshold fields (site: one f64 per vertex; bond: one per physical edge; non-finite = typed rejection) grown from explicit duplicate-free seeds through a binary-heap minimum-threshold frontier with ALL policies explicit — tie = `(total_cmp(threshold), lowest dense ID)` (site: vertex ID, bond: edge ID; `-0.0 < +0.0` distinct keys, pinned by gate), stop = `GrowAll` / `FirstSpanning{from,to}` (checked after planting and after every step; two satisfied seeds stop at step zero) / `Steps(n)`, trapping = optional Wilkinson-Willemsen `Trapping { exit }` (invadable iff connected to an un-invaded exit through the defending phase; permanent monotone retirement; invaded outlets stop providing escape). Site mode queues a site once carrying its first invader's edge; bond mode queues boundary edges and discards stale pops; parallel edges are independent entries; self-loops never enter; the invasion graph is a forest (one tree per seed, trees never merge) by construction. Every accepted operation drives the G1 substrate with invasion times as the caller-defined step coordinate (trajectory `[0, 1, 1, 2, 2, ...]`). Gates (roadmap §13 verbatim): hand-worked 7-vertex parallel+self-loop multigraph — sequence, live-frontier evolution tables, and parent tree by hand in BOTH modes, cross-checked against the scan reference; 576+ heap-vs-complete-frontier-scan comparisons (chain/square/random-multigraph/two-component x both modes x continuous + tie-quantized draws x single/double seeds x all stop policies) plus a trapping grid — order, edges, thresholds, parent, stop reason/step, all outcome scalars; forest invariants via independent union-find replay + one-tree-per-seed BFS partition; with unique thresholds non-trapping bond invasion = Prim growth order exactly (single seed: Prim from the seed AND terminal tree = unique MST via independent Kruskal; k seeds: Prim from a virtual super-root joined to the seeds by -inf edges, min-weight k-tree forest; `Steps(n)` = exact prefix); tie hand examples incl. the signed-zero pin; BFS-verified spanning steps; trapping hand traces incl. retired-site-consumes-parallel-entries; bitwise adapter draw contract (dense-order); C0 one-vector-per-realization history compression; typed rejections; reuse. Asserted 0-allocation probe on five topologies (site, bond, and trapping variants, 128 realizations each); Criterion `invasion_percolation` 1.3-16.2 M steps/s. Literature: Wilkinson-Willemsen (1983) predates arXiv (Tarjan precedent — no archivable PDF); the Prim/MST identity is verified against independent in-repo implementations and anchored by the archived Sweeney & Middleton PRE 88, 032129 (2013). Evidence in `CMC.rs/PERCOLATION.md` (G3) |
| P1 particle candidate-pair adapter | **Complete** | `percolation/continuum.rs` on `feat/cmc-percolation`: borrowed flat-coordinate `CandidatePairGraph` adapter plus `candidate_pairs_from_configuration` (zero-copy) and `candidate_pairs_with_cell_list` (validated prepared-list reuse) entry points over the particle module's untouched `OrthorhombicCell`/`CellList`; boundary-inclusive candidate pairs `(i, j)`, `i < j`, in a documented cell-major canonical order; output is an ephemeral `PairGraph` (owned pair list implementing `UndirectedGraphView` with derived incidences) that exists for one F3 analysis pass — the transient counterpart to the banned persistent per-snapshot CSR. Gates: exact all-pairs parity vs an independent `O(N^2)` minimum-image scan on random 2D/3D periodic configurations, dyadic exact-boundary inclusivity (direct and through images), periodic-image-only pairs, single/zero-particle cases (zero is a typed rejection), determinism across entry points, lattice-shift sequence identity plus generic-shift set invariance, F3 `analyze`/`analyze_with_labels` vs independent adjacency-BFS parity on every summary field and canonical label, stale-list rejection, and typed invalid-domain rejections. Criterion `continuum_pairs` data (peak ~1.0e7 pairs/s; low-coordination over-subdivision cost recorded honestly). Evidence in `CMC.rs/PERCOLATION.md` (P1). Independent review passed 2026-10-23 ([report](../../reviews/2026-10-23-gpr-entry-reviews.md)): a reference built on explicit image enumeration (independent of the production round() metric) confirmed exact boundaries, image-only pairs, anisotropic cells (per-axis stencil proof), the c = L/2 tie regime, and label-level F3 parity; all bench numbers reproduced exactly; follow-ups applied (axis-order gloss, borrowed metric, tie documentation) |
| P2 Boolean percolation | **Complete** | `percolation/boolean.rs` on `feat/cmc-percolation`, integrator-verified 2026-10-23: `Radii` (monodisperse/finite-polydisperse, typed validation, re-validated at every enumeration entry point) + `PointProcess` (binomial fixed count; Poisson via a documented chunked-Knuth count-draw contract, cap 1e6) sampling, plus the equilibrium **measurement** path — `overlap_pairs`/`overlap_pairs_from_configuration` over caller-supplied coordinates and the frozen `BooleanPercolationMC::measuring` adapter whose sweep is a documented no-op (connectivity is a measurement of accepted coordinates, never a second geometry sampler). Overlap enumeration composes P1 at `2 r_max` with the inclusive squared per-pair filter `d^2 <= (r_i + r_j)^2`; monodisperse radii pass through unfiltered (bit-identical criterion, unit-tested against the filtered path). Ephemeral `OverlapGraph` (P1 transient pattern) feeds the F3 analyzer; the adapter records the F4 scalar plan per realization (Poisson vacuum realizations record the all-zero summary), `thermalization_sweeps = 0`. Gates: exact all-pairs overlap parity vs an independent `O(N^2)` minimum-image scan (2D/3D polydisperse, three seeds), generic-shift set+summary invariance and exact-lattice-shift sequence identity, particle-permutation relabeling invariance, periodic-representative invariance, exact zero-radius limit (only coincident pairs), fully-overlapping limits (1D complete graph at cutoff exactly `L/2`; 2D single-component saturation within the minimum-image ceiling), Poisson count mean/variance CI gates direct and through `ActiveVertexCount`, frozen-configuration zero-variance remeasurement + zero-copy path parity, bitwise per-seed determinism, typed rejections, and the 2D disk covered-fraction bracket — spanning probability Wilson-95% below 1/2 at `eta_c/1.15` and above at `1.15 eta_c` (`eta_c = 1.12808737(6)`, Mertens & Moore PRE 86, 061109 (2012), archived; statistical bracket, never exact). Probe asserts zero adapter-added allocations per realization beyond the documented P1 transient; Criterion `boolean_percolation` data. Evidence in `CMC.rs/PERCOLATION.md` (P2) |
| R1 external bond-configuration analysis | **Complete** | Decoupled analysis contract delivered in-branch on `feat/cmc-percolation`, awaiting independent review: `ExternalBondConfiguration` borrows a caller-supplied bond mask (one `bool` per physical edge of any `UndirectedGraphView`) and implements `EdgeActivity` (mask) + `VertexActivity` (all-active), feeding the shared F3 `analyze`/`analyze_with_labels` and W2 `analyze_winding` directly; `analyze_bonds` is the one-call raw-mask convenience; `BondSample` is an owned sorted/deduplicated open-edge record with mask round-trips (R3-readiness — no sampler implemented here, by design). Gates: independent adjacency-BFS parity (summary/labels/queries) on chain/square/random-multigraph at p 0.1/0.5/0.9 across seeds through both entry points; typed length/empty-mask rejections with legal all-false semantics; exhaustive 2^E parity on a parallel-edge + self-loop multigraph (open self-loop = union no-op); 2x2-torus `analyze_winding` bridge vs `StaticConfiguration` with the same bonds over 64 draws; Bernoulli-vs-checkerboard sampler-decoupling demo through one API; Criterion bench on the five shared topologies shows the external path equal to the `StaticConfiguration` path within quick-mode noise. Evidence in `CMC.rs/PERCOLATION.md` (R1 section) |
| R2 q=1 random-cluster reduction | **Complete** | q = 1 reduction delivered in-branch on `feat/cmc-percolation`, integrator-verified 2026-10-23: `RandomClusterMC` (percolation/random_cluster.rs) with an explicit `ClusterQ` pinned to `q = 1` (any other value, NaN included, is the typed `UnsupportedQ` rejection citing this roadmap's separate-plan requirement for generic real-q samplers); sampling by delegation to the F2 `BondBernoulli` law — no duplicated sampling code, bitwise draw identity verified — and analysis through the R1 external-bond path (`analyze_bonds`), exposing each realization as a `BondSample` (the R3 contract surface). The exact statement validated: the q = 1 FK measure `p^{\|A\|} (1-p)^{E-\|A\|} q^{k(A)}/Z(p,q)` equals the product Bernoulli measure because the cluster factor is exactly `1^{k(A)} = 1` and `Z(p,1) = (p + (1-p))^E = 1`. Gates: exhaustive 2^E enumeration on a 5-edge chain-ish graph (square + pendant, k(A) not a function of \|A\|) and a 6-edge self-loop/parallel multigraph at uniform p 0/0.25/0.5/0.75/1 plus heterogeneous fields — per-configuration FK-vs-Bernoulli weight equality, Z = 1, and ⟨k⟩/⟨S_max⟩ distribution-identity under an independent flood fill cross-checked against the production analyzer; bitwise same-seed mask identity through `BondSample::to_mask`; multi-seed two-sample z-scores (default 16 seeds, `SCUTTLE_ZSCORE_SEEDS`) of GiantFraction/finite-cluster susceptibility vs the `StaticPercolationMC` bond path at matched p = 0.5 on the 8x8 open square, per-seed \|z\| < 4 and \|mean z\| < 2; Criterion realization-throughput data. Evidence in `CMC.rs/PERCOLATION.md` (R2 section) |
| R3 q=2 Edwards-Sokal / Swendsen-Wang exposure | **Complete** | Delivered in-branch on `feat/cmc-percolation`, independently reviewed 2026-10-23 (PASS): `EdwardsSokalMC` (percolation/edwards_sokal.rs) — an owned-`CsrLattice` Swendsen-Wang adapter over strictly ferromagnetic `Bond.weight > 0` couplings (typed `Coupling` rejection otherwise) whose `sweep()` performs the conditional FK bond draw `p_e = 1 - exp(-2 beta J_e)` on same-spin edges in a documented dense-order RNG contract (bitwise replayed, F2 endpoint convention) and exposes that ACTUAL draw as an R1 `BondSample` (frozen-spins mode pins the spin configuration for pure conditional draws); the following `measure()` analyzes that same sample through the shared F3 analyzer and F4 `ObservablePlan` with the `AllActiveVertices` vertex source. Gates: frozen-spin exhaustive 2^E subset enumeration vs the product Bernoulli law (\|z\| < 5 per subset, two 4-spin configs, per-subset `BondSample`-vs-`analyze_bonds` analysis parity); partition-function agreement `Z_Ising = exp(beta ΣJ_e) · Σ_s Σ_A p^{\|A\|}(1-p)^{E-\|A\|}` — double sum over spins AND bond subsets vs the random-cluster single sum with own-union-find `2^{k(A)}` vs independent brute-force Ising, to 1e-12 on a heterogeneous graph, the 4-ring, and a parallel+self-loop multigraph at beta 0.3/0.9; `P_FK(i ~ j) = <s_i s_j>` on a thermalized proper ES chain (boundary-crossing mean of the exposed draw and the same chain's spin correlator, each \|z\| < 5 vs the exact brute-force correlator and to each other under the exact conditional-mean stderr); worm-representation consistency at the observable level only (even/two-defect tanh-weight subgraph sums vs spin enumeration) with the never-relabel assertion documented in module and test headers. Edwards-Sokal (1988) and Swendsen-Wang (1987) pre-date arXiv — no archivable PDFs (Tarjan precedent); the exact double-sum identities are the anchors. Criterion quick data (draw+exposure ~50-70 M edges/s, full cycle 16-21 M edges/s on the five shared topologies). Evidence in `CMC.rs/PERCOLATION.md` (R3 section) |
| E2 bootstrap percolation | **Complete** | `BootstrapWorkspace` + `BootstrapPercolationMC` on the E1 primitives: one O(V + E) monotone activation-threshold closure, generic over `UndirectedGraphView` with seeds as any `VertexActivity`, homogeneous `r` or per-vertex `u32` thresholds, and TWO explicit clocks — `SynchronousRounds` (Holroyd's parallel rounds) and `QueueClosure` (order-independent FIFO monotone closure) — recorded as `BootstrapActiveCount` (ONE C0 vector sample per realization, fixed truncate/pad shape) + `BootstrapFinalActive`/`BootstrapActivatedTotal`/`BootstrapRounds`\|`BootstrapEvents` scalars. Gates: exhaustive 2^V seed masks on path/cycle/star/complete for every r in 0..=max_degree+1 AND seeded heterogeneous vectors vs an independent repeated-full-scan reference (~1920 closure-set + 672 heterogeneous + 528 trajectory production runs, both clocks); r=0 activates all / r>max_degree preserves seeds verbatim; FIFO = one-random-vertex-at-a-time = shuffled-sweep = synchronous closure on random multigraphs (test-RNG orders); hand-worked heterogeneous cascade; bitwise replay of the one-f64-per-vertex density draw contract; Holroyd (PTRF 125, 195 (2003), arXiv:math/0206132, archived PDF) finite-size metastability brackets — his own full-activation observable I(L,p) is 1.0000 ± 0.0013 at p+ = pi^2/(18 ln L) and <= 0.0025 at p+/4 for L = 32/64/128; the bisected half-majority density p_hat = 0.0692/0.0577/0.0480 decays monotonically (the 1/ln L law) with p_hat ln L = 0.240/0.240/0.233 strictly inside (lambda/4, lambda); the metastable window is sharper than a factor 1.5 (P(p_hat/1.5) <= 0.018, P(1.5 p_hat) >= 0.992); 0-allocation realization probe on five topologies (asserted); Criterion data. Evidence in `CMC.rs/PERCOLATION.md` (E2). Independent review passed 2026-10-23 ([report](../../reviews/2026-10-23-e-track-reviews.md)): hand traces matched production exactly; the Holroyd anchor verified against the archived paper's own slow-convergence note; doc counts corrected upward (720/15/~4000); the Criterion chain row's two densities labeled honestly; both clocks now allocation-asserted |
| W1 periodic embedding capability | **Complete** | Displacement data layer complete on `feat/cmc-percolation`: `PeriodicEmbedding`/`LatticeEmbedding` store one explicit integer cell displacement per directed incidence (physical edge, endpoint), recorded by the pbc chain/square/hypercubic and triangular builders at bond emission and never inferred from `BondType` or vertex IDs; honeycomb/kagome embeddings Not implemented (adjacency construction drops per-bond direction). Hand-derived 2x2 torus, geometry property, parallel-edge, OBC-rejection, plaquette-cocycle, and fundamental-winding gates pass; independent review passed (reviewer hand-recomputed the 2x2 triangular displacements and cocycle sums). Evidence in `CMC.rs/PERCOLATION.md` |
| W2 winding analyzer | **Complete** | Integer-potential union-find (`analyze_winding`/`WindingWorkspace`, checked `i32`, no silent wrap) over the W1 displacements reports per-component winding rank and echelon generators; gates green on `feat/cmc-percolation`: non-wrapping loop / single- and multi-axis / diagonal-generator hand cases, independent finite-cover lift-BFS parity, origin/gauge/edge-reversal invariances, OBC explicit errors, full 2x2 (2^8) and 3x3 (2^18) bond-torus enumeration vs reference (P(wrapping) = 187/256 and 186684/262144 at p = 1/2), 16-seed scheduler z-scores vs the exact 2x2 enumeration (verified at 64 seeds), `Wrapping`/`WindingRank` plan observables on `StaticPercolationMC::with_embedding`, 6 zero-allocation probe records, and 15 Criterion IDs (winding scan ~2.3-2.6x F3 per edge). Evidence in `CMC.rs/PERCOLATION.md`; Experimental until review passes; Independent review passed 2026-10-22 ([report](../../reviews/2026-10-22-stage-reviews.md)): the reviewer re-derived both torus wrapping fractions (187/256, 186684/262144 and the mean ranks) with two from-scratch brute forces, hand-verified the union-find potential and SL(2,Z) echelon algebra, and re-derived the i32 overflow bound; follow-ups applied (exact enumeration anchors asserted, benchmark ratio prose corrected) |
| Z1 Newman-Ziff parameter scans | **Complete** | Site/bond incremental activation scans with reusable `NewmanZiffWorkspace`, `ScanCurve`/`ScanPoint`, validated-`Probability` binomial reweighting (`ReweightWorkspace`), and independent-scan `ScanEstimator`; every occupation-count prefix matches a fresh F3 `analyze` across chain/square/cubic/triangular/honeycomb/kagome/random/multigraph, full 3x3 site/bond enumeration at p = 0.3/0.5/0.7 feeds |z| < 4.5 gates (256 default scans, verified at 64 seeds), direct Bernoulli CI, endpoints/reproducibility/empty cases, 10 zero-allocation probe records, and a Criterion 64-replica/17-p comparison showing 1.83-4.65x lower time than direct Bernoulli+F3 analysis. Core scalar API only; no scheduler adapter or C0 vector output. Evidence in `CMC.rs/PERCOLATION.md` and `CMC.rs/VALIDATION.md`; Experimental until review passes; Independent review passed 2026-10-22 ([report](../../reviews/2026-10-22-stage-reviews.md)): the reviewer verified the reweight recurrence against the archived Newman-Ziff paper Eq. (1) to 9.5e-14 up to n=10^6, closed the prefix-parity loop with a from-scratch flood fill on multigraphs, and reproduced the speedup (1.80-4.87x); follow-ups applied (saturating_add, excerpt note) |
| D/E/G/P/R tracks | **Not started — entry stages unblocked** | F0-F4 and Z1 prerequisites are complete; parallel entry stages D1, E1, G1, P1, and R1 may start in any order, subject to their dependency gates above |

PR #4 carries the F0-F4 final architecture plus the in-branch N1 delivery as its vehicle: the
provisional API it originally proposed is deleted in-branch (never merged, so
no compatibility debt), and every applicable validation gate was migrated to
the final API before the F4 review landed.

## M6: Coverage Audit — 2026-10-23

Final-state audit of the platform against this roadmap, performed on
`feat/cmc-percolation` after all 27 development stages landed.

### Domain matrix (CMC.rs/PERCOLATION.md)

| Status | Rows | Notes |
|---|---|---|
| Validated | 16 | Independent review passed; production on the `dev` merge |
| Integrator-verified | 5 | G2, G3, P2, R2, R3: all gates green and integrator-checked 2026-10-23 (hand computations, code audits, 47/47 tests); independent reviewer re-run recommended before the `dev` merge (checklist embedded in `docs/reviews/2026-10-23-final-round-reviews.md`) |
| Experimental | 3 | Deliberate: F1 topology API, F2 configuration/law substrate, and W1 embedding remain provisional-until-merge BY DESIGN (the never-merged-API rule) — all three passed independent review and everything above is built on them |
| Not implemented | 1 | Generic real-q random-cluster (Sweeny/Chayes-Machta) — explicitly a separate future plan per §15 |

### Cross-module evidence

- **Tests**: workspace 1322 passed, 0 failed (Carlo.rs 321+ incl. hdf5; CMC.rs suite + lib; QMC, MCMC, scamper all green). MPI suites green at np=2 (last run at the D3 round; carlo-framework nightly continues to monitor).
- **Benchmarks**: 28 Criterion/allocation bench targets across the workspace; every production family reports throughput + zero-allocation probes with asserted counts (the kcore/bootstrap/invasion/achlioptas probes assert; the historical prints are retained in earlier reports).
- **Literature provenance**: 11 bibliography rows, 9 arXiv PDFs archived under `docs/ref/` with Markdown conversions (`convert_to_md.py`, incremental); every physics claim in PERCOLATION.md cites its archived anchor with the agreement numbers recorded. Pre-arXiv classics (Tarjan 1972, Wilkinson-Willemsen 1983, Edwards-Sokal 1988, Swendsen-Wang 1987, Fortuin-Kasteleyn 1972) are noted where no PDF can exist and replaced by definitional/exhaustive anchors.
- **Reviews**: 13 reports under `docs/reviews/` — 22 stages carry independent reviewer verdicts; the final five carry the integrator-verification record with the honest weaker-independence caveat.
- **Code surface**: the percolation module is 16,064 lines across the layered architecture (topology → activity/law → analyzer → observables → adapters), with the process/continuum/bonds/random-cluster/edwards-sokal extensions composed on the shared substrate rather than forking it.

### Definition-of-done residue

- The five integrator-verified stages have not had a *fresh-eyes* reviewer re-run since their delivery commit; the embedded checklist exists for exactly this and is the one recommended pre-merge action.
- Adaptive residual-degree attack (roadmap §9 N3 tail), heterogeneous Newman-Ziff scans, honeycomb/kagome embeddings, generic real q, and continuum winding remain Not implemented by design and are the natural next-roadmap items.
- PR #4 remains the delivery vehicle: the provisional API deleted in-branch, all validation assets migrated to the final API (verified across F4/N1/N2/N3/D1-D3/E2-E3/R1-R3 gate suites).

**Milestone M6 status: complete.** The platform is ready for the
independent-review re-run of the final five stages and the `dev` merge.
