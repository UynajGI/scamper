# CMC.rs Percolation Status and Scientific Definitions

> F0 baseline established 2026-09-03. This document describes the current
> `feat/cmc-percolation` reference implementation and freezes vocabulary for the
> staged platform roadmap. It does not make the provisional API stable.

The normative design and phase order are the [platform philosophy](../docs/metaphysics/2026-09-03-percolation-platform-philosophy.md)
and [implementation roadmap](../docs/plans/2026-09-03-percolation-platform-roadmap.md).

## Status vocabulary

- **Supported:** an API with explicit semantics exists, but not every production
  validation gate necessarily passes.
- **Validated:** the family-specific semantic, exact/reference, statistical,
  invalid-domain, reproducibility, and performance gates pass.
- **Experimental:** usable reference behavior exists, but it is not the intended
  production API or has not passed all production gates.
- **Not implemented:** no supported claim; a generic graph routine producing a
  number does not change this status.

## Support matrix

| Family or capability | Status | Current scope or reason |
|---|---|---|
| Read-only undirected topology capability | Experimental | Independently reviewed F1 `GraphView`/`UndirectedGraphView` contracts use dense IDs and static dispatch; `CsrLattice` and validated zero-copy `BorrowedUndirectedCsr` implement them. The API remains provisional until the platform merges into `dev`. |
| Directed topology capability sketch | Experimental | Separate outgoing/incoming signatures exist only crate-private with compile fixtures; no directed API is re-exported before D1 fixes arc identity and storage contracts. |
| Static configuration and Bernoulli-law substrate | Experimental | F2 provisional API separates private activity storage from independent site/bond/mixed laws, validates probability domains and topology lengths, and samples allocation-free after preparation. The independently reviewed F3 analyzer/workspace and F4 observables build on this substrate; the layer remains Experimental until the platform merges into `dev`. |
| Static site, bond, and mixed Bernoulli percolation on arbitrary owned or borrowed undirected graphs | Validated | All family-specific gates pass: F3 component analysis with reusable workspace, canonical identities, and multi-query crossing; F4 `ObservablePlan` observables and owned `StaticPercolationMC` with every PR #4 validation asset migrated, the zero-allocation adapter gate, and 16-seed z-scores; independent review passed. Production-ready once the branch merges into `dev`. |
| Heterogeneous occupation probabilities | Validated | N1 promotes the F2 `ProbabilityField` substrate to the `StaticPercolationMC` adapter: `StaticLaw::site_heterogeneous`/`bond_heterogeneous`/`mixed_heterogeneous` own one validated `Probability` per vertex and/or per physical edge, with adapter-constructor domain validation and scheduler scheduling through `Run::from_parts()`. Gates passed: exact non-identically-distributed small-graph enumeration (independently reproduced), bitwise uniform and mixed-pure reductions, endpoint determinism, typed length/invalid-probability rejection, 16-seed scheduler z-scores for site/bond/mixed (enumeration and Poisson-binomial references), uniform-vs-heterogeneous Criterion comparison, and a 0-allocation heterogeneous adapter probe. Independent review passed; production on the `dev` merge. Network/domain facades over heterogeneous fields (random failure, degree observables) remain Not implemented (N2). |
| Network random failure and degree-conditioned observables | Not implemented | No network facade, probability field, or original-degree profile. |
| Targeted or adaptive attack and robustness curves | Not implemented | No removal process or vector-valued realization output. |
| Directed static reachability, weak/strong components | Not implemented | `CsrLattice` is undirected; directedness is not a mode flag. |
| Layered space-time directed percolation | Not implemented | No causal-layer process or survival/arrival profiles. |
| Bootstrap percolation | Not implemented | No activation-threshold closure process. |
| k-core percolation | Not implemented | No residual-degree peeling process. |
| Achlioptas processes | Not implemented | No candidate-edge process, rule, or process-time output. |
| Invasion percolation | Not implemented | No frozen thresholds or minimum-frontier growth. |
| Continuum/Boolean percolation | Not implemented | Particle geometry and cell lists exist elsewhere, but no continuum percolation adapter or validated overlap analyzer exists. |
| FK/random-cluster q=1 | Not implemented | Ordinary bond Bernoulli exists, but no external bond-configuration/FK contract or q=1 reduction validation exists. |
| FK q=2 Edwards-Sokal/Swendsen-Wang | Not implemented | SW exists, but its sampled FK bonds are not exposed for this analysis. The worm high-temperature even subgraph is not FK. |
| Generic real-q random-cluster model | Not implemented | No Sweeny or Chayes-Machta sampler. |
| Boundary crossing between explicit vertex sets | Validated | F3 `BoundaryQuery` is validated for owned/borrowed generic undirected views, supports several queries in one component pass, and defines empty/overlap/dedup semantics; F4 migrates the full crossing test assets. The square left/right-column and chain endpoint defaults live in the `StaticPercolationMC` parameter schema; all other parameter-built topologies require explicit site sets. |
| Periodic embedding capability (W1) | Experimental | Directed-incidence integer cell displacements exist as explicit recorded data: the `PeriodicEmbedding` capability (`topology/embedding.rs`) and the `LatticeEmbedding` table serve the pbc chain/square/hypercubic builders and the triangular builder, which record each displacement at bond emission and never infer it from `BondType` or vertex IDs. Independently reviewed; see the W1 section below. |
| Periodic wrapping/winding analysis (W2) | Not implemented | No winding analyzer, wrapping query, or winding-rank reporting over the W1 displacements. Crossing is not wrapping. |
| Newman-Ziff parameter scans | Not implemented | No incremental activation curve or binomial reweighting. |

At F0 no family was **Validated**. Since F4 passed independent review, the
static family and boundary crossing are the first Validated rows; N1 adds the
heterogeneous-probability family on the same production-conditioned basis.
Their production availability is conditioned on the branch merging into `dev`.

## Scientific definitions

Let one configuration induce tracked connected components
`C = {C_1, ..., C_K}` with sizes `s_i = |C_i|`. In site and mixed percolation,
tracked vertices are occupied sites. In bond percolation every vertex is
tracked, including isolated vertices. An edge connects in site mode when both
endpoints are occupied, in bond mode when the edge is occupied, and in mixed
mode only when the edge and both endpoints are occupied.

### Components and moments

- **Component size:** `s_i` is the number of tracked vertices in component
  `C_i`; edges are never counted as component size.
- **Largest component:** `S_max = max_i s_i`, with `S_max = 0` when `K = 0`.
  The current reference reports only its size, not identity.
- **Giant-component fraction:** `G = S_max / V_norm`. The normalization must be
  named. The platform default is total topology vertices `V_norm = |V|`;
  normalization by active/tracked vertices is a distinct observable. This is
  defined vocabulary, not a current PR #4 output.
- **Raw cluster moment:** `M_k = sum_i s_i^k` over the stated component set,
  with no normalization. PR #4 reports `M_2` as `SecondMoment` and includes the
  largest component.
- **Finite-cluster susceptibility:** choose one largest component per sample,
  breaking equal-size ties by lowest stable component/root ID, and remove it.
  Then `chi_f = (sum_{i != max} s_i^2) / (sum_{i != max} s_i)`. The exclusion
  happens before ensemble averaging. If the denominator is zero, the observable
  is undefined and must be represented explicitly as missing/invalid, not
  silently as zero or NaN in an advertised scalar result. Removing every tied
  largest component or dividing by total `|V|` are different named policies.
  PR #4 does not compute `chi_f` and its aggregate outputs cannot reconstruct it.

### Connectivity events

- **Boundary crossing:** for explicit vertex sets `A` and `B`, crossing is true
  iff one tracked connected component intersects both sets. Sets may overlap;
  an active vertex in `A intersect B` therefore satisfies the event. The event
  is defined by the supplied sets, not inferred geometry. Its ensemble mean is
  the crossing probability.
- **Wrapping/winding:** on a periodic embedded topology, lift each occupied
  path using explicit integer cell displacements. A component wraps when a
  cycle has nonzero accumulated displacement; winding records generators,
  axes, or rank rather than merely boundary-set connectivity. Vertex IDs,
  `BondType`, and periodic adjacency alone do not define winding. PR #4 does
  not implement wrapping.

### Network conditioning

- **Original degree** is degree in the immutable input graph before occupation
  or removal. An observable conditioned on original degree keeps this label for
  the complete realization or attack trajectory.
- **Residual degree** is degree in the currently active induced graph and may
  change during a process. Residual-degree conditioning or adaptive attack is a
  different process and must be named explicitly.
- Empty degree bins require an explicit missing/undefined representation; they
  are not zero-valued observations.

### Samples and lifecycle

A **realization** is one complete draw or process trajectory under a stated law.
An **ensemble sample** is one realization contributed to an estimator. For
static Bernoulli percolation these coincide: one `sweep()` redraws all relevant
occupancies independently and one following `measure()` analyzes that draw.
Samples are i.i.d.; set `thermalization_sweeps = 0`. For attacks, cascades,
growth, or directed space-time processes, all steps of one trajectory remain
one realization and must not be misreported as independent scalar samples.

## PR #4 behavioral baseline

PR [#4](https://github.com/UynajGI/scamper/pull/4) is open against `dev`; the
recorded implementation head is `f6031bf88bea05f61ecbc2eb669268712c65c866`.
It provided `PercolationMode`, public `Vec<bool>` fields in `OccupancyState`,
`cluster_stats(&CsrLattice, ...)`, `UnionFind`, and scheduler-ready
`PercolationMC`. These names, storage choices, constructor/parameter schema,
and observable names are provisional and are not compatibility commitments.
**F4 has since deleted this entire public API without a compatibility
wrapper** (per the roadmap section 7 migration rule for never-merged APIs) and
migrated every applicable validation asset to the final API; the description
below is the frozen behavioral record the migration was checked against.

Reference behavior:

- each sweep independently redraws uniform site/bond occupancy;
- `cluster_stats` allocates fresh union-find and boundary scratch on every call;
- physical undirected edges are scanned once, then tracked vertices are scanned;
- site/mixed omit closed sites; bond mode counts every isolated site;
- outputs are `MaxCluster`, raw `SecondMoment`, `NClusters`, `Spanning`, plus
  `Occupied` for pure modes or `OccupiedSites`/`OccupiedBonds` for mixed mode;
- crossing defaults are open square left/right columns and open chain endpoints;
  all other parameter-built topologies require explicit site sets.

The validation assets documented in [VALIDATION.md](VALIDATION.md) — exact
full enumeration and hand-derived 1D/2x2 forms, an independent flood-fill
reference, mixed-to-pure reductions, cross-family and random-graph parity,
scheduler tests, fixed-seed reproducibility, 16-seed z-score gates, and
ignored square/cubic critical-region checks — have all been migrated to the
F4 API (see the F4 section below).

The former PR #4 limitations — public storage layout, fresh per-analysis
allocations, no component identity, no giant fraction or finite-cluster
susceptibility — are resolved by F3/F4. Still true for the whole platform:
a uniform-only `FromParams` parameter path (heterogeneous laws enter via direct construction, N1 Validated),
owned `CsrLattice` coupling in the parameter path (borrowed graphs use the
core APIs directly), scalar results, and no wrapping, directed, process,
continuum, FK, attack, or Newman-Ziff implementation.

## Performance baseline procedure

`percolation_baseline` measures the F0 reference sample+analyze path, rebuilt on
the final API after F4: an i.i.d. redraw into a `Vec<bool>` reference
configuration plus the production `analyze` with a workspace allocated per call
(the historical PR #4 allocation profile). The 2026-09-03 F0 table below was
recorded with the original PR #4 implementation at commit `f6031bf`. It does not
enable, depend on, or link
`allocation-counter`; the allocation instrumentation lives in the separate
`percolation_allocations` target behind the `allocation-probe` feature. Cases
are named with topology, `V`, `E`, `p_site`, and `p_bond`; all modes use
`p_site = 0.5927`, `p_bond = 0.5`. Topologies are an open chain (`V=4096`),
open square (`64x64`), open cubic (`16^3`), seeded sparse Erdos-Renyi-like graph
(`V=4096`, `E=16384`), and seeded degree-proportional preferential-attachment
graph (`V=4096`, four distinct attachments per added vertex). The latter starts
with edge `(0,1)` and attachment pool `[0,1]`; setup asserts both `pool.len() =
2E` and, vertex by vertex, pool multiplicity equals graph degree.

```bash
# Criterion throughput only; allocation-counter is disabled and unlinked.
cargo bench -p cmc-rs --bench percolation_baseline --no-default-features -- --test
cargo bench -p cmc-rs --bench percolation_baseline --no-default-features

# Standalone allocation/storage probe; no Criterion timing occurs here.
cargo bench -p cmc-rs --bench percolation_allocations \
  --features allocation-probe -- --test
```

The standalone probe prints `VEC_BOOL_PROBE`, `PERCOLATION_ALLOC`, and
`TOPOLOGY_VIEW_ALLOC` records. The topology construction/scan probe uses a
separate measurement region from the F0 sample/analyze allocation baseline.
`site_logical_capacity` and `bond_logical_capacity` are the capacities reported
by `Vec<bool>` in logical elements; they are not byte counts by API contract.
`estimated_owned_storage_bytes` sums the owned CSR vector capacities plus
occupancy storage using the byte-per-element behavior measured for this exact
compiler. It excludes struct headers, allocator metadata, boundary query
vectors, RNG state, and analyzer scratch, so it is a documented estimate rather
than RSS or an allocator measurement of the complete case. Its per-vertex and
per-edge columns are two normalizations of the same estimate, not additive.

`allocation-counter` separately measures 128 post-warmup sample+analyze calls
in one thread and reports allocator-observed calls and bytes per sample plus the
peak live allocation count and bytes in the measurement region. There is no CI
performance threshold because allocator and timing measurements are
environment-sensitive.

### Recorded F0 run

Recorded on 2026-09-03 with:

- CPU: Intel Xeon Gold 6148 at 2.40 GHz; Linux 6.8.0-134-generic x86_64.
- Compiler: `rustc 1.98.0 (88d9e12ae 2026-08-18)`, LLVM 22.1.8;
  `cargo 1.98.0 (797e8a9bc 2026-08-05)`.
- Throughput command: `cargo bench -p cmc-rs --bench percolation_baseline
  --no-default-features`.
- Allocation command: `cargo bench -p cmc-rs --bench percolation_allocations
  --features allocation-probe -- --test`.
- Criterion short-baseline configuration: 10 samples, 100 ms warm-up, 250 ms
  measurement per case. Results below are point estimates; timing is not a CI
  gate.

The throughput target's normal dependency tree contained no `allocation-counter`,
and its final executable contained no `allocation_counter` marker. Thus the
following measurements do not include that crate's global allocator wrapper.

| Topology | V | E | Estimated owned storage bytes | Est. bytes/V | Est. bytes/E | Site samples/s | Bond samples/s | Mixed samples/s |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| chain | 4096 | 4095 | 303,079 | 73.994 | 74.012 | 14,790 | 14,823 | 11,653 |
| square | 4096 | 8064 | 565,128 | 137.971 | 70.080 | 9,027 | 7,544 | 7,190 |
| cubic | 4096 | 11,520 | 941,320 | 229.814 | 81.712 | 7,469 | 6,392 | 5,308 |
| sparse ER | 4096 | 16,384 | 1,101,832 | 269.002 | 67.250 | 4,243 | 4,203 | 3,852 |
| preferential attachment | 4096 | 16,374 | 1,101,374 | 268.890 | 67.264 | 6,515 | 5,553 | 4,552 |

## F1 read-only topology capability

The F1 API remains **Experimental and provisional** after independent review.
`GraphView` and
`UndirectedGraphView` expose dense `VertexId`/`EdgeId`, physical edge endpoints,
and incidence iteration through static dispatch. The public iterator contract is
only `Iterator`; current CSR implementations may provide stronger iterator
properties internally without imposing them on future mmap, generated, or
filtered views. `GraphView::vertex_ids()` is a default convenience, not a
performance requirement for every implementation.

IDs are created publicly only by `GraphView::vertex_id` and
`UndirectedGraphView::edge_id`. Their transparent newtypes guarantee a dense
numeric representation but carry no graph provenance. An ID must be used only
with the unchanged validated view that produced it. Cross-graph IDs or IDs kept
across structural mutation violate the caller contract and may panic at
infallible accessors, but cannot cause memory unsafety. This matters for
`CsrLattice`: its existing public fields remain mutable, so construction or an
explicit successful `CsrLattice::validate()` is the precondition for trait use.

`BorrowedUndirectedCsr` borrows offsets, neighbors, physical edge IDs, and AoS
`[[usize; 2]]` physical endpoints without copying. Its constructor accepts the
empty graph (`offsets = [0]`) and validates offsets, ranges, endpoint agreement,
and exactly two correctly oriented incidences per physical edge. Parallel edges
and self-loops are preserved; adjacency ordering is unrestricted.
Offsets/neighbors alone are insufficient because stable physical edge IDs cannot
be inferred for parallel edges or self-loops.

Construction uses temporary validation state of one byte per physical edge and
drops it before returning the zero-payload view. `try_reserve_exact` capacity or
allocation failures return `TopologyError::ValidationCapacity`; allocator
configurations that abort on OOM cannot be promised recoverable.

Directed out/in signatures remain a crate-private compile sketch. They are not
crate-root exports or compatibility commitments; D1 must define stable arc
identity before publishing them.

`topology_view` compares direct owned access, generic owned access, and generic
borrowed access on an open 128x128 square (`V=16,384`, `E=32,512`, `I=65,024`).
Run with:

```bash
cargo bench -p cmc-rs --bench topology_view --no-default-features -- --test
cargo bench -p cmc-rs --bench topology_view --no-default-features
```

Recorded on 2026-09-03 on the same Xeon Gold 6148 and rustc 1.98.0 toolchain as
F0, using 20 samples, 250 ms warm-up, and 1 s measurement per case. These are the
second complete run after the three-path benchmark was added:

| Scan | Direct owned | Generic owned | Owned delta | Generic borrowed | Borrowed delta |
|---|---:|---:|---:|---:|---:|
| Physical edges | 24.211 us | 23.472 us | -3.05% | 13.768 us | -43.13% |
| CSR incidences | 58.033 us | 56.896 us | -1.96% | 57.978 us | -0.09% |

The first complete run also stayed below the 5% regression gate: generic owned
was -0.08% for edge scan and +0.32% for incidence scan; generic borrowed was
-42.93% and +0.89%, respectively. Borrowed edge scans benefit from compact AoS
endpoint input compared with the wider `Bond` records; this is a data-layout
effect, not evidence that trait dispatch accelerates work. There is no CI timing
threshold.

The separate allocation probe constructs and validates the borrowed view outside
the scan measurement region, then repeats a generic edge+incidence scan 128
times. For its open 64x64 square (`V=4,096`, `E=8,064`, `I=16,128`), the 64-byte
view header is an x86_64 rustc 1.98 measurement, not a cross-platform guarantee.
Construction performed one 8,064-byte scratch allocation, exactly one byte per
physical edge. The view owns zero payload bytes. The four borrowed input arrays
total 419,848 bytes (`102.502 bytes/V`, `52.064 bytes/E`); each steady-state scan
performed **0 allocations**, allocated **0 bytes**, and had **0 peak live
allocations/bytes**. Run the probe with:

```bash
cargo bench -p cmc-rs --bench percolation_allocations \
  --features allocation-probe -- --test
```

## F2 activity configuration and Bernoulli laws

The provisional F2 API separates generation from storage. `StaticConfiguration`
owns private vertex and physical-edge masks; the mask type and storage accounting
are crate-private and are not returned by any public method. `VertexActivity`
and `EdgeActivity` accept dense `VertexId`/`EdgeId` values under the same
unchanged-view provenance contract as F1: IDs from another graph or after
structural mutation violate the caller contract and may panic, but cannot cause
memory unsafety. The `law` module is private; reviewed law types are selectively
re-exported. The intended crate-root F2 surface comprises `StaticConfiguration`,
`VertexActivity`, `EdgeActivity`, `Probability`, `ProbabilityField`,
`ProbabilityError`, `SamplingError`, `SiteBernoulli`, `BondBernoulli`, and
`MixedBernoulli`.

The internal mask has logical all, none, and dense states while retaining any
dense byte allocation across endpoint transitions. A new pure site law
materializes only its vertex mask and makes edges logically all-active; a new
pure bond law materializes only its edge mask and makes vertices logically
all-active. Mixed sampling materializes both. Every successful sample writes
complete semantics, so switching law types cannot expose stale state.
`StaticConfiguration::resize` changes both domains and clears logical activity
while retaining usable capacity. Sampling checks graph/configuration and all
heterogeneous-field lengths before any mutation or RNG consumption.

`Probability` validates finite `[0, 1]` values once. `ProbabilityField` supports
uniform, borrowed, and owned validated values without per-sample revalidation.
Uniform exact endpoint fields use logical all/none states and consume no RNG for
site, bond, or any mixed endpoint combination. Heterogeneous fields are sampled
through the dense element loop; exact endpoint elements skip RNG individually,
but no O(1) or whole-field no-RNG contract is claimed for heterogeneous data.
N1 has since promoted heterogeneous sampling through the `StaticPercolationMC`
adapter to the Validated heterogeneous family (see the N1 section); the law
module itself remains the Experimental F2 substrate it always was, and
network/domain facades over heterogeneous fields remain Not implemented (N2).

The F2 statistical gate aggregates activity across every requested seed rather
than taking a maximum over a growing set of per-seed scores. For each of six
separate domains (uniform site, uniform bond, uniform mixed vertex/edge, and
heterogeneous mixed vertex/edge), probabilities `p_i` give exact total moments

```text
mean     = n_seeds * samples_per_seed * sum_i p_i
variance = n_seeds * samples_per_seed * sum_i p_i * (1 - p_i).
```

Uniform fields reduce to the binomial formula; heterogeneous fields use the
exact Poisson-binomial mean and variance. Each aggregate requires `|z| < 4.5`.
For six approximately normal two-sided gates, the resulting family-wise false
alarm probability is about `6 * 2 * Phi(-4.5) = 4.1e-5`. Zero variance requires
an exact success count. There is deliberately no fixed per-seed maximum,
sign-fraction, or dispersion gate whose rejection probability changes with
`SCUTTLE_ZSCORE_SEEDS`. The targeted test passed at 1, 16, 64, and 4096 seeds;
the observed six-gate z ranges were `[-1.450, 0.935]`, `[-0.759, 1.417]`,
`[-0.568, 1.049]`, and `[-1.782, 0.815]`, respectively. The 4096-seed run took
37.00 seconds on the recorded Xeon host.

`percolation_laws` times only `sample` on an open 64x64 square (`V=4,096`,
`E=8,064`); graph, laws, probability fields, configuration, and RNG are prepared
outside Criterion iteration. The short run uses 10 samples, 100 ms warm-up, and
250 ms measurement per case:

```bash
cargo bench -p cmc-rs --bench percolation_laws --no-default-features -- --test
cargo bench -p cmc-rs --bench percolation_laws --no-default-features
cargo bench -p cmc-rs --bench percolation_allocations \
  --features allocation-probe -- --test
```

Recorded 2026-09-03 on the same Xeon Gold 6148 and rustc 1.98.0 toolchain as F0.
The 2026-09-03 target contained **14 Criterion IDs** (three site, three bond,
seven mixed uniform, one heterogeneous mixed); N1 later added heterogeneous
site/bond for 16 total — see the N1 section. Throughput denominators are `V` for site,
`E` for bond, and `V + E = 12,160` for every mixed case.

| Law | `p_vertex` | `p_edge` | Point time | Point throughput |
|---|---:|---:|---:|---:|
| Site uniform | 0.01 | all | 6.974 us | 587.30 M vertices/s |
| Site uniform | 0.50 | all | 7.053 us | 580.72 M vertices/s |
| Site uniform | 0.99 | all | 7.200 us | 568.87 M vertices/s |
| Bond uniform | all | 0.01 | 13.708 us | 588.27 M edges/s |
| Bond uniform | all | 0.50 | 13.637 us | 591.32 M edges/s |
| Bond uniform | all | 0.99 | 13.869 us | 581.44 M edges/s |
| Mixed uniform | 0.01 | 0.01 | 20.643 us | 589.07 M entities/s |
| Mixed uniform | 0.50 | 0.50 | 21.097 us | 576.39 M entities/s |
| Mixed uniform | 0.99 | 0.99 | 20.685 us | 587.86 M entities/s |
| Mixed uniform | 0.01 | 0.50 | 20.917 us | 581.34 M entities/s |
| Mixed uniform | 0.99 | 0.50 | 20.592 us | 590.53 M entities/s |
| Mixed uniform | 0.50 | 0.01 | 20.543 us | 591.94 M entities/s |
| Mixed uniform | 0.50 | 0.99 | 21.084 us | 576.73 M entities/s |
| Mixed heterogeneous | alternating 0.25/0.75 | repeating 0.2/0.6/0.6 | 32.144 us | 378.30 M entities/s |

The allocation probe warms a mixed dense configuration, switches through exact
`p=0/1`, then measures 128 iterations for five laws across five topologies (25
law/topology records). Site uniform, bond uniform, mixed uniform, mixed
heterogeneous, and repeated mixed endpoint switching all reported
**0 allocations/iteration**, **0 allocated bytes/iteration**, and zero peak live
allocations for chain, square, cubic, sparse-ER, and preferential-attachment
cases. An endpoint-switch iteration performs two samples, `(0,1)` then `(1,0)`,
after both dense buffers have been materialized. The logical worst-case dense
payload is `V + E` bytes: 8,191 (chain), 12,160 (square), 15,616 (cubic),
20,480 (sparse ER), and 20,470 (preferential attachment). These figures are the
expected byte payload for the current dense backend, not a measurement of
private `Vec` capacity or a public runtime introspection API. They describe the
worst-case reusable two-mask state; a fresh pure law does not materialize the
irrelevant mask. They exclude struct headers, RNG state, probability-field
storage, topology, and allocator metadata.

F2 passed independent review and its engineering gates, but its public API
remains **Experimental** until the foundation milestone stabilizes. The
independently reviewed F3 analyzer/workspace is described below; the F4
composition that turns this substrate into scheduler observables is described
in the F4 section.

## F3 undirected component analyzer

The F3 production path is generic over borrowed `UndirectedGraphView`,
`VertexActivity`, and `EdgeActivity` capabilities. It has no site/bond/mixed mode:
a tracked vertex is exactly an active vertex, and a physical edge is active only
when the edge and both endpoints are active. Thus site laws use all-active edges,
bond laws use all-active vertices, and mixed laws compose naturally. Physical
edges are scanned once in edge-ID order; parallel edges remain distinct and a
self-loop contributes one active edge.

`ComponentSummary` records active vertex/edge counts, component count, largest
component identity and size, and raw `M2`. Scientific identity is the canonical
minimum `VertexId` in a component, never the union-find implementation root.
Equal-size largest components choose the lowest canonical identity. Empty
activity and empty borrowed graphs return zero counts, no largest identity, and
zero `M2`.

`BoundaryQuery::new` checks raw indices against the graph, sorts and deduplicates
both sets, and retains their validated vertex domain. Empty sets are valid and
never cross. Overlap means an active vertex in the intersection crosses by
itself. `analyze` answers several queries without cloning their vectors;
`analyze_with_labels` explicitly requests per-vertex canonical labels, with
inactive vertices labelled `None`. `AnalysisResult<'_>` borrows outcomes and
optional labels from `ComponentWorkspace`, so the borrow checker prevents a
subsequent analysis from invalidating a live result.

Activity and query dimensions are checked before workspace mutation. Workspace
capacity failures are typed. Because preparation is a sequence of fallible
`Vec::try_reserve_exact` calls, successful earlier reservations may remain after
a later failure, but no partial scientific output is returned. Allocators that
abort on OOM remain outside the recoverability guarantee. Equal-size graph ID
domains have no runtime provenance token; using IDs/query semantics from another
graph of the same size remains the documented F1 caller-contract violation.

Validation retains the PR #4 analyzer and its existing tests unchanged. New F3
coverage uses a separate incidence/stack flood fill rather than production UF:
full 2x2 site/bond/mixed enumeration; 7 topology families x 3 laws x 128 seeded
configurations; owned/borrowed CSR parity; and self-loop, parallel-edge,
disconnected, isolate, empty, overlap, empty/duplicate/multiple-query, labels,
stable-identity, typed vertex/edge/query errors, and grow/shrink/labels/query
workspace reuse cases. Eight targeted integration tests and the internal
stamp-wrap unit test pass.

The standalone allocation probe prepares topology, query, dense activity, and
workspace before measuring 128 sample+analyze calls. All 15 topology/law records
(chain, square, cubic, sparse ER, power law x site/bond/mixed) report **0.000
allocations/sample** and **0.000 allocated bytes/sample**, compared with the F0
reference's 7 allocations and 86,016 bytes per sample. The feature-gated,
doc-hidden audit snapshot accounts for every owned workspace vector. At
`V=4,096`, one query, and labels prepared, allocator observation and summed
capacity both equal 196,609 bytes across six allocations:

| Workspace vector | Capacity bytes |
|---|---:|
| parent | 32,768 |
| component size | 32,768 |
| canonical identity | 32,768 |
| query root stamps | 32,768 |
| query outcomes | 1 |
| optional labels | 65,536 |
| **Total** | **196,609** |

This is `48.000 bytes/V` and, for this one-query probe, `196,609 bytes/query`.
The latter is total workspace normalized by query count, not marginal query
storage: each additional query outcome currently costs one byte while root
stamps remain vertex-sized. A no-label workload does not reserve label storage.
Struct headers and allocator metadata are excluded.

The Linux-only RSS target runs one case per fresh process and reads
`/proc/self/status`; it adds no dependency or unsafe code. From before topology
construction through 8 warmups plus 128 mixed sample+analyze calls, square
(`V=4,096`, `E=8,064`) measured `VmHWM 2,096 -> 2,796 KiB` (+700 KiB), while
sparse ER (`V=4,096`, `E=16,384`) measured `2,100 -> 3,572 KiB` (+1,472 KiB).
These are whole-process lifetime high-water deltas including topology,
configuration, RNG, runtime, and allocator effects; they are not workspace
payload or portable RSS guarantees.

`percolation_components` contains 63 short-run Criterion IDs: production
sample+analyze, production analyze-only, and equivalent independent flood-fill
sample+analyze for five topologies x three laws, plus three historical PR #4
square measurements and, since F4, fifteen `static-mc` sweep+measure records. Before timing each topology/law pair, 128 identical samples
assert equality of all six summary fields and query outcomes. The production and
equivalent-reference timed paths both use `TIMED_SAMPLE_SEED`, advance exactly
`TIMED_INITIAL_SAMPLES = 1` sample before Criterion, and advance exactly one
sample per logical iteration. A separate 129-state check (the initial sample plus
128 following samples) compares every vertex and physical-edge activity flag,
proving that corresponding iteration indices produce identical configuration
sequences. Criterion may request different total iteration counts for separate
benchmark functions, but each path starts from and advances through the same
indexed sequence. Labels are disabled on both timed paths. Criterion throughput
uses physical `E`, so it reports edges/s directly; samples/s below is
`1 / point_time`. Displayed edges/s is recomputed as `E * displayed samples/s`
and rounded to `0.01 M`, with an automated arithmetic check. Parameters are
`V=4,096`, `p_site=0.5927`, and `p_bond=0.5`.

| Topology | Law | Production samples/s | Production edges/s | Equivalent reference samples/s | Reference edges/s | Speedup |
|---|---|---:|---:|---:|---:|---:|
| chain | site | 13,493 | 55.25 M | 11,480 | 47.01 M | 1.175x |
| chain | bond | 14,799 | 60.60 M | 9,475 | 38.80 M | 1.562x |
| chain | mixed | 11,258 | 46.10 M | 8,319 | 34.07 M | 1.353x |
| square | site | 9,738 | 78.53 M | 4,250 | 34.27 M | 2.291x |
| square | bond | 8,683 | 70.02 M | 3,324 | 26.80 M | 2.612x |
| square | mixed | 7,577 | 61.10 M | 2,737 | 22.07 M | 2.769x |
| cubic | site | 7,608 | 87.64 M | 3,844 | 44.28 M | 1.979x |
| cubic | bond | 6,540 | 75.34 M | 2,963 | 34.13 M | 2.207x |
| cubic | mixed | 5,812 | 66.95 M | 2,260 | 26.04 M | 2.572x |
| sparse ER | site | 3,892 | 63.77 M | 2,935 | 48.09 M | 1.326x |
| sparse ER | bond | 4,220 | 69.14 M | 2,090 | 34.24 M | 2.019x |
| sparse ER | mixed | 4,026 | 65.96 M | 2,102 | 34.44 M | 1.915x |
| power law | site | 6,597 | 108.02 M | 3,122 | 51.12 M | 2.113x |
| power law | bond | 5,521 | 90.40 M | 2,098 | 34.35 M | 2.632x |
| power law | mixed | 4,806 | 78.69 M | 2,291 | 37.51 M | 2.098x |

Production analyze-only throughput for the same fixed configurations:

| Topology | Site samples/s / edges/s | Bond samples/s / edges/s | Mixed samples/s / edges/s |
|---|---:|---:|---:|
| chain | 17,440 / 71.42 M | 18,884 / 77.33 M | 15,536 / 63.62 M |
| square | 10,698 / 86.27 M | 10,187 / 82.15 M | 9,660 / 77.90 M |
| cubic | 8,237 / 94.89 M | 7,574 / 87.25 M | 7,255 / 83.58 M |
| sparse ER | 3,882 / 63.60 M | 4,780 / 78.32 M | 4,792 / 78.51 M |
| power law | 6,868 / 112.46 M | 6,625 / 108.48 M | 5,930 / 97.10 M |

The retained PR #4 square timing reproduces the historical per-call allocation
profile on the final analyzer; it is only a **historical non-equivalent migration
baseline**. It is not used to
claim speedup. The equivalent flood-fill table above is the F3 performance
comparison. Short-run point estimates are not CI thresholds. Run the evidence
with:

```bash
cargo test -p cmc-rs --test suite percolation_components --no-default-features
cargo bench -p cmc-rs --bench percolation_components --no-default-features -- --test
cargo bench -p cmc-rs --bench percolation_components --no-default-features
cargo bench -p cmc-rs --bench percolation_allocations \
  --features allocation-probe -- --test
SCUTTLE_RSS_CASE=square SCUTTLE_RSS_LAW=mixed \
  cargo bench -p cmc-rs --bench percolation_components_rss --no-default-features
SCUTTLE_RSS_CASE=sparse-er SCUTTLE_RSS_LAW=mixed \
  cargo bench -p cmc-rs --bench percolation_components_rss --no-default-features
```

F3 passed independent review and its engineering gates but remains
**Experimental**. F4 (below) adds the observable normalization,
finite-cluster susceptibility, and Carlo.rs adapter that turn this analyzer
into production-path observables; with the F4 review landed, the static family
is the platform's first Validated row.

#### `Vec<bool>` storage evidence

On this compiler, the reviewer concern is confirmed: these `Vec<bool>` values
are allocator-observed as byte-backed, despite `capacity()` being expressed in
logical `bool` elements. Each probe performed exactly one allocation:

| Logical len/capacity | Allocator-observed bytes |
|---:|---:|
| 1 | 1 |
| 63 | 63 |
| 64 | 64 |
| 65 | 65 |
| 4096 | 4096 |
| 8064 | 8064 |
| 11,520 | 11,520 |
| 16,384 | 16,384 |

This is empirical behavior of `rustc 1.98.0`, not a promised representation of
`Vec<bool>` across toolchains. The probe remains executable so future baselines
can detect a representation change instead of assuming bit packing or byte
storage from the type name.

For every topology and mode, the post-warmup sample+analyze measurement returned
**7.000 allocations/sample**, **86,016 allocated bytes/sample**, and a peak of
**7 live allocations / 86,016 live bytes** in the serial measurement region.
These are analyzer reference costs, not targets. Both bench targets' `--test`
commands completed successfully; Criterion exercised all 14 F2 law throughput
IDs listed above. Exact raw Criterion intervals remain in `target/criterion/` on
the measurement host.

## F4 scientific observables and production adapter

F4 passed independent review and all engineering gates. It composes
the F1-F3 foundation into the final scheduler surface and deletes the
provisional PR #4 public API in one migration (no compatibility wrapper, per
the roadmap rule for never-merged APIs). The provisional API PR #4 introduced
is superseded in place on this branch; PR #4 carries the final architecture as
its delivery vehicle. The static family is **Validated** and becomes
production-ready when the branch merges into `dev`.

### Observable definitions and Results semantics

`ObservablePlan` explicitly requests scalar observables and boundary queries.
Duplicate requests collapse (one sample per name per realization), and labels
are intentionally not part of the plan: every scalar follows from
`ComponentSummary`, whose canonical lowest-ID largest component realizes the
one-component tie policy. Each boundary query records `BoundaryCrossing0`,
`BoundaryCrossing1`, ... in plan order as a 0/1 indicator.

| Results name | Definition | Notes |
|---|---|---|
| `ActiveVertexCount` | occupied sites in the configuration | mixed mode counts sites only |
| `ActiveEdgeCount` | physical edges that are open under the law and have both endpoints occupied | site mode: edges between two occupied sites; bond mode: open bonds; mixed: both conditions |
| `ComponentCount` | tracked components among occupied sites | bond mode counts every isolate |
| `LargestSize` | `S_max`, `0` when nothing is occupied | |
| `GiantFraction` | `S_max / V`, denominator fixed to total topology vertices `V` | undefined and absent at `V = 0` |
| `RawSecondMoment` | `M2 = sum_i s_i^2`, no normalization, largest included | u128 accumulation |
| `FiniteClusterSusceptibility` | per sample, exclude one largest component (ties: canonical lowest identity), then `chi = (M2 - S_max^2) / (V_active - S_max)` | exclusion before averaging |
| `FiniteClusterSusceptibilityDefined` | 0/1 indicator recorded for every sample | 0 exactly when the chi denominator is zero |

Undefined samples never write the chi scalar: the `FiniteClusterSusceptibility`
row in `results.json` has an `n_bins` counting only defined samples, while the
indicator row counts every sample. This makes the defined fraction explicit
instead of silently emitting zero or NaN. The per-sample exclusion cannot be
reconstructed from separately aggregated `LargestSize`/`RawSecondMoment` means
(unit-tested counterexample against the ratio of means).

Moment observables convert u128 squares up to `V^2` into f64. Plans whose
topology could exceed the exact f64 integer range (`V^2 > 2^53`, i.e.
`V > 94,906,265`) are rejected at construction with
`ObservablePlanError::InexactMomentConversion`, for `RawSecondMoment` and for
the susceptibility numerator, so no silent precision loss is possible.

### Adapter and parameter schema

`StaticPercolationMC` composes an owned validated `CsrLattice`, the uniform
`StaticLaw` enum (site/bond/mixed, used only at this boundary),
`StaticConfiguration`, a prepared `ComponentWorkspace`, and the
`ObservablePlan`. `sweep()` only resamples occupancy; `measure()` analyzes
with the preallocated workspace and records exactly the planned scalars.
Configure `thermalization_sweeps = 0` (i.i.d. samples; the adapter cannot see
scheduler configuration to enforce it). The `FromParams` schema is:

- `mode`: `site` (default) / `bond` / `site-bond`; pure modes take `p`, mixed
  takes `p_site` + `p_bond`, and the wrong key for the mode is a typed error;
- lattice parameters reuse the standard builders (`chain`/`square`/`cubic`/
  `hypercubic`/...), `pbc` defaults to `false`;
- `observables`: comma-separated selection defaulting to the full plan;
- `spanning_from`/`spanning_to`: explicit comma-separated site lists, always
  both together; defaults are square left/right columns and chain endpoints,
  and any other parameter-built topology is rejected loudly (construct the
  adapter directly with a `BoundaryQuery` for arbitrary sets).

`validate_params` accepts exactly what `from_params` accepts. Borrowed graphs
and custom laws bypass the adapter through `analyze`/`analyze_with_labels` or
compose a runtime with `Run::from_parts()`.

### Validation

Every applicable PR #4 gate was migrated to the final API in
`tests/physics/percolation.rs` (14 tests, 2 `#[ignore]` long) and
`tests/physics/percolation_zscore.rs` (6 z-score tests): 2x2 exhaustive
enumeration vs hand-derived moments, site spanning polynomial, mixed closed
form `2 p_s^2 p_b - p_s^4 p_b^2`, mixed-to-pure reduction identities, chain
closed forms, union-find vs independent flood-fill across seven topology
families and the seeded random graph, crossing monotonicity (pure and mixed),
scheduler 2x2 exact moments and chain closed forms, fixed-seed bitwise
reproduction, `p = 0`/`p = 1` occupation extremes, the ignored 32x32 bond `p_c = 1/2` crossing and 16^3 cubic
critical-bracket long tests, and `BoundaryCrossing0`/`LargestSize` z-score
gates for site/bond/mixed (16 default seeds, `SCUTTLE_ZSCORE_SEEDS`-scalable,
verified at 64 seeds). New F4 unit coverage (10 tests in `observable.rs` +
`carlo.rs`): tie exclusion, interleaved duplicate collapse, no-active and
`V = 0`, single-giant-component undefined chi with indicator counts under the
scheduler, per-sample exclusion vs ratio-of-means counterexample, exact-f64
conversion limit, JSON schema, and the `FromParams` rejection matrix.

### Steady-state allocation gate

The `percolation_allocations` probe (feature `allocation-probe`) adds
`report_adapter_case`: it warms 8 full `sweep()+measure()` cycles, then
measures 128 more through `StaticPercolationMC` and a Carlo.rs `Context` with
a binsize strictly larger than warmup plus measurement, so no accumulator bin
completes inside the measured window. All **15 topology/law records** (chain,
square, cubic, sparse ER, power law x site/bond/mixed; 9 observables each)
report **0.000 allocations/sample**, **0.000 allocated bytes/sample**, and
zero peak live allocations/bytes:

```text
PERCOLATION_ADAPTER_ALLOC ... allocations_per_sample=0.000
allocated_bytes_per_sample=0.000 peak_live_allocations=0 peak_live_bytes=0
```

This gate required one behavior-preserving Carlo.rs change: with no
measurement namespace, `Context::measure*` previously allocated a fresh
`String` for every observable name on every call; `qualified_measurement_name`
now returns `Cow::Borrowed` in that case (`Cow::Owned(format!(...))` when a
parallel-tempering namespace is active, exactly as before). Recording nine
scalars per sample is otherwise allocation-free.

### Criterion data

`percolation_components` adds 15 `static-mc-sweep+measure` IDs (full adapter:
law sampling through the adapter, workspace analysis, and nine scalar
recordings) next to the F3 `sample+analyze` IDs. Recorded 2026-10-02 on the
same Xeon Gold 6148 / rustc 1.98.0 host, short-run configuration (10 samples,
100 ms warm-up, 250 ms measurement), **while the host carried an unrelated
sustained load of ~45/72 cores**, so absolute point estimates are inflated
relative to the 2026-09-03 F3 table by roughly 8-15% across all IDs —
including `analyze-only`, a path F4 does not touch and did not modify (the
F1-F3 core sources are byte-identical since the F3 record; the diff touches
only exports, the two new F4 files, benches, and tests). The load-independent
same-run comparison is therefore the meaningful F4 evidence: the adapter
composition costs the same as the bare core within contention noise (median
about -2%, scatter roughly +/-7%, one -23% power-law/site outlier where the
adapter measured faster).

| Topology | Law | Core sample+analyze (us) | Adapter sweep+measure (us) | Adapter samples/s | Adapter edges/s | Adapter vs core |
|---|---|---:|---:|---:|---:|---:|
| chain | site | 86 | 87 | 11,464 | 46.94 M | +1.9% |
| chain | bond | 85 | 79 | 12,695 | 51.99 M | -7.3% |
| chain | mixed | 100 | 101 | 9,940 | 40.71 M | +1.0% |
| square | site | 114 | 118 | 8,470 | 68.30 M | +3.5% |
| square | bond | 133 | 132 | 7,588 | 61.19 M | -1.2% |
| square | mixed | 155 | 145 | 6,891 | 55.57 M | -6.4% |
| cubic | site | 149 | 146 | 6,862 | 79.06 M | -2.4% |
| cubic | bond | 171 | 165 | 6,056 | 69.76 M | -3.2% |
| cubic | mixed | 191 | 188 | 5,308 | 61.15 M | -1.5% |
| sparse ER | site | 278 | 272 | 3,681 | 60.31 M | -2.4% |
| sparse ER | bond | 256 | 249 | 4,016 | 65.79 M | -2.8% |
| sparse ER | mixed | 274 | 265 | 3,774 | 61.83 M | -3.4% |
| power law | site | 218 | 168 | 5,966 | 97.68 M | -23.2% |
| power law | bond | 179 | 190 | 5,250 | 85.97 M | +6.5% |
| power law | mixed | 239 | 220 | 4,540 | 74.34 M | -8.0% |

For the F3 no-regression gate, the F3 core paths (`sample+analyze`,
`analyze-only`, `equivalent-reference`, laws) are unchanged code and their
same-run parity/sequence assertions still pass; the cross-day point-estimate
deltas under the current load are environmental, as evidenced by the equal
inflation of the untouched `analyze-only` IDs. Re-run on a quiet host to
reproduce absolute numbers:

```bash
cargo bench -p cmc-rs --bench percolation_components --no-default-features -- --test
cargo bench -p cmc-rs --bench percolation_components --no-default-features
cargo bench -p cmc-rs --bench percolation_allocations \
  --features allocation-probe -- --test
```

### Remaining F4 limitations

At F4 the adapter was uniform-only; N1 has since added the heterogeneous
`StaticLaw` variants on the same adapter (see the N1 section below), so the
uniform-only restriction no longer applies. The remaining limitations are
unchanged: the parameter path owns its lattice (borrowed graphs use
the core APIs); one boundary-query normalization per adapter plan; scalar
results only (vector observables need C0); no wrapping/winding (W2), no
parameter scans (Z1), no directed/process/continuum/FK/attack families; chi
tie conventions other than single canonical-largest exclusion are separate
named policies that do not exist.

## N1 heterogeneous Bernoulli laws

N1 promotes the F2 heterogeneous probability substrate to the production
adapter boundary. The analyzer, workspace, and observable plan are untouched
(the roadmap requires this); only the owned adapter surface grows.

### API

`StaticLaw` gains three owned heterogeneous variants:
`SiteHeterogeneous { probabilities: Vec<Probability> }`,
`BondHeterogeneous { probabilities: Vec<Probability> }`, and
`MixedHeterogeneous { site_probabilities, bond_probabilities }`. The
constructors `StaticLaw::site_heterogeneous` / `bond_heterogeneous` /
`mixed_heterogeneous` take `impl Into<Vec<Probability>>` (owned vectors or
borrowed slices) and cannot fail: every value is already a validated
`Probability`. Field lengths are checked against the adapter topology in
`StaticPercolationMC::new` — a mismatch returns the new typed
`StaticPercolationError::LawDomain(SamplingError::ProbabilityFieldLength {..})`,
reusing the F2 vocabulary — so a constructed adapter's `sweep()` can never hit
a domain error. `sweep()` borrows the fields (`ProbabilityField::Borrowed`),
so per-sweep sampling remains allocation-free. The `FromParams` schema stays
uniform-only (roadmap §9 N1 adds no probability-list input format);
heterogeneous adapters are scheduled through `Run::from_parts()` with a seeded
`Context`, exactly like other non-`FromParams` models. `StaticLaw::law()` now
returns a reference (the enum is no longer `Copy` because it owns vectors).

### Validation gates (roadmap §9 N1)

All gates live in `tests/physics/percolation_heterogeneous.rs` (7 tests) and
`tests/physics/percolation_heterogeneous_zscore.rs` (7 z-score tests), plus
3 adapter unit tests in `src/percolation/carlo.rs`:

- **Exact non-identically-distributed enumeration.** The 2x2 open square is
  enumerated with per-element weights `prod_i p_i^{x_i} (1-p_i)^{1-x_i}`
  (four distinct values per domain); scheduler statistics for
  `BoundaryCrossing0`, `LargestSize`, and `RawSecondMoment` match within
  |z| < 4 for site, bond, and mixed, over 4 seeds each at 200k sweeps. The
  weight products are checked to sum to one. Hand-derived heterogeneous chain
  closed forms `P = prod p_i` (site), `prod q_j` (bond), and the product of
  both (mixed) match through the full scheduler stack at L = 6 within |z| < 4.
  The 3x3 mixed enumeration reference of the z-score tests was independently
  reproduced during independent N1 review by a from-scratch Python
  flood-fill enumeration (different code and summation order) agreeing to
  ~1e-13-1e-14 relative difference.
- **Uniform field reduction, bitwise.** A heterogeneous field of identical
  `p` values consumes RNG in the same order as the uniform law (and consumes
  none at exact endpoints), so identical seeds reproduce identical `Results`
  means bit-for-bit across six observables — for site, bond, and mixed laws,
  at `p in {0.0, 0.3, 0.5, 1.0}`.
- **Mixed reduction identities, bitwise.** Heterogeneous mixed with an all-one
  edge field equals the heterogeneous site law, and with an all-one site field
  equals the heterogeneous bond law, bitwise at identical seeds.
- **Endpoint exactness.** All-0/1 fields produce deterministic configurations;
  hand-checked chain cases pin `ActiveVertexCount`, `ActiveEdgeCount`,
  `ComponentCount`, `LargestSize`, `RawSecondMoment`, and crossing exactly.
- **Length and invalid-probability rejection.** Wrong-length site/bond/mixed
  fields are typed `LawDomain` rejections at the adapter constructor (unit
  tests assert the exact payload); no unvalidated f64 can reach the
  heterogeneous constructors (`Probability::new` is the only entry).
- **Scheduler z-scores.** 16 default seeds (SCUTTLE_ZSCORE_SEEDS-scalable,
  verified at 64): `BoundaryCrossing0` and `LargestSize` against full
  heterogeneous enumeration on the 3x3 square for site (2^9), bond (2^12),
  and mixed (2^21, shared across its two tests), with per-seed |z| < 4,
  |z̄| < 1.5, and no one-sided bias; plus a Poisson-binomial activity gate
  (site `ActiveVertexCount` mean `sum p_i`, bond `ActiveEdgeCount` mean
  `sum q_j`, mixed `ActiveEdgeCount` mean `sum_j q_j p_u(j) p_v(j)`, |z| < 4
  per seed). The heterogeneous scheduler runs use the labeled seed base
  `0x4845_5445_0000` ("HETE").

### Uniform versus heterogeneous sampling benchmark

`percolation_laws` adds two Criterion IDs (heterogeneous site and bond next to
the existing heterogeneous mixed), for **16 IDs** total. Recorded 2026-10-08
on the same Xeon Gold 6148 / rustc 1.98.0 host family as the 2026-09-03 F2
record, short-run configuration (10 samples, 100 ms warm-up, 250 ms
measurement), open 64x64 square (`V=4,096`, `E=8,064`); heterogeneous fields
alternate mid-range values (0.25/0.75 per vertex, 0.2/0.6 per edge — no
endpoint skipping). Point estimates:

| Family | Uniform `p = 0.5` | Heterogeneous | Heterogeneous time delta |
|---|---:|---:|---:|
| Site | 7.009 us / 584.4 M vertices/s | 10.324 us / 396.7 M vertices/s | +47.3% |
| Bond | 13.543 us / 595.4 M edges/s | 20.453 us / 394.3 M edges/s | +51.0% |
| Mixed | 20.571 us / 591.1 M entities/s | 30.952 us / 392.9 M entities/s | +50.5% |

The ~50% cost is the per-element probability load and compare in the dense
heterogeneous loop versus the hoisted uniform scalar; the F2 record measured
the same effect (378.30 M entities/s for heterogeneous mixed on 2026-09-03).
These are point estimates, not CI thresholds. Run with:

```bash
cargo bench -p cmc-rs --bench percolation_laws --no-default-features -- --test
cargo bench -p cmc-rs --bench percolation_laws --no-default-features
```

### Heterogeneous adapter allocation gate

`percolation_allocations` adds `report_heterogeneous_adapter_case`: after 8
warmup `sweep()+measure()` cycles through `StaticPercolationMC` with
heterogeneous fields, 128 measured cycles report **0.000 allocations/sample,
0.000 allocated bytes/sample, and zero peak live allocations/bytes** for all
**15 topology/law records** (chain, square, cubic, sparse ER, power law x
heterogeneous site/bond/mixed; 9 observables each), exactly matching the
uniform adapter gate:

```text
PERCOLATION_ADAPTER_ALLOC ... law=heterogeneous-site|bond|mixed
allocations_per_sample=0.000 allocated_bytes_per_sample=0.000
peak_live_allocations=0 peak_live_bytes=0
```

### Remaining N1 limitations

The analyzer is untouched (heterogeneity enters only through activity
sampling, as the roadmap requires). The `FromParams` schema is uniform-only;
heterogeneous adapters are constructed directly and scheduled with
`Run::from_parts()`. Per-element probabilities must be validated upstream
(`Probability::new`); the adapter validates only field lengths. Network
facades over heterogeneous site fields (random failure, original-degree
observables — N2) and targeted attacks (N3) are Not implemented.

## W1 periodic embedding capability

Status: **Experimental**, independently reviewed (roadmap §16 W1; engineering gates
and tests below are green, independent review has not happened). W1 is the
data layer only — the W2 winding analyzer is Not implemented.

### API and semantics

- `topology::PeriodicEmbedding` is a read-only capability next to
  `GraphView`/`UndirectedGraphView`: `periodic_dimensions()` returns the
  embedding dimension `D >= 1`, and
  `edge_displacement(edge, from) -> Option<DirectedDisplacement>` returns the
  integer fundamental-cell displacement of one **directed incidence** — the
  pair (physical edge, endpoint) — or `None` when `from` is not an endpoint.
  The key is never the endpoint pair or a bond label, so parallel edges
  between the same pair keep independent displacements.
- `DirectedDisplacement<'_>` is a lazy signed view over the stored
  source-to-target vector; querying the same physical edge from the other
  endpoint negates it exactly, so antisymmetry
  (`d(v->u) = -d(u->v)`) is structural, not a validated copy. Accessors:
  `iter()`, `axis(a)`, `dimensions()`, `is_zero()`.
- `LatticeEmbedding` owns one vector per physical edge (flat, edge-major,
  `D` `i32` components, `16 + 4D` bytes per edge) plus the endpoint pair it
  was validated against. `try_from_edge_displacements` rejects, before
  storing anything: zero periodic dimensions, wrong table length, self-loop
  edges (the two directed incidences of a self-loop would need opposite
  displacements but share one `(edge, endpoint)` key), `i32::MIN` components
  (unnegatable), and unresolvable capacity. It does **not** re-verify cocycle
  consistency: that is a producer property, pinned for the built-in builders
  by the plaquette and winding tests below and exploited by the W2 analyzer
  rather than re-validated per construction.
- Open boundaries have no embedding at all: every `*_with_embedding` builder
  rejects `pbc = false` with `EmbeddingError::OpenBoundaries`. A non-periodic
  graph is represented by not having the capability, so a wrapping query on a
  non-embedded graph is an explicit consumer-side error, never a vacuous
  "no winding" answer.
- No-inference principle (the roadmap W1 constraint): displacements are
  recorded by the builder at the only stage where wrap information is
  unambiguous — while bonds are emitted. `BondType` and vertex IDs are never
  consulted for geometry: on a 2-wide torus the two parallel bonds between a
  pair of sites carry different displacements (`0` and `-1`), so endpoint
  geometry alone cannot identify them.
- `BorrowedUndirectedCsr` intentionally does not implement the capability: a
  CSR topology carries no displacement payload and post-hoc reconstruction
  would be the forbidden inference. An embedding is separate data supplied by
  whoever knows the geometry, beside the borrowed topology.

### Builder coverage

- `build_chain_with_embedding`, `build_square_with_embedding`,
  `build_hypercubic_with_embedding` (all `pbc = true`), via a shared
  `hypercubic_bonds` bond loop that tracks wraps; and
  `build_triangular_with_embedding` via `triangular_bonds`, whose diagonal
  bonds lift `(x, y)` to `(x + 1, y + 1)`. Each variant returns exactly the
  same lattice as its plain builder (pinned by test) plus the embedding.
- Honeycomb and kagome embeddings are **Not implemented**: those builders
  assemble symmetric adjacency lists, and `CsrLattice::from_adjacency` drops
  the per-bond construction direction while deduplicating physical edges. On
  small cells (honeycomb `Lx = 2` double horizontal bonds; kagome `Lx = 2`
  collapsed covering bonds) distinct displacements map to indistinguishable
  physical bonds, and recovering them would be guessing. Attaching embeddings
  requires restructuring those builders to emit bonds directly with recorded
  displacements — deliberately not done in W1.

### Validation gates

All in `CMC.rs` unit tests (`lattice::graph::tests`, `topology::embedding::tests`):

- **Hand-derived 2x2 torus:** every one of the 8 torus bonds is asserted in
  both directions against a hand-written displacement table; each of the four
  vertex pairs carries one in-cell and one wrapping bond.
- **Independent geometry re-derivation (property test):** for 8 dimension
  sets (1D L=2,3,5; 2x2, 3x3, 4x3; 3D 2x3x4 and 3x1x2), every bond's
  advancing axis and wrap flag are recomputed from the builder's documented
  mixed-radix coordinates and compared componentwise, including exact
  negation from the other endpoint and `is_zero`.
- **Parallel edges:** the L=2 ring's two parallel bonds carry `0` and `-1`;
  the 2x2 square pair `{0,1}` carries `(0,0)` and `(-1,0)`; the 2x2
  triangular double diagonals between sites 1 and 2 carry `(1,0)` and
  `(0,-1)`.
- **OBC semantics:** `pbc = false` is rejected with
  `EmbeddingError::OpenBoundaries` by the chain, square, and hypercubic
  embedding builders.
- **Cocycle consistency (Z^D embedding):** on the 3x3 square torus, all 9
  elementary plaquettes accumulate `(0,0)`; on the 3x3 triangular torus, all
  18 elementary triangles (both orientations) accumulate `(0,0)`.
- **Winding of fundamental cycles:** axis-ring walks on L=2/3/4 squares, a
  5-ring, and the 2x3x2 cubic torus accumulate exactly `+1` on the walked
  axis and `0` elsewhere, taking the bond emitted from each site (the L=2
  parallel-bond direction).
- **Table construction rejections:** zero dimensions, wrong length,
  self-loops, `i32::MIN`, non-endpoint queries, and generic-code acceptance
  of the `PeriodicEmbedding` trait; plus bit-identical lattice reproduction
  by every embedding builder.

### W2 preview (not implemented)

The intended consumer is one pass over physical edges:
`embedding.edge_displacement(edge, from)` next to
`UndirectedGraphView::edge_endpoints`, running union-find with integer
potentials; a closed component whose accumulated potential mismatches is a
winding generator. The `(edge, endpoint)` lookup key exists precisely so
that pass never needs `BondType` or endpoint-pair geometry. No part of W2 is
implemented in this stage.

### Limitations

No benchmarks were added: lookups are O(1) table reads with no allocation
(the storage note above is the capacity claim), and there is no analyzer to
bench. The cocycle property is a builder contract verified by tests, not
re-checked by `try_from_edge_displacements`. An embedding detached from its
lattice has no provenance guard — the same caller contract as the F1 view
IDs (both must stay unchanged together).
