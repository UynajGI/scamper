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
| Static site, bond, and mixed Bernoulli percolation on arbitrary owned undirected `CsrLattice` | Experimental | PR #4 reference implementation is extensively scientifically validated, but F1-F4 production layering, borrowed topology, reusable workspace, final observables, and stable API do not exist. |
| Heterogeneous occupation probabilities | Not implemented | Only one uniform probability per sampled entity kind exists. |
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
| Boundary crossing between explicit vertex sets | Experimental | Implemented and validated for the PR #4 reference path; square/chain have defaults, other topologies require explicit sets. |
| Periodic wrapping/winding | Not implemented | No incidence displacement/cocycle or winding analyzer. Crossing is not wrapping. |
| Newman-Ziff parameter scans | Not implemented | No incremental activation curve or binomial reweighting. |

No percolation family is marked **Validated** at F0. The existing static family
has strong validation assets, but the roadmap reserves production-ready status
until F4 and the uniform Definition of Done are complete.

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
It provides `PercolationMode`, public `Vec<bool>` fields in `OccupancyState`,
`cluster_stats(&CsrLattice, ...)`, `UnionFind`, and scheduler-ready
`PercolationMC`. These names, storage choices, constructor/parameter schema,
and observable names are provisional and are not compatibility commitments.

Reference behavior:

- each sweep independently redraws uniform site/bond occupancy;
- `cluster_stats` allocates fresh union-find and boundary scratch on every call;
- physical undirected edges are scanned once, then tracked vertices are scanned;
- site/mixed omit closed sites; bond mode counts every isolated site;
- outputs are `MaxCluster`, raw `SecondMoment`, `NClusters`, `Spanning`, plus
  `Occupied` for pure modes or `OccupiedSites`/`OccupiedBonds` for mixed mode;
- crossing defaults are open square left/right columns and open chain endpoints;
  all other parameter-built topologies require explicit site sets.

Validation assets retained for migration are documented in
[VALIDATION.md](VALIDATION.md): exact full enumeration and hand-derived 1D/2x2
forms, an independent flood-fill reference, mixed-to-pure reductions,
cross-family and random-graph parity, scheduler tests, fixed-seed
reproducibility, 16-seed z-score gates, and ignored square/cubic critical-region
checks.

Current limitations include owned `CsrLattice` coupling, public storage layout,
uniform probabilities only, fresh per-analysis allocations, scalar results,
one crossing query per analysis, no component identity/labels, no giant fraction
or finite-cluster susceptibility, no wrapping, and no directed, process,
continuum, FK, attack, or Newman-Ziff implementation.

## Performance baseline procedure

`percolation_baseline` measures the unmodified PR #4 `resample` plus
`cluster_stats` throughput path. It does not enable, depend on, or link
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

The standalone probe prints `VEC_BOOL_PROBE` and `PERCOLATION_ALLOC` records.
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
commands completed successfully; Criterion exercised all 15 throughput IDs.
Exact raw Criterion intervals remain in `target/criterion/` on the measurement
host.
