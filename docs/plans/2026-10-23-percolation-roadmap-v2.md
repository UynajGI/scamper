# CMC.rs Percolation Platform — Roadmap v2 (New Domains)

> Successor to `2026-09-03-percolation-platform-roadmap.md` (all 27 of its
> stages delivered; M1–M6 complete). The v1 philosophy document remains
> normative; this roadmap extends the same staged, literature-anchored
> method to the deferred tails and to percolation domains v1 never scoped.

**Status:** Draft for review
**Established:** 2026-10-23
**Base:** `feat/cmc-percolation` (PR #4 vehicle); Carlo.rs runtime + CMC.rs
percolation engine as merged.

## 0. Principles carried over, two added

1. Every v1 rule carries over (separate laws, distinct result types,
   independent references, definitional anchors where no paper exists,
   DoD per PR, Experimental-until-reviewed).
2. **New — anchor-first**: a stage may not enter implementation until its
   anchor PDF is archived in `docs/ref/` and the specific quantity to
   reproduce is named in the stage text (v1 did this de facto; v2 makes
   it blocking).
3. **New — residue-first**: the v1 audit residue (§M6 of v1) has priority
   over new domains; no new track starts while a same-track residue item
   remains open, unless the residue is explicitly re-planned here.

## 1. Dependency graph

```text
v1 residue (adaptive attack, hetNZ, honeycomb/kagome, continuum winding)
├── T1 adaptive residual-degree attack (N3 tail)
├── T2 heterogeneous Newman-Ziff reweighting (Z1 tail)
├── T3 honeycomb/kagome embeddings (W1 tail)
└── T4 continuum winding contract (W-track tail over P1)

F first-passage percolation
├── F1 edge-weight Dijkstra/heap layer on the topology substrate
│   └── F2 FPP observables + fluctuation exponents (ν, χ)
│       └── F3 shape theorem brackets / geodesic aging (research-grade)

L long-range percolation
├── L1 power-law edge-universe law (p ~ β/d^s on Z^d) on G1
│   ├── L2 chemical-distance scaling regimes (Biskup s-regimes)
│   └── L3 graph-diameter / small-world transitions

C correlated percolation
├── C1 power-law correlated activity fields (Schrenk generator)
│   └── C2 shifted thresholds + compactness (Pruessner tables)
│       └── C3 Weinrib–Halperin criterion map (dimension a vs d)

Q generic-q random cluster (the §15 v1 tail)
├── Q1 Sweeny edge-update sampler (Geyer's efficient variants)
│   ├── Q2 Chayes–Machta for integer q (bridges R3)
│   └── Q3 phase-diagram + latent heat across q (research-grade)
```

## 2. Residue track T (priority)

### T1: Adaptive residual-degree attack
The roadmap v1 §9 N3 explicitly deferred "adaptive residual-degree attack
(a different process)". Scope: removal order recomputed from **residual**
degree after each removal; one trajectory = one C0 vector sample on the
N3 adapter pattern. Gates: each step vs naive recomputation; monotone
active count; tie policy explicit; per-step cost model documented (lazy
priority-queue decrements vs full resort). **Anchor (archived
2026-10-22): Holme, Kim, Yoon, Han, *Attack vulnerability of complex
networks*, PRE 65, 056109 (2002), arXiv:cond-mat/0202410.** Quantity to
reproduce (statistical bracket, never exact claim): on BA scale-free
networks (their Fig. 7e configuration N=10^4, m0=5, m=3, ⟨k⟩≈6), the
adaptive recalculated-degree attack (their "RD") suppresses the
largest-component fraction S(N_rm/N) below the initial-degree attack
("ID", the existing N3 order) through the late stage, both crossing zero
earlier than random failure; RD ≈ ID in the early stage. Secondary
context: the N3 Cohen κ-crossover machinery applied to the adaptive
order. Their WS-model crossover (RD least harmful near N_rm/N ≈ 0.25)
is noted but out of scope (no WS generator on the platform).

### T2: Heterogeneous Newman-Ziff reweighting
Z1's binomial reweighting assumes uniform p. Scope: per-vertex /
per-edge probability fields in the incremental scan, with Poisson-binomial
(the N1 machinery) reweighting. Gates: prefix parity as v1; reweighted
estimates match direct heterogeneous Bernoulli (N1 path) within CI on
small graphs and at scale; benchmark advantage measured before merge.

### T3: Honeycomb/kagome embeddings
W1 v1 left these Not implemented (adjacency drops per-bond direction).
Scope: extend the bond emission of the two builders to carry directed
displacements; hand-derived plaquette/cocycle gates as W1's. Delivers the
winding observables on two standard lattice families.

### T4: Continuum winding contract
W-track v1 says continuum wrapping may reuse the W2 contract "only after
their own geometry and physical configurations are defined". Scope: a
displacement source for P1's PairGraph (minimum-image offsets already
exist in `OrthorhombicCell::displacement`) feeding `analyze_winding`;
validate on periodic disks spanning via unwrapped extent (P2's BFS
unfolding is the independent reference).

## 3. F track: first-passage percolation

**Anchor (archived): Auffinger–Damron–Hanson, *50 years of first passage
percolation*, AMS ULS 68 (2017), arXiv:1511.03262** — the survey's
Tables/appendices carry the standard exponents and bounds.

### F1: Weighted-shortest-path layer
Dijkstra with a binary heap over the topology substrate (generic over
`UndirectedGraphView`), reusing E1's `ProcessQueue`/counter discipline:
edge weights from `Bond.weight`, one-pass single-source, reusable
workspace, zero steady-state allocation. Gates: every dist vs Bellman–Ford
reference on small graphs (exact); hand cases incl. zero-weight edges and
parallel edges (min wins); typed non-finite/negative rejections.

### F2: FPP observables
Passage time to distance r, fluctuation χ(r) = Var(T_r), corridor
geometry on Z^2. One realization = one C0 vector sample (time profile
along the boundary). Gates: exact small-graph enumeration (all edge-weight
assignments on a 4x4 corridor at 2-level weights, T_r vs brute-force
Dijkstra); χ(r) growth bracket consistent with 1/2 ≤ β/(2ν) regime claims
(statistical brackets per v1 W2/E-track convention); the ADH survey's
collected bounds as context, never exact claims.

### F3: Shape theorem and geodesics (research-grade)
The asymptotic shape μ(x); geodesic counting/aging. Requires F2 at scale
plus provenance for weight laws. Gates defined at design time with a
domain physicist; anchor-first.

## 4. L track: long-range percolation

**Anchor (archived): Biskup, *On the scaling of the chemical distance in
long-range percolation*, Ann. Probab. 32, 2938 (2004), arXiv:math/0304418**
— regimes s ≤ d, d < s < 2d, s ≥ 2d with their distance scalings.

### L1: Power-law edge law on G1
Edges {(x,y): |x−y|^{-s} ≤ draw} over Z^d with explicit s and boundary
handling; one realization = one G1-driven incremental construction with
occupation-count = nothing (edge universe is infinite) — instead the
process coordinate is the **added-edge count within a box**. Gates:
degree-distribution moments vs the LRP measure on small boxes; edge
independence; bitwise draw contract.

### L2: Chemical-distance regimes
Graph distance on the LRP cluster vs |x−y| in each Biskup regime:
linear / stretched-exp / log scalings as statistical brackets on
1D boxes. C0 vector output = distance profiles.

### L3: Diameter transitions
Box graph diameter vs s across the regime boundaries.

## 5. C track: correlated percolation

**Anchors (archived): Schrenk et al., PRE 88, 052102 (2013),
arXiv:1309.2994 (generator + exponents); Pruessner et al., PRE 97,
062141 (2018), arXiv:1708.02296 (square/cubic thresholds, fractal
dimensions); Saberi review EPJ B 89, 253 (2016), arXiv:1404.5325.**

### C1: Correlated activity fields
Power-law correlated Gaussian copula → thresholded occupation (the
Schrenk generator: Fourier filtering with exponent a). A law over
`StaticConfiguration`, not a new analyzer — the F3 engine is reused
unchanged. Gates: two-point correlation function of the field matches the
prescribed power law within CI; a = 0 reduces bitwise to iid Bernoulli.

### C2: Shifted thresholds
p_c(a) on square/cubic vs Pruessner's tables (statistical brackets);
cluster compactness (fractal dimension drift) observable.

### C3: Weinrib–Halperin map
The extended Harris criterion a·ν_iid = 2 — verify the correlation-length
exponent flips (ν = 1/a) only above the criterion, brackets across (a,d).

## 6. Q track: generic-q random cluster

**Anchor (archived): Geyer, *Efficient simulation of the random-cluster
model*, J. Stat. Mech. (2015) P08010, arXiv:1307.6647** — Sweeny
implementations compared (single-bond, ACES/union-find variants).

### Q1: Sweeny sampler
Single-bond heat-bath updates at arbitrary real q on the shared F1
topology: add/remove edge with the q-cluster-weight ratio; acceptance via
the connected check (the union-find connectivity query cost is THE
bottleneck — Geyer's whole point). Gates: exhaustive 2^E enumeration on
small graphs against the exact FK weight `p^|A| (1-p)^{E-|A|} q^{k(A)}/Z`
(computed with own union-find k(A); Z by full sum); q=1 reduces bitwise
to R2; q=2 marginal statistics agree with R3's exposed draws; detailed
balance argued; acceptance/connectivity cost benchmarked against Geyer's
reported characteristics.

### Q2: Chayes–Machta for integer q
Activate/put-to-sleep + per-component color resampling; agrees with Q1
statistically at q=2 and with R3's ES chain. Bridges the two algorithm
families on one platform.

### Q3: Phase diagram across q (research-grade)
Potts latent heat / discontinuity line on Z^2 across q ∈ [1, 4+] using
Q1; the exact 2D Potts critical curve (duality) as the anchor — anchor
requires archiving the 2D FK/Potts duality reference at T-time
(anchor-first rule).

## 7. Milestones

| Milestone | Content | Gate to proceed |
|---|---|---|
| N1 (residue) | T1–T4 | All four delivered or explicitly re-planned |
| N2 (dynamics) | F1–F3 | F2 validated before F3 |
| N3 (long-range + correlated) | L1–L3, C1–C3 | L1 and C1 laws validated first |
| N4 (generic q) | Q1–Q3 | Q1 exhaustive gates before any curve output |
| N5 (audit) | v2 M6-equivalent | Matrix audit as v1's |

## 8. Stage status

All stages Not started. The first PR should be T2 (heterogeneous NZ) —
smallest scope, exercises the anchor-first rule end to end — or T1 if an
acceptable adaptive-attack anchor is identified first.

## 9. Archived anchors for v2 (already in docs/ref/)

| Stage family | Paper | File |
|---|---|---|
| F track | Auffinger–Damron–Hanson 2017 | `Auffinger_2017_50yrs_first_passage_percolation.pdf` |
| L track | Biskup 2004 | `Biskup_2003_long_range_chemical_distance.pdf` |
| Q track | Geyer 2015 | `Geyer_2015_sweeny_random_cluster.pdf` |
| C track | Schrenk 2013 / Pruessner 2018 / Saberi 2016 | three files under `docs/ref/` |
| T1 adaptive attack | Holme et al. 2002 (archived 2026-10-22) | `Holme_2002_attack_vulnerability.pdf` |

Ten anchor papers across v1+v2; every physics claim keeps citing its
archived source with the agreement numbers recorded in PERCOLATION.md.
