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
`UndirectedGraphView`): one-pass single-source, reusable workspace, zero
steady-state allocation. (Delivered with its own id-tied `BinaryHeap`
rather than E1's ProcessQueue — the id tie gives deterministic
equal-distance settling and 16-byte heap entries; documented in the
module rustdoc.) Edge weights are per-call slices over dense physical
edges from `Bond.weight`; hand cases incl. zero-weight edges and parallel
edges (min wins); typed non-finite/negative rejections. Gates: every
dist vs Bellman–Ford reference on small graphs (exact); hand cases incl.
zero-weight edges and parallel edges (min wins); typed non-finite/
negative rejections.

### F2: FPP observables
Passage time to distance r, fluctuation χ(r) = Var(T_r), corridor
geometry on Z^2. One realization = one C0 vector sample (time profile
along the boundary). Gates: exact small-graph enumeration (all edge-weight
assignments on a 4x4 corridor at 2-level weights, T_r vs brute-force
Dijkstra). **Anchor quantities (ADH survey, archived
`docs/ref/md/Auffinger_2017_50yrs_first_passage_percolation.md`): the
variance scaling Var T(0,x) ~ ||x||^{2χ} with the rigorous bound
0 <= χ(d) <= 1/2 (survey line ~1675, Alexander/Kesten) and the physicists'
prediction χ = 1/3 in d = 2 (line ~1571; Open Question 13, line ~1784);
the scaling relation χ = 2ξ − 1 (line ~158) as context.** Reproduce as
statistical brackets ONLY — e.g. the measured growth exponent of Var(T_r)
vs r on Z^2 corridors consistent with the (0, 1) open interval and never
claimed against 1/3 — plus the subadditive time constant μ = lim_r T_r/r
existing (bracket its convergence from above, T_r/r ≥ μ by subadditivity).

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
statistical brackets on 1D boxes. **Anchor quantities (Biskup, archived
`docs/ref/md/Biskup_2003_long_range_chemical_distance.md`): the main
theorem D(x,y) = (log|x−y|)^{Δ+o(1)} for s ∈ (d,2d) with Δ⁻¹ =
log₂(2d/s) (the sharp TYPICAL-distance statement — bracket this);
s = d gives D_N = Θ(log N / log log N) for the DIAMETER with the
typical distance "at most logarithmic" and the lower bound explicitly
missing (bracket the upper-envelope only); s < d is the stochastic-
dimension regime where Eq. (1.8) bounds the SUP/diameter distance by the
power law — the typical distance is smaller and no sharp statement is
proved (bracket the box DIAMETER, never the typical distance, there);
s ≥ 2d returns to linear chemical distance with a positive time constant
(the critical case s = 2d is explicitly not well understood — do NOT
bracket it).** On d = 1 the regime boundaries fall at s = 1 and s = 2:
bracket the (log r)^Δ typical distance inside (1,2) (Δ⁻¹ = log₂(2/s)
checkable per s), the diameter envelope at s ≤ 1, and linear-plus-
constant typical distance for s > 2. C0 vector output = distance
profiles.

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
prescribed power law within CI; the white-noise limit `a -> infinity`
reduces bitwise to iid Bernoulli. (`a = 0` is the constant-field limit
`C(r) = 1`, not iid — Pruessner's `a = infinity` row is the uncorrelated
reference; this corrects the earlier `a = 0` wording, which mixed the
`C(r) = (1+r^2)^{-a/2}` and Schrenk `H` conventions.)

### C2: Shifted thresholds
p_c(a) on square/cubic vs Pruessner's tables (statistical brackets).
**Anchor values (Pruessner et al., archived
`docs/ref/md/Pruessner_2018_correlated_thresholds.md`, Table ~1, 2D
square, Euclidean): a=inf 0.592746 (the uncorrelated reference), a=3
0.561406(4), a=2.5 0.556214(4), a=2 0.550143(5), a=1.75 0.546717(7),
a=1.5 0.54299(1), a=1.25 0.53895(2), a=1 0.53452(4), a=0.75 0.5296(1),
a=0.5 0.5239(3)** — bracket against the quoted uncertainties as
tolerance bands, never exact claims; cluster compactness (fractal
dimension drift) observable alongside.

### C3: Weinrib–Halperin map
The extended Harris criterion a·ν_iid = 2 — verify the correlation-length
exponent flips to **ν = 2/a** (Weinrib–Halperin long-range-correlated
value) only BELOW the criterion, brackets across (a,d). Pruessner's own
ν_a column shows exactly this on the 2D square lattice: ν = 4/3 (the
uncorrelated value) down to a = 1.5, then 2/a takes over — 8/5 at a=1.25,
2 at a=1, 8/3 at a=0.75, 4 at a=0.5 — i.e. the flip sits between a=1.5
and a=1.25, consistent with a·ν_iid = (3/2)·(4/3) = 2. Reproduce the
flip location and the two branches as statistical brackets.

## 6. Q track: generic-q random cluster

**Anchor (archived): Geyer, *Efficient simulation of the random-cluster
model*, J. Stat. Mech. (2015) P08010, arXiv:1307.6647** — Sweeny
implementations compared (single-bond, ACES/union-find variants).

### Q1: Sweeny sampler
Single-bond heat-bath updates at arbitrary real q on the shared F1
topology: add/remove edge with the q-cluster-weight ratio; acceptance via
the connected check (the union-find connectivity query cost is THE
bottleneck — Geyer's whole point). **Anchor formulas (Geyer, archived
`docs/ref/md/Geyer_2015_sweeny_random_cluster.md`): heat-bath acceptance
p_acc^HB(Δb, Δn) = q^{Δn} v^{Δb} / (1 + q^{Δn} v^{Δb}) (Eq. 4) with
Metropolis min(1, q^{Δn} v^{Δb}) (Eq. 5) as the comparison variant; the
acceptance ratio bound 1/2 <= pHB/pMH < 1, minimal 1/2 at criticality
v = sqrt(q) in the percolation limit q -> 1; Δn ∈ {0,−1} on insertion
(internal vs amalgamating) and {0,+1} on removal — the connectivity query
that decides Δn is the documented cost center.** Gates: exhaustive 2^E
enumeration on small graphs against the exact FK weight
`p^|A| (1-p)^{E-|A|} q^{k(A)}/Z` (computed with own union-find k(A); Z by
full sum); q=1 reduces bitwise to R2; q=2 marginal statistics agree with
R3's exposed draws; detailed balance argued; the pHB/pMH ratio bound
verified across (q, v); acceptance/connectivity cost benchmarked against
Geyer's reported characteristics.

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

All four T-track residue stages are delivered AND independently reviewed:
T1/T2 (2026-10-22) and T3/T4 (2026-10-23, reports in `docs/reviews/`)
are Validated. F1, F2, and L1 are delivered and independently reviewed
(PASS-with-follow-ups; follow-ups applied — L1's included a real
`coordinates()` axis-reversal fix for d >= 2, caught by the review) —
Validated. L2 is delivered, independently reviewed 2026-10-24
(no Blocker; the review's two Majors — a seed-fragile s=2.5 gate and a
pre-run `occupied_edge_count` bug — are fixed as follow-ups: the linear
gate is now paired per realization at R=24 with 95% CIs, and `prepare`
reserves adjacency capacity without resizing; the estimator, overflow,
and wording minors are also applied) — Validated. L3 is delivered with
its seven diameter-regime gates green (2026-10-24: the s>2 ratio/slope/
time-constant bracket, the sub-polynomial (1,2) bracket, the s=1
log-envelope bracket, the bounded s<1 case, the ungated s=2 critical
case, ladder replay parity, and typed rejections) and is awaiting
independent review. F3 is delivered with its six gates green (2026-10-24:
uniform l1-ball shapes, exact geodesic counts vs brute force and closed
forms, directional time constants converging from above with positive
paired shape inequalities, bitwise draw/determinism, an explicitly
EXPLORATORY aging statistic, and typed rejections) and its independent review PASSED
2026-10-24 (the Major overflow-edge miscount is fixed; minors recorded) —
Validated. L3's independent review PASSED 2026-10-24
(no Blocker/Major; the four minors — realization-count capacity panic,
overlapping row RNG streams, SE wording, and validation-order wording —
are fixed), so L3 is Validated. The C track is delivered as Experimental:
C1's field/adapter gates are green, C2's large-L bracket now contains Pruessner's a=1.5 threshold
(L=128/256/512, 64 replicas: p_c(inf)=0.5383, CI 0.5246-0.5519 vs
0.54299), and
C3's analytic map is exact while its statistical nu certification is
open. The Q track is delivered as Experimental: Q1 passes exact-FK
marginal and q=1 bitwise gates, Q2 agrees with Q1 at q=2 within noise, and
Q3 anchors the exact 2D self-dual curve with an explicitly exploratory
phase-gap scan. Independent review of C/Q remains; the N5 audit is recorded in
`docs/reviews/2026-10-24-v2-final-audit.md` (10 Validated + 6
Experimental, with C2's threshold-table limitation and the C3/Q2/Q3
research-grade items explicit).

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
