# T1 (adaptive residual-degree attack) — Independent Reviewer Report

> Independent review of commit `2a9165f` (branch `feat/cmc-percolation`), staged
> on top of `53b2019` (v1 platform complete). Reviewer did not write the code
> under review. All probes ran from `git archive` snapshots of `2a9165f` (post)
> and `53b2019` (pre) in `/tmp` — the working tree was never modified except for
> this report file. Reviewer harness: from-scratch reference implementation,
> from-scratch random multigraph generator, complete-graph key-sort oracle,
> N3 order dumps (pre vs post), and a heap-vs-naive timing probe.

**VERDICT: PASS WITH FOLLOW-UPS** (2 Minor, 6 Notes; no Blocker, no Major)

## 0. What was reviewed

- `CMC.rs/src/percolation/network/attack.rs` — `ResidualDegreeDescendingStableId`
  / `ResidualDegreeDescendingSeededKey`, `ResidualEntry`, `sample_residual_order`,
  workspace reservation.
- `CMC.rs/tests/physics/percolation_attack.rs` — T1 gates (parity, invariants,
  tie replay, Holme Fig. 7(e) bracket, BA builder).
- `CMC.rs/benches/percolation_attack.rs` (order group),
  `CMC.rs/benches/percolation_allocations.rs` (`RESIDUAL_ATTACK_ALLOC`).
- `CMC.rs/PERCOLATION.md` T1 row + T1 section.
- Spec: `docs/plans/2026-10-23-percolation-roadmap-v2.md` §2 T1.
- Anchor: `docs/ref/md/Holme_2002_attack_vulnerability.md` (Sec. IV E, Sec. VI,
  Fig. 7 caption/discussion, ~pages 6–9).
- Scope check: T2 (`7a9c5ef`) touched only `percolation_allocations.rs` on top;
  `attack.rs`, the test file, and `percolation_attack.rs` bench are byte-identical
  between `2a9165f` and the working tree.

## 1. Heap soundness, adversarial fuzz — REPRODUCED

Reviewer's own reference: full rescan per step, residuals recomputed from
adjacency (`alive` filter; self-loop = two incidences at its own vertex while
active), argmax under (residual desc, key asc, dense id asc) via explicit
three-way tie logic — no project code copied.

Graphs: **373 total** = 53 hand-built degenerates (single vertex; single
self-loop; **double parallel self-loop on one vertex** (degree 4 via 4 CSR
incidences); two vertices joined by 3 parallel edges; all-isolated graphs;
isolated+edge+self-loop mixes; stars K(1,1..8); complete K(2..10); chains and
cycles n=2..16) + 320 seeded random multigraphs (V uniform 1..=40, E uniform
0..=3V, endpoints uniform — self-loops and parallel edges allowed). Generator
feature audit over the 320 random graphs: **210 contain self-loops, 210 contain
parallel edges, 179 contain isolated vertices, 14 are edge-free, 10 are
single-vertex; max degree 16, mean 31.5 edges.**

Checks per graph: stable policy vs reference under identity keys with TWO seeds
(seed-independence asserted); seeded policy vs reference under the replayed draw
contract with THREE seeds (bitwise reproduction asserted); permutation validity.

Result (my harness, debug build, post snapshot):

```
[REVIEW FUZZ] 373 graphs (53 degenerates + random), 373 stable checks,
1119 seeded checks — all exact
```

**1492 exact full-order comparisons, zero mismatches.** No counterexample found.

## 2. Soundness argument — RE-DERIVED, HOLDS; no counterexample

Re-derivation (mine, from the code at `attack.rs:407-464`):

- **Invariant**: after n removals, every still-active vertex v has a *fresh*
  heap entry `(residual[v], key[v], v)`.
- **Base**: the initial loop pushes one entry per vertex 0..V
  (`attack.rs:416-425`), with `residual[v]` = original degree (self-loop = 2,
  matching `NetworkProfile`'s endpoint-counting convention).
- **Inductive step**: `residual[u]` is mutated in exactly ONE place
  (`attack.rs:455`), guarded by `residual_active[u]`, and every mutation is
  immediately followed by a re-push of the new value (`attack.rs:456-461`).
  So every decrement path pushes; no other path mutates residual degree.
  When the victim goes inactive it no longer needs an entry.
- **Pop loop**: pops are discarded only if `!residual_active[v]` or
  `entry.degree != residual[v]` — i.e. only entries that are provably NOT the
  vertex's fresh entry are discarded. A fresh entry is either popped as the
  heap maximum (becomes the victim) or stays in the heap. Hence the heap can
  never empty while an active vertex remains (the `expect` at `attack.rs:443`
  is unreachable), and no active vertex's fresh entry is ever silently lost.
- **Argmax**: `Ord` ranks (degree desc, key asc, **id asc**) — verified byte-level:
  `self.degree.cmp(&other.degree).then_with(|| other.key.cmp(&self.key))
  .then_with(|| other.vertex.cmp(&self.vertex))` — both comparisons reversed so
  the std max-heap pops the smaller key/id first. If the first fresh pop `e`
  were not the argmax, some active u's fresh entry would rank strictly above `e`
  yet still sit in the heap (fresh entries are never discarded, and `e` was the
  maximum of everything in the heap at that pop) — contradiction.
- **Answering the probe questions directly**: a vertex's fresh entry cannot be
  absent while the vertex is active (single guarded decrement site, all paths
  push; initial push covers all vertices). I additionally traced the two
  trickiest cases on paper — (a) parallel edges: each physical edge contributes
  one incidence in each endpoint's CSR row, so removing the first endpoint
  decrements the other exactly once per parallel edge; (b) self-loops: the two
  incidences live in the victim's own row and are skipped because the victim is
  already inactive — the self-loop correctly keeps counting 2 in the owner's
  residual while it is active and then dies with it. No underflow is possible:
  total decrements of u ≤ initial residual[u] by the incidence count.
- The empirical fuzz (item 1) includes the adversarial shapes (double
  self-loop, parallel-heavy) that would expose any hole here.

## 3. Tie semantics — REPRODUCED, bit-exact

- **Draw contract**: replayed in my own code — fresh
  `Xoshiro256PlusPlus::seed_from_u64(seed)`, then V draws of
  `rng.random::<u64>()` in ascending dense order — reproduces the production
  seeded-key order on all 1119 seeded fuzz cases.
- **Pin test**: on a complete graph K_n every step is fully tied, so the order
  IS the key sort. For (n, seed) in {(8,7), (11,99), (16,0xBADC0DE), (25,5)}
  the production order equals `vertices sorted by drawn key ascending`, an
  independent expectation computed without the reference. This proves the
  contract "one u64 per vertex, ascending dense order, before the first
  removal" bit-exactly, not just probabilistically.
- **Stable == naive with identity keys**: asserted on all 373 fuzz graphs
  (identity keys `0..V as u64`), plus seed-independence across two wildly
  different seeds.

## 4. N3 regression — REPRODUCED, bit-identical

- Structural: `git diff 53b2019 2a9165f -- attack.rs` shows the
  original-degree path is untouched — `sample_order` gained only an early
  return for residual policies; T1 adds fields (empty for N3 policies) and new
  code. The five pre-existing test bodies are character-identical (diff hunks
  are: module doc expansion, and a pure 471-line addition after the last N3
  test).
- Empirical: identical dump harness run in BOTH snapshots
  (`git archive` of each commit) printing 8 graphs x 2 original-degree policies
  x 3 seeds = 48 removal orders + 8 degree profiles (graphs: chain15, cycle15,
  star8, complete10, the multigraph fixture, isolated9, ER(30,60), dense
  ER(12,40)):
  `diff` of the two dumps: **identical (48/48 orders, 8/8 profiles).**
- The five pre-existing attack tests pass in the current tree (see item 7).

## 5. Holme anchor fidelity — PARTIALLY (2 Minor findings)

What the paper (Sec. VI, Fig. 7(e) discussion) actually claims for the BA
model: (i) "the differences among the removal procedures are not significant in
the early attack stage"; (ii) "as the removals proceed, the attack strategies
harm network more in the order, RB > RD > ID > IB". Fig. 7 caption: N = 1500
for panels (c)-(f); S normalized by the initial-network value; error bars from
~40 realizations. **There is no random-failure curve in Fig. 7 and Sec. VI
never quantifies random failure for BA.**

Test claims vs paper:

- (a) RD below ID, disjoint 95% CIs, f >= 0.15 — **consistent with (ii)**
  (restricted to RD vs ID; S lower = more harmful). Reproduced exactly; see
  item 7 numbers.
- (b) RD ~ ID, overlapping CIs, f <= 0.125 — **consistent with (i)**.
- (c) both attacks collapsed (upper CI < 0.05) by f=0.30 while random failure
  sits at ~0.69 — **true in the data and standard scale-free physics
  (Albert-Jeong-Barabási 2000, Ref. [20] of the paper), but it is NOT a Holme
  Fig. 7(e) claim.** PERCOLATION.md's T1 row groups it inside the "Holme
  Fig. 7(e) statistical bracket" parenthetical — loose attribution (Finding 3,
  Note). The T1 section's own "Readings" list words it as an observation
  without attributing it to Holme, which is fine.
- Late-stage separation reaching ~0.21 at f=0.25: verified (0.2189 - 0.0116 =
  0.2073).

BA builder vs Sec. IV E definition (paper: initial condition = **m0 vertices
and NO edges**; growth one vertex with m edges per step; attachment
P_u = k_u / sum_w k_w, iterated m times per growth step):

- The builder starts from a **complete graph on m0** (test
  `barabasi_albert_graph`, `tests/physics/percolation_attack.rs:784-816`;
  same in the bench). The doc comment says "the Holme et al. (2002) Sec. IV E
  construction: initial complete graph on m0 vertices" — **the paper says the
  opposite (no edges).** The complete core is the standard practical deviation
  (attachment probabilities are undefined when all degrees are 0) and is
  physically immaterial here (10 of 29,995 edges, 0.03%), but the comment
  misattributes it to the paper, and PERCOLATION.md does not flag the
  deviation (Finding 1, Minor).
- Endpoint-pool sampling gives exactly P_u = k_u / sum k_w; targets are drawn
  DISTINCT with redraws, and the pool is frozen during one vertex's m draws
  (degrees not updated between the m attachments of the same new vertex). The
  paper does not specify distinctness or per-edge updating; both choices are
  standard, with O(m / sum k) effect. Note only (Finding 6).
- N=10^4 vs paper's Fig. 7 N=1500: roadmap-mandated and explicitly
  acknowledged in PERCOLATION.md ("the paper itself uses N = 1500"). Honest.
- Edge count assert `10 + 3*(N-5)` = 29,995 verified by my own build.

## 6. Statistics deviation — ASSESSED CORRECT (with one Note)

The observable is a per-graph largest-component fraction S(f) in [0,1]
(continuous across graph realizations), averaged over R=24 independent
(graph, seed) pairs. The **empirical-SE interval mean ± 1.96 s/sqrt(R) is the
statistically appropriate** interval for the mean of an i.i.d. bounded
continuous observable; the Wilson interval models n Bernoulli trials of a
binary outcome and would be the wrong model here (it answers "what fraction of
graphs collapsed", not "where is the mean S"). Computing Wilson only as a
reported conservative companion is transparent and correct. The deviation note
in the test and PERCOLATION.md documents the reasoning accurately.

- **Design soundness**: one trajectory per graph with independent seeds across
  graphs samples the joint (graph, seed) product measure i.i.d. — the SE
  interval is valid for the combined graph+trajectory variance. Reusing one
  graph with many seeds as if independent would not be; they do not do that.
  For stable policies the trajectory is deterministic per graph, so the
  graph-to-graph variance is the whole story.
- **Disjoint-CI criterion**: requiring disjoint 95% CIs is *stricter* than a
  two-sample z test on the difference (threshold z(SE1+SE2) vs
  z*sqrt(SE1²+SE2²)), so the gate is conservative in the right direction.
- **z vs t (Finding 4, Note)**: at R=24 the normal quantile 1.96 is ~5.5%
  narrower than t(0.975, 23)=2.069. Recomputing the two marginal points with
  the t multiplier: f=0.15 stays disjoint (RD hi 0.7116 < ID lo 0.7151) and
  f=0.125 stays overlapping — every bracket decision is robust; still, a
  t-interval is the textbook choice at n=24.
- Multiple comparisons across 6 late grid points are uncorrected, but each
  disjointness event is roughly a 0.5%-level coincidence under the null and
  the direction is monotone/consistent — acceptable for a bracket.
- Wilson companion math checked by hand at p=0.9945, n=24 (half-width 0.0736)
  and matches the printed [0.8527, 0.9998]; the doc's "~0.14 half-width /
  two orders of magnitude" is a fair order-of-magnitude summary of what I
  measure across the grid (half-widths 0.07-0.19 vs empirical SE 0.0001-0.011).

## 7. Gates re-run — ALL PASS

Command correction first: **there is no `--test percolation_attack` cargo
target**; the attack tests live in the aggregated `tests/suite.rs` target
(`#[path = "physics/percolation_attack.rs"]`). Correct invocation (also the one
PERCOLATION.md documents elsewhere):

```
cargo test -p cmc-rs --test suite percolation_attack
-> 8 passed; 0 failed; 0 ignored   (16.2 s, debug)
```

The 8: 5 pre-existing N3 (`full_trajectory_matches_independent_bfs_every_step`,
`tie_policies_reproduce_and_differ_appropriately`,
`curve_coordinates_and_provenance_survive_results_json`,
`degree_cutoff_attack_matches_callaway_theory`,
`intentional_attack_bracket_matches_cohen_prediction`) + 3 new T1
(`adaptive_residual_order_matches_naive_rescan_every_step` — prints
"[T1 parity] exact full-order equality on 28 cases",
`adaptive_tie_policies_reproduce_and_differ_only_through_keys`,
`adaptive_attack_bracket_matches_holme_fig7e`).

Holme bracket numbers re-captured with `--nocapture` — **every row matches the
PERCOLATION.md table exactly** (deterministic seeds). Key rows (mine):

```
f=0.125: RD 0.7793 [0.7777,0.7808]  ID 0.7818 [0.7801,0.7835]  disjoint=false
f=0.150: RD 0.7092 [0.7069,0.7115]  ID 0.7174 [0.7153,0.7196]  disjoint=true
f=0.250: RD 0.0116 [0.0098,0.0133]  ID 0.2189 [0.2075,0.2304]  disjoint=true
f=0.300: RD 0.0014 [0.0013,0.0015]  ID 0.0072 [0.0060,0.0084]  fail 0.6905 [0.6900,0.6909]
```

So: disjoint for f in {0.15, 0.20, 0.25, 0.30, 0.35, 0.40}; overlapping for all
f <= 0.125; both attack upper CIs < 0.05 first at f=0.30 with random-failure
lower CI 0.6900. Wilson companion never disjoint (wilson_disjoint=false on all
14 rows) — the choice of interval is decision-relevant, and the chosen one is
the defensible one (item 6).

`cargo clippy -p cmc-rs --all-targets`: clean (no warnings, no errors).
`cargo test --workspace`: **exit code 0, all targets green** (full-suite wall
~10 min; the captured output tail alone shows the last 8 test targets with
102 passed / 0 failed, earlier targets truncated by the reviewer's tail
filter — the exit status covers the entire run).

## 8. Allocation & perf claims — REPRODUCED

Allocation probe (main working tree):

```
cargo bench -p cmc-rs --bench percolation_allocations --features allocation-probe -- --nocapture
RESIDUAL_ATTACK_ALLOC topology=chain     V=4096 E=4095  checkpoints=17 ... allocations_per_trajectory=0.000 allocated_bytes_per_trajectory=0.000 peak_live_allocations=0
RESIDUAL_ATTACK_ALLOC topology=square    V=4096 E=8064  ... 0.000 / 0.000
RESIDUAL_ATTACK_ALLOC topology=cubic     V=4096 E=11520 ... 0.000 / 0.000
RESIDUAL_ATTACK_ALLOC topology=sparse-er V=4096 E=16384 ... 0.000 / 0.000
RESIDUAL_ATTACK_ALLOC topology=power-law V=4096 E=16374 ... 0.000 / 0.000
```

0 allocations / 0 bytes on all five topologies x 128 trajectories; the probe's
hard `assert_eq!(info.count_total, 0)` passed. The reservation argument also
checks out analytically: pushes = V initial + exactly one per non-loop physical
edge (each non-loop edge decrements its later-active endpoint exactly once; a
self-loop never decrements), i.e. <= V + E pushes, and the constructor reserves
exactly `vertices + edge_count` (`attack.rs:286`).

Heap-vs-naive timing (my own release-mode probe, BA N=10^4, m0=5, m=3, same
builder): **heap 6.1 ms/trajectory vs naive full-rescan 785 ms/trajectory =
128x.** The committed Criterion claim is ~92x (5.84 vs 536.7 ms); my naive
reference allocates a per-step residual vector so it is slightly inflated, but
both numbers establish the same conclusion — the O((V+E) log V) vs
O(V(V+E)) cost model is realized at this size.

## 9. Docs / code quality — GOOD, with notes

- Rustdoc cost model (`AttackPolicy::ResidualDegreeDescendingStableId`):
  the push bound, the stale-skipping description, and the complexity statements
  are all exactly right (verified analytically in items 2 and 8). The
  "decrements exactly one endpoint once" sentence is precisely true as written.
- `ResidualEntry` doc ("ID only disambiguates a seeded-key collision") is
  accurate; for the stable policy key==id makes (degree, key) unique, and u64
  key collisions are ~V²/2^65 — the third comparator is effectively defensive
  dead code (harmless; Finding 8, Note).
- PERCOLATION.md T1 row: Experimental status is honest ("Experimental until
  independent review"); "28 cases" is arithmetically correct (7 graphs x [1
  stable seed + 3 seeded seeds]); the numbers table reproduces bit-for-bit.
  Two wording issues: Findings 1-3 below.
- Bench doc (`benches/percolation_attack.rs:15-17`) says the naive reference
  "duplicates the test-local reference ... kept in sync by the exact-parity
  gate there" — the parity gate covers the TEST-local copy, not the
  bench-local duplicate (which is additionally the identity-tie
  specialization, no key parameter). No automated link exists between the two
  copies; if the bench copy drifted the timing table would silently misstate
  the reference cost (Finding 5, Note).
- Clippy clean; no unsafe code in the diff; error paths typed via
  `NetworkPercolationError::Capacity`.

## Findings

**Blocker**: none.

**Major**: none.

**Minor**:

1. **BA builder doc misattributes the initial complete core to the paper**
   (`tests/physics/percolation_attack.rs:784-789`,
   `benches/percolation_attack.rs:78-82`): the comment says "the Holme et al.
   (2002) Sec. IV E construction: initial complete graph on m0 vertices", but
   Sec. IV E of the archived paper specifies "m0 vertices and **no edges**".
   The complete core is the standard, necessary deviation (P_u = k_u/sum k_w
   is undefined at all-zero degrees) and physically immaterial (10 of 29,995
   edges), but the comment asserts the paper says the opposite of what it
   says. Suggested fix: state the deviation explicitly ("paper specifies a
   bare m0 core; we use the standard complete seed so preferential attachment
   is defined"), and mention it in PERCOLATION.md's anchor paragraph.
2. **PERCOLATION.md presents the f<=0.125 overlap as a gate, but the assertion
   only pins the first grid point.** The test asserts overlapping CIs at
   f=0.005 only (`rd_lo0 <= id_hi0 && id_lo0 <= rd_hi0`); overlap through
   0.125 is the observed outcome of this seed set, not enforced. Either extend
   the assert over the early prefix or reword the row/section ("overlap
   observed through f=0.125, asserted at the first checkpoint").

**Notes**:

3. The random-failure contrast (S_fail = 0.69 at f=0.30) is grouped inside the
   "Holme Fig. 7(e) statistical bracket" in the PERCOLATION.md T1 row; Holme's
   Fig. 7 contains no random-failure curve. It is true, standard context
   (Albert et al. 2000), and the "Readings" bullet phrases it as an
   observation — but the row's parenthetical should not imply the paper
   quantifies it.
4. z=1.96 vs t(0.975,23)=2.069 at R=24: ~5.5% narrower intervals; all bracket
   decisions survive the t multiplier (f=0.15 disjoint margin 0.0034 under t),
   so no change required — a t-interval would be the textbook choice.
5. Bench-local naive copy is not actually gated against the test-local copy
   (see item 9); consider a one-line comment fix or a shared cfg(test) helper.
6. BA attachment freezes the endpoint pool during a new vertex's m draws and
   enforces distinct targets by redraw — both standard, unspecified by the
   paper; effect O(m / sum k). No action needed.
7. Wilson-companion phrasing "~0.14 half-width / two orders of magnitude" is
   approximate (measured 0.07-0.19 across the grid; ratio 30x-1900x depending
   on f) — acceptable as an order-of-magnitude statement.
8. The third `ResidualEntry` comparator (vertex id) is unreachable in practice
   (stable keys are ids; seeded u64 collisions ~V²/2^65) — fine as defensive
   determinism documentation.
9. `--test percolation_attack` does not exist as a cargo target (tests are
   aggregated in `suite`); PERCOLATION.md's own invocations use the correct
   `--test suite <filter>` form, so this is only a trap for people reading the
   stage tasking text.

## Recommendation

**PASS WITH FOLLOW-UPS.** The core claims all reproduce independently: the
lazy-heap adaptive order is exactly the naive full-rescan argmax order under
both tie policies on 373 adversarial graphs (1492 order comparisons, zero
mismatches); the seeded draw contract is pinned bit-exactly by a
complete-graph key-sort oracle; N3 is bit-identical to 53b2019 both structurally
and empirically; the allocation probe reports hard zeros; the ~92x perf claim
is the right magnitude (my measurement 128x); and the Holme bracket table
reproduces verbatim. The statistics deviation (empirical-SE intervals gated,
Wilson reported) is the correct modeling choice for a mean over graph
realizations. The two Minor findings are documentation-accuracy issues (BA
seed-core attribution; early-stage gate scope) that do not affect any
simulation result; the stage can be promoted from Experimental once they are
addressed.

## Reviewer probe provenance

- Snapshots: `/tmp/t1review/pre` (53b2019), `/tmp/t1review/post` (2a9165f),
  via `git archive` — repo working tree untouched.
- Harness files (reviewer-written, snapshot-local only — never in the repo):
  `tests/t1_review_fuzz.rs`, `tests/t1_review_n3_dump.rs`,
  `tests/t1_review_genstats.rs` under `/tmp/t1review/{pre,post}/CMC.rs/tests/`.
- N3 dump diff: `/tmp/t1review/n3_pre.txt` vs `/tmp/t1review/n3_post_norm.txt`
  (48 order lines + 8 profile lines, identical).
