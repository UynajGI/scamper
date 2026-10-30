# Independent Review — T4 Continuum Winding (commit `d663435`, branch `feat/cmc-percolation`)

Reviewer: independent agent; did not author the code. Scope:
`CMC.rs/src/percolation/continuum_winding.rs`,
`CMC.rs/tests/physics/continuum_winding.rs`, the `CONTINUUM_WINDING_ALLOC`
probe in `CMC.rs/benches/percolation_allocations.rs`, the PERCOLATION.md T4
row + section, against the stage spec (roadmap v2 §2 T4). All verification
below was performed with a from-scratch harness (own RNG — xoshiro256++ from
splitmix64; own minimum image — `rem_euclid` + compare, a third convention
besides production's reciprocal-multiply and the committed test's division
form; own offset generator; own BFS unfolding + exact integer rank; own
Poisson sampler — exponential inter-arrival sums) at `/tmp/t4-review`, built
as a path dependency on `cmc-rs`. No project code was modified; no git
operations performed.

## Verdict per claim

| # | Claim (where) | Verdict | Independent evidence |
|---|---|---|---|
| 1 | Offset convention `n_c = round((d_c − r_c)/L_c)` is the boundary-crossing count of the minimum-image connection (rustdoc §"integer cell offset") | **Reproduced** | Re-derived from scratch: both descriptions of the pair differ by an integer number of box lengths, and for in-box representatives `n_c = −round(r_c/L_c)` with sign = net boundaries crossed. Hand checks through the public constructor: (0.1, 9.9) in box 10 → `[-1, 0]`; same-cell pair → `[0, 0]`; duplicate positions (distance 0) → `[0, 0]`. Swept beyond the committed runs: **1,156 randomized configurations** (2D+3D, isotropic `[12,12]`/`[12,12,12]` and anisotropic `[10,7]`/`[9,11,13]`, sizes 13–500, cutoff scales 0.3/0.7/0.95/0.999999 of L_min/2, six seeds, four modes: plain, snapped to exact boundary values {0, L, L/2} incl. duplicate particles, nudged to `L.next_down()`, and re-wrapped representatives shifted by random integer multiples of ±3L per axis), **1,831,913 edges / 4,319,115 axis entries: exact per-edge equality, 0 disagreements**, 0 constructor rejections, and **0 exact minimum-image ties** — a tie (axis separation exactly L/2) has Euclidean distance ≥ L/2 > cutoff, so P1 can never emit one at a legal cutoff; empirically confirmed. |
| 2 | Exact negation through the lazy `DirectedDisplacement` (W1 antisymmetry), for every entry | **Reproduced** | Forward vs reverse query compared as i64 (overflow-safe) on all 4,319,115 axis entries: `f == −b` everywhere, plus the committed 1,222,184. Boundary: `n = ±2,000,000,000` accepted and negates exactly (`i32::MIN` = −2,147,483,648 is unreachable since `MAX_CELL_OFFSET = 2e9 < 2^31−1`; the `as i32` cast cannot saturate — claim verified arithmetically and by probe). Non-endpoint queries return `None`. |
| 3 | Gauge invariance: representative re-wrapping shifts incident vectors telescopically, cycle sums and winding facts unchanged (rustdoc "Representative (gauge) freedom") | **Reproduced** | 90 re-wrap trials (3 configs: 2D dense n=160, 3D n=160, anisotropic 2D n=100; each trial re-wraps a random ~30% subset by random ±1..3L per axis). 241,230 per-edge checks: every stored vector shifted by exactly the predicted telescoping amount `n + a_i − a_j`. All fundamental-cycle sums (my own union-find spanning-forest cycle basis) unchanged. Production `analyze_winding` facts identical every trial: component count, per-component rank, wraps verdict, and exact generator vectors. |
| 4 | Precondition `cutoff < L_c/2` per axis: rejected with the axis named; just-below works | **Reproduced** | 7/7 probes: the constructor walks axes in order and names the **first** offender — `[4,10,6]` cutoff 2.0/3.0/100 → axis 0; `[10,8,6]` cutoff 3.0/3.1 → axis 2 (isolated: below 5 and 4), 4.5 → axis 1 (isolated), 5.0 → axis 0. Cutoff `L_min/2 − 1e-12` constructs. P1 **accepts** cutoff exactly L/2 = 5.0 on L=10 and emits the pair; T4 rejects the same request with axis 0 — the "P1 accepts the tie, T4 deliberately does not" claim verified end to end. |
| 5 | Why the precondition is load-bearing: double-image tie makes winding physically ambiguous ("sub-ulp coordinate noise at the tie can flip verdicts") | **Reproduced** | Re-derived: cutoff < L/2 ⇒ unique closest image (two in-range images would be ≥ L_min apart but ≤ 2·r_max < L_min apart — contradiction); so the minimum image is the one physical connection and `n` its crossing count. Beyond it, demonstrated the ambiguity at the formula level: a hand-built triple at 0, 5+ε, 5−ε in box 10 with cutoff exactly L/2 gives cycle winding sum **0 (ε=0), −1 (ε=+1e-15), +1 (ε=−1e-15)** — the same physical configuration flipping verdicts under ~4-ulp coordinate noise, exactly as documented. The constructor makes such inputs unrepresentable by policy. |
| 6 | The physics gate: production rank/wraps agree with an independent BFS unfolding, component by component | **Reproduced, 0 disagreements** | My own BFS unfolding (own min-image, own non-tree generator extraction, exact **integer** row-elimination rank instead of the committed float Gaussian elimination) over **47 fresh Poisson Boolean configurations** with seed namespaces disjoint from every committed seed: 27 committed-density analogs (2D L=40 N̄=2500 η∈{1.5,0.7,1.128}; 3D L=15 N̄=4000 η∈{0.6,0.15,0.34}) — 85,648 particles, 10,554 components — plus a 20-config small-box hunt (L=10, N̄=6, η=0.55, cutoff 3.412 < 5 — all preconditions satisfied), 117 particles, 32 components incl. rank-2 double-winding small clusters. Production vs reference: size, canonical identity, rank, wraps — **0 disagreements on all 10,586 components**. Density sanity also reproduces: spanning regimes wrap in 6/6 (2D) and 5/5 (3D) configurations; sub regimes 0/4 and 0/3; threshold regimes mixed (5/5, 2/4), as expected near η_c. |
| 7 | Tolerance `1e-9`: fractional/non-finite/>2e9 ratios rejected loudly, edge and axis named; never silent | **Reproduced** (with a doc nit, F4) | Through the public constructor: NaN/±inf coordinates → `NonFiniteCoordinate{particle:1, axis:0}` (named); finite ±1e308 whose raw delta overflows → `NonIntegerCellOffset{edge:0, axis:0, ratio:NaN}`; cutoffs {0, −1, NaN, +inf} → `InvalidCutoff`; the committed 1e8-boxes case reproduces with ratio −100000000.99999999 (residual ≈ 1.0e-8, matching the documented "~1.5e-8" ballpark). Magnitude bound: k = ±2e9 accepted exactly (n = ∓2000000000, negation exact); k = 2000000001, k = i32::MAX, k = 3e9 all rejected naming edge 0 axis 0. Committed unit tests for `integer_cell_offset` (private) re-verified by reading; the public-API probes above reach the same paths. |
| 8 | Allocation probe: 0 steady-state allocations; construction a per-realization transient | **Reproduced exactly** | `cargo bench -p cmc-rs --features allocation-probe --bench percolation_allocations`: `CONTINUUM_WINDING_ALLOC topology=poisson-2d N=2048 E=1588 cutoff=1 measured_samples=128 construction_allocations_per_realization=5735 construction_allocator_bytes_per_realization=719104 workspace_total_capacity_bytes=196640 workspace_bytes_per_vertex=96.016 allocations_per_sample=0.000 allocated_bytes_per_sample=0.000` — every PERCOLATION.md number exact; the 0-allocation assertion passed. |
| 9 | Gates rerun | **Reproduced** | `cargo test -p cmc-rs --test suite continuum_winding`: **6/6 pass in 8.42 s (debug)**, **6/6 in 0.59 s (release)**; lib unit tests **8/8**. Committed census numbers from the release rerun: negation sweep **515,079 edges / 1,222,184 axis entries / 0 disagreements** (exact); Boolean census components **50,484** (exact), particles **259,068** (see F3), spanning wraps 20/20 per dim, sub wraps 0/20, zero spanning-without-wrapping events — all matching. My independent estimate of the sweep totals from expected Poisson pair counts (≈515.5k edges) is consistent. |
| 10 | "W2's analyzer is consumed unchanged"; "discharges the v1 W-track note" | **Reproduced** | `git show --stat d663435`: only `continuum_winding.rs` (new), `lib.rs`/`percolation/mod.rs` re-exports, tests, bench — `analyzer/winding.rs` untouched. The v1 note exists verbatim in `docs/plans/2026-09-03-percolation-platform-roadmap.md` §16 ("Continuum and FK wrapping may reuse this analysis contract only after their own geometry and physical configurations are defined"); T4 defines the geometry (the derivation + precondition) and validates on physical Boolean configurations — the discharge claim is accurate. PERCOLATION.md's P1-accepts-tie and adapter-magnitude-gate statements also check out (`percolation/continuum.rs` docs + `carlo.rs` magnitude gate). |
| 11 | "winding around axis c implies a BFS-unwrapped extent ≥ L_c" (rustdoc lines 74–77; test doc lines 34–37 and assertion lines 568–579; PERCOLATION.md lines 3237–3238 and T4 row) | **NOT reproduced — the implication is false** (Finding F1) | Deterministic counterexamples below, including the committed hand configuration itself. |

## Findings

### F1 — MAJOR (docs overclaim + gate assertion that is not a theorem)

The claimed one-way implication "winding ⇒ BFS-unfolded extent ≥ L" is
mathematically false, and the committed physics gate *asserts* it
(`tests/physics/continuum_winding.rs:568-579`), so the gate's validity rests
on the committed 80-configuration window never containing the falsifying
event — an empirical accident of density, not a property.

What is true: a cycle with winding `w` drifts by `w·L` when followed via
minimum-image steps (my derivation and the committed reference both confirm
it), so the *cycle's* min-image lift has range ≥ L. What is measured in the
test is the **BFS-tree unfolding of the component's vertices** — a different
object. For a minimal wrapping ring the tree unfolds vertices near their
wrapped representatives and the wrap is closed by the non-tree edge, so the
extent is strictly below L:

- **The committed hand ring itself** (particles at x = 1, 3, 5, 7, 9 in box
  10, `tests/physics/continuum_winding.rs:210`): production correctly reports
  rank 1, generator `[1, 0]` (my reference agrees), but the BFS-unfolded
  x-coordinates are {−3, −1, 1, 3, 5} (seeded at vertex 0, wrap edge (0,4)
  unfolded through the minimum image) — extent 8 < 10. Run through the
  Boolean-gate assertion, this committed configuration would FAIL it.
- Minimal 3-ring (0, 4, 8 in box 10, cutoff 4.5 < 5): rank 1, extent 6 < 10.
- **Inside the Poisson Boolean ensemble itself** (all preconditions
  satisfied): 20 fresh configurations at L=10, mean N=6, η=0.55 (cutoff
  3.412 < L/2) — **3 of 20** contain a component where production rank ≥ 1
  with a wound axis whose extent is below L (sizes 13/7/5, ranks 2/1/1,
  spans [false,false] each time). Any of these would fail the committed
  assertion. The committed window (η = 1.5/0.6 at N̄ = 2500/4000) passes only
  because its wrapping components are giant, with tree-detour extents ≫ L;
  near-threshold or small-box ensembles, where the spanning cluster *is* a
  thin ring, hit the event routinely.

Production code is not wrong — its winding verdicts are correct everywhere
(0 disagreements against my reference on all 10,586 components, including
the three events). The defects are (a) three false "theorem" statements
(rustdoc `continuum_winding.rs:74-77`; test doc
`tests/physics/continuum_winding.rs:34-37` "necessarily has extent ≥ L_c −
tol (asserted)"; PERCOLATION.md:3237-3238 and the T4 row's "every wound axis
has unwrapped extent `>= L` (the exact direction of the extent criterion)"),
and (b) a gate assertion encoding the false direction while the test file
explicitly (and correctly) declines to assert the *other* non-implication
(spanning-without-wrapping is counted, not asserted). Follow-up: state the
extent claim for the cycle lift (or as "giant-component heuristic"), and
count winding-without-extent events in the census instead of asserting them
absent — exactly the treatment already given to spanning-without-wrapping.

### F2 — no blockers, no production-code findings

Every production behavior tested — table derivation, negation, gauge
invariance, precondition rejections (axis named, first-offender), tolerance
and magnitude rejections, boundary values, unwrapped representatives,
anisotropic boxes, 2D/3D — reproduced exactly. No input I could construct
produced a silent wrong integer or a wrong verdict inside the sanctioned
regime (see F4 for the boundary analysis).

### F3 — MINOR (docs): census particle total off by 100

PERCOLATION.md line 42 (T4 row) and line 3230 state **259,168** particles
across the 80 Boolean configurations; the deterministic rerun totals
**259,068** (49,744 + 49,865 + 79,522 + 79,937). Components (50,484),
configurations (80), all four per-regime component counts, and every other
census figure reproduce exactly — this looks like a transcription slip
(168 vs 068).

### F4 — MINOR (docs): tolerance-onset claim is not a threshold

Rustdoc lines 41-46 say the tolerance "is only exceeded by inconsistent
representatives separated by ~`1e8` box lengths **or more**". Measured
through the public constructor (box 0.1, fractional part 0.05): separations
up to 3e7 boxes all construct; 1e8 rejects; but **1e9 boxes was ACCEPTED**
with the correct integer (−1000000001) — the fp cancellation happened to
cancel. Rejection onset (~1e7–1e8 boxes) is right, but it is
cancellation-luck, not monotone in separation; "or more" overstates. The
important safety property does hold, and I verified it independently: within
the accepted magnitude bound |ratio| ≤ 2e9, the cancellation error is ≤
~2e9·2⁻⁵²·few ≈ 1e-6 ≪ 0.5, so an *accepted* offset can never be silently
off by one; beyond 2e9 the magnitude rejection fires first.

### F5 — MINOR (docs): per-edge storage arithmetic

`continuum_winding.rs:89-91` and PERCOLATION.md:3206 claim "`24 + 4 D` bytes
per edge" for the two owned vectors. `Vec<[usize; 2]>` costs 16 bytes/edge
plus `Vec<i32>` costs 4·D/edge — i.e. **16 + 4D**. The claim overstates by
8 bytes/edge (the measured construction accounting in the bench is
unaffected — it measures whole-realization allocations).

### F6 — NOTE: tie behavior for hand-built pair lists

For a hand-built pair at exactly half-box axis separation (unreachable from
P1 at a legal cutoff — distance ≥ L/2 > cutoff; 0 ties in 1.83M swept
edges), the constructor returns a well-defined integral offset via
round-half-away (e.g. separation 5.5 in box 1 → n = −6): the quotient-level
tie is integral by convention, no error. Consistent with the documented
precondition; recorded because the three minimum-image conventions
(production reciprocal-multiply, committed test division, my rem_euclid)
would differ by ±1 at such a tie, which is precisely why they cannot arise
in the sanctioned path.

### F7 — NOTE: committed gate 1's "independent" recomputation is weakly independent

`tests/physics/continuum_winding.rs:158-164` recomputes `n` with a
division-form minimum image — algebraically identical to production's
reciprocal-multiply form, differing only at quotient ties (which cannot
occur; see F6). It is a real arithmetic cross-check but not a convention
cross-check; my `rem_euclid`-based generator (a genuinely different
arithmetic) also agrees on all entries, which closes the gap.

## Overall verdict: **PASS-with-follow-ups**

The deliverable — a continuum displacement source that lets the unchanged
W2 analyzer run on continuum configurations — is correct as far as I could
break it: 4.3M axis entries exact against a from-scratch third-convention
generator, gauge invariance exact over 90 re-wrapping trials, 10,586
components of winding verdicts agreeing with an independent BFS-unfolding
reference on 47 fresh configurations, all rejection paths typed and named,
committed gates and bench reproducing exactly in debug and release. All
findings are documentation accuracy or gate-window fragility, not production
behavior; the one Major (F1) requires doc/gate follow-up but no analyzer or
constructor change. Suggested promotion path: apply F1/F3/F4/F5 doc and
assertion fixes, then promote from Experimental alongside T3.

Harness: `/tmp/t4-review` (untracked, outside the repo), full run 5.8 s.
