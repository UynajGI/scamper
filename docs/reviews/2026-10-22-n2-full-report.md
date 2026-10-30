# N2 (Random Network Failure + Degree Observables) — Full Reviewer Report (verbatim)

> Archived from the independent reviewer agent, 2026-10-22. Reviewed commit
> `bb11fb3`. See `2026-10-22-stage-reviews.md` for follow-ups applied.

**VERDICT: PASS WITH NOTES**

No blockers or majors; two minor items (both applied as follow-ups) and four
informational notes.

## Roadmap §9 N2 gates

| Gate | Result | Evidence |
|---|---|---|
| Random failure reduces to the site law | PASS | `failure.rs` calls the identical `SiteBernoulli::new(keep).sample(...)` as `StaticLaw::site`; uniform draw loop is per-vertex dense-index order; test 1 verifies bitwise identity vs adapter and raw law over 32 sweeps. |
| Path/cycle/star/complete hand calculations | PASS | Exhaustive 2^V enumeration with an adjacency-BFS reference (no analyzer code). 16 seeds: max \|z\| = 2.02, all mean \|z\| < 0.6. |
| Graph relabeling invariance | PASS | Constant arrays bitwise-identical; ratios and GiantFraction within 4 combined sigma, 2000 sweeps. |
| Stable shapes / empty-bin semantics | PASS | Full `[0..=4]` shape; bins {0,2,3} exactly 0.0 with 0.0 stderr; no NaN path exists. |
| Literature anchors | PASS | Both reproduced from the archived papers, independently re-derived (below). |

## Independently computed numbers (no project code reused)

**Poisson(2.5) ideal:** ⟨k⟩=2.5, ⟨k²⟩=8.75, **κ = 3.5000 → q_c = 0.4000**;
empirical graph κ = 3.5025 → q_c = 0.3996 (deviation 0.0004). Formula
verified against Cohen et al. Eq. (6) (`1 − p_c = 1/(κ0 − 1)`).

**Callaway Eq. (5) fixed point** (reviewer's own iteration, tol 1e-16):

| q | reviewer u* | reviewer S | test theory | MC measured | deviation |
|---|---|---|---|---|---|
| 0.5 | 0.814315 | **0.1857** | 0.1861 | 0.1857(1) | 0.0000 |
| 0.65 | 0.574166 | **0.4258** | 0.4259 | 0.4261(1) | 0.0003 |
| 0.8 | 0.362550 | **0.6374** | 0.6373 | 0.6377(1) | 0.0003 |

q_c = 1/G1′(1) verified as the paper's Eq. (12); bracket points outside the
N^(−1/3) ≈ 0.02 rounding window.

**Star (center + 4 leaves), keep = 1/2 — exact fractions, reviewer's own
flood-fill enumeration + closed forms (both agree):**
- ActiveFractionByDegree[1] = **1/2 exactly** (= keep)
- GiantMembershipByDegree[4] (center) = **1/2**
- GiantMembershipByDegree[1] = **47/128 = 0.3671875** (tie-policy-sensitive)
- DegreeFirstMoment [1]=9/8, [4]=3/2; DegreeSecondMoment [1]=7/2, [4]=5;
  empty bins {0,2,3} = 0 exactly
- keep=0 → all observables 0; keep=1 → all 1
- **V³ gate**: floor((2^53)^(1/3)) = **208063** (208063³ = 9,007,091,372,906,047
  ≤ 2^53 < 208064³ = 9,007,221,244,166,144); constructor comparison exact.

## Gate runs

- `cargo test --release ... percolation_network`: **5/5 passed** (78.6 s);
  printed numbers match PERCOLATION.md.
- `SCUTTLE_ZSCORE_SEEDS=64`: **5/5 passed**; max \|z\| = 1.98, max \|mean z\|
  = 0.81 — matches the doc verbatim.

## Facade semantics audit

Delegation bitwise-identical at code level; tie policy shared with the
scalar observables by construction (same analyzer pass); empty bins guarded
(`count > 0`), `bins().len() ≥ 1` so no empty-slice path; V³ gate verified;
steady-state allocation-free by inspection (preallocated buffers, Cow
names, amortized bin flush excluded by the documented binsize discipline).
Self-loops count 2, matching CsrLattice's even-incidence CSR convention.
S(k) ratio-of-means convention correct for Cohen et al.'s conditional S(k);
the susceptibility non-comparison (Callaway H0′(1) vs platform chi) is a
definitional distinction, correctly respected.

## Issues found

1. **[minor, FIXED]** PERCOLATION.md bracket sentence transcribed 0.1857
   (the keep=0.5 figure) where the q_c+0.10 run gives 0.1851 ± 0.0001 vs
   theory 0.1855.
2. **[minor, FIXED]** `InexactDegreeMoment` rejection branch untested;
   follow-up adds the exact 208063/208064 boundary test (edgeless lattices
   make construction cheap).
3. **[note]** Relabeling gate is single-run-per-labeling — weakest N2
   statistical gate, but with real power on occupied bins and bitwise
   constant-array checks.
4. **[note]** Pooled mean-z gate diluted by trivially-zero empty-bin entries
   (per-observable \|z\| < 4 gate unaffected).
5. **[note]** `record_degree_observables` with mismatched scratch/profile
   panics on index bounds — acceptable for its documented use.
6. **[note]** 64-seed path5/keep=0.7 mean z = −0.81 across 14 strongly
   correlated observables — consistent with common-mode Monte Carlo
   fluctuation; no evidence of bias.

## Recommendation

**PASS WITH NOTES.** Physics right (both literature anchors reproduce to
≤ 0.0005 via independent derivation); enumeration reference independent and
hand-verified to exact fractions; facade genuinely delegates. N3 can safely
build on this substrate.
