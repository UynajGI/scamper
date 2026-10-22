# Z1 (Newman-Ziff Parameter Scans) — Full Reviewer Report (verbatim)

> Archived from the independent reviewer agent, 2026-10-22. See
> `2026-10-22-stage-reviews.md` for the summary and follow-ups applied.

**VERDICT: PASS**

## Roadmap §10 gates

| Gate | Result | Evidence |
|---|---|---|
| Every occupation-count step matches naive recomputation | **PASS** | `scan_prefixes_match_naive_analysis` re-runs a fresh F3 `analyze` at every prefix (all V+1 / E+1 points), field-by-field, all 8 topologies × both modes × 4 seeds; passed at default and 64 seeds. |
| Reweighted estimates agree with direct Bernoulli within CIs | **PASS** | Combined-CI test at p=0.55 (8192 direct vs 256 scans) plus exact 3x3-enumeration gates (2^9 site / 2^12 bond, p = 0.3/0.5/0.7, all 5 observables, \|z\| < 4.5, 256 scans); deterministic bond-ActiveEdgeCount exact path and sign-balance check present. |
| Benchmarks demonstrate an actual parameter-scan advantage | **PASS** | Reviewer re-ran Criterion: site 41.73 ms vs 203.23 ms direct = **4.87x**, bond 127.56 ms vs 229.83 ms = **1.80x** — reproducing the documented 4.65x / 1.83x within noise. |

## Independently computed checks

- **Reweight recurrence vs the archived paper's Eq. (1)**
  (`Q(p) = Σ_n C(N,n) pⁿ(1−p)^(N−N−n) Qₙ`): recurrence
  `B(k+1)/B(k) = (n−k)/(k+1) · p/(1−p)` and the downward mirror are exactly
  `fill_weights` (reweight.rs:96-106). Python mirror: recurrence weights vs
  log-space exact binomials (n up to 997, p ∈ {1e-6 … 1−1e-9}): worst
  normalized deviation **9.5e-14**; weighted mean vs Eq. (1): worst relative
  error **1.2e-14**; mode `floor((n+1)p)` attains the maximum for
  n ∈ {1..4096}; endpoints n=0 / p=0 / p=1 exact; at n = 10⁶ the recurrence
  matches log-space to ≤ 4.1e-15 relative (no log-space safeguard needed).
- **u128→f64 guard**: V = 94,906,265 → V² ≤ 2⁵³; V+1 exceeds — threshold
  exactly right.
- **Brute-force prefix parity (reviewer's own code)**: independent
  union-find + flood fill, every prefix of 500 random permutations on a
  4-vertex multigraph (self-loop + parallel edges), 500 on triangle+pendant,
  100 on 3x3 square, both modes: identical at every prefix. End-to-end NZ
  scan + direct binomial reweight over 4000 scans converges to exhaustive
  2ⁿ enumeration means (3x3, site & bond, p = 0.3/0.5/0.7) within Monte
  Carlo error (7e-5–4e-3 relative). Bond ActiveEdgeCount = 12p to 5e-14.
- **Benchmark fairness**: same observable and deliverable on both sides; the
  direct path includes a full F3 analyze per sample; identical graphs and
  RNG family; reweighted NZ per-sample sd = 0.32x direct at p=0.5, so the
  equal-replica comparison is conservative for NZ.
- **Doc claims vs live probe**: 10 allocation records all 0.000/sample;
  retained-capacity bytes reproduce exactly (chain 335,927/335,919; square
  562,104/593,848; cubic 759,096/818,488; sparse-ER 1,036,344/1,134,648).

## Issues found

1. **note** — reweight.rs:90 `(trials + 1)` could overflow-panic in debug at
   `usize::MAX` (siblings use saturating; consistency nit; applied in
   follow-up).
2. **note** — `ScanEstimator` one-pass variance (`sum_sq − sum²/n`) can lose
   precision when \|mean\| ≫ std; `.max(0.0)` guards negativity; fine for
   scan observables; Welford if the Z-track grows.
3. **note** — retained-capacity table in PERCOLATION.md shows 8 of 10 probe
   records without saying "excerpt" (doc nit; fixed in follow-up).
4. **note** — the parity oracle is F3's union-find (not literally a flood
   fill); mitigated because F3 itself is validated against an independent
   flood fill, the two union-finds share no code, and the reviewer's
   from-scratch flood fill closed the loop directly on the scan.

No blockers, no majors. RNG draw-order, numeric drift at N=10⁶, estimator
semantics, and off-by-one failure modes were each investigated and dismissed
with evidence.

## Recommendation

**PASS.** Faithful implementation of Newman-Ziff Eq. (1) with a numerically
sound mode-anchored recurrence; correct at every occupation count on all
eight topologies; statistical gates real and reproducible at higher power;
claimed speedup fair and reproduced.
