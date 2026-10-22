# D3 (Layered Space-Time Directed Percolation) — Full Reviewer Report (verbatim)

> Archived from the independent reviewer agent, 2026-10-22. Reviewed commit
> `91d07b5`; reviewer artifacts under `/tmp/d3trace/`. See
> `2026-10-22-stage-reviews.md` for follow-ups applied.

**VERDICT: PASS WITH NOTES**

No blockers, no majors. Core algorithm verified by hand trace, an
independent from-scratch model, and the explicit-DAG parity gate; physics
anchors genuine and convergent; memory/allocation claims hold.

## Per roadmap D3 requirement

All five requirements PASS (layers explicit; one sweep = one complete
realization; never T independent measurements — n_bins == sweeps gate;
frontier DP == explicit path enumeration via the D2 engine; two-frontier
bitset memory with 0 allocations/realization on five topologies).
"Arrival profiles" realized as activity/density + extinction time (no
first-arrival observable; recorded explicitly in the follow-up).

## Independent checks

- **From-scratch xoshiro256++ model** (SplitMix64 seeding identified
  empirically; f64 = (u64>>11)·2^-53 confirmed bit-exact): reproduces every
  production trace exactly — hand-checked 4-ring T=3 draw-by-draw; L=150
  3-word T=40 parity (334 draws = 2·Σactive, exactly one draw per directed
  incidence); L=70 Full trailing-bit mask; L=64; p=0; p=1 bipartite
  alternation; layers=1; bit-identical repeats.
- **Paper anchors verified in the archived md**: Table 1 bond DP
  p_c = 0.6447001(1) (with the DK p2,c = p(2-p) cross-check); Table 2
  delta = 0.159464(6), beta = 0.276486(8), nu_parallel = 1.733847(6).
- **Convergence**: release at 8000 realizations: delta_hat = 0.1604
  (survival) / 0.1590 (decay) vs 0.159464 — residuals +0.0009/-0.0005,
  no detectable systematic bias at this precision.
- **Statistical honesty**: propagated sigma(delta_hat) at N=2000 is 0.0145;
  seed spread (0.1659/0.1607/0.1487) consistent; the 0.05 gate is ~3.4
  sigma; finite-size contamination ruled out (t=1200 << L^z ~ 19,000).
- **Parity-gate independence confirmed**: both incidences APIs walk the
  same CSR slice in the same order.

## Issues found

- **M1 (minor, FIXED)** — statistical gates used one giant bin, printing
  vacuous "± 0.0000" error bars; now binned (sweeps/20) so printed
  uncertainties are real.
- **M2 (minor, FIXED)** — the docs' derived throughput (~6.7 M
  site-layers/s) contradicted the cited bench (~401 M/s in Criterion's
  element convention); corrected to the bench's own numbers.
- **N3 (note, FIXED)** — "O(L/8) words" was ceil(L/64) words (= L/8
  bytes); wording corrected.
- **N4 (note, RECORDED)** — no first-arrival-per-site observable; the
  roadmap wording's defensible reading is recorded explicitly.
- **N5 (note, FIXED)** — xi_parallel ~ 200 corrected to ~330 (pure
  power-law); conclusion unchanged.

## Recommendation

**Accept.** Physics engine bitwise-reproducible, provably equal to
explicit path enumeration, statistically anchored to independently
verified Hinrichsen values with seed-robust tolerances, allocation-free.
