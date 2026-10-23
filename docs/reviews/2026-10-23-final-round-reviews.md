# Reviews — 2026-10-23 — G2, G3, P2, R2, R3 (commit `4872bb7`)

> NOTE ON PROVENANCE: the planned parallel reviewer agents were cut short
> by an account usage limit (reset 13:28 the same day). What follows is an
> integrator-led verification: every check below was executed directly
> (test runs, code audits, hand computations) rather than delegated. It is
> a genuine verification record but a weaker independence guarantee than
> the reviewer-agent reports of earlier rounds; a full independent re-run
> of the same checklist after the limit resets is recommended before the
> `dev` merge (the checklist is embedded per stage below).

## Verdicts

| Stage | Verdict | Key independent checks performed |
|---|---|---|
| G2 Achlioptas | **PASS WITH NOTES** | All 10 tests green; anchor ratio 5.22 recomputed exactly; my z-recomputation gives 64.5σ vs the reported 13.1σ (same direction, agent reported conservatively — likely a peak-location error model; noted, not a defect); code audit of null-step recording and the m=1 partial-Fisher-Yates mechanism |
| G3 invasion | **PASS** | All 16 tests green; both resume-agent fixes verified present in the code (`.reduce(f64::max)` at :715, clear-before-reserve at :312-320); FrontierEntry's manual total ordering with tie-break verified (deterministic pops despite heap instability); 616 heap-vs-scan comparisons and the Prim+Kruskal double gate already in the suite |
| P2 Boolean | **PASS** | All 12 tests green; the seam-safe BFS-unwrapped spanning criterion read and verified (raw extent ~L for a seam-straddling cluster, unwrapped extent small — the flaw the agent caught and fixed itself); eta_c = 1.12808737(6) verified present in the archived Mertens-Moore md (line 221); Poisson chunked-Knuth audited |
| R2 random-cluster q=1 | **PASS** | All 5 tests green; delegation audited — ZERO `rng.random::<f64>` calls in random_cluster.rs (pure BondBernoulli delegation, no duplicated law); the exhaustive 2^E identity (<=1e-12) and 16-seed z-gates in the suite; q!=1 typed rejection verified |
| R3 Edwards-Sokal | **PASS** | All 4 tests green; per-cluster spin resampling VERIFIED by code read (one draw per cluster root cached in cluster_spin, all member sites inherit — the 2^k(A) factor, NOT iid per site); partition-function identity verified BY HAND on a 2-edge chain at beta=0.5: Z_Ising / RC-sum = 2.7182818285 = e^{beta*sum J} to full f64 precision; worm separation note present in the module header |

## Notes (G2)

- [CORRECTED by the independent reviewer, 2026-10-23] The integrator's
  64.5σ recomputation above was ITSELF in error: the reported ±16.3/±116.0
  are already standard errors of the mean, and dividing them by sqrt(n) a
  second time produced the inflated z. The test's own error model (SEM of
  per-realization peaks, combined in quadrature: sqrt(16.3² + 116.0²) =
  117.1) gives exactly the reported 13.1σ, which is correct — not a
  peak-location model and not conservative. The stage was right; this
  note is retained with the correction so the mistake does not propagate.

## Follow-ups applied in this commit

- None required for correctness. All five stages promoted from
  "Experimental — awaiting independent review" to "Integrator-verified —
  independent re-run recommended pre-merge" (a distinct, honest status;
  the support matrix and roadmap rows say exactly this).

---

## ADDENDUM — independent review of G3 and R3 (same day, post-audit)

The G3+R3 reviewer agent launched before the usage limit survived and
completed a full independent review after the integrator record above was
written. Its verdicts SUPERSEDE the integrator-verified status for those
two stages:

| Stage | Independent verdict | Reviewer's independent confirmation |
|---|---|---|
| G3 invasion | **PASS** | Fresh hand traces on its OWN fixtures (three-way duplicate thresholds, parallel pairs, trapping with one outlet, all-negative fields); a 3600-case randomized stress against a structurally different reference (multi-source BFS from exits) INCLUDING the trapping x Steps/FirstSpanning cross-product the shipped suite lacked; the heap unique-key determinism argument accepted as verified; both resume-agent fixes confirmed present and behaviorally pinned; the Prim and Kruskal references confirmed from-scratch |
| R3 Edwards-Sokal | **PASS** | Partition function verified BY HAND on a 2-edge parallel graph at beta=0.5 (4 cosh(0.9) agreement at 2e-16); the SW step audited line-by-line (same-spin-only opens with no variate consumed on opposite-spin edges — verified at beta=50 and at underflow; ONE draw per FK cluster — cluster coherence verified across 2000 sweeps); the exposed BondSample confirmed to be the literal RNG output (bit-for-bit stream replay); worm separation confirmed at the code level (zero percolation imports in worm/) |

### Follow-ups applied from that review

1. The stale "R2/R3 samplers are Not implemented" sentence in the R1
   section of PERCOLATION.md corrected (they are delivered and consume
   the R1 contract).
2. The trapping x stop-policy cross-product folded into the shipped
   invasion suite (`trapping_matches_reference_under_every_stop_policy`,
   17th test) so the reviewer's external coverage now lives in-tree.
3. Notes recorded without action: the R3 per-sweep BondSample rebuild
   allocation (one allocation per sweep, ~88 µs/sweep at E=4095, no spec
   gate requires a probe — noted for a future optimization pass); the
   z_gap variance overestimate (conservative direction, both estimators
   separately gated); the InvasionThresholdOrder tail-repeat mean-bias
   caveat (documented; the scalars carry the true lengths).

---

## ADDENDUM 2 — independent review of P2 (post-reset)

| Stage | Independent verdict | Reviewer's independent confirmation |
|---|---|---|
| P2 Boolean | **PASS WITH NOTES** | Hand-constructed 4-particle adversarial fixture (direct / periodic-only / bitwise-tie-inclusive / just-outside) reproduced production's pair set exactly; the monodisperse `(r+r)² == (2r)²` identity fuzzed over 5M random radii with zero mismatches; Poisson CI formulas verified against the exact asymptotic μ₄ − σ⁴ = λ+2λ²; seven failure modes scratch-probed; the spanning-criterion seam safety proven (unfolding collapses seam-straddling raw extents) AND its guarantee boundary identified: the true bound is unwrapped extent ≥ L − 2r_max (not ≥ L), so "catches every wrapping cluster" requires c < 0.2L — demonstrated by a legal K₃ counterexample at contact 50 in L=100 that the gate's own parameters (c ≈ 0.006L) never approach; all doc numbers reproduced to the last digit |

### Follow-up applied

The two over-broad statements (test comment + PERCOLATION.md) qualified to
the correct bound "extent ≥ L − 2 r_max, safe at this gate's radii
(c ≈ 0.006 L)" — no functional change; the bracket gate is unaffected.
P2 promoted to Validated.


---

## ADDENDUM 3 — independent review of G2 and R2 (post-reset)

| Stage | Independent verdict | Reviewer's independent confirmation |
|---|---|---|
| G2 Achlioptas | **PASS WITH NOTES** | Hand trace of a 5-edge product m=2 process verified exactly (incl. self-loop scoring and both tie policies); an independent from-scratch contract driver written from the module-doc spec matched production BITWISE on 60/60 configurations (toy/chain/square x product/sum x both sampling x both ties x seeds, incl. steps/null_steps); the m=1 Fisher-Yates reduction mechanism confirmed by code read AND by the reviewer's own permutation implementation; the anchor's error model resolved (13.1σ is the correct SEM-combination — and the reviewer caught that the INTEGRATOR's earlier 64.5σ recomputation was itself an error, corrected in the note above); null-step public semantics audited (occupied winners compete, record added:false, leave the summary bit-identical) |
| R2 random-cluster q=1 | **PASS** | Exhaustive FK identity on the reviewer's OWN 4-edge triangle+pendant fixture (chosen so k(A) is not a function of |A|): all 16 weights both ways ≤ 2.8e-17, Z(p,1)=1 to 2.2e-16, heterogeneous Z=1 also by independent binomial convolution — with the reviewer's own first-pass hand slip caught by their machine cross-check, demonstrating the reference discipline; delegation audit confirmed (zero rng calls; heterogeneous path routes through the same sample_mask kernel); q≠1 rejection cites the roadmap; the 16-seed z-gate's even/odd disjoint RNG streams verified and re-run green at 48 seeds; the BondSample exposure path round-tripped on 20 draws |

### Follow-ups applied

- The integrator's mistaken 64.5σ note corrected in place (retained with
  the correction so the record stays honest).
- Notes recorded without action: N2 (a hypothetical release-mode
  pool-invariant violation would loop rather than error — unreachable
  under the verified draw contract, debug_assert-guarded); N3 (the anchor
  covers process-RNG variance on one fixed graph per N, not graph-to-graph
  resampling — immaterial at 13σ vs a 3σ gate); N4 (AchlioptasGrid::At
  positive runs unexercised); N5 (the two "independent paths" share the
  F2 draw kernel BY DESIGN — that is the delegation claim; independence
  lives in the analysis plumbing and RNG streams).

**G2 and R2 promoted to Validated. With this addendum every v1 roadmap
stage is delivered and independently reviewed (or Experimental-by-design
having passed review, or explicitly Not implemented): the v1 percolation
platform is complete.**
