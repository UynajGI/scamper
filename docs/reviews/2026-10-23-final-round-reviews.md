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

- The 13.1σ claim vs my 64.5σ recomputation from the same reported
  numbers: both reject "no growth" overwhelmingly; the agent likely
  propagated peak-location uncertainty rather than naive per-run stderr.
  The gate passes under either model. Recorded for transparency.

## Follow-ups applied in this commit

- None required for correctness. All five stages promoted from
  "Experimental — awaiting independent review" to "Integrator-verified —
  independent re-run recommended pre-merge" (a distinct, honest status;
  the support matrix and roadmap rows say exactly this).
