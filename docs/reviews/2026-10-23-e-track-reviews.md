# Independent Reviews — 2026-10-23 — E2 and E3 (commit `5a21962`)

Both stages were reviewed in parallel by independent reviewer agents (the
code was authored by parallel development agents; reviewers verified from
scratch). Full verbatim reports follow the summary table.

| Stage | Verdict | Reviewer's independent confirmation | Follow-ups applied |
|---|---|---|---|
| E2 bootstrap | **PASS WITH NOTES** | Hand traces (5-cycle + chord, r∈{0,1,2}) matched production exactly including per-round counts; clock audit (no within-round contamination representable; queue monotone); E1 stamp/queue reuse audited incl. a V=5/V=3 alternation probe; Holroyd anchor verified in the archived md including the paper's own slow-convergence note (ASA's 0.245±0.015 vs λ=0.548), p_hat·ln L recomputed to 0.2399/0.2398/0.2331 strictly inside (λ/4, λ); edge cases (parallel edges per incidence, self-loops never self-count, typed rejections) probed in a scratch crate | Combination counts corrected upward (720 mask×r; 15 vectors; ~4000 runs — docs had understated by omitting K5); the Criterion chain row's twofold queue deficit was recorded at density 0.5927 while the committed bench uses 0.1 — both numbers now labeled honestly; allocation probe extended to assert the QueueClosure clock (was asserted for SynchronousRounds only); `is_active` pre-first-run behavior documented |
| E3 k-core | **PASS WITH NOTES** | Hand traces (K4 + pendant path, k∈{2,3,4} incl. the mid-round-crossing cascade) matched exactly; FIFO↔synchronous-rounds equivalence re-proved and verified on 3200 random multigraph runs vs a self-written simultaneous full-scan reference (per-vertex deletion rounds, ties included); own two-method fixed-point solver reproduces λ₃=3.35091887 and ψ₃=0.451667 (measured 0.45309±0.00134, 1.06σ); Janson-Luczak Theorem 1.1 verified in the archived md including its G(n,m) extension; 400-realization C0 shape attack (round counts genuinely varying) passed | The kcore allocation probe only printed its zero result — now asserts (matching the bootstrap probe); an in-repo multigraph gate added (self-loop/parallel-edge hand cases + 120-round random multigraph fuzz vs an independent adjacency peeler; the reviewer's external fuzz had closed this gap but the repo now carries it); `residual_degree` panic bound documented |

Both stages promoted from "Experimental — awaiting independent review" to
reviewed status; production availability remains conditioned on the `dev`
merge.

## E2 verbatim reviewer report

> **VERDICT: PASS WITH NOTES.** No blockers, no majors. The physics core is
> correct (verified independently by hand traces, adversarial probes, and
> reruns). All findings are documentation-accuracy items, mostly in the
> conservative direction (docs understate actual test coverage).
>
> All nine roadmap/spec gates PASS (exhaustive 2^V masks × r∈0..=deg+1 on
> P5/C5/star(4)/K4/K5 vs independent full-scan references, both clocks; r=0
> activates all; r>max_degree preserves seeds; order independence tested
> across FIFO/one-at-a-time/shuffled/synchronous on 7 graphs incl.
> multigraphs; explicit seed/threshold/clock API; E1 constraint upheld;
> Holroyd anchor reproduced bit-for-bit; 0-allocation probe with honest
> warmup discipline).
>
> Independent checks: hand traces matched exactly (closure {0} for r=2
> single seed; traj [2,3] for seeds {0,1}; [1,4,5] rounds / [1,2,3,4,5]
> events for r=1); a discriminating probe confirmed the synchronous clock
> has no within-round contamination ([1,4,5,6] not [1,4,6]); stamp
> wraparound, queue reset, counter re-zero verified; the archived paper
> itself documents the slow convergence (finite-size estimate 0.245±0.015
> vs λ=0.548, simulations up to L=28,800), so the test's statistical-
> bracket framing is the paper's own state of affairs; p_hat·ln L =
> 0.2399/0.2398/0.2331 recomputed and strictly inside (λ/4, λ).
>
> Issues (all applied): doc combination counts understated (528→720, 12→15
> vectors, ~2600→~4000 runs); the Criterion chain row mixed a 0.5927-density
> calibration with the committed 0.1 bench (both now labeled); the
> allocation probe asserted only SynchronousRounds (QueueClosure now
> asserted); is_active pre-first-run semantics undocumented (now
> documented). Notes: heterogeneous per-round trajectories not reference-
> checked in-suite (reviewer closed externally — passes); torus vs open
> square labeled honestly.

## E3 verbatim reviewer report

> **VERDICT: PASS WITH NOTES.** No blockers, no majors. Two minor
> enforcement gaps and five notes. Every roadmap E3 gate is present,
> physically correct, and reproduced.
>
> Independent checks: hand traces on K4 + pendant path (k=2/3/4, the last
> exercising mid-round crossing) matched production exactly (summary,
> in_core, deletion_round, frozen residuals, deletion order); multigraph
> hand cases (three parallel edges sustain a 3-core; self-loop = 2
> incidences, its edge survives in core_edge_count) match the documented
> convention; the FIFO↔synchronous equivalence was re-proved and verified
> on 3200 runs over 400 random multigraphs against a self-written
> simultaneous full-scan reference; order independence verified over 90
> distinct deletion orders; own two-method solver: λ₃ = 3.35091887 at
> µ* = 1.79328213, ψ₃(µ₃(3.500919)) = 0.451667 — measured 0.45309±0.00134
> (1.06σ); the archived paper's Theorem 1.1 states exactly the used system
> and explicitly covers G(n,m); KCoreMC's draw parity with
> StaticPercolationMC is structural (32-sweep bitwise test passed); a
> 400-realization shape attack confirmed C0 shape stability with genuinely
> varying round counts and exact pad semantics.
>
> Issues (applied): the kcore allocation probe printed without asserting
> (now asserts); no in-repo multigraph corpus (added: hand cases + 120-
> round fuzz vs an independent adjacency peeler); residual_degree panic
> bound undocumented (documented). Notes: "0.00000±0.00000" is zero sample
> variance, not emptiness proof (n stated, largest-core-0 also reported);
> the log's "Janson-Luczak/PSW: 3.350919" is the formula evaluation, not a
> printed paper digit (accurate as an identity check); subset initial
> activity only via the core function (documented, deliberate); measure()
> before sweep records zeros (consistent with siblings).
