# v3 physics-correctness review

Date: 2026-10-30. Scope: the physics-bearing claims added by the v3 stack
(A5e, X2b, X4a) plus the research records R2/R3. Refactors A3c/A4c/A6a are
also checked for golden/reference-order conservation. This is an internal
review with independent numerical controls; it is not a substitute for the
independent per-stage reviews required by the roadmap.

## 1. Verdict

| Stage | Verdict | Basis |
|---|---|---|
| A3c C0 writer/recording | Pass | Pure refactor; C0 contents and golden output unchanged; `FixedSeries` padding tests. |
| A4c policy traits | Pass | Closed enums stay the defaults; external reimplementations reproduce built-in trajectories bitwise; a deliberately different rule changes the reference order. |
| A5e FSS bracket | Pass as a finite-size calibration | Generic-toolkit white-noise crossing 0.5995 at L=64/128 (v2 0.5987, exact 0.592746); seed spread at L=32/48/64 is ~0.02. Not a precision measurement. |
| X2b multiplex bracket | Pass with a strengthened realism gate | Cascade semantics match the two-layer mutual-component definition; the 0.2-fraction crossing approaches Buldyrev's 2.4554 from below as N grows; realized mean degree now asserted. |
| X4a contact process | **Corrected** | The event implementation is correct (independent first-reaction simulator agrees), but the original 0.535 crossing was presented as a threshold. It is a finite-time transient: T=100 extinguishes N<=2000 at lambda=0.55. Tests/docs now state only a finite-time activity bracket and include a negative control. |
| R2 two-mode q-ladder | Pass as exploratory only | q=8 separation exceeds q=2 across three independent seeds; no latent-heat or FSS certification is claimed. |
| R1/R3 | Records only | Experimental status and negative results unchanged; no new physics claim. |

## 2. Independent checks

### X4a — contact process

- **Independent algorithm**: a test-local first-reaction simulator (one
  exponential clock per possible event, minimum fires) agrees with the
  library direct kernel:
  - lambda=0.20: direct 0/6, first-reaction 0/6 survivors;
  - lambda=0.90: direct 4/6, first-reaction 6/6 survivors (both finite-time
    active).
- **Horizon sensitivity**: on 3-regular configuration-model graphs,
  lambda=0.55 at N=1000/2000 has survival 1.0 at T=40 but 0.0 at T=100;
  lambda=0.50 also extinguishes by T=100. The committed negative control
  (N=200, lambda=0.55, T=100) observes 0/4 survivors.
- **Conclusion**: the 0.05-density crossing 0.535 at T=15 is a finite-time
  observable, not a quasi-stationary or thermodynamic threshold. The tree
  branching value 0.5 remains only a reference line. The contact-process
  event rates and draw contract are correct; the limitation is statistical
  horizon, not the implementation.

### X2b — two-layer cascade

- The cascade definition and the independent reference-map tests are
  unchanged and pass.
- Committed 10-seed ladder (N=1024, c step 0.1): 0.2-fraction crossing
  c=2.359; a 20-seed review ladder at N=512/1024/2048 gives
  2.361/2.426/2.414. The crossing approaches Buldyrev's c*=2.4554 from
  below, consistent with a finite-size first-order transition.
- The configuration-model layers realize the requested mean degree (the
  review initially printed half the value by averaging `E/N` instead of
  `2E/N`; the committed test now asserts `(E_a+E_b)/N = c` within 0.25).

### A5e — generic FSS

- Fast L=32/48/64 ladder, three review seeds: crossing means
  0.5984/0.6113/0.6179 against the exact 0.592746; the committed
  per-crossing tolerance 0.06 / mean tolerance 0.03 covers this spread.
- Nightly L=64/128 single-seed calibration: 0.5995, within the committed
  0.02 bracket. This is a deterministic calibration, not an ensemble CI.

### R2 — two-mode q-ladder

- Three seeds at L=12: q=8 separation 3.19–3.75, q=2 2.61–2.76; the
  asserted q=8 > q=2 and q=8 > q=1 ordering is stable. The estimator is
  still a definitional 2-means proxy, not a latent-heat measurement; Q3
  remains Experimental.

## 3. Required changes applied

1. X4a test renamed/scoped to finite-time activity; threshold claim removed;
   independent first-reaction reference and a long-horizon negative control
   added.
2. X2b test now asserts the realized mean degree of the generated layers and
   uses 10 seeds on a finer c grid.
3. Roadmap, `PERCOLATION.md`, and the v3 audit text updated to match the
   corrected X4a/X2b claims.

## 4. Merge recommendation

PR #4 remains open. After this correction commit, the physics claims in the
PR are scoped honestly: A5e is a calibration, X2b is a finite-size bracket,
X4a is a finite-time activity bracket, and R1–R3 remain Experimental. The
independent stage reviews and the X4b temporal-network work remain open.
No automatic merge is performed.
