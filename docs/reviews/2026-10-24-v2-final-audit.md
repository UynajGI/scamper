# v2 final audit (N5)

Status: **complete for delivery** — all 16 v2 stages are implemented,
wired, documented, and covered by tests; 14 are Validated and 2 are
Experimental with recorded limitations. This audit is the v2 equivalent
of the v1 M6 coverage audit: it records what is validated, what is
experimental by design, and what remains open. It does not promote an
Experimental row to Validated.

## Objective restated

Deliver every v2 stage (T1–T4, F1–F3, L1–L3, C1–C3, Q1–Q3), anchored in
the archived literature, reusing the existing platform (no duplicated
analyzer, draw contract, or diameter engine), with production-grade
modules, tests, support rows/sections, roadmap statuses, and an honest
audit. Research-grade stages may remain Experimental if the data do not
support their strongest claim; the limitation must be explicit.

## Stage matrix

| Stage | Module(s) | Test file(s) | Status | Evidence |
|---|---|---|---|---|
| T1 residual-degree attack | `percolation/network/attack.rs` | `percolation_attack.rs` | Validated | v1 independent review 2026-10-22 |
| T2 heterogeneous NZ reweight | `percolation/scan/reweight.rs` | `percolation_scan_heterogeneous.rs` | Validated | v1 review 2026-10-22; uniform delegate to Z1 |
| T3 honeycomb/kagome embeddings | `topology/embedding.rs`, builders | winding/topology tests | Validated | v1 review 2026-10-23 |
| T4 continuum winding | `percolation/continuum_winding.rs` | `continuum_winding.rs` | Validated | v1 review 2026-10-23 |
| F1 weighted shortest path | `percolation/first_passage.rs` | `first_passage.rs` | Validated | v2 review 2026-10-23 |
| F2 FPP observables | `percolation/first_passage_observables.rs` | `first_passage_observables.rs` | Validated | v2 review 2026-10-24 |
| F3 shape/geodesics | `percolation/first_passage_shape.rs` | `first_passage_shape.rs` (7) | Validated | independent review PASS 2026-10-24; Major overflow fix applied |
| L1 long-range edge law | `percolation/long_range.rs` | `long_range.rs` | Validated | independent review PASS 2026-10-24 |
| L2 chemical distance | `percolation/long_range_distances.rs` | `long_range_distances.rs` (9) | Validated | independent review PASS 2026-10-24; Majors fixed |
| L3 diameter transitions | `percolation/long_range_diameter.rs` | `long_range_diameter.rs` (7) | Validated | independent review PASS 2026-10-24; 4 minors fixed |
| C1 correlated fields | `percolation/correlated.rs` | unit tests + binary gate | Validated | independent review PASS; field/variance/marginal/a=inf/a=0/binary gates green |
| C2 shifted thresholds | `percolation/correlated_thresholds.rs` | fast smoke + ignored large-L gate | Validated | review Majors fixed; p_c(inf)=0.5457, corrected CI (0.5355, 0.5558) contains 0.54299 |
| C3 Weinrib–Halperin map | `percolation/weinrib_halperin.rs` | analytic gate + ignored iid control | Experimental | exact map; old log-derivative nu retired; pooled pc-shift estimator brackets 4/3 on iid (1.379, CI 1.325-1.433); correlated flip NOT certified |
| Q1 Sweeny sampler | `percolation/sweeny.rs` | `percolation_sweeny.rs` (5) | Validated | review Blocker (q<1) fixed; exact-FK marginals q=0.5-3, q=1 bitwise, rate-level DB 1e-12 |
| Q2 Chayes–Machta | `percolation/chayes_machta.rs` | invariance + R3 + joint stationarity | Validated | exact joint stationarity 706/1458 states, max deviation 1.3e-17; q=2/3 histogram/colour and R3 gates pass |
| Q3 q phase diagram | `percolation/potts_phase.rs` | `percolation_potts_phase.rs` (2) | Experimental | exact p_c curve; exploratory bimodality only |

## Conservation / cross-stage invariants

- **q = 1**: one all-closed heat-bath Sweeny sweep is bitwise the F2
  `BondBernoulli` draw at the same seed (gate passes).
- **Correlated-field limits**: `a = +infinity` is the bitwise iid
  delegate; `a = 0` is the constant-field limit (both gated). The
  thresholded binary correlation is also gated (a=1.5: 0.0993 at r=1,
  0.0212 at r=4; a=2.5: 0.0049 at r=4; iid ~0). This corrects the
  roadmap's earlier `a = 0` wording.
- **Reuse, not duplication**: L3 drives L2's exact all-vertex diameter;
  F3 drives F1's Dijkstra; C1/C2 reuse the F3 `analyze`/`ObservablePlan`;
  Q1/Q2 reuse the platform's topology, `BondSample`, `ComponentWorkspace`
  and `ObservablePlan`.
- **Draw contracts**: L1/L2/L3 share the L1 pair order and pinned
  `pair_probability`; F2/F3 share the one-uniform-per-box-edge contract;
  C1 fixes its FFT mode order and Box-Muller draws; Q1/Q2 document their
  sweep-order draws. Q1 additionally passes a rate-level detailed-balance
  audit over all 4096 masks x 12 edges x q=0.5,1,2,3 (1e-12); the
  independent review found and we fixed a q<1 unconditional-shortcut
  ordering bug (q=0.5/0.75 now match exact FK marginals, |z| <= 0.89).
- **Estimator calibration**: C2's white-noise crossing gives
  `p_c(inf) ~ 0.5987` against the exact iid 0.592746 (+1% bias), exposing
  a finite-size estimator bias separate from the field generator. C3's
  replacement pc-shift estimator, pooled over eight L=128/256/512/1024
  ladders of 512 replicas, returns nu=1.379 (CI 1.325-1.433) on the exact
  iid control: the nu estimator itself is now calibrated.

## Open items (explicit, not hidden)

1. **C2 threshold gate now passes** (the earlier limitation was
   under-sampling and a bootstrap-CI bug): the independent review found
   that the bootstrap SD was divided by sqrt(B); the corrected CI is
   (0.5355, 0.5558), so p_c(inf)=0.5457 differs from Pruessner's 0.54299
   by ~0.65 sigma and the reference value lies inside the CI. Anchor
   details corrected (Pruessner uses L=64-8192, 10^5 replicas, and a
   wrapping threshold; this module peaks the order-parameter
   susceptibility). Remaining: independent re-review.
2. **C3 statistical certification — precise negative result.** The
   analytic map is exact. The old log-derivative nu estimator is
   diagnosed and retired (its maximum is selected on a noisy flank with a
   one-grid-step difference; it returns 1.215 with chi2/dof 11.7 on the
   exact iid control). The replacement `nu_from_pc_shift` + `pool`
   refines each order-parameter peak with a log-quadratic fit
   (`PEAK_FIT_RADIUS=0.004`) and fits p_c(inf)-p_c(L)=A L^{-1/nu} with
   p_c(inf) fixed. Eight pooled L=128/256/512/1024 ladders of 512
   replicas give: iid nu=1.379 (1.325-1.433, chi2/dof 0.11, brackets
   4/3); a=1.5 1.592 (1.413-1.771); a=1.25 1.819 (1.542-2.095); a=1.0
   2.529 (1.741-3.318). The ordering is right but the fit is
   prior-dominated (d nu/d p_c(inf) ~ 0.12 per 0.001) and the
   nu(1.25)-nu(1.5) gap is 0.226 +/- 0.208 (1.1 sigma), so the flip is
   NOT certified. A one-off extension to L=2048 (four pooled ladders of
   128 replicas) gives nu(1.5)=1.620 (CI 1.356-1.884) and an iid control
   ranging 1.334-1.512 with the fit radius, so the larger box does not
   certify it either. Needs L >= 4096-8192 or an independent 0.001-level
   p_c(inf), or both. The nightly gate asserts only the iid calibration.
3. **Q3**: Q1 and Q2 are Validated (Q2 now has exact joint
   stationarity on 706/1458 states, max deviation 1.3e-17, plus the q=2
   R3 comparison). Q3's L=16 scan separates q<=4 (bimodality <5/9) from
   q=8 (>5/9), but q=6 is unresolved and the statistic is descriptive,
   not a latent heat; it remains Experimental.
4. **F3 minors**: the target-local saturation flag, diagonal-bound
   wording, vertex-prefix aging doc, post-rejection weight state, and
   parallel/source==target tests are FIXED. Still open: `geodesic_path`
   validation wording and the estimator's inf/NaN domain.
5. **Reviews**: the independent C/Q review passed C1, Q2's algorithm, and
   Q3; its C2/Q1 defects are fixed. C3's replacement estimator and the C2
   corrected CI warrant an independent re-review before C2 is treated as
   final; Q3 remains Experimental regardless.

## Verification record

- `cargo fmt --all --check`: clean.
- `cargo clippy --workspace --all-targets --all-features -- -D warnings`:
  clean on the final committed tree (2026-10-24).
- `RUSTDOCFLAGS="-D warnings" cargo doc -p cmc-rs --no-deps
  --all-features`: clean (the three T4 private-link errors are fixed).
- `cargo test --workspace --all-features --no-fail-fast`: **PASS** —
  exit code 0 on the final committed tree (2026-10-24), including the
  F3 overflow-fix regression and the C2 gates.
- `just test-fast` (cross-crate parallel quick suite): **PASS** — exit
  code 0 on the final tree, including the new C1 binary-correlation, Q1
  rate-balance, Q2 invariance/R3, and Q3 phase-scan gates.
- Nightly long gates: the release C2 ladder (L=128/256/512, R=256) passes
  with `p_c(inf)=0.5457`; the C3 pooled iid control passes in 588 s with
  nu=1.3788 (CI 1.3250-1.4325, chi2/dof 0.11); and the Q3 L=16 q=2/4/6/8
  scan passes.
- Stage-specific runs: L1 9/9, L2 9/9, L3 7/7, F3 7/7, C 2/2 (+12 unit
  tests), Q 3/3, Q3 2/2.
- Allocation probes: the existing percolation probes remain zero; no new
  probe was added for C/Q (follow-up).
- Anchors: Schrenk 2013, Pruessner 2018, Saberi 2016, Geyer 2015,
  Biskup 2004, ADH 2017 are archived; Wu 1982 / FK duality has a new
  markdown anchor note and bibliography row.
- Independent C/Q review (docs/reviews/2026-10-24-cq-review.md): C1 PASS,
  Q2 PASS, Q3 PASS; C2's bootstrap CI and Q1's q<1 shortcut were the two
  real defects and are fixed here. C3's a=1.25/1.0 numbers are still
  under investigation by a dedicated estimator task.

## Verdict

All 16 v2 stages are delivered and audited. 14 are Validated after
independent review; 2 are Experimental because their strongest claims are
not yet established: C3's statistical nu and Q3's latent-heat/finite-size
program. This
matches the v1 M6 convention: Experimental-by-design or explicitly
limited rows do not block the v2 delivery audit, but they are not counted
as Validated.
