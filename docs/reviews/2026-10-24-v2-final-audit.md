# v2 final audit (N5)

Status: **complete for delivery** — all 16 v2 stages are implemented,
wired, documented, and covered by tests; 10 are Validated and 6 are
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
| C1 correlated fields | `percolation/correlated.rs` | unit tests + adapter parity | Experimental | field/variance/marginal/a=inf/a=0 gates green; binary correlation and independent review open |
| C2 shifted thresholds | `percolation/correlated_thresholds.rs` | fast smoke + ignored large-L bracket | Experimental | large-L bracket contains Pruessner 0.54299; review pending |
| C3 Weinrib–Halperin map | `percolation/weinrib_halperin.rs` | analytic branch gate | Experimental | exact map/table branches; statistical nu certification open |
| Q1 Sweeny sampler | `percolation/sweeny.rs` | `percolation_sweeny.rs` (3) | Experimental | exact 2^E FK marginals z<=1.36, q=1 bitwise, fuzz; review pending |
| Q2 Chayes–Machta | `percolation/chayes_machta.rs` | q=2 agreement | Experimental | agrees with Q1 within 4 sigma; detailed balance/R3 comparison open |
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
- **Time constant**: white-noise calibration of the C2 estimator gives
  `p_c(inf) ~ 0.580` against the exact iid 0.592746, exposing a
  finite-size estimator bias separate from the field generator.

## Open items (explicit, not hidden)

1. **C2 threshold gate now passes** (the earlier limitation was
   under-sampling and a bootstrap-CI bug): the independent review found
   that the bootstrap SD was divided by sqrt(B); the corrected CI is
   (0.5355, 0.5558), so p_c(inf)=0.5457 differs from Pruessner's 0.54299
   by ~0.65 sigma and the reference value lies inside the CI. Anchor
   details corrected (Pruessner uses L=64-8192, 10^5 replicas, and a
   wrapping threshold; this module peaks the order-parameter
   susceptibility). Remaining: independent re-review. C3's nu estimate at
   a=1.5 is 0.934 (CI 0.516-1.352), still not certifying the flip.
2. **C3 statistical certification**: the analytic map and branch values
   are exact. A release L=128/256/512, R=256 diagnostic gives p_c =
   0.5430/0.5317/0.5393 at a=1.5/1.25/1.0 (a=1.5 matches 0.54299) but
   log-derivative nu = 0.934/1.020/0.804 with wide CIs, systematically
   below the predicted 4/3/8/5/2. The estimator is not reliable at these
   sizes, so the flip is not certified; larger L or a better estimator is
   needed.
3. **Q2/Q3**: Q1 has the rate-level detailed-balance audit; Q2 now has
   q=2/3 exact-histogram and colour-marginal invariance plus the q=2 R3
   Edwards–Sokal comparison. A full joint detailed-balance proof remains
   optional/open. Q3's L=16 scan separates q<=4 (bimodality <5/9) from
   q=8 (>5/9), but q=6 is unresolved and the statistic is descriptive, not
   a latent heat.
4. **F3 minors**: the target-local saturation flag, diagonal-bound
   wording, vertex-prefix aging doc, post-rejection weight state, and
   parallel/source==target tests are FIXED. Still open: `geodesic_path`
   validation wording and the estimator's inf/NaN domain.
5. **Reviews**: C1/C2/C3/Q1/Q2/Q3 still need independent review before
   any Validated promotion.

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
  with `p_c(inf)=0.5457`, and the Q3 L=16 q=2/4/6/8 scan passes.
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

All 16 v2 stages are delivered and audited. 10 are Validated after
independent review; 6 are Experimental because their strongest claims
(C2 table reproduction, C3 statistical nu, Q2 detailed balance, Q3
latent heat, C1 binary correlation, and independent review) are not yet
established. This matches the v1 M6 convention: Experimental-by-design
or explicitly limited rows do not block the v2 delivery audit, but they
are not counted as Validated.
