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
| C2 shifted thresholds | `percolation/correlated_thresholds.rs` | machinery/limitation gate | Experimental | does NOT reproduce Pruessner Table 1 (see open items) |
| C3 Weinrib–Halperin map | `percolation/weinrib_halperin.rs` | analytic branch gate | Experimental | exact map/table branches; statistical nu certification open |
| Q1 Sweeny sampler | `percolation/sweeny.rs` | `percolation_sweeny.rs` (3) | Experimental | exact 2^E FK marginals z<=1.36, q=1 bitwise, fuzz; review pending |
| Q2 Chayes–Machta | `percolation/chayes_machta.rs` | q=2 agreement | Experimental | agrees with Q1 within 4 sigma; detailed balance/R3 comparison open |
| Q3 q phase diagram | `percolation/potts_phase.rs` | `percolation_potts_phase.rs` (2) | Experimental | exact p_c curve; exploratory bimodality only |

## Conservation / cross-stage invariants

- **q = 1**: one all-closed heat-bath Sweeny sweep is bitwise the F2
  `BondBernoulli` draw at the same seed (gate passes).
- **Correlated-field limits**: `a = +infinity` is the bitwise iid
  delegate; `a = 0` is the constant-field limit (both gated). This
  corrects the roadmap's earlier `a = 0` wording.
- **Reuse, not duplication**: L3 drives L2's exact all-vertex diameter;
  F3 drives F1's Dijkstra; C1/C2 reuse the F3 `analyze`/`ObservablePlan`;
  Q1/Q2 reuse the platform's topology, `BondSample`, `ComponentWorkspace`
  and `ObservablePlan`.
- **Draw contracts**: L1/L2/L3 share the L1 pair order and pinned
  `pair_probability`; F2/F3 share the one-uniform-per-box-edge contract;
  C1 fixes its FFT mode order and Box-Muller draws; Q1/Q2 document their
  sweep-order draws.
- **Time constant**: white-noise calibration of the C2 estimator gives
  `p_c(inf) ~ 0.580` against the exact iid 0.592746, exposing a
  finite-size estimator bias separate from the field generator.

## Open items (explicit, not hidden)

1. **C2 threshold table**: at a=1.5 with 32 replicas, the susceptibility
   peak over L=64/128/256 converges to 0.508–0.513, while
   largest-fraction crossings give 0.558/0.572; Pruessner's 0.54299 lies
   between them. The threshold depends on the full discrete FFM
   correlation function, not only the exponent. Follow-up: implement
   Pruessner's exact discrete generator or calibrate a crossing estimator.
2. **C3 statistical certification**: the analytic map and branch values
   are exact, but the finite-size `nu(a)` estimate is not yet precise
   enough to certify the a=1.5/a=1.25 flip.
3. **Q2/Q3**: Q2 needs a full detailed-balance audit and the q=2 R3
   comparison; Q3 needs larger L and a two-mode latent-heat estimator.
   The Q3 bimodality values (0.362 at q=2, 0.505 at q=8, L=8) are
   descriptive only.
4. **F3 minors**: global saturation flag, diagonal bound wording, aging
   vertex-vs-edge comparison, untested parallel/source==target cases,
   `realize_with_weights` post-rejection state, and `geodesic_path` doc.
5. **Reviews**: C1/C2/C3/Q1/Q2/Q3 still need independent review before
   any Validated promotion.

## Verification record

- `cargo fmt --all --check`: clean.
- `cargo clippy --workspace --all-targets --all-features -- -D warnings`:
  clean.
- `cargo test --workspace --all-features --no-fail-fast`: **PASS**
  (recorded below).
- Stage-specific runs: L1 9/9, L2 9/9, L3 7/7, F3 7/7, C 2/2 (+12 unit
  tests), Q 3/3, Q3 2/2.
- Allocation probes: the existing percolation probes remain zero; no new
  probe was added for C/Q (follow-up).
- Anchors: Schrenk 2013, Pruessner 2018, Saberi 2016, Geyer 2015,
  Biskup 2004, ADH 2017 are archived; Wu 1982 / FK duality has a new
  markdown anchor note and bibliography row.

## Verdict

All 16 v2 stages are delivered and audited. 10 are Validated after
independent review; 6 are Experimental because their strongest claims
(C2 table reproduction, C3 statistical nu, Q2 detailed balance, Q3
latent heat, C1 binary correlation, and independent review) are not yet
established. This matches the v1 M6 convention: Experimental-by-design
or explicitly limited rows do not block the v2 delivery audit, but they
are not counted as Validated.
