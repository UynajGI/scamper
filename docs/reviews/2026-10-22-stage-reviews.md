# Independent Review — 2026-10-22 — `feat/cmc-percolation`

Four stages were reviewed in parallel by independent reviewer agents that
did not author the code. Each reviewer ran the gates themselves, re-derived
the headline numbers from scratch (their own scripts, no project code), and
audited the implementations against the roadmap specs. Full reviewer reports
are archived verbatim in this directory; verdicts and follow-up actions are
summarized below.

| Stage | Commit under review | Verdict | Reviewer's independent confirmation | Follow-ups applied |
|---|---|---|---|---|
| W2 winding analyzer | `8357965^..00b6a5b` | **PASS WITH NOTES** | 2x2 P(wrapping) = 187/256 and 3x3 = 186684/262144 re-derived by two from-scratch methods (DSU with unbounded cell offsets; lifted BFS); union-find potential algebra and SL(2,Z) echelon folding hand-verified; i32 checked arithmetic confirmed; overflow bound re-derived | Asserted 187/256, 252/256, 186684/262144, 259658/262144 as exact regression anchors (previously printed only); fixed benchmark ratio prose (~2.0-2.6x / ~2.8-3.5x) |
| Z1 Newman-Ziff scans | `8357965` | **PASS** | Reweight recurrence verified against the archived Newman-Ziff paper Eq. (1): weights match log-space binomials to 9.5e-14 up to n = 10^6; from-scratch flood-fill prefix parity on multigraphs; speedup reproduced (1.80-4.87x vs claimed 1.83-4.65x) | `trials + 1` -> saturating_add for consistency; capacity table marked as an 8-of-10 excerpt |
| C0 vector estimates | `cf701b8` | **PASS WITH NOTES** | Hand-derived per-element statistics bit-exact vs scalar accumulators; legacy JSON parsed; old code proven to silently zip-truncate in release (probe), new code errors; multispin migration preserves the same physical quantity | **Fixed merged stderr** (pre-existing flaw, review issue 1): pooled-bin variance `sum_i n_i s_i^2 / N^2` replaces the `/n_runs^2` form that underestimated merged error bars by sqrt(n_runs); regression test added. **Fixed per-element tau**: now estimated from completed bins only, exactly matching the scalar convention (outlier partial bins no longer inflate tau); parity test added |
| N2 random failure + degree observables | `bb11fb3` | **PASS WITH NOTES** | Cohen q_c re-derived (kappa = 3.5000 -> 0.4000 ideal; empirical 0.3996); Callaway Eq. (5) fixed point re-solved independently (S = 0.1857/0.4258/0.6374 vs MC 0.1857/0.4261/0.6377); star-graph degree observables hand-enumerated to exact fractions (GiantMembershipByDegree[1] = 47/128, tie-policy-sensitive); V^3 gate threshold recomputed (208063 boundary) | PERCOLATION.md bracket transcription corrected (S = 0.1851, not 0.1857); V^3-gate rejection branch now tested at the exact boundary |

N3 (targeted degree attack + robustness curves) was reviewed later the same
day after its delivery commit `5f2dd2a`:

| Stage | Commit under review | Verdict | Reviewer's independent confirmation | Follow-ups applied |
|---|---|---|---|---|
| N3 targeted attacks | `5f2dd2a` | **PASS WITH NOTES** | End-to-end rebuild of both anchor graphs in a throwaway crate; independent Callaway Eq. (14) solver reproduces all eight cutoff values digit-for-digit; Cohen f_c values reproduced exactly; Fisher-Yates uniformity verified by exhaustive enumeration; path6 attack hand-simulated | Corrected the "0.02-0.03" misquote of Cohen et al. (their continuum value at alpha=2.5 is ~0.056; the range belongs to alpha~3), aligned the roadmap bracket wording (0.3 f_c), removed a stray pdftotext dump, documented the config-model repair bias |

The network edge-list input was reviewed after its delivery commit
`539a16f`:

| Stage | Commit under review | Verdict | Reviewer's independent confirmation | Follow-ups applied |
|---|---|---|---|---|
| Network edge-list input | `539a16f` | **PASS WITH NOTES** | Throwaway-crate parser attack (CRLF, Unicode whitespace, overflow, token counts — all typed or correct); CSR semantics verified empirically (self-loop double incidence, parallel edges distinct, file-order emission); benchmark re-run reproduced 5.5 M / 1.26 M edges/s | `vertex_count = 0` panic replaced by a typed `InvalidVertexCount` error; duplicate detection now covers repeated self-loops; pre-declared review claim replaced by this archived report; doc comments aligned |

D1 (directed topology) was reviewed after its delivery commit `c2ea61f`:

| Stage | Commit under review | Verdict | Reviewer's independent confirmation | Follow-ups applied |
|---|---|---|---|---|
| D1 directed topology | `c2ea61f` | **PASS WITH NOTES** | Builder hand-traced on a parallel-arc + self-loop multigraph (exact arrays predicted and confirmed); every validator attack predicted in writing then confirmed empirically in debug and release; trait separation compile-checked; the arc-table SIMD explanation verified from struct-layout evidence (Arc 16B vs Bond 32B) | Checked arithmetic for a magic-value declared count that could overflow a validation error path; stale F1 sentence updated; per-element figure corrected to the recorded 0.72 ns; default()/validate allocation behavior documented |

All reviewed stages (W2, Z1, C0, N2, N3, network input, D1) are promoted
from "Experimental — awaiting independent review" to reviewed status; their production availability remains conditioned
on the branch merging into `dev`, per the platform's standing convention.

## The one substantive bug the reviews caught

`Results::merge` (scalar, complex, and — after C0 — array observables)
computed the combined standard error as `sqrt(sum_i s_i^2 / R^2 / N)` with R
runs and N total bins. The variance of the pooled-bin weighted mean is
`sum_i n_i s_i^2 / N^2`, so equal-size runs got error bars `sqrt(R)` too
small — an overconfidence that would have propagated into every merged
percolation curve the platform produces. The formula pre-dates C0 (C0's
compatibility promise was honored); the reviewer flagged it as major because
C0 extends it to the curve outputs N3/D3/E2/E3/G2/G3 depend on. Fixed in
this commit with an exact-value regression test; MPI/PT merged results now
report correct (larger) error bars.

## Reviewer notes accepted without code change

- W2: the adapter's `measure()` panics loudly (never wraps) on adversarial
  embedding tables that pass the static magnitude gate — documented in
  `with_embedding`; builder embeddings cannot trigger it.
- W2: `PeriodicEmbedding` is effectively sealed (`DirectedDisplacement::new`
  is crate-private), so the adversarial surface is internal-only.
- C0: shape-mismatched array observables vanish from merged results with a
  `tracing::warn` (documented); autocorr floor conventions differ across
  pre-existing estimators (merge.rs tau=0-for-uncorrelated vs accumulate-side
  floor 1.0) — pre-existing, not worsened.
- C0: `finalize_split` intentionally stops reporting array observables as
  collapsed scalars; in-repo consumers migrated, external JSON consumers of
  collapsed values would see the change (documented in VALIDATION/README).
- N2: the relabeling gate is single-run-per-labeling (real but the weakest
  N2 statistical gate); the pooled mean-z gate is diluted by trivially-zero
  empty bins (the per-observable |z| < 4 gate is unaffected).
- N2: `record_degree_observables` with mismatched scratch/profile panics on
  index bounds rather than erroring — acceptable for its documented use.

## Reproduction

```bash
cargo test -p carlo-rs --test suite                                   # 310 pass
cargo test -p carlo-rs --features hdf5 --test suite                   # 321 pass
cargo test -p cmc-rs --test suite -- percolation_winding              # 16 pass, anchors asserted
cargo test --release -p cmc-rs --test suite percolation_network       # 6 pass (incl. V^3 boundary)
just mpi-test 2                                                       # 5 MPI suites pass
```
