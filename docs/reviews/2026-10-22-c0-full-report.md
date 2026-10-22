# C0 (Carlo.rs Vector Estimates + Provenance) — Full Reviewer Report (verbatim)

> Archived from the independent reviewer agent, 2026-10-22. Reviewed commit
> `cf701b8`; probes ran from `git archive` snapshots. See
> `2026-10-22-stage-reviews.md` for follow-ups applied.

**VERDICT: PASS WITH NOTES**

## Test gates (run by reviewer)

| Gate | Result |
|---|---|
| `cargo test -p carlo-rs --test suite` | **308 passed, 0 failed** |
| `cargo test -p carlo-rs --features hdf5 --test suite` | **319 passed, 0 failed** |
| `cargo test -p cmc-rs multispin --release` | **3 passed** (9.23 s) |
| Reviewer's 7-test adversarial probe, NEW code, `--release` | 7/7 passed |
| Same probe, OLD code (`cf701b8~1`), `--release` | proved silent truncation |
| Same probe, OLD code, debug | panicked on `debug_assert` (expected) |

## Roadmap §8 exit gate

- Per-element independent mean/stderr/autocorr: **PASS** — hand-derived
  2-element/5-sample/binsize-2 case matched to 1e-12 and bit-exact vs
  independent scalar accumulators.
- JSON round-trip + old fixtures: **PASS** — reviewer's hand-written
  pre-C0-shaped JSON (with `complex_observables`, without
  `array_observables`/`attributes`) parsed with empty arrays/attributes;
  ArrayEstimate round-trip bit-exact.
- ResultTools/readback/merge preserve shapes: **PASS**.
- Scalar byte/semantic compatibility: **PASS** — scalar-only `to_json`
  byte-identical (modulo volatile version/timestamp); `merge_weighted` and
  `variance_ratio_autocorr` are verbatim factor-outs of pre-C0 code
  (git-diff verified); two intentional documented changes (complex stderr
  weighting; empty samples no longer auto-register).

## Key independent checks

- **Hand-derived statistics**: element 0 [1,3,5,7,9] → bins [2,6]+partial
  [9] → mean 17/3, stderr √37/3; element 1 [10,12,20,16,30] → [11,18,30] →
  mean 59/3, stderr √277/3 — matched the code exactly and bit-for-bit vs
  scalar accumulators (`to_bits()` equality).
- **merge_weighted derivation (the finding)**: for R independent equal runs
  (n bins, variance s²) the code produced merged stderr² = s²/(R²n) while
  the variance of the n-weighted pooled-bin mean is s²/(Rn) — an
  underestimate by √R. Verified numerically: two 2-bin runs [0,2],[0,2] →
  merged stderr 0.500 vs pooled-bin reference 0.577. Pre-existing
  (character-for-character the pre-C0 scalar merge), so the C0 compatibility
  claim was honored — but C0 propagates it to array elements, i.e. to the
  merged curve outputs downstream stages produce. **Fixed in the follow-up
  commit** (pooled-bin variance `Σ nᵢ sᵢ² / N²`), with a regression test.
- **Shape-mismatch fix in release**: old code probe — `add_array(&[1,2])`
  then `add_array(&[3])` silently accepted, bins became `[3, 0]`; new code —
  `CarloError::ShapeMismatch { expected_len: 2, got_len: 1 }`, accumulator
  untouched. All in-workspace callers handle the Result with named expects.
- **Multispin migration**: `get_array("Energy_replica")` + element-mean
  average equals the old collapsed replicas-and-time ⟨E⟩ — same physical
  quantity.
- **e2e + HDF5**: Run-level curve stays element-preserving through finalize
  and JSON (means 4.5·scale exact); checkpoint round-trip restores arrays
  including the partial bin (hand-checked 7/3 and 70/3).

## Issues found

1. **Major (pre-existing, propagated; FIXED in follow-up)** — merged stderr
   underestimates by √R (derivation above).
2. **Minor (FIXED in follow-up)** — scalar `finalize()` computes tau over
   completed bins only while `finalize_array()` included the normalized
   partial bin: reviewer probe 4 complete + partial → array tau 1.104 vs
   scalar 1.0; parity claim was only true where both floor at 1.0. Now
   aligned (completed bins only) with an outlier-partial-bin parity test.
3. **Minor** — shape-mismatch merge skip is observable only via
   `tracing::warn` and absence from output (documented; accepted).
4. **Minor (doc, pre-existing)** — autocorr floor conventions differ across
   estimators (merge.rs τ=0-for-uncorrelated vs accumulate-side floor 1.0).
5. **Note (intentional, documented)** — `finalize_split` removes collapsed
   array scalars from `observables`; external consumers of collapsed values
   see the change.
6. **Note (pre-existing)** — mixing `add()` then `add_array()` on one name
   silently re-initializes the accumulator.
7. **Note** — partial bin counts as a full bin in n_bins weighting
   (pre-existing semantics); VALIDATION said "18 new tests" where the diff
   adds 20; benchmark numbers not independently re-run.

## Recommendation

**Approve C0 (PASS WITH NOTES).** The stage does exactly what the roadmap
asks; compat promises real at byte level; the shape-mismatch fix genuine in
release. File the merged-stderr underestimation before merged curves are
used for scientific error bars — **applied as the follow-up fix**.
