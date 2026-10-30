# D1 (Directed Topology) — Full Reviewer Report (verbatim)

> Archived from the independent reviewer agent, 2026-10-22. Reviewed commit
> `c2ea61f`; reviewer artifacts under `/tmp/d1_review/`. See
> `2026-10-22-stage-reviews.md` for follow-ups applied.

**VERDICT: PASS WITH NOTES**

## Per roadmap §11 D1 bullet

| Bullet | Result |
|---|---|
| owned directed CSR views | PASS — builder hand-traced and empirically confirmed |
| borrowed directed CSR views | PASS — zero-copy, one-time validation, parity |
| separate outgoing/incoming capabilities | PASS — sibling traits; compile-checked absence of `in_arcs` on the out-only view |
| arc multiplicity validated | PASS — 0 and >=2 caught on both sides; saturating_add cannot mask |
| endpoints validated | PASS — table-level and per-incidence sense checks |
| optional transpose arrays | PASS — out-only needs none; full view validates against the offsets-derived count |

## Independent checks

- **Builder trace**: V=4, arcs (0->1, 0->1, 2->2), vertex 3 isolated —
  predicted and confirmed exact out/in offsets, neighbors, arc IDs (input
  order preserved; self-loop once per side; isolated vertex empty both
  sides). Endpoint validation precedes allocation (no panics, typed only).
- **Validator attacks (predictions written first, then confirmed in debug
  and release via a /tmp binary)**: last-offset excess, boundary endpoint,
  neighbor-is-source mismatch, transpose-missing arc, transpose-at-source,
  empty offsets, huge arc counts — every prediction matched the exact typed
  variant; self-loop double-listing caught at count==2.
- **Bench arithmetic**: 21.8 us / 48832 = 0.446 ns; the +17% arc-table
  inversion is the direct baseline accelerating (Arc = 16 bytes vs Bond =
  32 bytes; direct fold ~0.38 ns/arc, sub-scalar = SIMD), not the generic
  path degrading — every compact-layout generic scan lands at 0.42-0.50
  ns/element on both F1 and D1.
- **Trait design**: dense #[repr(transparent)] ArcId; no dyn, no
  directed-flag anywhere; arc_ids() mirrors vertex_ids() including its
  caveat; out-only callers never forced to fake a transpose.

## Issues found

1. **MINOR (FIXED)** — `validate_transpose_structure` computed
   `vertex_count + 1` unchecked: a corrupt magic-value public field
   (`usize::MAX`) panicked in debug on the error path. Now
   `saturating_add` in both validators.
2. **MINOR (FIXED)** — stale F1-section sentence still described the
   directed signatures as a crate-private sketch.
3. **MINOR (FIXED)** — cited F1 per-element 0.77 ns matched no recorded
   figure (recorded 0.72 ns; same-day rerun 0.89 ns); corrected to the
   recorded value. Claim direction robust under every measurement.
4. **NOTE (DOCUMENTED)** — `DirectedCsr::default()` fails validate (typed);
   empty graphs via `from_arcs(0, vec![])`; documented.
5. **NOTE (DOCUMENTED)** — `validate` allocates temporary scratch (endpoints
   copy per side + one byte per arc), F1-parity cold-path behavior;
   documented.
6. **NOTE (FIXED)** — "12 typed rejections" clarified as 14 assertions
   across 12 distinct variants.
7. **NOTE** — generic-vs-generic F1/D1 comparison slightly understates F1
   (F1's loop pays an Option+expect per element); the like-for-like
   borrowed-vs-borrowed band is 0.42-0.50 ns/element on both.

## Recommendation

**Merge.** Follow-ups applied were small and non-blocking; no gate or
correctness change was required.
