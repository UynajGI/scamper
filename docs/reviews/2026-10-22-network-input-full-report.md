# Network Edge-List Input — Full Reviewer Report (verbatim)

> Archived from the independent reviewer agent, 2026-10-22. Reviewed commit
> `539a16f`; reviewer's throwaway parser-attack crate lived under /tmp.
> See `2026-10-22-stage-reviews.md` for follow-ups applied.

**VERDICT: PASS WITH NOTES**

## Per roadmap §9 bullet

| Bullet | Result |
|---|---|
| edge-list input builds an owned CSR | PASS — parity test (V/E, degree profile, full-activity summary) |
| external CSR arrays borrow zero-copy | PASS — truly borrowed (no clone of incidence arrays); owned vs borrowed analysis equal |
| parsing reports line/endpoint/duplicate errors explicitly | PASS (with the vertex_count=0 panic caveat, fixed in follow-up) |
| no binary format before million-edge justification | PASS — benchmark genuine and reproduced |

## Independent checks

- **Parser attack** (throwaway crate against the published API): CRLF,
  leading `+`, zero padding, tabs, Unicode whitespace all parse; `#`
  mid-line and extra tokens are typed Parse errors; integer overflow is a
  typed Parse error (no panic); whitespace-only lines skipped; duplicates
  detected across normalized order (`2 1` then `1 2`); file-order emission
  verified (`(5,6),(1,2),(0,9)` in, same order out).
- **CSR semantics (empirical)**: self-loop = two incidences (multigraph
  degree 5 from two self-loops + one edge); three parallel edges stay
  distinct; `from_edges` cannot reject anything the parser accepted except
  `n_sites == 0` (the panic path, fixed).
- **Benchmark re-run**: er16k 2.97 ms / 5.51 M edges/s (doc 3.0 ms /
  ~5.4 M); er1m 796 ms / 1.26 M edges/s / 13,110,148 bytes (doc 814 ms /
  ~1.23 M / 13.1 MB) — doc numbers slightly conservative, honest.

## Issues found

1. **MAJOR (FIXED)** — `vertex_count = 0` with empty or comment-only input
   panicked inside `CsrLattice::from_edges`; now a typed
   `NetworkInputError::InvalidVertexCount`.
2. **MAJOR (FIXED)** — the PERCOLATION.md matrix row pre-declared
   "Independent review passed" before any review existed; the claim is now
   substantiated by this archived report and the follow-up commit.
3. **MINOR (FIXED)** — self-loops were exempt from duplicate detection;
   Reject now covers repeated self-loops too (a repeated edge is a repeated
   edge), with a regression test.
4. **MINOR (FIXED)** — stale "construction-level reason" doc comment
   (variant now exists); test header claimed N3 coverage where only N2 is
   exercised.
5. **NOTE** — a UTF-8 BOM before a leading `#` is a typed Parse error
   (acceptable; now documented).
6. **NOTE** — text at ~16 MB/s means ~80 s at 100M edges; a future binary
   format may be justified beyond the roadmap's million-edge bar (outside
   current scope, correctly not introduced).

## Recommendation

**PASS WITH NOTES.** Core claims hold under adversarial testing; both
majors were small, well-localized fixes applied in the follow-up commit.
