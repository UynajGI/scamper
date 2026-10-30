# D2 (Static Directed Analysis) — Full Reviewer Report (verbatim)

> Archived from the independent reviewer agent, 2026-10-22. Reviewed commit
> `e6120f8`; reviewer artifacts under `/tmp/d2review/`. See
> `2026-10-22-stage-reviews.md` for follow-ups applied.

**VERDICT: PASS WITH NOTES**

No correctness defect found despite adversarial attack; all deviations are
documentation/hygiene level.

## Per roadmap D2 gate

| Gate | Result |
|---|---|
| Three distinct result types | PASS |
| Tarjan vs Kosaraju | PASS — plus the reviewer's own third implementation (per-vertex BFS mutual reachability, own LCG) agreed on 300 fresh multigraphs |
| Weak vs undirected flood-fill | PASS — plus 300 additional multigraphs incl. parallel arcs/self-loops, all summary fields |
| Chains/trees closed forms | PASS — exhaustive corpus subsumes all oriented 3/4-vertex trees |
| SCC vs transpose | PASS — 40 random rounds partition-invariant |
| Performance | PASS — bench reproduces docs within noise; 0.000 allocations/round |

## Independent checks

- **Tarjan hand-trace** (0->1, 1->2, 2->0, 1->3): frame-by-frame
  index/lowlink evolution; parent propagation confirmed to use the CHILD's
  lowlink (not index); pop-loop arithmetic underflow-impossible (break at
  member == vertex fires at position >= 0); labels assigned BEFORE stack
  truncation (the historical bug verified fixed); no unreachable expects.
- **Union-find canonical identity**: traced merge chains; canonical
  maintained at roots, stale non-root values never read after find_root.
- **Stamp semantics**: proven sound across grow cycles and the u32
  rollover (cleared every grown slot); EMPIRICALLY confirmed by running
  >2^32 queries through the rollover on one workspace; interleaved
  big/small graph queries on one reused workspace correct (a stale-state
  path the repo tests did not cover).
- **Additional attacks**: V=0/1 through all three APIs; 1,000,000-vertex
  path (iterative DFS, no overflow, all-singleton labels); 300k-vertex
  ring (single SCC); equal-size tie to lowest identity; complete digraph
  with self-loops (acyclic=false).
- **Exhaustive gate is STRONGER than documented**: the corpus enumerates
  all ordered pairs including (v,v), so 2^9 = 512 and 2^16 = 65,536
  digraphs (docs said 64/4096 — understated; fixed in follow-up). Ground
  truth is independent Floyd-Warshall closure; Kosaraju shares no code.
- **O(V+A) / zero-alloc honesty**: adjacency cleared+refilled with
  retained capacity = rebuild without allocation; probe confirms.

## Issues found

1. **MINOR (FIXED)** — exhaustive counts understated (64/4096 vs actual
   512/65,536 with self-loops) in roadmap, PERCOLATION.md, and the commit
   message; "simple digraph" wording wrong for a self-loop corpus.
2. **MINOR (FIXED)** — unused import in the allocation-probe target made
   `clippy -D warnings` fail there (the target has required-features, so
   `--all-targets` skipped it); the pre-existing `0xa11_a` byte grouping
   in the same file fixed as well.
3. **MINOR (FIXED)** — bench header comment said A=48832; the grid builds
   A=40,576.
4. **NOTE (FIXED)** — Tarjan's u32 index bound (UNSET = u32::MAX) is
   implicit; now guarded by a debug_assert with an explanatory comment.
5. **NOTE (FIXED)** — the reachability queue grew via plain
   `Vec::with_capacity` (abort on OOM) while the visited array used typed
   capacity errors; now consistent (try_reserve_exact).
6. **NOTE** — the exhaustive corpus tests each arc subset in one canonical
   CSR order; arc-order independence rests on the algorithm plus the
   random corpus. Acceptable. The random corpus really does contain
   self-loops (188) and parallel arcs (553 instances) for the fixed seed,
   though nothing asserts it.

## Recommendation

**Accept D2 after the three minor follow-ups** (applied). The algorithmic
core survived every attack; the only failures were prose and lint.
