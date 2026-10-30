# N3 (Targeted Degree Attack + Robustness Curves) — Full Reviewer Report (verbatim)

> Archived from the independent reviewer agent, 2026-10-22. Reviewed commit
> `5f2dd2a`; reviewer artifacts under `/tmp/n3review/`. See
> `2026-10-22-stage-reviews.md` for follow-ups applied.

**VERDICT: PASS WITH NOTES**

## Per roadmap §9 N3 gate

| Gate | Result | Evidence |
|---|---|---|
| Each step matches naive recomputation | PASS | Exact (non-statistical) equality vs independent adjacency BFS on star5/path6/er12, both policies, every checkpoint incl. susceptibility |
| Active count strictly monotone | PASS | Survivor mask recomputed per step |
| Ties reproduce for fixed seeds | PASS | Stable order seed-independent, equals hand-derived path8 order; seeded shuffle bitwise for fixed seed, reshuffles for new seed; degree-group structure preserved |
| Coordinates + provenance | PASS | `RobustnessRemovedFraction[i] == i·5/40` exactly; policy attribute round-trips Results JSON |
| Literature anchors | PASS | All eight Callaway cutoff values and both Cohen f_c values reproduced digit-for-digit by the reviewer's independent solver on independently rebuilt graphs |

## Independent verification highlights

- **Callaway Eq. (14) from-scratch solver**: closed-form two-delta checks to
  1.1e-16; ideal Poisson(2.5) values distinct from empirical-histogram
  values (confirms the test uses each graph's own bins); end-to-end rebuild
  of both anchor graphs (V=60k ER with z=2.5 exactly, kappa0=3.494; power
  law alpha=2.5, m=2 with 0 self-loops, degree sums preserved) reproduces
  K=7..4 -> 0.8851/0.8627/0.7954/0.5690 and K=80..10 ->
  0.9898/0.9694/0.9175/0.7707 identically. Boundary argument verified: the
  first sum_{k>K} n_k entries of any degree-descending order are exactly
  the k>K vertices, so each checkpoint realizes q_k = 1{k<=K} exactly.
  The test's F1 normalization (full original z, per the paper's Eq. 2 as
  printed) is correct.
- **Path6 hand simulation**: removals [1,2,3] give S = 4/6, 3/6, 2/6 and
  chi = 1 throughout; intact/empty checkpoints undefined per the F4
  indicator convention.
- **Fisher-Yates**: forward form with absolute indices is uniform
  (verified by exhaustive choice-sequence enumeration for n = 2, 3, 4);
  RNG consumption depends only on tie-group sizes; StableId consumes no
  RNG. The historical double-offset bug is gone.
- **Cohen bracket**: the reviewer's from-scratch discrete implementation of
  Eqs. (7)/(9)/(2) + kappa = 2 reproduces f_c = 0.0211 and 0.2428 exactly;
  the binomial-thinning identity kappa = p_tilde + (1-p_tilde)kappa_0
  verified analytically.
- **Allocation probe re-run**: 0.000 allocations/bytes per trajectory on
  all five topologies.

## Issues found

1. **MINOR (doc accuracy; FIXED)** — the "inside their predicted 0.02-0.03"
   characterization misquotes Cohen et al.: their continuum equations give
   p_c ~ 0.056 at alpha = 2.5 (0.02-0.03 belongs to alpha ~ 3). Docs and
   the test comment now state our discrete-empirical estimate explicitly.
2. **MINOR (doc inconsistency; FIXED)** — roadmap row said "survives
   0.5 f_c"; the gate uses 0.3 f_c.
3. **MINOR (hygiene; FIXED)** — stray `docs/ref/txt_cohen2001.txt` removed.
4. **NOTE** — with StableId the 16-seed loop is decorative (deterministic
   trajectories); the printed z is Delta-S/0.006, honestly documented.
5. **NOTE** — self-loop repair partner exclusion mildly biases the builder
   ensemble vs a uniform loopless configuration model (documented in the
   builder; realized graphs have 0 residual self-loops and exact degree
   sums).
6. **NOTE** — `walk_trajectory` trusts its order argument (debug_assert
   only); the `measure()` path is always self-consistent.
7. **NOTE** — grid sorting/dedup/validation, ER G(n,m) uniformity, the
   mirrored-index fix alignment, and the V^2 exactness gate all verified.

## Recommendation

**Accept N3 (PASS WITH NOTES).** All core claims survive adversarial
re-derivation, including an independent end-to-end rebuild of both anchor
graphs. Follow-ups applied were documentation-honesty and hygiene fixes
only; no gate or code change was required.
