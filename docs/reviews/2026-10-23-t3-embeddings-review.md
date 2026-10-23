# T3 independent review — honeycomb/kagome periodic embeddings (commit 9211ea5)

Reviewer: independent agent (not the stage author). Reviewed at branch
`feat/cmc-percolation`; code at HEAD `cf65fd4` is byte-identical to `9211ea5`
for both in-scope code files (only PERCOLATION.md / roadmap docs landed after,
in `501fa10`). No code was modified during this review; all verification ran
out-of-repo (scratch crate at `/tmp/t3review` with a path dependency on the
real `cmc-rs`) or via read-only git.

Scope reviewed: `CMC.rs/src/lattice/graph.rs` (honeycomb_bonds / kagome_bonds /
`build_honeycomb_with_embedding` / `build_kagome_with_embedding` + unit tests),
`CMC.rs/tests/physics/percolation_winding.rs` (two T3 tests),
`CMC.rs/PERCOLATION.md` (T3 row + section, W1 note edit),
`docs/plans/2026-10-23-percolation-roadmap-v2.md` §2 T3.

**Overall verdict: PASS** (two Minor documentation follow-ups, no code defects
found; every engineering claim reproduced independently, most extended beyond
the shipped gates).

## Verdict per claim

| # | Claim under review | Verdict |
|---|---|---|
| 1 | Shared bond-emission restructure; with-embedding builders return exactly the plain lattices (pinned 10+10 sizes) | **Reproduced** (+ extended: 12 honeycomb, 13 kagome sizes) |
| 2 | Hand-derived 2x2 displacement tables (6 honeycomb + 24 kagome bonds, both directions) | **Reproduced** (derived from rustdoc alone before reading the test tables; 30/30 rows match test AND production) |
| 3 | Per-bond displacement geometry follows the documented cell arithmetic (independent re-derivation) | **Reproduced** (+ extended: from-scratch generator, exact per-bond + emission-order match over 25 sizes incl. degenerate/asymmetric) |
| 4 | Plain builders bit-identical to the pre-T3 builders (`2a9165f`) | **Reproduced** (verbatim pre-T3 builders run side-by-side; all CSR arrays + bond kinds + weight bits identical over 25 sizes) |
| 5 | Honeycomb `Lx=2`: `ly` doubled horizontal pairs, each carrying `(0,0)` then `(-1,0)`; wrapping count `ly + lx/2` | **Reproduced** (all closed forms exact at every swept size) |
| 6 | Kagome has NO parallel bonds at any `Lx,Ly >= 2`; W1 "collapsed covering bonds" is a reconstructability collapse, not multiplicity | **Reproduced** (0 repeated endpoint pairs at 13 sizes; algebraic proof; reinterpretation sound — see §4) |
| 7 | Fully occupied tori wind with rank 2 and det ±1 generators; cutting wrapping bonds never winds | **Reproduced** (+ extended: my own lift-BFS on 5 fresh sizes per family; adversarial rank-1/rank-0 configs) |
| 8 | PERCOLATION.md T3 row/section, W1 note edit, gates list honesty; vertex numbering reproducible from rustdoc | **Mostly reproduced** — two Minor wording findings (§Findings) |

## 1. Hand derivation from the rustdoc (before reading the tests)

Method: I derived both tables using ONLY the `build_honeycomb` / `build_kagome`
rustdoc (site numbering + neighbour families), the displacement convention in
the `build_honeycomb_with_embedding` rustdoc ("a `+1` step along an axis from
coordinate `c` crosses one fundamental-cell boundary iff `c + 1` is the axis
length, a `-1` step crosses one iff `c == 0`"), and the once-per-bond
smaller-endpoint emission, then compared against (a) the tests' asserted
tables and (b) production output. All 30 rows matched both.

### Honeycomb 2x2 (lx = 2, ly = 2; sites `(x,y) -> y*2+x`)

- site 0 = (0,0), even x: right → (1,0)=1, disp (0,0) [x+1=1 ≠ 2]; left → (1,0)=1, disp (−1,0) [x=0]; matched → (1,1)=3, disp (0,0). All targets > 0 → 3 rows.
- site 1 = (1,0), odd x: right → 0 (+1,0) skipped; left → 0 (0,0) skipped; matched → (0,1)=2, disp (0,−1) [y=0]. → 1 row.
- site 2 = (0,1), even x: right → 3, (0,0); left → 3, (−1,0); matched → (1,0)=1 skipped. → 2 rows.
- site 3 = (1,1), odd x: all candidates target ≤ 3 → none.

My table (edge, source→target, displacement):

```
e0: 0→1 (0,0)   e1: 0→1 (-1,0)   e2: 0→3 (0,0)
e3: 1→2 (0,-1)  e4: 2→3 (0,0)    e5: 2→3 (-1,0)
```

Matches `honeycomb_2x2_displacements_match_hand_derivation` 6/6 (graph.rs:1354-1361),
including the doubled pairs {0,1} and {2,3} each carrying (0,0) before (−1,0),
and matches production. 6 edges on 4 distinct pairs; 3 wrapping bonds
(= ly + lx/2 = 2 + 1).

### Kagome 2x2 (lx = 2, ly = 2; `S + 3(x + 2y)`: cells (0,0)=0-2, (1,0)=3-5, (0,1)=6-8, (1,1)=9-11)

Applying the documented neighbour lists (A: B,C,B(x−1,y),C(x,y−1);
B: A,C,A(x+1,y),C(x−1,y); C: A,B,A(x,y+1),B(x+1,y)) with the ±1-boundary
displacement rule and smaller-endpoint emission:

```
e0: 0→1 (0,0)    A(0,0)-B(0,0)      e12: 4→5 (0,0)   B(1,0)-C(1,0)
e1: 0→2 (0,0)    A(0,0)-C(0,0)      e13: 5→9 (0,0)   C(1,0)-A(1,1)
e2: 0→4 (-1,0)   A(0,0)-B(1,0) wrap e14: 6→7 (0,0)   A(0,1)-B(0,1)
e3: 0→8 (0,-1)   A(0,0)-C(0,1) wrap e15: 6→8 (0,0)   A(0,1)-C(0,1)
e4: 1→2 (0,0)    B(0,0)-C(0,0)      e16: 6→10 (-1,0) A(0,1)-B(1,1) wrap
e5: 1→3 (0,0)    B(0,0)-A(1,0)      e17: 7→8 (0,0)   B(0,1)-C(0,1)
e6: 1→5 (-1,0)   B(0,0)-C(1,0) wrap e18: 7→9 (0,0)   B(0,1)-A(1,1)
e7: 2→4 (0,0)    C(0,0)-B(1,0)      e19: 7→11 (-1,0) B(0,1)-C(1,1) wrap
e8: 2→6 (0,0)    C(0,0)-A(0,1)      e20: 8→10 (0,0)  C(0,1)-B(1,1)
e9: 3→4 (0,0)    A(1,0)-B(1,0)      e21: 9→10 (0,0)  A(1,1)-B(1,1)
e10: 3→5 (0,0)   A(1,0)-C(1,0)      e22: 9→11 (0,0)  A(1,1)-C(1,1)
e11: 3→11 (0,-1) A(1,0)-C(1,1) wrap e23: 10→11 (0,0) B(1,1)-C(1,1)
```

Matches `kagome_2x2_displacements_match_hand_derivation` 24/24 (graph.rs:1392-1417)
and production. 24 edges on 24 distinct endpoint pairs (no parallel bonds);
6 wrapping bonds = 2·A–B left + 2·B–C left (x) + 2·A down (y) = `2*ly + lx`.

The tests' 4x2 (full 12-bond), 3x2 (all 7 wrapping), and 4x3 spot tables were
also cross-checked against production output at assertion level: 12/12, 7/7
(verified exhaustive — no wrapping bond of 3x2 is missing from the list), and
3+3/3+3 respectively.

## 2. Independent geometry checker (from-scratch generator, size sweep)

I wrote my own displacement generators in JS sharing nothing with the builder
code: lifted coordinates `(x+δx, y+δy)` with displacement computed by
**Euclidean mod + integer division** (`d = (lifted − mod(lifted, L)) / L`) —
a genuinely different arithmetic path from the builder's boundary-flag
comparisons (`x+1 == lx`, `x == 0`) — plus my own emission-order code
(ascending source, stable sort by target, family order preserved). Sizes
swept (superset of every size appearing in any T3 gate):

- honeycomb (12): 2x2, 2x3, 2x4, 2x5, 2x7, 4x2, 4x3, 6x2, 6x3, 6x5, 8x4, 10x3
- kagome (13): 2x2, 2x3, 2x4, 2x5, 3x2, 3x3, 3x5, 4x2, 4x3, 5x2, 5x4, 6x3, 7x2

(includes asymmetric shapes and both Lx=2 and Ly=2 degenerates; production
tables extracted from the real crate via a path-dependency binary.)

**Result: exact per-bond equality with production in emission order at every
one of the 25 sizes** — source, target, dx, dy for every edge; additionally
bond kinds (honey_x/honey_y/kagome) and weight bits (all exactly 1.0) match my
family-based prediction. No counterexample found. Closed forms all exact:

- honeycomb: edges = 3·lx·ly/2 everywhere; doubled pairs = ly at lx=2 (2x2:2,
  2x3:3, 2x4:4, 2x5:5, 2x7:7), 0 at lx ∈ {4,6,8,10}; every doubled pair's two
  rows are [(0,0), (−1,0)] with the non-wrapping bond first; wrapping bonds =
  ly + lx/2 (e.g. 2x5: 6, 4x3: 5, 6x5: 8, 10x3: 8).
- kagome: edges = 6·lx·ly everywhere; distinct pairs = edges (zero parallel
  bonds) at every size; wrapping bonds = 2·ly + lx (e.g. 2x2: 6, 3x2: 7,
  5x4: 13, 7x2: 11).

## 3. Plain-builder parity vs pre-T3 (commit 2a9165f)

I extracted the pre-T3 `build_honeycomb` / `build_kagome` **verbatim** from
`git show 2a9165f:CMC.rs/src/lattice/graph.rs` into the scratch crate (they
run on the current crate's public `CsrLattice::from_adjacency`, which I
verified byte-identical to the 2a9165f version), built both old and new for
all 25 sweep sizes, and compared every public field:
`offsets`, `neighbors`, `edge_ids`, `edges` (endpoints + kind + weight bits),
`n_sites`, `n_bonds`.

**Result: bit-identical at all 25 sizes, both families.** The in-repo parity
test (`honeycomb_and_kagome_embedding_builders_reproduce_the_plain_lattices`,
graph.rs:1281-1342) does compare the right things (full `edges` incl. kinds —
`Bond` derives `PartialEq` over weight too — plus all CSR arrays, 10+10 sizes)
but only pins embedding-variant vs plain-variant, both of which share the new
emission; the historical parity is argued (BTreeMap order equivalence), not
pinned in-repo. My out-of-repo check closes that gap empirically. (Downstream
consumers `percolation_scan.rs`, `percolation_components.rs`,
`usage_exact.rs` also exercise the refactored plain builders and are green —
extra regression evidence.)

## 4. Kagome no-parallel-bonds and the W1 reinterpretation

**Enumeration:** zero repeated endpoint pairs at all 13 swept sizes
(distinct pairs == edge count everywhere).

**Proof (why it holds at every size ≥ 2):** a site would need two of its
documented neighbours to coincide. For A: B(x,y)=B(x−1,y) or C(x,y)=C(x,y−1)
would require lx | 1 or ly | 1 — impossible for L ≥ 2 (A vs B, A vs C are
different sublattices); similarly for B (A(x,y)=A(x+1,y), C(x,y)=C(x−1,y)) and
C (A(x,y)=A(x,y+1), B(x,y)=B(x+1,y)). So no site ever has a repeated
neighbour, hence no parallel bonds.

**Why W1 expected collapse, and whether the reinterpretation is sound:** at
`Lx = 2` there are two A–B cross-family bonds per row:
{A(0,y), B(1,y)} — A's left step, displacement (−1,0) (wraps the x=0 boundary)
and {A(1,y), B(0,y)} — B's right step, displacement (0,0) (no wrap; 0+1 ≠ 2).
My 2x2 table shows them as distinct endpoint pairs (e2: 0→4 vs e5: 1→3), so
there is **no multiplicity collapse** — the pre-T3 lattice data kept the two
bonds perfectly distinguishable. The genuine obstruction is different: from
endpoint cells alone, {A(0,y), B(1,y)} has B in cell 1 ≡ −1 (mod 2), so the
covering-lattice displacement (+1,0) and (−1,0) are both valid images and
equidistant — the wrap cannot be attributed without knowing which family
generated the bond. (For Lx ≥ 3 the nearest image is unambiguous, e.g.
{A(0,y), B(lx−1,y)} → (−1,0).) This is exactly the "reconstructability
collapse, not bond multiplicity" reading in the T3 section and the rustdoc
(graph.rs:749-754), and the per-bond emission-time attribution is the correct
fix. **The reinterpretation is sound.** The one residue is a wording leftover
in the W1 note (Finding M1 below).

## 5. Winding — own lift-BFS reference

I implemented my own bounded-cover lift BFS (radius 4, per-seed components,
offsets at non-zero-cell returns to the seed) on **my own** generated bond
tables, with exact integer rank (2x2-minor gcd for the subgroup index; for
rank 2, index = gcd of all minors, so index 1 ⇔ subgroup is all of Z^2 ⇔ a
generator pair has |det| = 1).

**Full occupancy, fresh sizes (not in any T3 test list):**

| family | sizes | result |
|---|---|---|
| honeycomb | 6x2, 8x3, 4x4, 10x3, 2x5 | 1 component each, all sites connected, rank 2, index 1 (= Z^2, det ±1) |
| kagome | 6x2, 4x4, 7x3, 3x5, 5x2 | 1 component each, all sites connected, rank 2, index 1 |

This independently reproduces the T3 winding claim (the repo tests cover
honeycomb 2x2/2x3/4x2/4x3/6x4 and kagome 2x2/2x3/3x2/4x3/5x4 against the
analyzer; mine covers 10 disjoint fresh sizes against my own reference).

**Adversarial partial configurations (honeycomb 8x3, kagome 7x3, plus one
nasty case):**

| configuration | honeycomb 8x3 | kagome 7x3 |
|---|---|---|
| cut every nonzero-displacement bond | 1 component, rank 0, no wrap | 1 component, rank 0, no wrap |
| keep only dx = 0 bonds (kill x-wraps) | rank 1 | rank 1 |
| keep only dy = 0 bonds (kill y-wraps) | rank 1 | rank 1 |

Extra adversarial: honeycomb 2x5 keeping **only** the wrapping bonds (all
(0,0) bonds deleted, so each doubled pair keeps exactly its (−1,0) half):
4 components (largest 4/10 sites), every component rank 0 — consistent with
the geometry (the kept graph is a chain of rungs with no closed winding walk).
All outcomes match the expected rank structure; no anomaly.

## 6. Gates re-run

- `cargo test -p cmc-rs --test suite percolation_winding` — **18 passed, 0
  failed** (includes `honeycomb_and_kagome_full_tori_wind_both_axes` and
  `honeycomb_and_kagome_without_wrapping_bonds_do_not_wind`).
- `cargo test -p cmc-rs --lib` — **214 passed, 0 failed**; the six T3 unit
  tests run by name under `lattice::graph` (19 passed in that filter).
- `cargo test -p cmc-rs` (full crate) — **214 lib + 488 suite + 1 doc-test
  passed, 0 failed** (17 suite tests feature-ignored, pre-existing).
- `cargo test --workspace` — **all test targets passed, 0 failed** (Carlo.rs,
  QMC.rs, MCMC.rs, scamper, cmc-rs; re-run with full result-line capture after
  a first background run's log was truncated by `tail`; the captured result
  lines show every target `ok`).

Also verified out-of-repo: plan-level `Wrapping` on the plain honeycomb/kagome
builders is rejected loudly (`StaticPercolationError::WindingWithoutEmbedding`
for honeycomb 4x3 and kagome 3x3) — the T3 Limitations claim holds.

## 7. Docs honesty

- **T3 row (PERCOLATION.md:51)**: every enumerated gate matches the shipped
  tests exactly (10+10 parity sizes, 6+24 hand-table bonds, 4x2/3x2/4x3 spot
  checks, 8+9 re-derivation sizes, winding size lists, closed forms) — verified
  list-for-list against graph.rs:1281-1804 and percolation_winding.rs:989-1099.
- **T3 section small-cell numbers**: 2x2 has 2 doubled pairs (6 edges on 4
  distinct pairs) ✓; 2x3 has 3 ✓; Lx ≥ 4 doubles nothing ✓; kagome 2x2 has
  exactly 6 wrapping bonds with the stated A–B/B–C/A-down composition ✓ (my
  derivation: e2, e16 A–B; e6, e19 B–C; e3, e11 A-down); `2*ly + lx` at every
  size ✓.
- **W1 note edit**: the deferral-removal rewrite is accurate except one clause
  (Finding M1).
- **Vertex numbering reproducibility**: demonstrated empirically — my
  generators were written from the rustdoc text alone and reproduced
  production exactly, including emission row order. Caveat: the full
  emission-order contract (ascending (source,target); within a site by target;
  stable family order) is documented on the *private* `honeycomb_bonds` /
  `kagome_bonds` rustdoc (graph.rs:672-682, 774-778) and in PERCOLATION.md; the
  public `build_honeycomb_with_embedding` doc surfaces only the doubled-pair
  part (graph.rs:659-661). Note N2.
- **Roadmap §2 T3 letter** asks for "hand-derived plaquette/cocycle gates as
  W1's"; T3 delivered hand-derived full displacement tables + per-bond
  independent re-derivation + lift-BFS winding parity instead. The
  substitution is mathematically at least as strong: each displacement is a
  difference of the covering-lattice potential (floor division by L), so
  per-bond agreement with cell arithmetic implies zero sum around every
  contractible cycle; my 25-size per-bond sweep covers it. Deviation from the
  spec letter is not called out anywhere — Note N1.

## Findings

**Blocker**: none.

**Major**: none.

**Minor**:

- **M1 — W1 note keeps an inaccurate clause for kagome.**
  PERCOLATION.md:957-964 (the edited W1 builder-coverage bullet) still says
  the pre-T3 obstruction was that "distinct displacements mapped to
  indistinguishable physical bonds" for both families. That is true for
  honeycomb Lx=2 (two parallel bonds between one pair) but false for kagome —
  by T3's own (verified) analysis kagome has no parallel bonds at any size;
  the bonds stayed endpoint-distinguishable and what was unrecoverable is the
  covering-image attribution (±1 ambiguity at Lx=2). The T3 section states
  this correctly ("reconstructability collapse, not a multiplicity one",
  PERCOLATION.md:3313-3319); the W1 bullet should be aligned (e.g. "…for
  honeycomb, distinct displacements mapped to indistinguishable parallel
  bonds; for kagome, which covering image a bond uses became unattributable").
  Documentation-only; no behavioral impact.
- **M2 — "bit-identical … (pinned by test)" conflates two different parity
  claims.** PERCOLATION.md:3282-3286 (and the commit message) say the
  refactored plain builders "stayed bit-identical (pinned by test)". The
  in-repo test pins embedding-variant == plain-variant (both use the new
  emission); parity against the *historical* `2a9165f` builders is only
  argued via the BTreeMap-order equivalence, not pinned by any in-repo test.
  This review verified historical parity externally (verbatim pre-T3 builders,
  25 sizes, all fields identical), so the substance is true — but a reader
  could take "pinned by test" as covering the historical claim. Consider
  either a comment-level clarification or (nicer) keeping a frozen copy of the
  old construction as a test reference.

**Notes**:

- **N1** — Roadmap §2 T3 (roadmap:84-88) specifies "hand-derived
  plaquette/cocycle gates as W1's"; T3 delivered a different (mathematically
  subsuming, arguably stronger) gate set. The deviation deserves a sentence
  in the roadmap status or the T3 section.
- **N2** — The emission-order contract that fixes displacement-table row
  order is documented on private helper rustdoc; public docs expose only the
  doubled-pair ordering. A researcher using rustdoc HTML alone (private items
  not rendered) can reproduce every bond and displacement but not the general
  row order; PERCOLATION.md documents it, so this is a polish note.
- **N3** — graph.rs:1387-1391 comment inside
  `kagome_2x2_displacements_match_hand_derivation` is garbled: "the x = 0 left
  steps of A and B (six sites -> but only the pairs whose smaller endpoint
  emits them)". At 2x2 all six wrapping family steps do produce emitted bonds
  (none is dropped); "six sites" should read six family steps. Assertions are
  correct; comment-only.
- **N4** — Positive: routing the builders through
  `LatticeEmbedding::try_from_edge_displacements` on every construction means
  the W1 validation (length/self-loop/i32::MIN/capacity) is exercised by
  every T3 test for free, exactly as gate 4 claims.

## Conclusion

Every claim in the T3 commit message, PERCOLATION.md row/section, and rustdoc
contracts reproduced under independent verification, most under strictly
stronger conditions than the shipped gates (25-size per-bond exact sweeps vs
8/9; 10 fresh winding sizes with a from-scratch lift-BFS; verbatim pre-T3
builder parity). The small-cell degeneracy analysis — honeycomb doubled pairs
with (0,0)-before-(−1,0) emission order, kagome's no-parallel-bonds fact and
the reconstructability reinterpretation of W1's note — is correct. All test
gates green, full workspace green. The only follow-ups are documentation
wording (M1, M2, N1-N3); no code changes required.

**Verdict: PASS.**
