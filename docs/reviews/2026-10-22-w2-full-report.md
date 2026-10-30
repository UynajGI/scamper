# W2 (Winding/Wrapping Analyzer) — Full Reviewer Report (verbatim)

> Archived from the independent reviewer agent, 2026-10-22. Reviewer scripts
> lived outside the repository. See `2026-10-22-stage-reviews.md` for the
> summary and the follow-ups applied.

**VERDICT: PASS WITH NOTES**

## Roadmap §16 W2 gates

| Gate | Result | Evidence |
|---|---|---|
| Non-wrapping loop | PASS | `interior_plaquette_is_contractible` (rank 0, lift-BFS parity) |
| Single- and multi-axis wrapping | PASS | `single_axis_ring` (generator exactly `[1,0]`), `double_axis_wrap` (rank 2, `[1,0]`,`[0,1]`), `diagonal_wrap` (rank 1, generator `[1,1]`), cubic rank 3 |
| Independent finite-cover lift-BFS reference | PASS | `lift_reference` (5^D cover, f64 Gaussian rank); used in every hand case and **per-configuration** across all 2^8 and 2^18 masks |
| Origin-translation / representative / edge-reversal invariance | PASS | 8 random samples: gauge shift `d+φ(t)−φ(s)`, full table negation, cyclic relabeling; lattices compared via test-side HNF (`assert_same_lattice`) |
| Open/non-embedded wrapping → explicit error | PASS | `WindingWithoutEmbedding` (adapter), param-path rejection, `EmbeddingError::OpenBoundaries` (builders), honeycomb rejection (carlo.rs:1072 unit test). No silent-false path found |
| Union-find with integer potentials, reports generators/rank not Boolean | PASS | `analyze_winding` verified by hand-trace and independent enumeration |

Gates run by reviewer: `cargo test ... percolation_winding` → **16/16 pass**;
`SCUTTLE_ZSCORE_SEEDS=64` z-score → **pass**; `cargo fmt --all --check` →
clean; `cargo clippy -p cmc-rs --all-targets -- -D warnings` → clean.

## Independent re-derivation (reviewer's own constructions, no repo code)

Two from-scratch methods — (1) DSU with unbounded integer cell offsets,
(2) lifted BFS on a 5x5 cover — with the reviewer's own bond emission and
wrap sign convention:

| Torus | Method 1 (DSU) | Method 2 (lift BFS) | Claimed | Match |
|---|---|---|---|---|
| 2x2 P(wrapping) | 187/256 | 187/256 | 187/256 | YES |
| 2x2 mean largest rank | 252/256 | 252/256 | 252/256 | YES |
| 3x3 P(wrapping) | 186684/262144 | 186684/262144 | 186684/262144 | YES |
| 3x3 mean largest rank | 259658/262144 | 259658/262144 | 259658/262144 | YES |

(The reviewer's first DSU draft returned 179/256 due to a bug in the
reviewer's own script — stale implicit entries after re-rooting — found and
fixed; both methods then agreed with the repository.)

## Union-find / potential math audit (reviewer findings)

- Path compression (`find_root_with_potential`, winding.rs:604): potentials
  recomputed as suffix sums from the root — correct; hand-traced through the
  diagonal (1,1) 4-cycle.
- Union potential (winding.rs:425-471): both attach directions verified —
  new-root potential `±(disp + pot_R − pot_L)`.
- Cycle generators: same-root mismatch `vector − pot_L` is exactly the
  fundamental-cycle sum; fundamental cycles generate the winding subgroup.
- Echelon basis (`basis_insert`, winding.rs:660): SL(2,Z) unimodular folding
  `[[x,y],[−incoming/g, leading/g]]` preserves the Z-lattice; hand-verified
  (2,0)+(3,0) → [1,0]. Rank ≤ D by distinct pivots; terminates.
- Z-rank = R-rank for subgroups of Z^D, so the hypothesized Z-vs-R rank
  failure mode cannot occur.
- i32 safety: every potential-path operation `checked_*`; basis combinations
  in checked i64 then `i32::try_from`; the adapter static gate
  `(2V+2)·max_disp` provably covers all potential sums (bound re-derived).
- Tie policy identical to F3 (lowest canonical identity); parity-tested over
  32 random samples.

## Issues found

1. **Minor** — enumeration anchors printed, not asserted (fixed in follow-up:
   now asserted exactly).
2. **Minor** — benchmark ratio ranges overstated in prose (fixed in
   follow-up).
3. **Note** — adapter `measure()` panics loudly on adversarial tables
   passing the static gate; documented in `with_embedding`, disclosed,
   acceptable.
4. **Note** — `PeriodicEmbedding` effectively sealed
   (`DirectedDisplacement::new` is `pub(crate)`); adversarial surface
   internal-only.
5. **Note** — analyzer surfaces short embedding tables as typed
   `MissingDisplacement`; adapter adds the explicit `EmbeddingMismatch`
   check.

## Recommendation

**PASS WITH NOTES.** Physics and arithmetic correct; headline enumeration
evidence confirmed by two independently written brute-force implementations.
After the two minor follow-ups (applied), W2 can be promoted out of
Experimental.
