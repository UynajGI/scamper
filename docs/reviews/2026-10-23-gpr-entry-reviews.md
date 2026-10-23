# Independent Reviews — 2026-10-23 — G1, P1, R1 (commit `bd04574`)

All three M5 entry stages were reviewed in parallel by independent
reviewer agents (code authored by parallel development agents). Verdicts:
**G1 PASS, P1 PASS, R1 PASS** — no blockers, no majors anywhere; all
findings were documentation polish or API ergonomics. Full verbatim
reports follow.

| Stage | Reviewer's independent confirmation | Follow-ups applied |
|---|---|---|
| G1 incremental union | Hand trace (6-vertex multigraph, 12 interleaved ops incl. edge-before-endpoints, self-loop-inside-merge, second-parallel-edge no-union) exact on all six fields; 3,384-op drift audit vs own LCG + flood fill — zero second-moment drift; u128 overflow proven impossible (M2 ≤ V² < 2^128); trajectory budget proven (each op accepted at most once between resets); **complete Achlioptas product-rule and invasion drivers written against the public substrate, compile and run unchanged** — the without-forcing-a-law contract verified by construction | Stale process/mod.rs header updated (E track → E+G tracks); same-(V,E)-graph positional-reinterpretation sentence added |
| P1 candidate pairs | Reference built on **explicit image enumeration** (independent of the production round() metric), closing the gap the repo test itself acknowledged; exact-boundary, image-only, anisotropic (per-axis stencil proof), D=1/3, c = L/2 tie regime, boundary particles, stale-list rejection — all confirmed; all 12 bench pair counts and throughput reproduce exactly; three of the reviewer's own hand-expectations were wrong and production right each time | "axis 0 slowest" gloss corrected to fastest (strides[0]=1); the duplicated two-line norm fold replaced by the particle module's own `distance_squared`; c = L/2 minimum-image tie convention documented |
| R1 external bonds | Hand trace on a 5-vertex multigraph (self-loops + parallels) exact incl. canonical labels; decoupling audit clean (grep hits doc-comments only; the winding bridge genuinely monomorphizes the external type); BondSample sorted invariant proven at every constructor; genericity verified over `BorrowedUndirectedCsr`; the bench's bit-identity assertion proven non-vacuous (a law-side RNG change would abort it) | Doc overclaim scoped to dimension mismatch (same-dimension cross-graph use is positional, the library-wide ID contract); `BondSample` now implements `EdgeActivity` directly (O(log k) membership) so an R3 SW adapter hands its draw to the analyzers without a dense-mask materialization; bench wording softened |

All three stages promoted from "Experimental — awaiting independent
review" to reviewed status; production availability remains conditioned on
the `dev` merge.

## G1 verbatim reviewer report

> **VERDICT: PASS.** The substrate is exactly what the spec asks for and
> nothing more. All aggregates are exact under adversarial interleavings
> (3,384 ops, zero drift); the budget proof is airtight (each vertex/edge
> accepted at most once between resets); the u128 second-moment width is
> provably sufficient (M2 ≤ V² < 2^128). The `note_largest` tie policy is
> provably correct online (an equal-size merge never destroys a max-size
> component). G2/G3 drivers (a complete Achlioptas product-rule process
> and an invasion frontier process) compile and run against the public
> substrate unchanged — including graceful EdgeAlreadyActive handling.
> Issues: stale process/mod.rs header (fixed); same-size-graph positional
> acceptance follows the documented crate-wide no-provenance ID convention
> (now restated in the module doc); the bounds-check errors are
> defense-in-depth, not publicly exercisable.

## P1 verbatim reviewer report

> **VERDICT: PASS.** All spec items and documented numbers reproduce
> exactly (N=10000 test record bit-identical; all 12 bench pair counts
> exact; 1.05-1.08e7 pairs/s at 16384/z=12). My reference used explicit
> image enumeration — genuinely independent of the production round()
> metric — and confirmed exact boundaries (inclusive d² ≤ c², no
> epsilon), image-only pairs, anisotropic cells (per-axis stencil
> argument verified including the exact L−w boundary and minimum
> count=2), D=1/3, the c = L/2 tie regime (pair set well-defined), cell
> boundary particles, stale CellList rejection, and full canonical-label
> F3 parity. Adversarial anecdote in the implementation's favor: three of
> my own hand-expectations were wrong (missed image pairs, wrap-pair
> geometry) and production was right each time. Issues: an inverted
> "axis 0 slowest" gloss (contract itself correct — fixed); a duplicated
> two-line norm fold replaced by the module's own distance_squared
> (applied); the c = L/2 tie convention now documented; prepared-list
> savings slightly overstated (validate is O(N) per call) — informational.

## R1 verbatim reviewer report

> **VERDICT: PASS.** Spec fully met: zero sampler coupling (grep hits are
> doc comments only), both shared analyzers (F3 and W2) consume the
> external type directly, semantics match F3 to the last detail under an
> independent multigraph hand-trace (self-loop = counted no-op union;
> parallel edges counted separately; canonical labels exact), and every
> rejection path is typed. BondSample's sorted invariant holds at every
> constructor so the O(log k) membership claim is valid; the winding
> bridge monomorphizes the external type with no StaticConfiguration
> copy; genericity verified over BorrowedUndirectedCsr; the bench's
> bit-identity assertion is non-vacuous. Issues: a doc overclaim (same-
> dimension cross-graph analysis is positional reinterpretation, not a
> typed rejection — scoped); BondSample lacked an EdgeActivity impl
> (added — an R3 SW adapter now feeds the analyzers without a dense
> mask); write_mask reuses a MaskLength message that reads oddly for an
> output buffer (cosmetic, left); the draw-contract wording softened to
> reference the runtime assertion rather than a nonexistent law doc.
