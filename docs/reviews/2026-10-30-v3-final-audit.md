# Percolation platform v3 — final audit

Date: 2026-10-30. Scope: the working-tree implementation of the
[v3 roadmap](../plans/2026-10-25-percolation-roadmap-v3.md) on branch
`feat/cmc-percolation`, audited against the roadmap's uniform definition of
done (§6). This is an **internal** audit: independent per-stage reviews for
the stages added in this pass are still pending.

## 1. Verdict by stage

| Stage | Verdict | Evidence / residual gap |
|---|---|---|
| A1 observables | Done (pre-existing) | A1a/A1b gates recorded in `PERCOLATION.md`; unchanged in this pass. |
| A2 laws / views | Done (pre-existing) | A2a/A2b gates recorded; unchanged. |
| A3 adapter/estimator/C0 | **In progress** | A3a `Moments`, A3b generic static adapter, A3c `FixedSeries` + shared record helpers. Bond-dynamics, continuum, and scan adapters still own their structs; the ≤3-adapter target is not met. |
| A4 policy traits | **In progress** | A4a `AttackOrder`, A4b `ScanObservableSpec`, A4c `AchlioptasScore`/`AchlioptasTie`, `InvasionStopRule`/`InvasionTrappingRule`, `BootstrapThresholdRule`/`BootstrapSeedRule`, `CoreThresholdRule`. Structural modes (site/bond, sampling, clock) and checkpoint grids remain closed by design. Directed incidence is first-class through D1. |
| A5 FSS toolkit | **In progress** | A5a–d tools plus A5e generic-toolkit reproduction of the v2 white-noise C2 crossing calibration (0.5995 vs v2 0.5987, exact 0.592746). The correlated C2 peak bracket and C3 `nu` bracket remain on the v2 pc-shift estimator; certification is R1. |
| A6 spec/provenance | **In progress** | A6a `PercolationSpec`, Params round-trip, bitwise scheduler replay, and versioned Results provenance for the independent-occupation family. Non-static families and custom-observable specs are not covered. |
| X1 network I/O | Done (pre-existing) | Independent review still pending per the roadmap row. |
| X2 multiplex | **In progress** | X2a analysis, X2b explicit initial activity + `MultiplexCascadeMC` + Buldyrev bracket (10-seed c=2.359; review 20-seed N=512/1024/2048: 2.361/2.426/2.414, approaching 2.4554 from below). A layer-aware `SampleView` and >2 layers remain. |
| X3 hypergraph | Done (pre-existing) | Independent review still pending. |
| X4 temporal/directed | **In progress** | X4a feature-gated contact process on the Gillespie substrate, deterministic ensembles, finite-time activity bracket (λ=0.535 at T=15), independent first-reaction reference, and a T=100 negative control showing the crossing is not a threshold. Temporal-network input and a scheduler adapter (X4b) remain. |
| R1 C3 certification | Done — Experimental | Negative result recorded: iid control calibrated; correlated branch flip not certified at L≤2048 (1.5-vs-1.25 gap 1.1σ). |
| R2 Q3 phase diagram | Done — Experimental | Two-mode estimator + q-ladder delivered; L=12 separation orders q=1/2/4/8 as 2.63/2.74/2.87/3.63 but is not a certified latent heat. |
| R3 F3 large-scale | Done — Experimental | F3 remains Validated at its recorded scale; aging exploratory; no large-scale ADH bracket certified. |

## 2. What changed in this pass

- **A3c** (`percolation/series.rs`, `percolation/recording.rs`): `FixedSeries`
  owns the fixed-shape C0 buffer with explicit repeat-last/fixed-tail padding
  and in-place rewriting; `record_component_analysis` /
  `record_bond_analysis` replace seven hand-written analyze/label/record
  tails. Invasion, bootstrap, k-core, and long-range C0 contents are
  bitwise unchanged.
- **A4c** (`process/{achlioptas,bootstrap,invasion,kcore}.rs`): closed
  policy enums became default type parameters behind static-dispatch traits;
  every built-in reference order is preserved and externally reimplemented
  policies reproduce it bitwise.
- **A5e** (`tests/physics/fss_brackets.rs`): an independent iid
  `SiteBernoulli` ladder feeding `BinderEstimator`, `Moments`, `FssCurve`,
  and `crossing` reproduces the v2 white-noise calibration; the L=64/128
  gate is `#[ignore]` for nightly release.
- **A6a** (`percolation/spec.rs`, `tests/physics/percolation_spec.rs`):
  explicit `PercolationSpec`, lossless Params round-trip for 15 preset
  shapes, preset parity, bitwise scheduler replay, and the
  `cmc.percolation.provenance` Results attribute.
- **X2b** (`network/multiplex.rs`, `tests/physics/percolation_multiplex_bracket.rs`):
  `mutual_component_from`, `MultiplexCascadeMC`, and the configuration-model
  Poisson ER bracket.
- **X4a** (`percolation/temporal/`, feature `temporal`,
  `tests/physics/percolation_temporal.rs`): contact process with direct
  reference-order events and survival/density ensembles. The physics review
  found the original 0.535 threshold wording overclaimed a T=15 transient;
  the test/doc now records only a finite-time activity bracket, adds an
  independent first-reaction reference, and includes a T=100 negative
  control.
- **R2** (`potts_phase.rs`, `tests/physics/percolation_potts_phase.rs`):
  raw density samples, deterministic two-mode estimator, exploratory
  q-ladder.
- Documentation reconciled: roadmap progress notes/status rows,
  `PERCOLATION.md` rows.

## 3. Gate results (this working tree)

- `cargo fmt --check`: clean.
- `cargo clippy -p cmc-rs --all-features --all-targets`: clean.
- `cargo test -p cmc-rs --all-features --lib`: 303 passed.
- `cargo test -p cmc-rs --all-features --doc`: 8 passed, 1 ignored.
- `cargo test -p cmc-rs --all-features --test suite --no-fail-fast`:
  580 passed, 0 failed, 21 ignored.
- `cargo check -p cmc-rs` (default features, no `temporal`/`io`): clean.

Golden conservation: no built-in RNG draw order, Results column name,
observable indicator, or C0 schema changed. The new traits only add
extension points; the enum presets are the defaults. External policy tests
assert bitwise equivalence with the presets and deliberately different
results for different rules. The A6 spec path replays the scheduler
preset bitwise at the same seed (timing-only `_ll_` columns excluded).

Allocation: `FixedSeries` and the recording helpers allocate zero in steady
state; the old process buffers are replaced in place. `MultiplexCascadeMC`
allocates its cascade working vectors per sweep and documents it; the
contact process recomputes the full event catalog per event (the direct
reference order) and documents it.

## 4. Residual risks and open items

1. **Preset collapse (A3)** is incomplete: bond-dynamics, continuum, and
   scan adapters still duplicate construction/validation plumbing. The
   roadmap's ≤3-adapter target is unmet.
2. **A4 structural modes and grids** remain closed enums (site/bond,
   sampling, synchronous/queue clock, checkpoint grids). They select
   different algorithms rather than parameterizing one; opening them needs a
   separate design.
3. **A6 schema** covers only the independent-occupation family and built-in
   observables. Heterogeneous laws are programmatic-only for Params
   round-trip; custom observables and non-static processes are out of scope.
4. **X2** has no layer-aware `SampleView` and supports exactly two layers.
5. **X4** has no temporal-network input and no scheduler adapter; the
   contact process is static-topology only.
6. **R1/R2/R3** remain Experimental with recorded negative results. The C3
   flip and the Q3 latent-heat ordering are not certified; neither is a
   large-scale ADH shape bracket.
7. **Independent review** is missing for the stages added in this pass
   (A3c, A4c, A5e, A6a, X2b, X4a). The roadmap forbids merging a stage
   before its review report exists; this audit does not substitute for one.
8. The nightly `#[ignore]` statistical suite (including the new FSS
   calibration) has not been run as a whole in this pass; the ignored FSS
   calibration was verified individually in release.

## 4b. Physics review

[2026-10-30-v3-physics-review.md](2026-10-30-v3-physics-review.md) records
the independent checks and the X4a correction. Verdict: A3c/A4c/A6a are
golden/reference-order refactors; A5e is a finite-size calibration; X2b is
a finite-size bracket with a realized-mean-degree gate; X4a is a finite-time
activity bracket, not a threshold; R2 is exploratory only.

## 5. Merge readiness

Conditional. The code is green under the full quick suite and clippy, the
physics review corrections are applied, and the additive v3 contracts hold,
but the roadmap's M2/M3/M4/M5 gates are not fully closed (A3 preset
collapse, A4 grids, A6 breadth, X4b, independent reviews). No automatic
merge is performed. Recommended merge unit: the A3c + A4c + A5e + A6a + X2b + X4a
stack in per-stage commits, each with its review report, and a follow-up
track for the open items above.
