# Percolation platform roadmap v3 — pluggable research framework

Status: **Proposed** (2026-10-25). This roadmap turns the validated v2
stage collection into a general, customizable, research-enabling library
in the spirit of CMC/QMC/MCMC. It is architecture-first: generic
extension points are delivered and proven before new domains are added.

Predecessors: [v1 platform roadmap](2026-09-03-percolation-platform-roadmap.md),
[v2 roadmap](2026-10-23-percolation-roadmap-v2.md). Living evidence
ledger: [PERCOLATION.md](../../CMC.rs/PERCOLATION.md). v2 audit:
[2026-10-24-v2-final-audit.md](../reviews/2026-10-24-v2-final-audit.md).

## 0. Context

v2 delivered all 16 stages (14 Validated, 2 Experimental: C3's
Weinrib–Halperin flip and Q3's phase diagram). The physics is
validation-strong, but the **upper API is stage-shaped**:

- 15 concrete `MonteCarlo` adapters repeat lattice/configuration/
  workspace/plan validation and recording (~5 copies of the constructor
  sequence, 8 copies of C0 padding).
- 5 estimator types repeat count/sum/sum_sq logic.
- Core extension points are closed enums: `StaticObservable`,
  `StaticLaw`, `EdgeWeightLaw`, `AttackPolicy`, `ScanObservable`, and the
  process policy enums.
- Users can tune parameters and implement custom graph/activity types,
  but cannot add a law, process, observable, or estimator without editing
  CMC.rs.

The goal of v3 is **one generic engine, traits at every extension point,
and stage facades reduced to presets**. This is the percolation analogue
of CMC's `Hamiltonian` + `Algorithm<H>` + `ObservableSet`, QMC's
`QmcKernel` + `PositiveOperatorModel`, and MCMC's `LogDensity` +
`TransitionKernel<T>`.

## 1. Principles and constraints

1. **Additive evolution.** New traits/objects arrive first; legacy
   enums, adapters, and `ObservablePlan` remain as deprecated wrappers
   for one minor release. Removal only at the next major.
2. **Golden conservation.** No silent change to RNG draw order, Results
   column names, observable indicators, C0 JSON/HDF5 schema, or numeric
   semantics. Golden fixtures pin all of them.
3. **Preserve the reusable core.** `GraphView`/`UndirectedGraphView`,
   `VertexActivity`/`EdgeActivity`, generic `analyze`, `ComponentWorkspace`,
   `IncrementalUnionWorkspace`, RNG-free cores, C0 arrays, and
   `Run::from_parts` stay the foundation.
4. **Anchor-first.** Every scientific stage archives a PDF/markdown
   anchor and names the reproduced quantity, or records a definitional
   anchor. Statistical claims are finite-size brackets, never exact
   equalities.
5. **No GPU in v3.** GPU architecture remains issue #2's independent
   design track; percolation v3 is CPU-first.
6. **No new facade without a generic extension point.** A stage that can
   be expressed as a preset must not introduce another concrete adapter.
7. **Research honesty.** A stage that cannot certify its strongest claim
   stays Experimental with a precise negative result.

## 2. Dependency graph

```text
A1 Observable trait + open plan
 └─ A2 SampleView + OccupationLaw/WeightLaw
     └─ A3 generic adapter + estimator + C0 writer
         ├─ A4 policy traits (process/scan/attack, directed)
         ├─ A5 generic FSS toolkit
         └─ A6 generic specification/builder + provenance
X1 network I/O ── X2 multiplex/interdependent ── X3 hypergraph/simplicial
A2 + A4 + X1 ── X4 temporal/directed dynamics
A5 ── R1 C3 certification, R2 Q3 phase diagram, R3 F3 large-scale research
```

## 3. Architecture track A

Architecture stages have **definitional anchors**: bitwise/numerical
parity with the v2 golden outputs and an external-crate extension test,
not a physics paper.

### A1 — Observable trait and open plan

Problem: `StaticObservable` is a closed enum with a `pub(crate)` match;
new observables require crate edits.

Deliverable: `ScalarObservable` + `ObservableSpec` traits (name, required
graph/activity capabilities, compute from a `SampleView`), an open
`ObservablePlan` builder, built-ins reimplemented as zero-cost impls, and
a `StaticObservable` shim (deprecation attribute at the v4 release
boundary, once adapters no longer construct the enum internally). The
trait is named `ScalarObservable` because the crate root already exports
the Hamiltonian-generic `Observable` from `observables/`. `dyn` only for
optional third-party observables; built-in hot paths stay statically
dispatched.

Gates: every v2 scalar/array golden output equal under the same seeds;
an external test crate defines and runs a custom observable without CMC
edits; plan validation/rejection ordering preserved; allocation parity;
rustdoc example.

Deps: v2. Risks: object safety vs performance; plan ordering/duplicate
name semantics.

### A2 — SampleView, OccupationLaw, WeightLaw

Problem: `StaticConfiguration`, `ExternalBondConfiguration`,
`BondSample`, and the process-specific states do not share a view;
`StaticLaw` and `EdgeWeightLaw` are closed dispatch points.

Deliverable: `SampleView` constructible by third-party adapters
(counts, activity, labels, boundary queries); `OccupationLaw` trait for
site/edge/correlated occupation; `WeightLaw` trait for FPP edge weights.
Existing laws become blanket/built-in impls. Continuum and long-range
draws discover their edges during sampling, so they stay processes under
A4 rather than pretending to be occupation laws.

Gates: bitwise RNG replay for every existing law and draw contract;
typed rejection matrix; golden Results; zero steady-state allocation;
borrow/lifetime examples for a borrowed custom view.

Deps: A1. Risks: lifetime ergonomics; coherence with heterogeneous
fields; long-range laws whose edge universe is not materialized.

### A3 — Generic adapter, estimator, C0 writer

Problem: 15 adapters and 5 estimators duplicate plumbing.

Deliverable: one `PercolationAdapter<Law, Plan, Workspace>` with a single
`MonteCarlo` impl; a generic `Estimator` abstraction (count/sum/sum_sq,
mean/SEM) with reachability/defined-sample policies; a generic C0 vector
writer replacing the repeated pad/truncate code. Stage adapters become
thin presets (target: ≤3 adapter structs).

Gates: all v2 adapter/parameter tests pass; Results schema and C0
JSON/HDF5 round-trip unchanged; allocation parity; a custom law + custom
observable composed through the generic adapter.

Deps: A2. Risks: `FromParams` uniform-only constraint; column naming;
special cases (winding, i.i.d. sampling, sequential processes).

Progress note (2026-10-25): A3a adds `Moments` (count/sum/sum_sq,
unbiased variance, SEM) used by the passage, shape, chemical-distance, and
scan estimators. A3b makes the static adapter generic:
`PercolationAdapter<L = StaticLaw>` keeps `StaticPercolationMC` as a
compatibility alias, and a custom law composes with a custom observable
through `Run::from_parts`. Custom laws use the F2 `SamplingError`
vabulary (`Error = SamplingError`) so construction-time rejections keep
their typed `LawDomain` variant; laws with richer error contracts (the
correlated field) keep dedicated adapters for now. The generic C0 writer
and the remaining preset collapse (bond dynamics, continuum, scans) are
outstanding.

### A4 — Policy traits for processes, scans, and attacks

Problem: `AttackPolicy`, `AttackGrid`, `ScanObservable`, and the
Achlioptas/invasion/bootstrap/k-core policy enums are closed.

Deliverable: policy/grid/rule traits with built-in impls; a directed
incidence capability promoted to a first-class view trait; reference-order
guarantees documented for each policy.

Gates: reference-order equivalence for every built-in policy; directed
reuse tests; typed policy rejections; no changes to recorded coordinates
or provenance.

Deps: A3. Risks: semantic drift in tie-breaking; directed provenance
rules.

### A5 — Generic finite-size-scaling toolkit

Problem: C2/C3 estimators are stage-specific; crossing/Binder/
correlation-length/FSS fits are not reusable.

Deliverable: `percolation/fss/` — generic ladder runner, crossing
estimator, Binder cumulant, second-moment correlation length, collapse
fit, and bootstrap intervals. C2/C3 become instances.

Gates: reproduce the v2 crossing calibration and brackets; exact
small-lattice enumeration; Criterion target; documented finite-size
caveats.

Deps: A3. Risks: fit robustness; naming collisions with stage-specific
estimators.

### A6 — Generic specification, builder, and provenance

Problem: only `StaticPercolationMC` has `FromParams`; most stages require
bespoke construction and manually recorded provenance.

Deliverable: a `PercolationSpec`/builder covering topology + law +
process + plan + estimator, with `FromParams` presets and
`Run::from_parts` for non-schema laws; provenance attributes
(law/process/seed/estimator versions) recorded in Results.

Gates: round-trip configuration for all v2 presets; JSON provenance;
typed rejection matrix; no hidden defaults.

Deps: A3. Risks: schema creep; compatibility with existing Results
metadata.

## 4. Domain track X

### X1 — Network I/O and schema

Exactly **two** versioned formats: a text edge list (with weighted and
optional multiplex columns) and a JSON multiplex schema. No format zoo,
no silent coercion.

Anchor: definitional (format specification and round-trip identity).
Gates: round-trip and reference-graph equality; malformed-input
rejections with line/field provenance; benchmark for million-edge files;
feature-gated (`io`).

Deps: A2. Risks: format creep; encoding/versioning.

### X2 — Multiplex and interdependent networks

Layer-aware `SampleView`, mutual-component/cascade analysis, and preset
adapters for two or more layers.

Anchor (archive and verify at T-time): Buldyrev et al., *Catastrophic
cascade of failures in interdependent networks*, Nature 464, 1025 (2010);
Gao et al., *Networks formed from interdependent networks*, Nat. Phys. 8,
40 (2012). Gates: exact two-layer enumeration; anchor bracket; independent
review; documented cascade semantics (mutual vs. single-layer giant).

Deps: X1. Risks: under-specified cascade rules; node/layer identity.

### X3 — Hypergraph and simplicial percolation

Incidence view + hypergraph component analysis + a documented canonical
component-identity rule (vertex-min identity does not generalize
directly).

Anchor (archive and verify at T-time): Cooley–Kang–Koch high-order
percolation; Bianconi–Ziff random simplicial complex. Gates: incidence
hand cases; exact small-hypergraph enumeration; anchor bracket;
independent review.

Deps: X1. Risks: identity/tie rules; memory layout for incidence views.

### X4 — Temporal and directed dynamic percolation

Contact-process/directed dynamic percolation driven by the existing
`TransitionRate` substrate; temporal-network input.

Anchor: Hinrichsen, Adv. Phys. 49, 815 (2000) (archived) for absorbing-
state/directed percolation; a temporal-network reference (e.g. Holme–
Saramäki, Phys. Rep. 519, 97 (2012)) archived at T-time. Gates:
analytic/mean-field bracket; independent Monte Carlo reference;
statistical-only non-equilibrium caveats; independent review.

Deps: A2, A4, X1. Risks: non-equilibrium validation is weaker than
equilibrium exact enumeration.

## 5. Research track R

### R1 — C3 Weinrib–Halperin certification

Current state: calibrated pc-shift estimator; iid control ν=1.379 (CI
1.325–1.433); correlated rows ordered but prior-dominated and the
1.5-vs-1.25 gap is 1.1σ at L≤2048.

Outcome: promote C3 only with L≥4096–8192 or an independent 0.001-level
`p_c(inf)` (Binder cumulant/very-large-L crossing); otherwise keep
Experimental with the negative result. Deps: A5.

### R2 — Q3 phase diagram

Current state: exact self-dual curve; exploratory bimodality
(0.363/0.400/0.519/0.654 at q=2/4/6/8, L=16).

Outcome: latent-heat/two-mode estimator, q-ladder FSS, exact q≤4
cross-checks; promote only if certified. Anchor: Wu 1982 / FK duality
([note](../ref/md/Potts_FK_duality_note.md)). Deps: A5.

### R3 — F3 large-scale shape and aging

Current state: directional time constants, exact geodesic counts, and an
explicitly exploratory aging statistic.

Outcome: ADH-anchored large-scale brackets for the shape/time constant
and geodesic statistics; research-grade, may remain Experimental. Anchor:
Auffinger–Damron–Hanson 2017 (archived). Deps: A5.

## 6. Uniform definition of done (every stage)

1. Scope and explicit non-targets documented.
2. Archived anchor with named quantity, or a recorded definitional anchor.
3. Rustdoc contract and a runnable example.
4. Exact/reference gate on small systems.
5. Statistical gate with CI, or a precise negative result.
6. Bitwise draw-contract replay where randomness exists.
7. Typed invalid-domain rejection matrix.
8. Asserted zero-allocation/byte probe (or a documented reason).
9. Criterion target (or a documented reason).
10. `PERCOLATION.md` row/section, roadmap status, test registration,
    independent review report, and status transition.

## 7. Milestones (gate to proceed)

| Milestone | Content | Gate |
|---|---|---|
| M1 | A1–A2 | External-crate observable/law test; zero golden drift; bitwise law replay |
| M2 | A3–A4, A6 | Adapter collapse; Results/C0 schema parity; allocation parity; config round-trip |
| M3 | A5, X1 | FSS bracket reproduction; I/O round-trip and malformed-input matrix |
| M4 | X2–X4 | Exact two-layer/hypergraph gates; anchor brackets; independent reviews |
| M5 | R1–R3 | Research outcomes recorded; final v3 audit `docs/reviews/2026-*-v3-final-audit.md` |

## 8. PR sequence

Stacked, one stage per PR (≤500 LOC where practical):
`A1a` trait+shim → `A1b` migrate built-ins → `A2a` laws → `A2b`
SampleView → `A3a` estimator → `A3b` adapter collapse → `A3c` C0 →
`A6` config/provenance → `A4` → `A5` → `X1` → `X2`/`X3` (parallel) →
`X4` → `R1`/`R2`/`R3` → final audit.

No stage merges before its review report exists. Each stage updates its
roadmap status, `PERCOLATION.md`, and the relevant anchor/bibliography.

## 9. Compatibility and deprecation

- New traits are additive; legacy enums/adapters remain as deprecated
  wrappers for one minor release.
- Golden fixtures pin RNG draw order, column names, indicators, and C0
  schemas; a stage that changes them is rejected.
- New domains are feature-gated (`io`, `multiplex`, `hypergraph`,
  `temporal`); default features do not change behavior.
- `thermalization_sweeps = 0` remains the i.i.d. convention.

## 10. Documentation, validation, and CI

- The roadmap (this file), `PERCOLATION.md`, and `VALIDATION.md` are
  reconciled at every stage; stale status text is a release blocker.
- Heavy statistical tests are `#[ignore]` and run in nightly release;
  z-score tests honor `SCUTTLE_ZSCORE_SEEDS`.
- Criterion targets and allocation probes are added for every new
  algorithm family; the missing L/C/Q/F3 coverage from v2 is closed as
  part of A3/A5 and the X stages.
- Reviews follow `docs/reviews/YYYY-MM-DD-<stage>-review.md` with
  verdict, claim table, independent methods, findings, follow-ups.

## 11. Status

| Stage | Status |
|---|---|
| A1 Observable trait/open plan | Done (A1a+A1b: `ScalarObservable` + `SampleView`, capability-driven label/winding/moment requests, built-in `ObservableSpec` impls, `StaticObservable` shim; external-crate, golden, and 0-allocation gates green) |
| A2 SampleView/OccupationLaw/WeightLaw | Done (A2a: `OccupationLaw` for Bernoulli/correlated laws, `WeightLaw` for FPP with bitwise replay; A2b: public `SampleView::new`/`vacuum`/`from_analysis`, workspace-lifetime analyzer accessors, and checked public activity setters; external-crate law, golden, and rustdoc gates green) |
| A3 generic adapter/estimator/C0 | In progress (A3a: shared `Moments` accumulator replaces four duplicated count/sum/sum_sq estimators; A3b: `PercolationAdapter<L>` with `StaticPercolationMC` preset and external custom-law composition; C0 writer and remaining preset collapse pending) |
| A4 policy traits | Not started |
| A5 FSS toolkit | Not started |
| A6 specification/provenance | Not started |
| X1 network I/O | Done (feature-gated `io`: weighted edge-list v1 text + multiplex JSON v1 over one validated document; round-trip, reference-graph, rejection, and Criterion gates green; independent review pending) |
| X2 multiplex/interdependent | Not started |
| X3 hypergraph/simplicial | Not started |
| X4 temporal/directed dynamics | Not started |
| R1 C3 certification | Not started |
| R2 Q3 phase diagram | Not started |
| R3 F3 large-scale research | Not started |

## 12. Non-goals

- GPU execution/design (issue #2 remains separate).
- Distributed/MPI percolation and parallel union-find.
- Python bindings, visualization, or arbitrary I/O formats.
- Rewriting Validated physics.
- New stage-specific facades without a proven generic extension point.

## 13. Risks and open questions

- Static dispatch vs trait objects for observables and laws.
- `SampleView` lifetime/borrow design across owned and borrowed samples.
- Semantics of multiplex cascades and hypergraph component identity.
- Research stages may not certify and must stay Experimental.
- Deprecation window and semver-major timing.
- CI capacity for the growing nightly statistical suite.

## 14. Critical files

- Architecture: `CMC.rs/src/percolation/{observable.rs,carlo.rs,law,analyzer,workspace}`,
  `CMC.rs/src/topology/view.rs`.
- Domains: new `percolation/{io,multiplex,hypergraph,temporal}` modules and
  `topology/incidence.rs`.
- Validation: `CMC.rs/tests/suite.rs`, `CMC.rs/benches/percolation_allocations.rs`,
  `CMC.rs/Cargo.toml`, `.github/workflows/nightly.yml`, `justfile`.
- Documentation: this file, `CMC.rs/PERCOLATION.md`, `CMC.rs/VALIDATION.md`,
  `docs/ref/BIBLIOGRAPHY.md`, `docs/reviews/`.
