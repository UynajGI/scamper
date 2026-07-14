# QMC.rs architecture

## Workspace boundary

```text
Carlo.rs: scheduling, RNG streams, warmup/measurement, bins, backends, results
CMC.rs:   classical Monte Carlo (unchanged)
QMC.rs:   quantum representations, model catalogs, updates, estimators
```

No Carlo.rs or CMC.rs source is modified by this delivery.

## Shared QMC foundation

`src/algorithm.rs` defines `QmcKernel<C,R>` and `UpdateSchedule`. A measured
sweep always performs a fixed amount of work. Runtime adaptation is confined
to warmup.

## Lattice pipeline

```text
CsrGraph
  -> LocalHilbertSpace
  -> PositiveOperatorModel / sparse OperatorTerm catalog
  -> LatticeConfiguration + WorldlineIndex
  -> diagonal insertion/removal
  -> local detailed-balance directed loops
  -> estimators
  -> LatticeSpinQmc (Carlo.rs adapter)
```

### `graph.rs`

`CsrGraph` stores a unique typed/weighted edge table plus CSR neighbor rows.
Algorithms never infer dimensionality or shape. `from_csr` mirrors the data
layout used by CMC without creating a crate dependency between the two physics
packages.

### `local_space.rs`

`LocalHilbertSpace` is the finite local-basis contract. `SpinSpace` is the
production implementation and supports site-resolved arbitrary `S`. Statistics
are explicit. Fermions are reserved but rejected by the positive engine until
a sign-aware backend exists.

### `lattice/model.rs`

A physical spin Hamiltonian is compiled into positive sparse matrix elements of
`K=C-H`. `OperatorTerm` is local-space agnostic; `SpinModelBuilder` supplies
spin algebra and automatic shifts. The model compiler:

1. validates couplings;
2. solves a graph-wide Marshall `Z2` gauge;
3. rejects incompatible/frustrated signs;
4. emits diagonal and off-diagonal vertex kinds;
5. builds a local scattering table;
6. builds an importance proposal distribution over terms.

The engine has no model-name branches.

### `lattice/configuration.rs`

A configuration is a product-basis state at time zero plus an unsorted packed
vector of `(tau, term, kind)` insertions. `WorldlineIndex` sorts endpoint events
per site and builds periodic leg links. Topology is immutable during loop
updates, so the same index can be reused across a loop block.

### `lattice/scattering.rs`

An extended local state is `(kind, entrance_leg, delta)` with `delta=±1` for a
raising/lowering discontinuity. Compatible states form a local graph.
Symmetric path flows enforce

```text
W_a P(a -> b) = W_b P(b -> a).
```

The default residual-flow policy reduces bounce; a symmetric-proposal
Metropolis policy is retained as a reference. Bounce preserves the worm charge
for arbitrary Spin-S. The global engine samples simple non-self-intersecting
loops; a self-intersection is reversal-symmetrically rejected instead of using
spin-1/2-only flip assumptions.

### `lattice/updates.rs`

Diagonal proposals choose add/remove with a fixed probability one half,
including at expansion order zero. A zero-order removal is a null proposal.
This boundary rule is required by the acceptance ratios

```text
R_add    = beta W / ((n_diag + 1) q_term)
R_remove = n_diag q_term / (beta W).
```

Directed loops modify matrix-element kinds, follow periodic worldline links,
and close on the original discontinuity. Journaling provides exact rollback on
self-intersection, incompatible closure or a safety limit.

### `lattice/observables.rs`

Implemented raw estimators include uniform/gauge-staggered magnetization,
moments and static susceptibilities, expansion-order energy, vertex orders,
edge `Sz Sz` correlation and update diagnostics.

## Spin-boson pipeline

`src/impurity` remains the continuous-time retarded-interaction wormhole
backend. It is a sibling to `lattice`, not a special case inside the lattice
engine, because its two-time bath vertices and bath proposals have different
configuration semantics.

## Extension path

- More spin models: emit new positive `OperatorTerm` catalogs.
- Bosons with a finite cutoff: implement `LocalHilbertSpace` and positive terms.
- Fermions: add a sibling sign-aware/determinant backend; do not reuse the
  positive engine without parity handling.
- Longer-range interactions: add weighted edges; no algorithm change.
- Multi-site operators: extend `TermLocation`/catalog construction while keeping
  the same leg-link and scattering concepts.
- Projector/VMC: add sibling modules implementing Carlo-facing adapters and
  shared QMC scheduling traits.
