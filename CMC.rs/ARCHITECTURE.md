# CMC.rs architecture after the graph-core refactor

## Layering

```text
Carlo.rs
  Context / Scheduler / Backend / Measurements / Results / Parallel Tempering
      |
      v
ClassicalMC<H, A, O>
  H: Hamiltonian
  A: Algorithm<H>
  O: ObservableSet<H>
      |
      +-- System: configuration + physical energy + beta
      +-- CsrLattice: physical edges + CSR incidences
```

`ClassicalMC` implements Carlo.rs `MonteCarlo`, `FromParams`, and `ParallelTemperingCompatible`. It derives the update phase from Carlo.rs sweep counters, so adaptive kernels freeze before the first measurement sweep.

## Graph representation

`CsrLattice` contains:

- `edges: Vec<Bond>`: each physical undirected bond exactly once;
- `neighbors` and `edge_ids`: one incidence at each endpoint;
- `offsets`: CSR row boundaries;
- `n_bonds`: compatibility name for directed incidence count;
- `n_edges()`: physical edge count.

This supports arbitrary dimensions, irregular graphs, weighted/disordered couplings, parallel bonds and self-loops without an implicit divide by two.

## Model hierarchy

- `Hamiltonian`: general energy contract, directly implementable by multi-site/factor models.
- `PairInteraction`: optimized convenience contract for onsite + pair models; blanket-implements `Hamiltonian`.
- `Initializable`: hot/cold initialization independent of update proposals.
- `Proposable`: model-level symmetric Metropolis proposal.
- `ClusterModel`: complete Wolff/SW bond and transformation policy, with no panic defaults.
- `HeatBathable` / `ContinuousHeatBathable`: exact conditional samplers.
- `LocalFieldModel`: exact microcanonical over-relaxation capability.
- `Measurable`: model-native scalar order parameter.

XY and Heisenberg are aliases of `ONModel<2>` and `ONModel<3>`. `ONModel<D>` supports arbitrary spin dimension with const generics.

## Update kernels

All high-frequency temporary arrays are retained as algorithm workspaces.

- `MetropolisCore<S>`: random-order Metropolis-Hastings, log-domain acceptance, optional exact energy checks.
- `WolffCore`: one cluster with model-defined auxiliary and endpoint-specific bond probabilities.
- `SWCore`: physical-edge union-find and independent per-root transformations.
- `HeatBathCore`: exact discrete site conditionals without neighbor-vector allocation.
- `ContinuousHeatBathCore`: exact XY/Heisenberg conditionals.
- `MicrocanonicalCore`: local-field reflection with exact energy repair.
- `HybridCore<A,B>`: static composition without trait-object overhead.

## Observable flow

An `Observable` declares explicit raw moments (`E2`, `M2`, `M4`, etc.). `ObservableSet` records directly into Carlo.rs `Context`; `ClassicalMC` no longer identifies observables by matching their names.

## Compatibility boundary

The public `ClassicalMC<Model, Algorithm>` shape, built-in model names, `System` flat storage, CSR neighbor access, Carlo.rs scheduling and result keys are preserved. The changed model/cluster traits are intentionally stricter to make invalid algorithms impossible or explicit.
