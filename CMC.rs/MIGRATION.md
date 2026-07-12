# Migration notes

## Hamiltonian implementations

For normal onsite/pair models, replace a direct `Hamiltonian` implementation with `PairInteraction`:

- keep `spin_dim()` and `coupling()`;
- replace `local_energy()` with `bond_energy(left, right, bond)`;
- optionally add `onsite_energy(site, spin)`.

The blanket implementation supplies exact total and one-site delta energies using physical edges.

For genuine multi-site interactions, implement `Hamiltonian` directly. `compute_total_energy` and `local_energy` are explicit so the framework does not assume pairwise double counting.

## ClusterModel

The previous scalar/vector methods (`fk_bond_probability`, `opposite_spin`, `reflect`, and panic defaults) were replaced by:

- `wolff_auxiliary`;
- `sw_bond_auxiliary`;
- `sw_cluster_auxiliary`;
- `cluster_bond_probability`;
- `transform_cluster_spin`.

This supplies the endpoint spins and physical `Bond`, which is necessary for weighted graphs and correct O(N) embedded-cluster probabilities.

## ProposalStrategy

`propose` now returns `ProposedSpin`, including `log_reverse_over_forward`. Strategies receive `record_result` for every attempt and `finish_sweep(adaptation_enabled)` after each sweep.

The historical `OPSSStrategy` name is retained, but its implementation is now a symmetric adaptive random plane rotation. Adaptation is disabled in the Carlo.rs measurement phase.

## Algorithms

Custom algorithms implement `sweep_with_phase`. The compatibility `sweep` method remains and executes a frozen measurement-phase kernel.

## Observables

Custom observable collections can be used with:

```rust
ClassicalMC<MyModel, MyAlgorithm, MyObservableSet>
```

Raw moments are declared by `MomentSpec`; no special observable names are inspected by the wrapper.

## Lattice construction

Unknown `lattice_type` values now return `CarloError::InvalidConfig`. The accepted built-ins are `chain`, `square`, `cubic`, `hypercubic`, `triangular`, `honeycomb`, and `kagome`.

Use `CsrLattice::from_edges` or `try_from_edges` for arbitrary graphs. Existing code reading `neighbors`, `offsets`, or `n_bonds` continues to work.
