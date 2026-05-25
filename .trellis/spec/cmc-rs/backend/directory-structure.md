# CMC.rs Directory Structure

## Module layout

```
CMC.rs/src/
├── lib.rs           — pub mod + re-exports
├── lattice.rs       — CsrLattice, BondType, builders (chain, square, hypercubic, triangular, honeycomb, kagome)
├── system.rs        — System { lattice, spins, energy, beta }
├── hamiltonian.rs   — Hamiltonian, ClusterModel, Proposable, Measurable, HeatBathable traits
├── models.rs        — IsingModel, PottsModel, XYModel, HeisenbergModel
├── algorithm.rs     — Algorithm<H> trait + MetropolisCore<S>, WolffCore, SWCore, HeatBathCore
├── proposal.rs      — ProposalStrategy<H> trait + StandardStrategy, OPSSStrategy
├── classical_mc.rs  — ClassicalMC<H, A> + FromHamiltonianParams + JSON checkpoint
├── observables.rs   — Observable<H> trait + TotalEnergy, Magnetization, DefaultObservableSet
├── postprocess.rs   — Derived observables: susceptibility(), specific_heat(), binder_cumulant()
└── multi_spin.rs    — MultiSpinIsing (64-replica bit-packed) + MonteCarlo + FromParams + PT
```

## Layer dependency graph

```
classical_mc ──► algorithm ──► proposal
    │               │
    ▼               ▼
  system ◄──────── hamiltonian ◄── models
    │                   ▲
    ▼                   │
  lattice          postprocess (read-only)
    ▲
    │
multi_spin (standalone, impl MonteCarlo directly)
```

- `lattice` is pure data, no dependencies
- `system` depends on `lattice`
- `model` depends on `lattice` (for energy signatures)
- `algorithm` depends on `model` + `system`
- `proposal` depends on `model` + `system`
- `classical_mc` composes everything → impl MonteCarlo + FromParams

## Layer responsibilities

### Lattice (`lattice.rs`)
- `BondType` enum — direction labels (ChainX, SquareX/Y/Z, CubicX/Y/Z)
- `Neighbor { target: usize, bond_type: BondType }` — adjacency entry
- `Lattice { sites: Vec<Vec<Neighbor>>, n_sites, n_bonds }` — adjacency list
- Builders: `build_chain(n, pbc)`, `build_square(w, h, pbc)`, `build_hypercubic(dims, bond_types, pbc)`

### System (`system.rs`)
- `System { lattice, spins, energy }` — **all fields `pub`**
- `spins` is flattened: `spins[site*spin_dim .. (site+1)*spin_dim]`
- `energy` is the running total, incrementally updated by algorithms
- Helpers: `spin_at(site, spin_dim) -> &[f64]`, `spin_at_mut(site, spin_dim) -> &mut [f64]`

### Model (`model.rs`)
- `Model` trait — stateless physics computations
- Required methods: `spin_dim()`, `coupling()`, `beta()`, `local_energy(spins, lattice, site, proposed)`, `propose(rng)`, `magnetization(spins)`, `random_cluster_spin(rng)`, `opposite_spin(spin, rng)`
- Default methods: `fk_bond_probability()` = `1 - exp(-2βJ)`, `compute_total_energy()`, `normalize_spin()`, `random_spin()`
- Concrete models: `IsingModel(j, beta)`, `PottsModel(j, beta, q)`, `XYModel(j, beta)`, `HeisenbergModel(j, beta)`

### Algorithm (`algorithm.rs`)
- `Algorithm<M: Model>` trait — single method: `sweep(&mut self, system: &mut System, model: &M, rng: &mut impl Rng)`
- `MetropolisCore<S>` — generic over `ProposalStrategy<M>`, default `S = StandardStrategy`
- `WolffCore` — BFS cluster growth with FK bonds
- `SWCore` — union-find percolation with cluster flip
- Algorithm **must** update `system.energy` directly — no return value

### Proposal (`proposal.rs`)
- `ProposalStrategy<M: Model>` trait — `propose(&mut self, model, system, site, rng)` + `adapt_after_sweep(&mut self, model)`
- `StandardStrategy` — delegates to `model.propose()`
- `OPSSStrategy` — over-relaxation with adaptive sigma; works for scalar (Ising) and vector (XY, Heisenberg) spins

### ClassicalMC (`classical_mc.rs`)
- `ClassicalMC<H: Hamiltonian + Measurable, A: Algorithm<H>> { system, model, algorithm, observables }`
- `impl MonteCarlo` — sweep → `algorithm.sweep(system, model, &mut ctx.rng)`, measure → iterates observables + records E²/M²/M⁴ moments
- `impl FromParams` — parses L/Lx/Ly/Lz/pbc for lattice, delegates model parsing to `FromHamiltonianParams`
- `FromHamiltonianParams` trait — one impl per model type
- JSON checkpoint: `save_snapshot() -> Json`, `load_snapshot(&Json) -> Result<()>`

### Derived Observables (`postprocess.rs`)
- Pure functions taking `&carlo_rs::Results`: `susceptibility()`, `specific_heat()`, `binder_cumulant()`
- Depend only on Carlo.rs — no CMC internal imports
- Compute derived quantities from E²/M²/M⁴ moments (recorded in ClassicalMC::measure)

## Adding a new model

1. Add struct to `model.rs`
2. Implement `Model` trait
3. Add `FromModelParams` impl in `classical_mc.rs`
4. Add unit tests (energy formulas, magnetization, FK probability)
5. Add re-export in `lib.rs`

## Adding a new algorithm

1. Add struct to `algorithm.rs`
2. Implement `Algorithm<M>` — must update `system.spins` AND `system.energy`
3. Add test verifying convergence to known ground state
4. Add re-export in `lib.rs`
