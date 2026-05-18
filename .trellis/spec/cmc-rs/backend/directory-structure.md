# CMC.rs Directory Structure

## Module layout

```
CMC.rs/src/
├── lib.rs           — pub mod + re-exports (Algorithm, Model, System, ClassicalMC, etc.)
├── lattice.rs       — Lattice, BondType, Neighbor, builders
├── system.rs        — System { lattice, spins, energy }
├── model.rs         — Model trait + IsingModel, PottsModel, XYModel, HeisenbergModel
├── algorithm.rs     — Algorithm<M> trait + MetropolisCore<S>, WolffCore, SWCore
├── proposal.rs      — ProposalStrategy<M> trait + StandardStrategy, OPSSStrategy
└── classical_mc.rs  — ClassicalMC<M, A> impl MonteCarlo + FromParams
```

## Layer dependency graph

```
classical_mc ──► algorithm ──► proposal
    │               │
    ▼               ▼
  system ◄──────── model
    │
    ▼
  lattice
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
- `ClassicalMC<M: Model, A: Algorithm<M>> { system, model, algorithm }`
- `impl MonteCarlo` — sweep → `algorithm.sweep(system, model, &mut ctx.rng)`, measure → `ctx.measure("Energy", energy)` + `ctx.measure("Magnetization", mag)`
- `impl FromParams` — parses L/Lx/Ly/Lz/pbc for lattice, delegates model parsing to `FromModelParams`
- `FromModelParams` trait — one impl per model type

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
