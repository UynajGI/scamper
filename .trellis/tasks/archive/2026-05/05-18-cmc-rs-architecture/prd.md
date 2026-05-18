# brainstorm: CMC.rs architecture design

## Goal

Rewrite CMC.rs with clean separation of concerns: lattice topology, spin state, physics model, algorithm, proposal strategy — all composable into a `MonteCarlo` impl that runs through Carlo.rs Scheduler.

## Final Architecture

```
src/
├── lib.rs              — pub mod + re-exports
├── lattice.rs          — Lattice, BondType, builders (chain, square, hypercubic)
├── system.rs           — System { lattice, spins, energy } — mutable state
├── model.rs            — Model trait + Ising, Potts, XY, Heisenberg impls
├── algorithm.rs        — Algorithm trait + MetropolisCore, WolffCore, SWCore
├── proposal.rs         — ProposalStrategy trait + StandardStrategy, OPSSStrategy
└── classical_mc.rs     — ClassicalMC<M, A> impl MonteCarlo + FromParams
```

### Layer responsibilities

```
ClassicalMC<M: Model, A: Algorithm<M>>     ← impl MonteCarlo + FromParams
  ├── system: System                        ← mutable state (spins, energy)
  ├── model: M                              ← stateless physics formulas
  └── algorithm: A                          ← sweep logic, mutates system

Algorithm.sweep(&mut self, system, model, rng)
  ├── reads model.local_energy() / model.propose()
  ├── reads system.spins, system.lattice
  └── writes system.spins, system.energy

Model trait:
  ├── spin_dim(), coupling(), beta()
  ├── local_energy(spins, lattice, site, proposed) → f64
  ├── propose(rng) → Vec<f64>
  └── fk_bond_probability()

ProposalStrategy<M: Model> trait (independent, pluggable into Metropolis):
  ├── propose(&mut self, model, system, site, rng) → Vec<f64>
  └── adapt_after_sweep(&mut self, model)
```

### Key decisions

| # | Decision | Rationale |
|---|----------|-----------|
| D1 | Split System / Model / Algorithm | Clear ownership: state vs formulas vs update logic |
| D2 | System: pub fields (spins, energy, lattice) | Fine-grained access for algorithm flexibility |
| D3 | Algorithm trait: `sweep(&mut self, system, &M, rng)` only | Minimal trait, algorithm directly updates system.energy |
| D4 | Model trait: local_energy + propose + fk_bond_probability | Physics defined once, algorithms are model-agnostic |
| D5 | ProposalStrategy kept independent | OPSS has its own state (sigma), orthogonal to Model physics |
| D6 | ClassicalMC<M, A>: Model + Algorithm both generic | Each impls FromParams, ClassicalMC composes them |
| D7 | Energy tracking: Algorithm mutates system.energy directly | Cluster algorithms compute energy deltas incrementally during sweep |

### Carlo.rs integration

```rust
impl<M: Model + FromParams, A: Algorithm<M> + FromParams> MonteCarlo for ClassicalMC<M, A> {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        self.algorithm.sweep(&mut self.system, &self.model, &mut ctx.rng);
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        ctx.measure("Energy", self.system.energy);
        ctx.measure("Magnetization", self.model.magnetization(&self.system));
    }

    fn name(&self) -> &'static str {
        self.algorithm.name()
    }
}
```

### Type alias sugar

```rust
type IsingMetropolis = ClassicalMC<IsingModel, MetropolisCore>;
type IsingWolff     = ClassicalMC<IsingModel, WolffCore>;
type PottsSW         = ClassicalMC<PottsModel, SWCore>;
```

## Requirements

* All 5 layers (lattice, system, model, algorithm, proposal) are independent crates/modules
* `ClassicalMC<M, A>` implements `MonteCarlo` + `FromParams`
* Users can compose manually for custom behavior (ignore ClassicalMC)
* Support algorithm × model combinations: Metropolis×5, Wolff×3 discrete (Ising/Potts/XY), SW×3 discrete
* Whole rewrite — existing CMC.rs code is not preserved

## Acceptance Criteria

* [ ] `ClassicalMC<IsingModel, MetropolisCore>` runs end-to-end through `Scheduler.run_one()`
* [ ] Onsager validation test passes through Carlo.rs pipeline
* [ ] All 3 algorithms work with Ising model
* [ ] `FromParams` per model and per algorithm
* [ ] Proposal strategy pluggable for MetropolisCore

## Definition of Done

* Tests: unit (energy formulas, bond counting) + integration (Scheduler.run_one)
* Lint / clippy / test green
* Onsager exact solution validated: energy per site ± 1e-3, specific heat peak near Tc=2.269
* Carlo.rs lib.rs docs already updated (this session)

## Out of Scope

* New algorithms beyond Metropolis/Wolff/SW
* New models beyond Ising/Potts/XY/Heisenberg
* MPI-specific optimizations
* HDF5 checkpointing (Carlo.rs default hooks suffice)
* Carlo.rs lattice.rs modification (CMC.rs uses its own lattice module)

## Implementation Plan

1. **`lattice.rs`** — Lattice, BondType, builders — mostly from existing CMC, cleaned up
2. **`system.rs`** — System with pub fields
3. **`model.rs`** — Model trait + IsingModel (verify with Onsager first)
4. **`algorithm.rs`** — Algorithm trait + MetropolisCore
5. **`classical_mc.rs`** — ClassicalMC<M, A> impl MonteCarlo + FromParams
6. **Verify** — Ising + Metropolis through Scheduler.run_one, check energy converges
7. **`proposal.rs`** — ProposalStrategy + Standard/OPSS
8. **Expand** — WolffCore, SWCore
9. **Expand** — PottsModel, XYModel, HeisenbergModel
10. **Validate** — Onsager test through Carlo.rs pipeline

## Technical Notes

* Carlo.rs docs in `Carlo.rs/src/lib.rs` updated this session with full API reference
* `LatticeParams` in Carlo.rs is 2D-only and separate from CMC's lattice module
* `ctx.rng` is pub field — algorithms access RNG directly from context
* `FromParams::from_params` receives RNG for random initialization (Scheduler passes `&mut ctx.rng`)
* CMC.rs Cargo.toml already depends on `carlo-rs`, no changes needed
