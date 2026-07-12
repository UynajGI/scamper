# CMC.rs

CMC.rs is the classical lattice Monte Carlo toolbox in the Scuttle workspace. It is intentionally built **on Carlo.rs**:

- Carlo.rs owns RNG contexts, scheduling, thermalization/measurement phases, backends, binning/results and parallel tempering.
- CMC.rs owns graph topology, physical models, update kernels and classical observables.

## Built-in composition

```rust
use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use cmc_rs::{ClassicalMC, IsingModel, MetropolisCore};

type IsingMC = ClassicalMC<IsingModel, MetropolisCore>;

let mut params = Params::new();
params.set("Lx", 32);
params.set("Ly", 32);
params.set("beta", 0.44);
params.set("J", 1.0);
params.set("initial_state", "hot");

let results = Scheduler::new(RayonBackend::new(1), RunConfig::default())
    .run_one::<IsingMC>(&params);
```

The same wrapper works with `WolffCore`, `SWCore`, `HeatBathCore`, `ContinuousHeatBathCore`, `MicrocanonicalCore`, and statically composed `HybridCore` updates when the model implements the corresponding capability trait.

## Arbitrary weighted graph

`CsrLattice` is not restricted to square lattices. It stores each physical undirected edge exactly once and keeps CSR incidences for fast local access:

```rust
use cmc_rs::{Bond, BondType, CsrLattice};

let graph = CsrLattice::from_edges(
    3,
    vec![
        Bond::new(0, 1, BondType::Generic, 1.0),
        Bond::new(1, 2, BondType::Generic, 0.7),
        Bond::new(0, 2, BondType::Generic, -0.2),
    ],
);
```

Parallel edges and self-loops are represented explicitly. `neighbors`, `offsets`, and the historical directed-incidence `n_bonds` field remain available for compatibility; `n_edges()` is the physical bond count.

## Implementing a model

Most pair models only implement `PairInteraction`, then opt into the algorithms they support:

```rust
use cmc_rs::{Bond, PairInteraction};

struct MyModel { j: f64 }

impl PairInteraction for MyModel {
    fn spin_dim(&self) -> usize { 1 }
    fn coupling(&self) -> f64 { self.j }

    fn bond_energy(&self, left: &[f64], right: &[f64], bond: &Bond) -> f64 {
        -self.j * bond.weight * left[0] * right[0]
    }
}
```

Models with genuine multi-site/factor interactions can implement `Hamiltonian` directly, including exact `local_energy`, `delta_energy` and `compute_total_energy` behavior.

## Correctness changes in this refactor

- O(N) Wolff/SW activation probabilities depend on both endpoint projections: `1-exp(-2 beta J w (s_i·r)(s_j·r))` for equal projected signs.
- Every Swendsen-Wang cluster receives an independent state/reflection decision.
- Cluster updates recompute physical energy after the batch transformation, removing order-dependent cache drift.
- Metropolis uses the full Hastings correction and log-domain acceptance.
- Adaptive proposals receive every accept/reject result and adapt only while Carlo.rs is thermalizing.
- Hamiltonians return physical energy only; beta is applied once by the canonical algorithm.
- Unknown lattice names and invalid dimensions are configuration errors rather than silent fallbacks.

See [`ARCHITECTURE.md`](ARCHITECTURE.md) and [`MIGRATION.md`](MIGRATION.md).
