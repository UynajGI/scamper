# CMC.rs

CMC.rs is Scuttle's classical Monte Carlo sampling layer. It is intentionally built **on Carlo.rs** rather than duplicating execution infrastructure:

- **Carlo.rs** owns RNG contexts, explicit run phases, scheduling, backends, accumulation/results, checkpoint orchestration and parallel tempering.
- **CMC.rs** owns classical configurations, physical energy models, transactional trial moves, target ensembles, update kernels and observables.

This revision is a foundation refactor. It keeps the existing lattice-spin functionality and public `ClassicalMC<Model, Algorithm>` composition, but moves its updates onto reusable sampling primitives. Particle systems, Wang-Landau and worm sectors are deliberately not included yet.

## Existing user entry point

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

The same wrapper continues to support `WolffCore`, `SWCore`, `HeatBathCore`, `ContinuousHeatBathCore`, `MicrocanonicalCore`, and statically composed `HybridCore` updates when the model implements the required capability trait.

## Module organisation (Phase 1)

Source code is organised into four subdirectories plus three top-level adapter modules:

| Directory | Purpose |
|-----------|---------|
| `core/` | Move types, caches, `TrialEvaluator`, `Ensemble`, `AcceptanceRule`, visit schedules |
| `lattice/` | `CsrLattice` graph, `System` state, `Hamiltonian` traits, built-in models, `ProposalStrategy` |
| `algorithms/` | `Algorithm<H>` trait, 6 kernels (Metropolis, Wolff, SW, heat bath, microcanonical, hybrid) |
| `observables/` | `Observable<H>`, `DefaultObservableSet`, energy, magnetisation, correlation |
| Top-level | `classical_mc.rs` (Carlo.rs adapter), `multi_spin.rs`, `postprocess.rs` |

The public API is re-exported flat from `lib.rs` — user code sees no change.

## Sampling foundation

The local Metropolis path now has five independent parts:

```text
ProposalStrategy
    -> ProposedMove<M>
    -> TrialEvaluator<Model, M>::evaluate_trial (no accepted-state mutation)
    -> Ensemble<Delta>::log_weight_ratio
    -> AcceptanceRule<Delta>::log_acceptance
    -> commit_trial only when accepted
```

Important public building blocks are:

- `ProposedMove<M>`: move payload plus the Hastings proposal-density correction;
- `TrialEvaluator<Model, Movement>`: state-specific trial evaluation and atomic cache commit;
- `ThermodynamicDelta`: physical extensive-variable changes;
- `Ensemble<Delta>` and `CanonicalEnsemble`: target-weight policy;
- `AcceptanceRule<D>` and `MetropolisHastingsAcceptance`: log-acceptance formula (extensible to Barker, rejection-free);
- `EnergyPatch` / `BatchEnergyPatch`: reusable cache workspaces;
- `SiteSpinMove` / `BatchSpinMove`: the current lattice-spin move backend;
- `VisitSchedule` / `SiteOrder`: reusable site traversal without per-sweep allocation.

A custom move can use the generic Metropolis-Hastings driver directly:

```rust,ignore
let outcome = cmc_rs::metropolis_hastings_step(
    &mut state,
    &model,
    &proposal,
    &ensemble,
    &MetropolisHastingsAcceptance,
    &mut patch,
    &mut rng,
);
```

Evaluation is transactional: a rejected move never mutates the accepted configuration or its cache. An accepted move commits configuration and cache changes once.

## Batch updates and energy caches

`Hamiltonian::batch_delta_energy` is the general multi-site path. Its default implementation uses a reusable scratch configuration and one exact energy evaluation, so direct multi-body Hamiltonians remain correct.

`PairInteraction` overrides this with an affected-edge implementation:

- changed onsite terms are visited once;
- each affected physical edge is visited once using generation stamps;
- parallel edges and self-loops are handled explicitly;
- no full graph recomputation is needed for Wolff clusters;
- Swendsen-Wang uses the same atomic batch contract.

An optional exact audit remains available on Metropolis:

```rust
use cmc_rs::MetropolisCore;

let algorithm = MetropolisCore::new().with_energy_check_interval(1_000);
```

## Explicit Carlo.rs lifecycle

`carlo_rs::Context` now carries a `RunPhase`:

```text
Initialization -> Thermalization -> Measurement -> Finished
```

Schedulers and `Run` call `MonteCarlo::on_phase_start` / `on_phase_end` at boundaries. CMC receives the mapped phase and permits proposal adaptation only during `Thermalization`; production kernels are frozen before the first measurement sweep. Legacy checkpoint phase is inferred from stored counters when needed.

For future convergence-driven samplers, Carlo.rs also provides `AdaptiveRunControl` and `Scheduler::run_controlled`. This revision only supplies the lifecycle/control protocol; CMC does not yet add Wang-Landau or another adaptive ensemble.

## Arbitrary weighted graph

`CsrLattice` stores each physical undirected edge exactly once and keeps CSR incidences for local access:

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

Parallel edges and self-loops are represented explicitly. Historical `neighbors`, `offsets`, and directed-incidence `n_bonds` fields remain available; `n_edges()` is the physical bond count.

## Implementing a pair model

Most onsite/pair models only implement `PairInteraction`, then opt into the algorithms they support:

```rust
use cmc_rs::{Bond, PairInteraction};

struct MyModel {
    j: f64,
}

impl PairInteraction for MyModel {
    fn spin_dim(&self) -> usize {
        1
    }

    fn coupling(&self) -> f64 {
        self.j
    }

    fn bond_energy(&self, left: &[f64], right: &[f64], bond: &Bond) -> f64 {
        -self.j * bond.weight * left[0] * right[0]
    }
}
```

Models with genuine multi-site/factor interactions can implement `Hamiltonian` directly and inherit the correct scratch-backed batch path.

See [`ARCHITECTURE.md`](ARCHITECTURE.md), [`MIGRATION.md`](MIGRATION.md), and [`VALIDATION_REPORT.md`](VALIDATION_REPORT.md).
