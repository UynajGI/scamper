# CMC.rs

CMC.rs is Scamper's classical Monte Carlo sampling layer. It is intentionally built **on Carlo.rs** rather than duplicating execution infrastructure:

- **Carlo.rs** owns RNG contexts, explicit run phases, scheduling, backends, accumulation/results, checkpoint orchestration and parallel tempering.
- **CMC.rs** owns classical configurations, physical energy models, transactional trial moves, target ensembles, update kernels and observables.

This revision keeps the existing lattice-spin functionality and public `ClassicalMC<Model, Algorithm>` composition, while adding continuous-system backends (periodic Lennard-Jones NVT/NPT/μVT particles, rigid molecules) and generalized-ensemble methods (Wang-Landau DOS estimation, multicanonical, umbrella sampling with canonical reweighting), a persistent classical worm framework with a ferromagnetic Ising high-temperature graph backend, and explicit classical-dynamics kernels (Kawasaki, Gillespie/BKL and hard-sphere event-chain Monte Carlo).

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

## Module organisation (Phases 1–2)

Source code is organised into five subdirectories plus three top-level adapter modules:

| Directory | Purpose |
|-----------|---------|
| `core/` | Move types, caches, `TrialEvaluator`, `Ensemble`, `AcceptanceRule`, visit schedules |
| `lattice/` | `CsrLattice` graph, `System` state, `Hamiltonian` traits, built-in models, `ProposalStrategy` |
| `algorithms/` | `Algorithm<H>` trait, 6 kernels (Metropolis, Wolff, SW, heat bath, microcanonical, hybrid) |
| `observables/` | `Observable<H>`, `DefaultObservableSet`, energy, magnetisation, correlation |
| `particle/` | Periodic cells, AoS coordinates, pair potentials, packed cell lists, translations and NVT/NPT/μVT adapters, rigid molecules with an optional dipolar external field |
| `generalized/` | Wang-Landau, frozen biases, DOS/histograms, exact enumeration and reweighting |
| `worm/` | Persistent physical/worm sectors, generic local driver and Ising graph representation |
| `dynamics/` | Kawasaki exchange, direct Gillespie, Fenwick BKL/n-fold way and hard-sphere event chains |
| `percolation/` | i.i.d. site/bond occupancy sampling, union-find cluster statistics, spanning-set crossing |
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


## Continuous Lennard-Jones NVT

The continuous path reuses exactly the same sampling transaction as the lattice path:

```text
TranslateParticle
    -> ProposedMove<ParticleTranslation<D>>
    -> ParticleSystem::evaluate_trial (read-only accepted state)
    -> CanonicalEnsemble
    -> MetropolisHastingsAcceptance
    -> ParticleSystem::commit_trial (accepted only)
```

A scheduler-ready monatomic simulation is available as `LennardJonesNvt<D>`:

```rust,ignore
use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use cmc_rs::LennardJonesNvt;

let mut params = Params::new();
params.set("n_particles", 108usize);
params.set("density", 0.8);
params.set("beta", 1.0);
params.set("cutoff", 2.5);
params.set("max_displacement", 0.1);

let results = Scheduler::new(RayonBackend::new(1), RunConfig::default())
    .run_one::<LennardJonesNvt<3>>(&params);
```

`OrthorhombicCell<D>` supports periodic minimum-image geometry in two, three, or other const-generic dimensions. `LennardJones` supports truncated, shifted-potential and shifted-force cutoffs plus Lorentz-Berthelot species mixing. `CellList<D>` stores packed particle buckets and applies accepted membership changes in O(1), without rebuilding on each trial. Translation scale adaptation is restricted to thermalization and frozen before production measurements.

The transaction and parameter contract follows the same five-part sampling foundation described above.

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

Debug/test builds audit automatically every 1024 completed sweeps. Optimized builds can enable the same cross-kernel policy with the `cache-audit` feature. Audits detect rather than repair lattice energy drift, particle cell-list inconsistencies, and generalized-ensemble macrostate cache mismatches.

## Explicit Carlo.rs lifecycle

`carlo_rs::Context` now carries a `RunPhase`:

```text
Initialization -> Thermalization -> Measurement -> Finished
```

Schedulers and `Run` call `MonteCarlo::on_phase_start` / `on_phase_end` at boundaries. CMC receives the mapped phase and permits proposal adaptation only during `Thermalization`; production kernels are frozen before the first measurement sweep. Legacy checkpoint phase is inferred from stored counters when needed.

For convergence-driven samplers, Carlo.rs provides `AdaptiveRunControl`, `Scheduler::run_controlled` and `Scheduler::run_controlled_with_state`. CMC.rs adds `WangLandauRunControl` and `IsingWangLandau` as adaptive ensemble implementations, with `WangLandauCore` for user-supplied axes and `EnergyBiasCore` for frozen umbrella/multicanonical production.

## Classical persistent worm

`IsingGraphWormMC` samples the ferromagnetic Ising high-temperature graph representation in explicit physical and two-defect worm sectors. Open, close and local head moves include their complete Hastings factors and are accepted in log space. The chain may remain open across sweep and checkpoint boundaries; optional endpoint-pair samples provide two-point correlation estimators.

**Multi-component lattices are supported.** The high-temperature graph ensemble factorizes over connected components, so `IsingGraphWormMC::from_lattice` (and the `IsingGraphWormEnsemble` it wraps) decomposes any disconnected lattice — isolated sites included — into per-component sub-lattices and runs one independent two-defect worm per component on a domain-separated derived stream (`RngStreamKey`; one salt per component per sweep from the shared context stream, so a checkpointed run replays exactly). Observables combine additively; total energy is measured when every component is physical, which preserves the product ensemble. The raw `IsingGraphWormModel` + `WormKernel` pair remains restricted to connected lattices: its single defect pair would silently freeze the other components, so `IsingGraphWormModel::new` rejects disconnected input loudly for direct users. Multi-defect / multi-leg worm algorithms are not implemented — the two-defect worm per component is the algorithm.

```rust,ignore
use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use cmc_rs::{Bond, BondType, CsrLattice, IsingGraphWormEnsemble, WormConfig};
use carlo_rs::RngStreamKey; // domain-separated per-component streams

// Connected lattice through the parameter/scheduler path:
let mut params = Params::new();
params.set("lattice_type", "square");
params.set("Lx", 16);
params.set("Ly", 16);
params.set("beta", 0.44);
params.set("worm_updates_per_sweep", 512);

let results = Scheduler::new(RayonBackend::new(1), RunConfig::default())
    .run_one::<cmc_rs::IsingGraphWormMC>(&params);

// Arbitrary (possibly disconnected) lattice directly:
let two_rings = CsrLattice::from_edges(8, vec![
    Bond::new(0, 1, BondType::Generic, 1.0),
    Bond::new(1, 2, BondType::Generic, 1.0),
    Bond::new(2, 3, BondType::Generic, 1.0),
    Bond::new(3, 0, BondType::Generic, 1.0),
    Bond::new(4, 5, BondType::Generic, 1.0),
    Bond::new(5, 6, BondType::Generic, 1.0),
    Bond::new(6, 7, BondType::Generic, 1.0),
    Bond::new(7, 4, BondType::Generic, 1.0),
]);
let ensemble = IsingGraphWormEnsemble::new(two_rings, 0.44, 1.0, WormConfig::default()).unwrap();
assert_eq!(ensemble.n_components(), 2);
```

The reusable `WormModel`/`WormKernel` boundary is intended for future integer-current, dimer and loop-gas representations without pretending that their defect constraints are identical.

## Rigid molecules and external fields

`MolecularMetropolisCore` moves whole rigid molecules (topology-grouped translations and plane rotations) with the same transactional Metropolis-Hastings path as the particle kernels. An optional one-body `DipolarExternalField` couples per-atom point charges to a uniform field: the molecular dipole is measured through minimum-image displacements (wrap-invariant), every trial's acceptance includes the `-E·μ` energy change, and the one-body energy is available through `external_field_energy`. Non-neutral molecules are rejected loudly, because a net charge would couple to the wrapped absolute position and break periodicity. `ParticleSystem::energy` remains the pair energy.

```rust,ignore
use cmc_rs::{DipolarExternalField, MolecularMetropolisCore};

let field = DipolarExternalField::new([2.0, 0.0], vec![0.5, -0.5])?;
let kernel = MolecularMetropolisCore::new(topology, 0.3, 0.4)?
    .with_external_field(field)?; // validates neutrality + charge-table size
```

## Classical dynamics and event time

Stage 6 adds three distinct dynamic paths:

- `KawasakiCore` for magnetization-conserving canonical spin exchange;
- `GillespieKernel` and `BklIsingKernel` for continuous-time rejection-free events;
- `HardSphereEventChain<D>` for lifted rejection-free hard-sphere chains.

Carlo.rs now records sweeps, attempts, accepted/executed moves and event time as separate clocks. `KineticIsingBklMC` advances fixed event-time observation windows, while event-chain lifted distance remains a separate geometric quantity.

## Site and bond percolation

`percolation/` samples ordinary percolation on any `CsrLattice` as i.i.d.
configurations rather than a Markov chain: every sweep redraws occupancy
(sites or bonds open independently with probability `p`), every measurement
runs union-find over the occupied subgraph. Set `thermalization_sweeps = 0`;
there is nothing to equilibrate.

```rust
use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use cmc_rs::PercolationMC;

let mut params = Params::new();
params.set("lattice_type", "square");
params.set("Lx", 32);
params.set("Ly", 32);
params.set("mode", "bond");   // or "site" (default)
params.set("p", 0.5);
let config = RunConfig {
    thermalization_sweeps: 0,
    measurement_sweeps: 100_000,
    binsize: 100,
    ..Default::default()
};
let results = Scheduler::new(RayonBackend::new(1), config)
    .run_one::<PercolationMC>(&params);
```

Measured observables: `Occupied`, `MaxCluster`, `SecondMoment`
(`sum(s_i^2)`), `NClusters` and `Spanning` (mean = crossing probability).
Crossing is tested between `spanning_from`/`spanning_to` site sets;
square lattices default to the left vs. right column, chains to the two end
sites, and arbitrary graphs take explicit comma-separated site lists.
`cluster_stats` and `UnionFind` are public for direct, RNG-free analysis of
fixed configurations.

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

## Validated domain

Per-solver validated domains, exact references, analytic limits and known
limitations live in [VALIDATION.md](VALIDATION.md). Highlights: q-state
Potts (q = 3, 4) is validated against full enumeration for the heat-bath,
Swendsen-Wang and Wolff kernels; the rigid-molecule solver reproduces the
Langevin/von Mises free-rotor answers under a dipolar external field; and
the Wang–Landau estimator terminates loudly (`UnreachableBins`) when the
configured visited fraction exceeds the physically reachable set. Invalid
input is rejected with errors across all solvers — never silently accepted
(see the input-validation audit in VALIDATION.md).

## Roadmap — algorithms and models not yet included

Everything below is **not implemented and not validated** today; it is
recorded here so the validated domain stays unambiguous (see VALIDATION.md
for what *is* covered). None of these block the production status of the
existing solvers — they are the next application frontiers.

### Sampling methods (the biggest holes)

- **Gibbs ensemble (Panagiotopoulos)** — two boxes with particle and volume
  exchange; the standard tool for vapour–liquid phase equilibrium. We have
  NPT and μVT but not the combination, so phase diagrams require workarounds.
- **Configurational-bias MC (CBMC) / recoil growth** — the standard route to
  chain-molecule and polymer insertion; without it long chains are effectively
  impossible to insert into dense phases.
- **Nested sampling** — complementary to Wang–Landau/multicanonical; directly
  compresses phase-space volume, strong for first-order transitions and
  density of states.
- **Transition-matrix Monte Carlo (TMMC) / broad histogram** — a more stably
  converging alternative to Wang–Landau.
- **Invaded cluster / probability-changing cluster** — automatic
  critical-point location.
- **Luijten–Blöte long-range weighted cluster** — without it, long-range
  models fall back to single-spin flips.

### Dynamics and irreversible methods

- **Momentum HMC / Langevin / Brownian-dynamics integrators** — molecular
  systems are pure MC moves today; the HMC in MCMC.rs is statistical-posterior
  HMC, not a physical momentum coupling.
- **Geometric cluster algorithm (Dress–Krauth)** — global reflection moves
  for hard disks/polygons; complements event chain, which covers hard
  spheres only.
- **Creutz demon / Q2R microcanonical dynamics** — a microcanonical family
  distinct from the existing over-relaxation.

### Model surface

- **Clock models Z_q (q = 5, 6)** — discrete XY with ESR-type topological
  transitions; Potts permutation symmetry cannot stand in for cyclic
  symmetry.
- **Anisotropic O(N)** — XXZ/easy-axis/easy-plane, single-ion terms,
  Dzyaloshinskii–Moriya; `ONModel` currently carries an isotropic `j` only.
- **Long-range interactions** — dipolar, 1/r^σ, Ewald summation. A gap on
  both sides: no long-range lattice bonds, and the particle potentials are
  Lennard-Jones (three cutoff treatments) and hard sphere — no Coulomb/Ewald.
- **Edwards–Anderson spin glass** — ±J/Gaussian random bonds as a first-class
  citizen with parallel-tempering coupling; the frustrated triangle is
  validated, random-bond distributions are not.
- **Vertex/ice models** — 6/8-vertex, F-model, spin ice, with the matching
  loop/cluster algorithms; the worm currently serves the Ising HT graph only.
- **Close-packed dimers** — the dimer model and a dimer worm.
- **Diluted Ising, random-bond Potts** — workhorses for first-order-transition
  studies.
- **Anisotropic hard particles** — ellipsoids, Gay–Berne, polygonal disks.
- **Lattice polymers / bond-fluctuation models** — practical only once CBMC
  exists.

Suggested priority if expansion continues: Gibbs ensemble → clock +
anisotropic O(N) → Ewald long-range → CBMC — these map onto the four largest
application fronts (phase equilibrium, magnetism, charged/dipolar systems,
polymers) and each reuses the existing solver skeletons and validation
framework.
