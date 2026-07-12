# CMC.rs sampling-core architecture

## Scope of this revision

This revision changes infrastructure, not the modeled domains. The existing lattice-spin models and algorithms are rebuilt on a reusable trial/ensemble/cache foundation. No particle backend, density-of-states sampler, Wang-Landau controller or worm configuration is added in this package.

## Layering

```text
Carlo.rs
  Context / RunPhase / Scheduler / Run / Backend / Measurements / PT
      |
      v
ClassicalMC<H, A, O>
  H: Hamiltonian
  A: Algorithm<H>
  O: ObservableSet<H>
      |
      +-- System: accepted lattice-spin configuration + energy cache + beta
      +-- CsrLattice: physical edges + CSR incidences
      +-- sampling core:
            ProposedMove
            TrialEvaluator
            ThermodynamicDelta
            Ensemble
            reusable Patch
```

`ClassicalMC` remains the adapter implementing Carlo.rs `MonteCarlo`, `FromParams`, and `ParallelTemperingCompatible`.

## Carlo.rs lifecycle extension

Carlo.rs exposes an explicit `RunPhase`:

- `Initialization`;
- `Thermalization`;
- `Measurement`;
- `Finished`.

`Context::phase()` is the source of truth. `Scheduler`, incremental `Run`, and MPI parallel tempering enter phases explicitly and invoke the default-compatible `MonteCarlo::on_phase_start` / `on_phase_end` hooks. `ClassicalMC` maps that lifecycle to its source-compatible two-variant `SimulationPhase`, removing counter inference while preserving custom algorithm implementations.

`Context::advance_sweep` reaches the fixed-run thermalization boundary after exactly the configured number of warmup sweeps. While an explicit `Thermalization` phase is active, it remains authoritative even if that fixed counter is exceeded. Serialized legacy checkpoints without a phase field are normalized from their counters; HDF5 restoration reconstructs the phase and re-enters its start hook.

`AdaptiveRunControl<MC>` is the optional algorithm-driven path. `Scheduler::run_controlled` asks it for one of four decisions after each completed sweep: continue adaptation, begin production, continue production, or stop. The existing fixed-count `run_one` and `run_parallel` paths are unchanged.

## Transactional trial path

```text
proposal mechanics             state backend                 target policy              acceptance rule
------------------             -------------                 -------------              ---------------
ProposedMove<Movement>  ->  TrialEvaluator::evaluate  ->  Ensemble::log_weight_ratio  ->  AcceptanceRule::log_acceptance
       |                         |                                  |                           |
       | Hastings correction     | Delta + reusable Patch           | log π(new)-log π(old)      | MH, Barker, …
       +-------------------------+-------------------------------+--+---------------------------+
                                                                   |
                                                          accept/reject in log domain
                                                                   |
                                                       accepted -> commit_trial
                                                       rejected -> discard patch
```

The accepted state is never modified during evaluation. This prevents rollback bugs and lets a future backend patch multiple caches atomically.

### Current generic types

- `ProposedMove<M>` stores `M` and `ln q(old|new)-ln q(new|old)`.
- `TrialEvaluator<Model, Movement>` defines `Delta` and reusable `Patch` types.
- `ThermodynamicDelta` currently drives canonical lattice updates and reserves energy, particle-count, volume and Jacobian changes for later backends.
- `Ensemble<D>` maps a physical delta to a log target-weight ratio.
- `AcceptanceRule<D>` converts ensemble ratio + proposal asymmetry into log acceptance probability. `MetropolisHastingsAcceptance` is the only implementation shipped in Phase 1.
- `CanonicalEnsemble` applies beta exactly once.

## Lattice-spin backend

The compatibility backend implements the generic contract with:

- `SiteSpinMove` + `EnergyPatch`;
- `BatchSpinMove` + `BatchEnergyPatch`;
- `System` as `TrialEvaluator` for both movement types.

`System::energy` is always a cached physical energy, never beta-weighted.

### General Hamiltonians

A direct `Hamiltonian` implementation receives a correct default `batch_delta_energy`: the proposed batch is materialized in a reusable scratch vector and total energy is evaluated once. This is intentionally conservative and supports arbitrary multi-site interactions.

### Pair interactions

The `PairInteraction` blanket implementation uses generation-stamped workspaces to evaluate only affected terms:

```text
O(number of changed sites + number of incident physical edges)
```

Each physical edge is counted once, including parallel bonds and self-loops. Wolff therefore avoids the previous unconditional full-graph energy recomputation. SW changes all sites in the usual case, but still uses the same atomic cache contract.

## Graph representation

`CsrLattice` contains:

- `edges: Vec<Bond>`: each physical undirected bond once;
- `neighbors` and `edge_ids`: incidences at endpoints;
- `offsets`: CSR row boundaries;
- `n_bonds`: compatibility name for incidence count;
- `n_edges()`: physical bond count.

It supports irregular graphs, arbitrary dimension, weighted/disordered interactions, parallel bonds and self-loops without an implicit divide-by-two convention.

## Source module layout (Phase 1)

| Directory | Contents |
|-----------|----------|
| `core/` | `move.rs`, `cache.rs`, `trial.rs`, `ensemble.rs`, `acceptance.rs`, `visit.rs` |
| `lattice/` | `graph.rs`, `state.rs`, `interaction.rs`, `models.rs`, `proposal.rs` |
| `algorithms/` | `common.rs` (trait + phase), 6 kernel files |
| `observables/` | `energy.rs`, `magnetization.rs`, `correlation.rs`, `common.rs` |
| Top-level | `classical_mc.rs`, `multi_spin.rs`, `postprocess.rs` |

All public types are re-exported flat from `lib.rs`.

## Update kernels on the foundation

- `MetropolisCore<S>` uses the generic MH driver, a reusable energy patch and cached site order.
- `WolffCore` uses generation-stamped membership and an atomic `BatchSpinMove`.
- `SWCore` uses physical-edge union-find, independent per-root transformations and an atomic batch commit.
- `HeatBathCore` and `ContinuousHeatBathCore` use the site trial evaluator to update the energy cache after exact conditional sampling.
- `MicrocanonicalCore` uses the same site evaluator after local-field reflection, retaining numerical cache consistency.
- `HybridCore<A, B>` composes existing kernels statically without trait-object dispatch.

## Visit schedules

`VisitSchedule` currently provides:

- `RandomPermutation` (default), with storage reused across sweeps;
- `Sequential`, restoring identity order if the workspace was previously shuffled.

This is the extension point for checkerboard or graph-color schedules later; those are not included in this revision.

## Observable flow

An `Observable` declares explicit moments. `ObservableSet` writes into Carlo.rs `Context`; no string matching is used to decide whether `E2`, `M2`, or `M4` should be recorded.

## Compatibility boundary

Preserved:

- `ClassicalMC<Model, Algorithm>` user composition;
- built-in model and algorithm names;
- `System` flat spin storage;
- CSR neighbor access and physical-edge graph;
- existing parameter recipes, result keys, snapshots and PT adapter;
- direct `Algorithm::sweep`, which uses the frozen measurement phase.

New extension surface:

- `RunPhase` and lifecycle hooks in Carlo.rs;
- `ProposedMove`, `TrialEvaluator`, `Ensemble`, `ThermodynamicDelta`;
- site/batch move and patch types;
- reusable visit schedules.
