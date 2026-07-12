# Migration notes

## Existing simulations

Normal users of:

```rust
ClassicalMC<Model, Algorithm>
```

need no scheduling or call-site changes. Carlo.rs still drives the run and CMC still exposes the existing built-in kernels.

## Carlo.rs lifecycle

`Context` now exposes `phase() -> RunPhase`. Custom schedulers should call `enter_phase` and the model lifecycle hooks at boundaries:

```rust,ignore
ctx.enter_phase(RunPhase::Thermalization);
mc.on_phase_start(RunPhase::Thermalization, &mut ctx);
// warmup sweeps
mc.on_phase_end(RunPhase::Thermalization, &mut ctx);

ctx.enter_phase(RunPhase::Measurement);
mc.on_phase_start(RunPhase::Measurement, &mut ctx);
```

The hook methods have default empty implementations, so existing `MonteCarlo` implementations remain source-compatible.

Algorithms with a convergence-defined warmup can implement `AdaptiveRunControl<MC>` and use `Scheduler::run_controlled`. Fixed-count simulations continue to use `run_one` / `run_parallel` unchanged.

CMC's two-variant `SimulationPhase` remains source-compatible. `ClassicalMC` maps Carlo.rs `RunPhase` into it; scheduled update sweeps receive `Thermalization` or `Measurement`, while any out-of-contract initialization/finished sweep uses the frozen measurement kernel.

## Custom algorithms

The `Algorithm` interface remains `sweep_with_phase`. Direct `sweep` still executes a frozen production/measurement kernel.

A custom local Metropolis-like algorithm can now reuse:

- `ProposedMove<M>`;
- `TrialEvaluator<Model, M>`;
- `Ensemble<Delta>`;
- `metropolis_hastings_step`.

The generic driver handles the Hastings correction and log-domain acceptance. The state backend owns evaluation and atomic cache commit.

## Custom state/move backends

Implement `TrialEvaluator<Model, Movement>` for the accepted state:

```rust,ignore
impl TrialEvaluator<MyModel, MyMove> for MyState {
    type Delta = MyDelta;
    type Patch = MyReusablePatch;

    fn evaluate_trial(
        &self,
        model: &MyModel,
        movement: &MyMove,
        patch: &mut Self::Patch,
    ) -> Self::Delta {
        // Do not mutate accepted state.
    }

    fn commit_trial(&mut self, movement: &MyMove, patch: &Self::Patch) {
        // Apply configuration and cache changes together.
    }
}
```

The owning kernel constructs the patch once and should reuse it between trials; the generic trait does not require a particular constructor.

## Hamiltonian implementations

For onsite/pair models, implement `PairInteraction`:

- keep `spin_dim()` and `coupling()`;
- implement `bond_energy(left, right, bond)`;
- optionally implement `onsite_energy(site, spin)`.

The blanket implementation supplies total energy, one-site delta and affected-edge batch delta.

For genuine multi-site interactions, implement `Hamiltonian` directly. The default `batch_delta_energy` uses reusable scratch storage and an exact total-energy difference. Override it only when a proven incremental implementation is available.

The compatibility `beta` argument remains in `local_energy` / `compute_total_energy`, but returned values must be physical energy. Ensemble policy applies beta exactly once.

## Proposal strategies

`ProposalStrategy::propose` returns `ProposedSpin`, including `log_reverse_over_forward`. It receives every decision through `record_result` and a per-sweep `finish_sweep(adaptation_enabled)` call.

`OPSSStrategy` adapts only in `RunPhase::Thermalization`. It is frozen in measurement, initialization and finished phases.

## Batch moves

Use `BatchSpinMove` for an atomic set of distinct site replacements. Duplicate sites are rejected during evaluation. `System` updates all changed spins and its energy cache in one commit.

For pair models, the batch evaluator counts every affected physical edge exactly once. For direct Hamiltonians, it falls back to exact scratch recomputation.

## Visit order

Local kernels default to `VisitSchedule::RandomPermutation`. Cache-oriented sequential traversal is selectable:

```rust
use cmc_rs::{MetropolisCore, VisitSchedule};

let algorithm = MetropolisCore::new().with_visit_schedule(VisitSchedule::Sequential);
```

## Model and lattice traits retained from the prior refactor

- `ClusterModel` uses explicit Wolff/SW auxiliary and endpoint-specific bond probability methods.
- `Initializable`, `Proposable`, `HeatBathable`, `ContinuousHeatBathable`, `LocalFieldModel`, and `Measurable` remain capability traits.
- `CsrLattice::from_edges` remains the arbitrary weighted multigraph entry point.
- Unknown `lattice_type` values remain configuration errors.
