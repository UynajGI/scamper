# CLAUDE.md

Scuttle — Monte Carlo framework with four Rust crates: Carlo.rs (core), QMC.rs (quantum), CMC.rs (classical), MCMC.rs (statistical).


Hooks via **lefthook** (`lefthook.yml` + `.lefthook/` scripts). pre-commit runs `fmt --check` → `cargo check` → `clippy` → `typos` on staged `.rs` (affected crates only); commit-msg enforces Conventional Commits; post-commit runs `codegraph sync` (non-blocking); pre-push runs `cargo deny` only (test moved to CI for fast pushes — run `just test` manually before pushing). Install: `lefthook install` or `just hooks`. Skip: `LEFTHOOK=0 git commit`. Lint policy (incl. `unsafe_code = "deny"`) is codified in `[workspace.lints]` (`Cargo.toml`) — no `RUSTFLAGS` needed; the sole exception is `Carlo.rs/src/backend/mpi.rs` (`#![allow(unsafe_code)]` for MPI FFI). CI (`.github/workflows/ci.yml`) runs fmt + clippy + test + deny (`--all-features`) as parallel jobs.

## Workspace

| Crate | Role | Description |
|-------|------|-------------|
| Carlo.rs | Core framework | `MonteCarlo` trait, `Scheduler`, `Context`, `Measurements`, `Merge`, `Backend` |
| QMC.rs | Quantum MC | General continuous-time lattice QMC (`LatticeSpinQmc` implements `MonteCarlo` + `FromParams`), impurity wormhole QMC
| CMC.rs | Classical MC | Lattice + particle + generalized ensembles + classical worm (Ising HT graph) + classical dynamics (Kawasaki, Gillespie/BKL, hard-sphere event-chain). `ClassicalMC` wrapper, `IsingGraphWormMC`, `KineticIsingMC`, `HardSphereEventChainMC` Carlo.rs adapters |
| MCMC.rs | Statistical MC | Euclidean-state kernels (RW/component/slice/Gibbs/composed, StaticHMC, NUTS), dual-averaging step-size tuning, windowed metric adaptation, unit/diagonal/dense metrics, constrained transforms, replica exchange, multi-chain diagnostics (R-hat, ESS, E-BFMI), Carlo.rs adapter |

## Carlo.rs Architecture

MonteCarlo trait → Scheduler.run_one() → Results flow:

| Module              | File                    | Purpose                                                            |
| ------------------- | ----------------------- | ------------------------------------------------------------------ |
| `MonteCarlo` trait  | `monte_carlo.rs`        | Core: `sweep(ctx)`, `measure(ctx)`, `Rng` type, lifecycle hooks    |
| `FromParams` trait  | `monte_carlo.rs`        | Construct model from `Params` dict                                 |
| `Context`           | `context.rs`            | RNG, measurements, sweep counter, `RunPhase`, checkpoint state     |
| `Run`               | `run.rs`                | Single run lifecycle, `step()`, `from_parts()` (no-`FromParams` path), `finalize_with_mc()` (recover MC after run), checkpoint/restart |
| `Scheduler`         | `scheduler.rs`          | Thermalization → measurement loop, `run_one`/`run_parallel`/`run_controlled`/`run_controlled_with_state` |
| `Backend`           | `backend/`              | `RayonBackend` (threads), `MpiBackend` (MPI)                       |
| `Measurements`      | `measurements.rs`       | Binned `Accumulator`, complex observables                          |
| `Merge`             | `merge.rs`              | Rebinning, autocorr time, `merge_results`, `merge_task_results`    |
| `Evaluable`         | `evaluable.rs`          | Jackknife resampling, `Evaluator`, `MultiplexEvaluator`            |
| `ResultTools`       | `output/resulttools.rs` | `dataframe()`, `measurement_from_obs()` — read-back `results.json` |
| `ParallelTempering` | `parallel_tempering.rs` | PT MC with chain scheduling                                        |
| `RunPhase`          | `phase.rs`              | Explicit lifecycle: Initialization→Thermalization→Measurement→Finished |
| `AdaptiveRunControl`| `run_control.rs`        | Algorithm-driven phase transitions, `RunDecision`                  |
| `CLI`               | `cli.rs`                | `carlo run/status/merge/delete`                                    |
| `Job`               | `job/`                  | `JobInfo`, `TaskInfo`, `TaskMaker`, progress tracking              |

## CMC.rs Architecture

Orthogonal traits instead of a monolithic trait — each concern is a separate trait.
Directories group related modules; public API is re-exported flat from `lib.rs`.

| Directory | Files | Purpose |
|-----------|-------|---------|
| `core/` | `move.rs`, `cache.rs`, `trial.rs`, `ensemble.rs`, `acceptance.rs`, `visit.rs` | Move types, incremental patches, `TrialEvaluator`, `MetropolisHastingsAcceptance`, visit schedules |
| `lattice/` | `graph.rs`, `state.rs`, `interaction.rs`, `models.rs`, `proposal.rs` | `CsrLattice` + builders, `System`, `Hamiltonian` + capability traits, built-in models, `ProposalStrategy` |
| `algorithms/` | `metropolis.rs`, `wolff.rs`, `swendsen_wang.rs`, `heat_bath.rs`, `microcanonical.rs`, `hybrid.rs`, `common.rs` | `Algorithm<H>` trait + 6 kernels, `SimulationPhase`, `checked_probability` |
| `observables/` | `energy.rs`, `magnetization.rs`, `correlation.rs`, `common.rs` | `Observable<H>`, `DefaultObservableSet`, `TotalEnergy`, `Magnetization`, `compute_correlation_1d` |
| `particle/` | `potential.rs`, `cell.rs`, `cell_list.rs`, `configuration.rs`, `state.rs`, `movement.rs`, `algorithm.rs`, `mc.rs`, `error.rs`, `batch.rs`, `grand.rs`, `mixture.rs`, `molecule.rs`, `volume.rs` | `PairPotential` trait, `LennardJones` 12-6 (Lorentz-Berthelot mixing, 3 cutoff modes), `OrthorhombicCell<D>`, packed `CellList`, `ParticleSystem` transactional evaluate/commit, `ParticleMetropolisCore` (NVT), `ParticleNptMetropolisCore` (NPT), `ParticleGrandCanonicalCore` (μVT), `MolecularMetropolisCore` (rigid molecules), `MoveMixture`, `LennardJonesNvt/Npt/MuVt` Carlo.rs adapters |
| `generalized/` | `axis.rs`, `bias.rs`, `histogram.rs`, `macrostate.rs`, `multicanonical.rs`, `wang_landau.rs`, `exact.rs`, `reweight.rs`, `error.rs` | `MacrostateAxis` (Binned/Discrete), `LogDensityOfStates`, `Histogram`, `LogBias` (Fixed/HarmonicUmbrella/Multicanonical), `EnergyBiasCore` (frozen production), `WangLandauState` + `WangLandauCore` (adaptive DOS estimation, flat-histogram + 1/t refinement, JSON checkpoint), `IsingWangLandau` (Carlo.rs adapter), `WangLandauRunControl` (AdaptiveRunControl), `canonical_reweight` (log-sum-exp), exact Ising DOS enumeration |
| `worm/` | `error.rs`, `model.rs`, `state.rs`, `kernel.rs`, `ising.rs`, `mc.rs` | Generic `WormKernel<M>` persisting open/close/step/bounce transitions, `WormState` physical/worm sector container, `IsingGraphWormModel` (HT graph, tanh weights), versioned JSON checkpoint `cmc-rs-ising-worm-v1`, `EndpointPairHistogram`, exact 2ⁿ enumeration (≤24 edges) |
| `dynamics/` | `error.rs`, `rate.rs`, `ising.rs`, `gillespie.rs`, `event_chain.rs`, `mc.rs` | `TransitionRate` trait, `KineticIsingModel` (Glauber/Kawasaki), `GillespieKernel` (rejection-free BKL), `EventChainKernel` (hard-sphere event-chain with lift), `KineticIsingMC`, `HardSphereEventChainMC` Carlo.rs adapters, exact event-time tracking |
| Top-level | `classical_mc.rs`, `multi_spin.rs`, `postprocess.rs` | `ClassicalMC` Carlo.rs adapter, `MultiSpinIsing`, derived observables |

Key patterns:
- `ClassicalMC<IsingModel, MetropolisCore>` → `Scheduler.run_one()` → `Results`
- Models are stateless — temperature (β) lives in `System`, not in model structs
- Algorithm trait: `sweep(&mut self, system, model, rng)` — directly mutates system.energy
- CSR lattice: cache-friendly neighbor iteration via `lattice.neighbors(site)` returning `&[usize]`
- `SmallVec<[f64; 3]>` for spin proposals (stack-allocated, no heap for spin_dim ≤ 3)
- Heat-bath (Glauber): discrete (Ising/Potts → `HeatBathable`) + continuous (XY/Heisenberg → `ContinuousHeatBathable`, vMF/Best-Fisher sampling)
- Over-relaxation: `MicrocanonicalCore` reflects spins across local field (ΔE=0), no acceptance
- Derived observables measured post-run from E²/M²/M⁴ moments stored in `Results`
- `lattice_type` param: `"chain"`, `"square"`, `"triangular"`, `"honeycomb"`, `"kagome"`
- Users can ignore `ClassicalMC` and compose manually for custom behavior
- Particle: `ParticleSystem` implements `TrialEvaluator` for `ParticleTranslation`, `ParticleBatchMove`, `GrandCanonicalMove`, and `IsotropicVolumeChange`. Evaluate never mutates accepted state; commit updates position + cell-list + energy atomically. `ParticleMC<D, P, A>` is generic over dimension, potential, and algorithm. Pre-built: `LennardJonesNvt<D>`, `LennardJonesNpt<D>`, `LennardJonesMuVt<D>`. `CanonicalParticleKernel` marker trait gates parallel-tempering eligibility (NVT/molecule kernels only, not NPT/μVT). `MolecularMetropolisCore<D>` composes `RigidMoleculeTranslation` + `RigidMoleculeRotation` via `MoveMixture`, with `MoleculeTopology` for atom grouping and `TorsionRotation` for local dihedral moves. `LogVolumeScale` adapts ln(V) random-walk step size toward target acceptance. `InsertDeleteParticle` proposes reversible insert/delete with species weights and particle-count bounds; `validate_potential()` checks species/potential compatibility upfront. `GrandCanonical` and `IsothermalIsobaric` ensembles drive acceptance via `ThermodynamicDelta`.
- Generalized ensembles: `MacrostateAxis` maps scalars to stable bins. `EnergyBiasCore<A, B>` runs frozen umbrella/multicanonical production on the lattice `TrialEvaluator` path. `WangLandauCore<A>` does adaptive DOS estimation: `Discovery → Adaptation (flat-histogram + optional 1/t) → FrozenProduction → Finished`. `WangLandauRunControl` plugs into `Scheduler::run_controlled_with_state` to return the estimator after scheduling. `CanonicalLatticeKernel` marker gates lattice PT to canonical kernels only (generalized-ensemble kernels never implement it). `canonical_reweight` uses log-sum-exp for stable canonical reconstruction from ln g(E). `IsingWangLandau` is the scheduler-ready reference with exact axis enumeration (≤24 sites).
- Classical worm: `WormKernel<M>` drives open/close/step/bounce transitions in an extended physical/worm sector space. `WormModel` trait owns defects, proposals, `log_reverse_over_forward` Hastings factors, weight deltas and transactional patches. Open: `ln η + ln P_close + ln N`; close: exact negative; step: `log_weight_ratio + log_reverse_over_forward + ln(1-P_close·1_new_coincident) - ln(1-P_close·1_old_coincident)`. `EndpointPairHistogram` estimates two-point correlations via `count(tail,head)/count(tail,tail)`. `IsingGraphWormMC` implements `MonteCarlo`+`FromParams` with versioned `cmc-rs-ising-worm-v1` JSON checkpoint. Zero-field, non-negative `J`, self-loop rejection.
- Classical dynamics: `GillespieKernel` uses BKL rejection-free sampling over `TransitionRate` trait. `KineticIsingModel` supports Glauber (spin-flip) and Kawasaki (spin-exchange conserving magnetization). `EventChainKernel<D>` performs hard-sphere event-chain Monte Carlo with exact lift-at-contact (zero-gap = contact, not overlap). Explicit event-time tracking for non-equilibrium observables.
- Model validation: `PairInteraction::validate_spin` (Ising ±1, Potts [0,q), O(N) unit-norm). Self-loops excluded from local field / heat-bath conditionals (constant energy contributes zero to spin-flip ΔE). `local_energy` uses edge-ID set for self-loop dedup, robust to CSR incidence reordering. `ClassicalMC` snapshot loading validates all spins.

## QMC.rs Architecture

Boundary: Carlo.rs owns runtime (scheduling, RNG, measurements, I/O); QMC.rs owns physics (representations, kernels, estimators).

Two production backends:

### Lattice QMC (primary)

| Module | File | Purpose |
|--------|------|---------|
| `graph` | `graph.rs` | `CsrGraph` — typed/weighted CSR adjacency + unique edge table; builders (chain, square, hypercubic, adjacency, CSR) |
| `local_space` | `local_space.rs` | `LocalHilbertSpace` trait, `SpinSpace` — site-resolved arbitrary `S` with exact integer `2S` algebra |
| `lattice::model` | `lattice/model.rs` | `SpinModelBuilder` → `SpinLatticeModel` — sparse positive `OperatorTerm` catalog, Marshall `Z2` gauge, Heisenberg/XY/XXZ/XYZ/tfim/generic |
| `lattice::vertex` | `lattice/vertex.rs` | `VertexKind` (positive local matrix elements), `Vertex` (sampled insertion), `Event` (worldline endpoint) |
| `lattice::configuration` | `lattice/configuration.rs` | `LatticeConfiguration` — product state + unsorted vertex vector; `WorldlineIndex` — time-ordered leg links |
| `lattice::scattering` | `lattice/scattering.rs` | `ScatteringTable` — exact local-detailed-balance directed-loop routing (LowBounce + Metropolis policies) |
| `lattice::updates` | `lattice/updates.rs` | `ContinuousLatticeEngine<M>` — diagonal add/remove + directed-loop blocks; journal-based rollback |
| `lattice::observables` | `lattice/observables.rs` | Magnetization, staggered magnetization, susceptibility, energy, vertex orders, edge SzSz correlation |
| `lattice::mc` | `lattice/mc.rs` | `LatticeSpinQmc` — Carlo.rs adapter (`MonteCarlo` + `FromParams`); warmup schedule adaptation |

### Spin-Boson QMC

| Module | File | Purpose |
|--------|------|---------|
| `impurity` | `impurity/` | Continuous-time retarded-interaction wormhole QMC (quantum impurity, bath samplers, scattering table, diagonal/loop updates, observables) |

### Shared

| Module | File | Purpose |
|--------|------|---------|
| `algorithm` | `algorithm.rs` | `QmcKernel<C,R>` trait, `UpdateSchedule` — reusable kernel contract |

Key patterns:
- `LatticeSpinQmc` wraps `ContinuousLatticeEngine<SpinLatticeModel>` and implements `MonteCarlo` + `FromParams` — drop-in for `Scheduler.run_one()`
- Models are compiled from physical couplings into positive `K=C-H` operator catalogs; the engine has no model-name branches
- `CsrGraph` is data, not an algorithm — `from_csr` mirrors CMC's layout without crate coupling
- Marshall gauge: BFS solves bipartite `Z2` phases; rejects frustrated/non-stoquastic models
- `rand` 0.10: use `RngExt` for `random()`/`random_range()`, `Rng` for trait bounds (`R: Rng + ?Sized`)
- `UpdateSchedule` (diagonal proposals, directed loops, max loop steps) fixed during measurement; adaptation allowed during thermalization
- `from_sparse()` on `OperatorTerm` is the extension point for future bosonic/fermionic catalogs

## Features

- `hdf5` — checkpoint/measurement files (needs `libhdf5-dev`)
- `mpi` — distributed backend (needs `libopenmpi-dev`)

## Key Patterns

1. Define struct with config state
2. `impl MonteCarlo` — `type Rng`, `fn sweep()`, optional `fn measure()`
3. `impl FromParams` — construct from `Params` (or use `Run::from_parts()` for closure/sampler-based models that can't implement `FromParams`)
4. Call `ctx.measure("Name", value)` in `sweep()` or `measure()`
5. Optional: `impl MonteCarloCheckpoint` for HDF5 save/load

### MCMC.rs Architecture

Statistical inference layer with Euclidean-state and Hamiltonian kernels:

| Module | Directory | Purpose |
|--------|-----------|---------|
| `LogDensity` / `DifferentiableLogDensity` | `target.rs` | Value-only and combined value/gradient target contracts; closure adapters for both |
| `ChainState` | `state/` | Position + cached log density + synchronized accepted-state gradient + PT exchange |
| `TransitionKernel<T>` | `kernel/` | Target-typed trait — RW/component/slice/Gibbs, static composition, `StaticHmc<M>` and `Nuts<M>` |
| `metric` | `metric/` | Unit, diagonal and dense inverse-mass geometries; momentum, velocity and kinetic energy |
| `LeapfrogIntegrator` | `integrator.rs` | Private `PhasePoint` integration workspace and invalid-trajectory reporting |
| `adaptation` | `adaptation/` | RW covariance/scale adaptation plus HMC dual averaging and windowed metric adaptation |
| `Bijector` / `DifferentiableBijector` | `transform/` | Value transforms, gradient pullback and log-Jacobian gradients |
| `TransformedTarget` | `transform/target.rs` | Constrained target wrapper implementing value-only or differentiable unconstrained density |
| `tempering` | `tempering.rs` | Rayon local transitions + alternating generic neighboring exchange |
| `MemoryTrace` | `trace/` | Contiguous posterior storage, thinning, divergence/energy/tree-depth/depth-limit flags, JSON/HDF5 |
| `diagnostics` | `diagnostics/` | Rank-normalized R-hat, bulk/tail ESS, MCSE, per-chain E-BFMI, depth-limit-hit totals |
| `ChainCheckpoint` | `checkpoint.rs` | Serde envelope for kernel workspace, warmup, RNG, state and trace |
| `McmcSampler` | `carlo_adapter.rs` | Carlo lifecycle adapter including HMC scalar diagnostics |

Key patterns:
- `McmcSampler<T, K, Tr>` implements `MonteCarlo`; use `Run::from_parts()` to avoid `FromParams`
- `TransitionKernel<T>` remains target-typed, enabling model-specific Gibbs kernels and compile-time HMC gradient capability
- `StaticHmc<M>` integrates only private workspace and atomically commits position/log-density/gradient after acceptance
- Metrics store inverse mass `G=M^-1`; dense momentum sampling solves `L^T p=z` for `G=L L^T`
- HMC divergence is an invalid trajectory or absolute energy error above the configured threshold
- `HmcWarmup` uses dual averaging plus fast/slow/fast windows; metric changes reset dual averaging and terminal warmup retunes step size
- Configured HMC warmup length must match the runner/scheduler warmup transition count; incomplete warmup is rejected
- Differentiable transforms provide analytic pullback and Jacobian gradients for Positive/Interval/Ordered/Simplex/Product
- `TransitionReport.subtransitions` counts elementary transitions in composed kernels; `merge()` aggregates across children
- Traces never mix chain IDs; multi-chain uses Rayon with deterministic per-chain seeds
- Component-wise, slice, Gibbs and HMC kernels preserve accepted-state validity on proposal/integration errors
- Replica exchange keeps targets/kernels/RNGs/traces fixed to ladder slots and swaps only synchronized states
- `Nuts<M>` doubles a binary Hamiltonian trajectory until U-turn, divergence or depth exhaustion; multinomial candidate selection in the log domain; optional `StepSizeSearch` (pre-warmup step-size tuning); shares dual averaging and windowed metric adaptation with `StaticHmc<M>`; metric-aware U-turn via `displacement_dot_velocity`; gradient-check validation available
- `TransitionReport` includes `energy`, `acceptance_statistic`, `tree_depth`, `max_tree_depth_reached` for dynamic-HMC diagnostics; `merge()` uses weighted averaging for `acceptance_statistic` and OR for `max_tree_depth_reached`
- Default-feature MCMC.rs v0.5 validation passes fmt, check, Clippy and all tests on Rust 1.90.0
- Optional HDF5 validation is blocked: `hdf5-sys 0.8.1` does not recognize the installed HDF5 1.14.5 header format


## Behavioral Guidelines

**Tradeoff:** These guidelines bias toward caution over speed. For trivial tasks, use judgment.

## 1. Think Before Coding

**Don't assume. Don't hide confusion. Surface tradeoffs.**

Before implementing:
- State your assumptions explicitly. If uncertain, ask.
- If multiple interpretations exist, present them - don't pick silently.
- If a simpler approach exists, say so. Push back when warranted.
- If something is unclear, stop. Name what's confusing. Ask.

## 2. Simplicity First

**Minimum code that solves the problem. Nothing speculative.**

- No features beyond what was asked.
- No abstractions for single-use code.
- No "flexibility" or "configurability" that wasn't requested.
- No error handling for impossible scenarios.
- If you write 200 lines and it could be 50, rewrite it.

Ask yourself: "Would a senior engineer say this is overcomplicated?" If yes, simplify.

## 3. Surgical Changes

**Touch only what you must. Clean up only your own mess.**

When editing existing code:
- Don't "improve" adjacent code, comments, or formatting.
- Don't refactor things that aren't broken.
- Match existing style, even if you'd do it differently.
- If you notice unrelated dead code, mention it - don't delete it.

When your changes create orphans:
- Remove imports/variables/functions that YOUR changes made unused.
- Don't remove pre-existing dead code unless asked.

The test: Every changed line should trace directly to the user's request.

## 4. Goal-Driven Execution

**Define success criteria. Loop until verified.**

Transform tasks into verifiable goals:
- "Add validation" → "Write tests for invalid inputs, then make them pass"
- "Fix the bug" → "Write a test that reproduces it, then make it pass"
- "Refactor X" → "Ensure tests pass before and after"

For multi-step tasks, state a brief plan:
```
1. [Step] → verify: [check]
2. [Step] → verify: [check]
3. [Step] → verify: [check]
```

Strong success criteria let you loop independently. Weak criteria ("make it work") require constant clarification.

---

**These guidelines are working if:** fewer unnecessary changes in diffs, fewer rewrites due to overcomplication, and clarifying questions come before implementation rather than after mistakes.

## CodeGraph

This project has a CodeGraph index (`.codegraph/`). ALWAYS use
`codegraph explore` (CLI) or `codegraph_explore` (MCP) BEFORE Grep/Glob/Read
to explore the codebase. One call returns verbatim source + call paths +
blast-radius summary — faster and cheaper than reading files yourself.

Syncs automatically on post-commit via lefthook.
