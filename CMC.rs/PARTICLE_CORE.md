# CMC Particle Core

The `cmc_rs::particle` module implements continuous-system Monte Carlo:
periodic orthorhombic cells, AoS particle coordinates, species-aware pair
potentials, Lennard-Jones 12-6 cutoffs, packed cell lists, and three
Metropolis-Hastings ensembles (NVT, NPT, grand-canonical μVT) plus
rigid-molecule collective moves.

## Ensembles and scheduler-ready types

| Ensemble | Target | Scheduler type |
|----------|--------|----------------|
| `CanonicalEnsemble` | NVT: `π ∝ exp(-βE)` | `LennardJonesNvt<D>` |
| `IsothermalIsobaric` | NPT: `π ∝ exp[-β(E+PV)]` | `LennardJonesNpt<D>` |
| `GrandCanonical` | μVT: `π ∝ exp(-βE + log_activity · N)` | `LennardJonesMuVt<D>` |
| (NVT) | Rigid molecules | `MolecularMetropolisCore<D>` (via `ParticleMC`) |

All three pre-built Lennard-Jones types implement `FromParams` and `MonteCarlo`.

## Core transaction (NVT single-particle)

```text
TranslateParticle
    -> ProposedMove<ParticleTranslation<D>>
    -> ParticleSystem::evaluate_trial
    -> CanonicalEnsemble
    -> MetropolisHastingsAcceptance
    -> ParticleSystem::commit_trial (accepted only)
```

The same `TrialEvaluator` + `Ensemble` + `AcceptanceRule` pipeline drives NPT
volume changes (`IsotropicVolumeChange`), grand-canonical insert/delete
(`GrandCanonicalMove`), and multi-atom batch moves (`ParticleBatchMove`).

`evaluate_trial` never mutates accepted positions, physical energy, or cell-list
membership.

## NPT

`IsotropicVolumeChange` with `LogVolumeScale` adaptive proposals. Volume commit
rebuilds the CellList from scratch. The minimum-image cutoff guard rejects moves
where `cutoff > 0.5 * min_cell_length`.

## Grand-canonical

`InsertDeleteParticle` with configurable species weights, branch probabilities,
and particle-count bounds. `ParticleGrandCanonicalCore<D>` runs per-particle
translations plus exchange attempts.

## Rigid molecules

`MoleculeTopology` maps atoms to molecules. `MolecularMetropolisCore<D>` visits
each molecule once per sweep, selecting translation or rotation via a
`MoveMixture<K>`. `TorsionRotation` provides local dihedral-angle moves in 3D.

## Lennard-Jones cutoff modes

- `CutoffTreatment::Truncated`
- `CutoffTreatment::ShiftedPotential`
- `CutoffTreatment::ShiftedForce`

Multi-species parameters use Lorentz-Berthelot mixing. Exact overlaps and
non-finite separations are represented as an infinite trial barrier and are
rejected without modifying accepted state.

## Scheduler-ready NVT example

```rust,ignore
let mut params = Params::new();
params.set("n_particles", 108usize);
params.set("density", 0.8);
params.set("beta", 1.0);
params.set("cutoff", 2.5);
params.set("max_displacement", 0.1);
let results = Scheduler::new(RayonBackend::new(1), RunConfig::default())
    .run_one::<LennardJonesNvt<3>>(&params);
```

NPT adds `pressure`, `max_log_volume_change`. μVT adds `log_activity` (or
`chemical_potential` + `thermal_wavelength`), `maximum_particles`,
`exchange_attempts`.

## Validation and benchmark targets

```bash
cargo test -p cmc-rs --test particle_core_test
cargo test -p cmc-rs --test particle_metropolis_test
cargo test -p cmc-rs --test particle_stage3_test
cargo bench -p cmc-rs --bench particle_bench
```

The integration suite checks read-only rejection, atomic accepted cache patches,
minimum-image boundary crossing, cell-list completeness versus brute force,
fixed-seed reproducibility, thermalization-only adaptation, the two-particle
canonical energy mean against deterministic relative-coordinate quadrature,
rigid-rotation bond-length preservation, torsion axial/radial invariance, NPT
volume-delta correctness, grand-canonical insertion/deletion cache consistency,
and ideal-gas Poisson particle-number mean.
