# CMC Particle Core v1

The `cmc_rs::particle` module implements the stage-2 continuous-system path:
periodic orthorhombic cells, AoS particle coordinates, species-aware pair
potentials, Lennard-Jones 12-6 cutoffs, packed cell lists, transactional
single-particle translation, and canonical NVT Metropolis-Hastings sampling.

## Core transaction

```text
TranslateParticle
    -> ProposedMove<ParticleTranslation<D>>
    -> ParticleSystem::evaluate_trial
    -> CanonicalEnsemble
    -> MetropolisHastingsAcceptance
    -> ParticleSystem::commit_trial (accepted only)
```

`evaluate_trial` never mutates accepted positions, physical energy, or cell-list
membership. The reusable `ParticleEnergyPatch` stores the local energy change,
old/new cell indices, and candidate-neighbor scratch. `commit_trial` updates the
coordinate, packed membership, and energy exactly once.

## Lennard-Jones cutoff modes

- `CutoffTreatment::Truncated`
- `CutoffTreatment::ShiftedPotential`
- `CutoffTreatment::ShiftedForce`

Multi-species parameters use Lorentz-Berthelot mixing. Exact overlaps and
non-finite separations are represented as an infinite trial barrier and are
rejected without modifying accepted state.

## Scheduler-ready type

```rust,ignore
use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use cmc_rs::LennardJonesNvt;

let mut params = Params::new();
params.set("n_particles", 108usize);
params.set("density", 0.8);
params.set("beta", 1.0);
params.set("sigma", 1.0);
params.set("epsilon", 1.0);
params.set("cutoff", 2.5);
params.set("max_displacement", 0.1);

let results = Scheduler::new(RayonBackend::new(1), RunConfig::default())
    .run_one::<LennardJonesNvt<3>>(&params);
```

Optional parameters include `box_length` or axis-specific `Lx`, `Ly`, `Lz`,
`cutoff_treatment`, warmup adaptation settings, and `energy_check_interval`.
The proposal scale adapts only during thermalization and is frozen for
measurement.


## Validation and benchmark targets

```bash
cargo test -p cmc-rs --test particle_core_test
cargo test -p cmc-rs --test particle_metropolis_test
cargo bench -p cmc-rs --bench particle_bench
```

The integration suite checks read-only rejection, atomic accepted cache patches,
minimum-image boundary crossing, cell-list completeness versus brute force,
fixed-seed reproducibility, thermalization-only adaptation, and the two-particle
canonical energy mean against deterministic relative-coordinate quadrature.
