# QMC.rs architecture

## Responsibility boundary

`Carlo.rs` is the runtime framework. It owns parameter containers, deterministic
RNG seeding, scheduler/backends, thermalization and measurement phases,
binning, result reduction, and optional HDF5/MPI infrastructure.

`QMC.rs` is the physics and representation layer. It owns Hamiltonian/model
catalogs, worldline or operator-string configurations, update kernels,
representation invariants, and raw estimators.

The boundary is intentionally narrow:

```text
Carlo.rs Scheduler
    -> SpinBosonQmc::from_params
    -> SpinBosonQmc::sweep / measure
        -> WormholeEngine
            -> WormholeConfiguration
            -> SpinBosonModel
            -> Bath
```

Future lattice QMC, variational Monte Carlo, projector QMC, and other methods
should add sibling modules and implement the same Carlo-facing pattern. They
must not be inserted into the impurity-specific `spin_boson` module.

## Generic QMC foundation

`src/algorithm.rs` defines:

- `QmcKernel<C, R>`: one representation-specific update engine;
- `UpdateSchedule`: a fixed amount of work per measured sweep.

The schedule may adapt during thermalization, but is frozen after
thermalization. This prevents state-dependent sweep length from biasing
sweep-sampled observables.

## Spin-boson module

```text
spin_boson/
├── bath.rs           normalized frequency/time proposals
├── configuration.rs  retarded vertices and circular worldline index
├── error.rs          construction/runtime errors
├── mc.rs             Carlo.rs adapter and Params conventions
├── model.rs          JC/XXZ/XYZ/rotated-spin-boson catalogs
├── observables.rs    diagonal worldline estimators
├── scattering.rs     generic local detailed-balance scattering
├── updates.rs        diagonal and wormhole directed-loop updates
└── vertex.rs         four-leg vertex primitives
```

### Configuration

A sampled retarded vertex is

```text
(interaction, kind, omega, tau_a, tau_b)
```

and has four spin legs:

```text
A_in, A_out, B_in, B_out.
```

`WorldlineIndex` sorts both endpoints in imaginary time and creates the
periodic links joining an outgoing leg to the next incoming leg. Retarded
endpoint pairing is implicit in the two endpoint legs of the same vertex.

### Bath proposal

The diagonal update samples the normalized factor

```text
calJ(omega) P(omega, delta_tau),
P(omega,tau) = omega D(omega,tau).
```

It therefore cancels from the Metropolis ratio. Implemented shapes are:

- one mode;
- sharp-cutoff power law;
- arbitrary positive tabulated spectral mass proportional to `J(omega)/omega`.

Directed JC vertices retain the sampled orientation. Hermitian coordinate
couplings choose either orientation with probability one half.

### Model catalog

A model does not implement Monte Carlo control flow. It provides independent
interaction channels, positive local vertex kinds, diagonal seed lookup, and a
local scattering table.

Implemented impurity catalogs:

- Jaynes-Cummings: directed `S_+(tau_a) S_-(tau_b)`;
- XXZ: exchange flips plus longitudinal diagonal interaction;
- XYZ: exchange and pair flips;
- original spin-boson/single-mode Rabi after rotation into the bath-coupling
  basis.

The current engine is deliberately restricted to sign-free positive catalogs.
Spatially nonlocal propagators with negative matrix elements and genuinely
complex spectral matrices require a different sign/phase treatment.

### Updates

A sweep contains a fixed number of:

1. diagonal add/remove proposals;
2. closed directed loops.

Insertion draws an interaction, `tau_a`, `omega`, and `delta_tau`, then selects
the diagonal local kind compatible with the current worldline spins. The
acceptance ratios are

```text
R_add    = beta W_v / ((n_diag + 1) p_interaction)
R_remove = n_diag p_interaction / (beta W_v).
```

The loop starts on a random leg. At each vertex it flips the entrance and a
sampled exit leg. An exit on the other endpoint is a wormhole jump. The loop
closes when it returns to its starting discontinuity.

`ScatteringTable` builds symmetric path weights on the graph of compatible
extended local states. The default residual-flow solver strongly reduces
bounce while preserving exact local detailed balance for every positive
catalog; a symmetric-proposal Metropolis table remains available as a reference
and fallback. Analytic or linear-programmed minimum-bounce policies can be
added without changing the engine API.

### Observables

Implemented raw estimators:

- time-averaged `sigma_z` and `S_z`;
- second/fourth powers of the time-averaged magnetization;
- longitudinal static susceptibility sample;
- half-period longitudinal correlation;
- total, diagonal, and off-diagonal expansion orders;
- shifted interaction-expansion energy `-n/beta`;
- update acceptance, loop length, bounce, and wormhole diagnostics.

Off-diagonal improved loop estimators and bath reconstruction estimators are
separate follow-on estimator modules; they are not required for correctness of
the sampler.

## Parameter conventions

All explicit `lambda*` parameters are normalized retarded vertex couplings and
take precedence over inferred couplings.

For a power-law bath:

```text
lambda_l = 2 alpha_l omega_c / s.
```

For a single mode whose coupling multiplies `S_l`:

```text
lambda_l = g_l^2 / omega0.
```

For the common standard Rabi convention

```text
H_int = g_sigma sigma_z (a + a^dagger),
```

use `g_sigma`; the rotated spin coupling is `2 g_sigma`, hence

```text
lambda = 4 g_sigma^2 / omega0.
```

Using `g` for the Rabi model instead means that `g` multiplies `S_z` before
rotation.

## Validation strategy

The source includes tests for:

- bath sample support;
- vertex-catalog completeness;
- scattering-row normalization;
- local detailed balance;
- circular worldline continuity and link involution;
- mixed diagonal/wormhole updates;
- construction and execution of every impurity catalog through Carlo.rs.

Production validation should additionally compare single-mode models against
exact diagonalization and benchmark integrated autocorrelation time against the
specialized Winter cluster backend.
