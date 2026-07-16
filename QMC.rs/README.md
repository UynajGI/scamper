# QMC.rs

QMC.rs is the physics layer of the Scuttle Monte Carlo workspace. Carlo.rs
provides scheduling, deterministic RNG setup, thermalization/measurement
phases, accumulation, parallel backends, and result analysis. QMC.rs provides
representations, model catalogs, update kernels, invariants, and estimators.

The module boundary and extension rules are documented in the workspace
`CLAUDE.md`.

## Continuous-time lattice QMC

`qmc_rs::lattice` is a general continuous-time interaction-expansion
directed-loop engine for sign-problem-free quantum spin models on arbitrary
CSR adjacency graphs with arbitrary quantum spin `S`.

### Architecture

```text
CsrGraph → LocalHilbertSpace → PositiveOperatorModel (sparse K=C-H catalog)
  → LatticeConfiguration + WorldlineIndex
  → diagonal add/remove + low-bounce directed loops
  → estimators → LatticeSpinQmc (Carlo.rs adapter)
```

### Implemented models

Set parameter `model` to:

- `heisenberg` — isotropic XXZ with `J`
- `xy` — `J_xy` only
- `xxz` — `J_xy` + `J_z`
- `xyz` — `J_x` + `J_y` + `J_z`
- `tfim` — transverse-field Ising (`J_z` + `h_x`)

### Topologies

Set parameter `topology` to: `chain`, `square`, `hypercubic`, `edges` (edge list), `adjacency`.

### Carlo.rs usage

```rust
use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use qmc_rs::LatticeSpinQmc;

let mut params = Params::new();
params.set("beta", 8.0);
params.set("model", "heisenberg");
params.set("topology", "square");
params.set("Lx", 4);
params.set("Ly", 4);
params.set("spin", 0.5);
params.set("J", 1.0);

let results = Scheduler::new(RayonBackend::new(1), RunConfig::default())
    .run_one::<LatticeSpinQmc>(&params);
```

A complete runnable example is in
[`examples/lattice_continuous.rs`](examples/lattice_continuous.rs).

### Common parameters

```text
beta                    required inverse temperature
model                   heisenberg | xy | xxz | xyz | tfim
topology                chain | square | hypercubic | edges | adjacency
spin / two_s            quantum spin magnitude
J / J_xy / J_z / J_x…   exchange couplings
h_x / h_z               transverse/longitudinal fields
D                       single-ion anisotropy
gauge                   auto (Marshall Z2) | identity
scattering              low_bounce | metropolis
adaptive_schedule       warmup-only work adaptation
diagonal_proposals      per measured sweep
directed_loops          per measured sweep
```

### Validated domain and known limitations

The lattice solver has been validated against exact diagonalization and
analytic limits for **S = 1/2** sign-problem-free models on bipartite
graphs (Heisenberg, XXZ, XYZ, TFIM).

**Spin-S > 1/2 caveat:** the directed-loop bounce fallback (used when
the scattering table has no residual flow for a given entrance channel)
preserves the worm charge by reflecting it. For S = 1/2 this is always
a legal spin-flip, but for S > 1/2 it can create an illegal
raising/lowering operation at the linked leg. Until per-level scattering
tables replace the bounce fallback, **S > 1/2 results should not be
trusted** without independent verification.

**Sign-problem-free requirement:** the solver auto-detects and rejects
frustrated (non-stoquastic) models via the Marshall Z2 gauge solver.
Fermionic statistics are reserved but not implemented.

## Continuous-time impurity wormhole QMC

`qmc_rs::impurity` implements a generic retarded-interaction directed-loop
engine for one spin-1/2 impurity. Quadratic bosons are integrated out and
represented by two-time four-leg vertices. Diagonal updates sample
`(interaction, omega, tau, tau')`; directed loops convert diagonal and
spin-flip vertices and can traverse the nonlocal endpoint connection—the
wormhole move.

### Basis rotation for the Rabi model

To make the retarded-interaction matrix elements sign-free, the wormhole
solver samples in a **rotated basis** for the `rabi`/`rotated_impurity`
model: `σz_sampled = σx_physical`, `σx_sampled = -σz_physical`, `σy` fixed
(see `BasisTransform::rotated_rabi()`). The observable `MagnetizationSigmaZ`
reported by the wormhole is therefore the **physical ⟨σx⟩**, not ⟨σz⟩.

The occupation solver (`OccupationWorldlineQmc`) does **not** rotate: it
samples in the physical σz basis. Its `OccupationSigmaZ` is the physical
⟨σz⟩, and `OccupationSigmaX` is computed from the transfer matrix.

Direct observable comparison between the two solvers requires accounting
for this basis difference. See `tests/impurity/cross_solver.rs` for the
full convention reconciliation notes.

### Implemented impurity models

Set parameter `model` to:

- `jc` / `jaynes_cummings`;
- `rw_crw` / `weber`;
- `xxz`;
- `xyz`;
- `rabi` / `impurity` / `rotated_impurity`.

### Implemented bath proposals

Set parameter `bath` to:

- `single`;
- `powerlaw`;
- `tabulated`.

### Carlo.rs usage

```rust
use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use qmc_rs::ImpurityQmc;

let mut params = Params::new();
params.set("model", "jc");
params.set("bath", "single");
params.set("beta", 8.0);
params.set("omega0", 1.0);
params.set("g", 0.35);
params.set("h_z", 0.4);

let results = Scheduler::new(RayonBackend::new(1), RunConfig::default())
    .run_one::<ImpurityQmc>(&params);
```

A complete runnable example is in
[`examples/impurity_wormhole.rs`](examples/impurity_wormhole.rs).

## Coupling conventions

An explicit `lambda`, `lambda_xy`, `lambda_x`, `lambda_y`, or `lambda_z`
always takes precedence.

For the sharp-cutoff power-law convention:

```text
lambda_l = 2 alpha_l omega_c / s.
```

For one mode with `g_l S_l (a + a^dagger)`:

```text
lambda_l = g_l^2 / omega0.
```

For standard Rabi coupling `g_sigma sigma_z (a + a^dagger)`, use
`g_sigma`, giving

```text
lambda = 4 g_sigma^2 / omega0.
```

Using `g` in the Rabi model instead means that `g` multiplies `S_z`.

The JC retarded flip weight is `lambda`; the XXZ exchange-vertex weight is
`lambda_xy / 2`; XYZ exchange and pair-flip weights are respectively
`(lambda_x + lambda_y)/4` and `abs(lambda_x - lambda_y)/4`.

## Common parameters

```text
beta                    required inverse temperature
model                   jc | rw_crw | xxz | xyz | rabi
bath                    single | powerlaw | tabulated
C                       optional positive diagonal shift
h_z                     diagonal field in the sampled basis
diagonal_proposals      fixed proposals per measured sweep
directed_loops          fixed closed loops per measured sweep
max_loop_steps_factor   loop safety limit multiplier
adaptive_schedule       warmup-only work adaptation
adaptation_interval     warmup adaptation window
correlation_samples     random origins for C_z(beta/2)
validate_each_sweep     expensive invariant checks
```

Bath-specific parameters:

```text
single:     omega0
powerlaw:   s, omega_c
tabulated:  bath_omegas="...", bath_weights="...", explicit lambda*
```

## Verification commands

With a Rust toolchain installed at the workspace root:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --all-features
cargo test --workspace --all-features
cargo run -p qmc-rs --example lattice_continuous
cargo run -p qmc-rs --example impurity_wormhole
```
