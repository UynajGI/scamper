# QMC.rs

QMC.rs is the physics layer of the Scuttle Monte Carlo workspace. Carlo.rs
provides scheduling, deterministic RNG setup, thermalization/measurement
phases, accumulation, parallel backends, and result analysis. QMC.rs provides
representations, model catalogs, update kernels, invariants, and estimators.

The module boundary and extension rules are documented in
[`ARCHITECTURE.md`](ARCHITECTURE.md).

## Continuous-time spin-boson wormhole QMC

`qmc_rs::spin_boson` implements a generic retarded-interaction directed-loop
engine for one spin-1/2 impurity. Quadratic bosons are integrated out and
represented by two-time four-leg vertices. Diagonal updates sample
`(interaction, omega, tau, tau')`; directed loops convert diagonal and
spin-flip vertices and can traverse the nonlocal endpoint connection—the
wormhole move.

### Implemented impurity models

Set parameter `model` to:

- `jc` / `jaynes_cummings`;
- `xxz`;
- `xyz`;
- `rabi` / `spin_boson` / `rotated_spin_boson`.

### Implemented bath proposals

Set parameter `bath` to:

- `single`;
- `powerlaw`;
- `tabulated`.

### Carlo.rs usage

```rust
use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use qmc_rs::SpinBosonQmc;

let mut params = Params::new();
params.set("model", "jc");
params.set("bath", "single");
params.set("beta", 8.0);
params.set("omega0", 1.0);
params.set("g", 0.35);
params.set("h_z", 0.4);

let results = Scheduler::new(RayonBackend::new(1), RunConfig::default())
    .run_one::<SpinBosonQmc>(&params);
```

A complete runnable example is in
[`examples/spin_boson_wormhole.rs`](examples/spin_boson_wormhole.rs).

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
model                   jc | xxz | xyz | rabi
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
cargo run -p qmc-rs --example spin_boson_wormhole
```
