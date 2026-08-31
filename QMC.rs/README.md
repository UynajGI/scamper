# QMC.rs

QMC.rs is the physics layer of the Scamper Monte Carlo workspace. Carlo.rs
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

## Variational family (L0–L3)

Continuum variational QMC lives in the `variational/` module family
(`qmc_rs::variational`, flat re-exports at the crate root), hosted by the
same Carlo.rs `sweep`/`measure` contract as the lattice solvers: the walker
population is solver-internal state and one sweep is one epoch of
single-particle Metropolis over all walkers, on per-walker RNG streams
derived through `RngStreamKey` (walker index in the replica field).

What exists (layer L0):

- `WaveFunction` trait with hand-derived analytic `∇ ln|ψ|`, `∇² ln|ψ|`
  and parameter gradients, plus allocation-free incremental
  `delta_log`/`commit_move` single-particle updates (flat particle-
  interleaved 3-D configurations, open box — no PBC/minimum image yet).
- Three ansätze: `GaussianTrap` (one-body harmonic-trap Gaussian),
  `McMillanJastrow` (`Π exp(−½(b/r)⁵)`) and `HarmonicJastrow`
  (`Π exp(−a r_ij²)`), combinable via the log-additive `Product`
  combinator (the L1 Slater–Jastrow seam).
- `ContinuumHamiltonian` (harmonic trap + Lennard-Jones or harmonic pair
  terms), the `local_energy` estimator, and `VmcKernel<W>` — a
  Metropolis-in-`|ψ_T|²` population kernel implementing `MonteCarlo`
  (drive it with `Run::from_parts`; `FromParams` is deliberately not
  implemented for the generic ansatz type). Versioned JSON checkpoints
  (`qmc-rs-vmc-v1`) cover walkers + parameters and reject unknown or
  corrupted snapshots loudly.

What exists (layer L1, fermionic family behind the same trait):

- `SlaterDeterminant`: two spin blocks of contracted cartesian Gaussians
  (exact HO shells 0..=2 via `harmonic_trap(ω, n_shells)` — 2/8/20
  electrons), dense LU via `nalgebra` (mature-crate policy), hand-derived
  `Tr(D⁻¹∇D)` gradient rows and Hessian-chain Laplacians.
- Single-particle fast path: the Sherman–Morrison column identity gives
  the Metropolis ratio as one O(N) dot product against the cached inverse,
  with an O(N²) rank-1 inverse update on accept and per-sweep `rebuild`
  re-anchoring (the K-rebuild policy).
- `Backflow` (Kwon–Ceperley–Martin quasiparticle displacement,
  electron-gas shape preset): full Jacobian/Laplacian chains; `λ = 0`
  reproduces the plain determinant bit-exactly. With backflow active
  `delta_log` is a full recompute, as the literature prescribes.
- Slater–Jastrow composes through the unmodified `Product` combinator.

Validated domain (exact statements, see `tests/variational/`):
`GaussianTrap` at `α = ω/2` gives `E_L ≡ 3Nω/2` to machine precision (zero
variance) through the full Metropolis pipeline; `HarmonicJastrow` is the
exact ground state (`E₀ = 3aN(N−1)`) of the pair-harmonic trap
`k = 4a²N`; the exact-HO-shell `SlaterDeterminant` is zero-variance at
every closed shell (E₀ = 3ω/18ω/60ω for 2/8/20 electrons) through the LU
path and through the kernel; Sherman–Morrison updates match fresh LU
(ratios and log-dets ≤ 1e-12, entrywise inverse at the conditioning
floor); `delta_log` matches full recomputes at the scale-aware machine
floor `16ε(1+|ln|ψ||)`; every gradient (incl. GTO exponents/coefficients
and the backflow scale) agrees with central finite differences; `λ = 0`
backflow is bit-exact; same-seed runs are bit-identical; the He-4-like
confined McMillan droplet respects the rigorous bound
`E ≥ −ε·N(N−1)/2` and the Slater–Jastrow fermion droplet respects
`E ≥ E₀(H₀)` under a repulsive pair, both with multi-seed z-score
consistency.

What exists (layer L2, parameter optimizers as an outer loop):

- `BlockStats`: block-moment accumulator (energy, variance, centered
  force `S_k = ⟨Ȯ_k E_L⟩`, SR metric `G_kl = ⟨Ȯ_k Ȯ_l⟩`, three-point
  moment `T_kl` for the linear method), fed per sweep by
  `VmcKernel::collect_block_stats`; importance-weighted pushes make
  deterministic quadrature statistics possible for theorem-level tests.
- `StochasticReconfiguration` (Sorella natural gradient,
  `(G + λ·diag G) Δp = −ε S` on `nalgebra` Cholesky) and `LinearMethod`
  (Umrigar–Nightingale generalized eigenproblem on the linearized
  displacement basis, symmetric reduction via Cholesky +
  `SymmetricEigen`), both with a diagonal trust-region shift
  (escalate-and-retry on rejected steps) and patience-based convergence
  on the natural force norm; `VmcKernel::update_wave_function_params`
  applies/reverts parameter updates with walker re-anchoring and
  rollback on out-of-domain steps.
- Validated: SR and LM converge to the closed-form optimum
  `α* = ω/2` of the one-particle Gaussian trap against quadrature-exact
  statistics; the LM predicted energy never exceeds the current block
  energy (subspace convexity); the force vanishes for an exact state
  (zero-variance principle on the statistics layer), deterministically
  and through the Metropolis kernel; trust-region escalation/relaxation
  and input rejection; SR improves the two-parameter
  Gaussian×McMillan LJ droplet beyond noise with physical parameters
  throughout.
- `VarianceMinimization` (third L2 entry point): Umrigar-style
  correlated-sampling variance minimization. Reference configurations
  carry their sampling-measure density (`ReferenceSample`), the
  reweighted two-pass local-energy variance is the objective, and the
  search is `argmin`'s Nelder–Mead simplex (mature-crate policy; the
  objective's parameter gradient would need third-order chains the
  `WaveFunction` trait deliberately does not carry, hence derivative-free).
  Out-of-domain candidates cost a finite penalty the simplex walks
  around. Validated: the objective reproduces the closed form
  `Var(α) = c(α)²·3/(8α²)` on uniform-grid samples (the importance base
  weight makes the reweighting exact quadrature at any candidate) and
  vanishes at the exact state; Nelder–Mead converges to `α* = ω/2`
  deterministically and, on kernel-sampled configurations from a poor
  start, reduces the variance by more than an order of magnitude.

What exists (layer L3, diffusion Monte Carlo):

- `DmcKernel`: the walker population as solver-internal state (one
  `step()` = one imaginary-time step over the whole population).
  Drift-diffusion moves `R' = R + τb + √τ·χ` with drift `b = ∇ln|ψ_T|`,
  Metropolis-accepted against the exact Green-function ratio
  `ln A = 2Δln|ψ| + (|D_fwd|² − |D_bwd|²)/2τ` (backward displacement
  `R−R'−τb(R')` — sign load-bearing, see VALIDATION.md); branching
  `⌊g+u⌋` with `g = exp(−τ(E_L−E_T))` and a population-safety cap;
  classic `E_T = E_ref − ln(N/N_target)/τ` feedback on an EMA reference
  energy; forward-walking pure estimators through per-walker lineage
  rings and a pending-measurement ring (a measurement is credited once
  per surviving descendant line — descendant weighting by construction).
  Versioned checkpoints (`qmc-rs-dmc-v1`) serialize the full state
  including lineage for bit-identical replay.
- Validated: DMC converges to the **derived** exact ground state
  `E₀ = (3/2)(ω + √(ω²+2k))` of the two-particle trap-plus-repulsive-pair
  system from a deliberately approximate nodeless trial state (CM/relative
  separation, fully derived in the test), strictly closer than VMC at the
  same ψ_T; for the exact ψ_T the entire machinery reduces to an identity
  (mixed = pure = E₀ at 1e-10); the population-control bias shrinks from
  N_w = 8 to 64 at equal walker-step budgets; same-seed bit-identical
  runs, checkpoint round-trip/replay, and loud rejection of corrupt or
  mismatched snapshots and invalid constructor inputs.
- Validated domain: **nodeless (bosonic) trial states**. The
  fixed-node constraint for fermionic determinants — sign tracking plus
  node-crossing rejection, which needs a signed-amplitude extension of
  the `WaveFunction` trait — is the named follow-up.

Not there yet: fixed-node enforcement, reptation (L4), NQS/t-VMC (L5,
mature-crate autodiff decision).
Architecture and layer plan: [`research/vmc/DESIGN.md`](../research/vmc/DESIGN.md).

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

## Roadmap — algorithm families not yet included

Current coverage is four solver families (lattice directed-loop, impurity
wormhole, cavity-QED occupation, longitudinal cluster), all living in the
sign-free positive-weight world of spins and bosonic occupation; fermionic
statistics are reserved at the `LocalHilbertSpace` boundary and rejected
loudly. Everything below is **not implemented and not validated**; it is
recorded so the validated domain stays unambiguous. None of it blocks the
production status of the existing solvers.

### Determinantal / fermionic family (the largest hole)

- **DQMC** (discrete-time auxiliary field), **CT-INT** (interaction
  expansion), **CT-AUX** (auxiliary-field continuous time), **CT-HYB**
  (hybridization expansion — the standard DMFT impurity solver; the existing
  wormhole handles spin-boson retarded self-interactions, not Anderson
  impurities), **Hirsch–Fye**. With **average-sign measurement** and
  parallel tempering in determinant space. This family decides whether
  correlated-electron physics is reachable at all.

### Bosonic family

- **Lattice bosons** (Bose–Hubbard worldline + worm) — the lattice side has
  `SpinSpace` only; occupation bases exist only inside the cavity-QED
  impurity solver.
- **Continuous-space PIMC** (permutation cycles), **superfluid density /
  one-body density matrix** worm estimators, **PIGS**.

### Discrete-time SSE and improved estimators

- A **discrete SSE** sibling of the continuous-time directed-loop (reuses
  the K=C−H operator catalog and scattering tables) carrying **loop
  improved estimators** (zero-variance χ, C(τ)) and the **dynamic structure
  factor S(q,ω)** spectral representation. The lattice solver currently
  measures only nearest-neighbor Sz correlations (the long-standing
  QMC-P2.2 deferral).

### Variational family (inside QMC.rs, `variational/` module family) (L0 implemented 2026-08-19)

- **VMC** (Jastrow / Slater–Jastrow / backflow / Pfaffian trial states) with
  the **optimization machinery** — stochastic reconfiguration, natural
  gradient, linear method — that is half the method; **NQS** (neural quantum
  states) as modern ansätze; **t-VMC** for real-time dynamics. Optimizer and
  autodiff machinery adopted from mature crates (`nalgebra`, `argmin`,
  candle-class AD at the NQS layer); custom code restricted to physics.
- Architectural note: the family fits the existing `sweep`/`measure`
  contract — the walker population lives as solver-internal state (one
  sweep = one imaginary-time step over all walkers; per-walker RNG streams
  derive from the existing `RngStreamKey` domain separation). What is
  missing is conveniences, not permission: weighted-observable conventions,
  descendant-tracking measurement windows for pure estimators, and
  intra-population parallel dispatch that composes with chain-level
  backends. The genuinely new machinery is wavefunction-gradient
  evaluation and the optimizer outer loop.

### Configuration-space projection family

- **DMC** (fixed-node), **GFMC**, **reptation QMC**, **AFQMC** (phaseless —
  straddles auxiliary-field sampling and a variational trial-state
  constraint). Complementary to the worldline solvers: ground-state focus,
  continuum and lattice, fermions via fixed-node rather than the sign
  problem.

### Sign-problem machinery beyond Marshall

- **Majorana-representation QMC**, **fermion-bag**, **PT-symmetric basis
  optimization** — the means to turn currently-rejected frustrated models
  into computable ones.

### Diagrammatic family

- **DiagMC** (bold/skeleton expansions), **CDet** (Fermi polaron),
  electron-phonon **CT-INT**.

### Non-equilibrium and impurity extensions

- **Keldysh real-time QMC**, **multi-impurity retarded interactions**
  (RKKY), **non-Gaussian/anharmonic baths** — `Bath` today is
  single-mode/power-law/tabulated, all Gaussian.

### Suggested priority (by reuse)

Discrete SSE + improved estimators → lattice-boson worldline + worm →
determinantal family (the `ParticleStatistics::Fermion` boundary and its
rejection logic were designed for exactly this extension) → AFQMC and
real-time. The variational family lives inside QMC.rs as the
`variational/` module family (see `research/vmc/DESIGN.md`); its
wavefunction-gradient machinery is new ground but shares the Carlo.rs
hosting and RNG-stream conventions.
