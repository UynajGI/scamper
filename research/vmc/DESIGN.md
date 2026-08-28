# Variational & Projector QMC Family: Architecture Design

> Status: design for the variational family **inside QMC.rs** (package
> `qmc-rs`, module family `QMC.rs/src/variational/`) — not a separate crate.
> Grounded in the corrected hosting analysis: walker populations live as
> solver-internal state inside the Carlo.rs `sweep`/`measure` contract
> (one sweep = one epoch / one imaginary-time step over the population).
> Learning corpus: `research/vmc/papers/` (classic papers, MinerU-converted).

## 1. Design goals

1. **Low technical debt** — layer seams chosen so no later feature (DMC, NQS)
   forces a rewrite of earlier layers. The `WaveFunction` trait is the
   firewall: all physics in ansätze, all statistics in kernels, all
   adaptation in optimizers. No layer reaches across.
2. **High performance** — zero heap allocation on the hot path, structure-of-
   arrays configurations, incremental O(N)/O(N²) local updates instead of
   O(N²)/O(N³) full re-evaluations, static dispatch (no dyn in the inner
   loop), reproducible per-walker RNG streams enabling future intra-
   population parallelism.
3. **Validation-first culture** — every layer ships with the house physics
   standard (machine-precision identities, exact references, z-scores,
   variational bounds; the 8 maturity criteria from `MATURITY_ASSESSMENT.md`).
4. **Carlo.rs hosting from day 1** — VMC/DMC kernels implement
   `MonteCarlo` with populations as internal state; `RngStreamKey`
   domain-separated streams per walker.

## 2. Layer plan (low → high)

| Layer | Content | Key acceptance criteria |
|-------|---------|------------------------|
| **L0** | Continuum VMC: McMillan-style Jastrow ψ_T, analytic gradients, Metropolis walker population | HO with Gaussian ψ_T: E_L exact constant at machine precision (zero variance); virial identity; Δlnψ ≡ full recompute bit-exact; E_VMC ≥ E0 statistically; He (literature window) |
| **L1** | Slater–Jastrow + backflow; determinant machinery | Sherman-Morrison rank-1 update ≡ full inverse (1e-12); non-interacting closed-shell kinetic exact; param-gradient vs central finite differences (1e-6) |
| **L2** | Optimizers: stochastic reconfiguration / natural gradient, linear method, variance minimization (argmin, correlated sampling) | SR reproduces the analytic optimum on the deterministic Gaussian toy (quadrature-exact block statistics); linear method converges on the same toy and its predicted energy ≤ current block energy (subspace convexity theorem) — step-count dominance over SR is a stiff multi-parameter phenomenon and is NOT gated on a one-parameter toy; trust-region monotonicity (reject/escalate, keep-parameters); force ≡ 0 for an exact state; variance minimization's reweighted objective ≡ closed form and ≡ 0 at the exact state (importance base weight), Nelder–Mead converges deterministically and ≥10× variance reduction on kernel-sampled configurations |
| **L3** | DMC: drift-diffuse + branching + population control; descendant-weighted pure estimators (nodeless implemented; fixed-node enforcement — signed amplitudes + node-crossing rejection — is the named follow-up) | Derived separable anchor `E₀ = (3/2)(ω + √(ω²+2k))` (2-particle trap + repulsive pair): nodeless DMC converges to it from an approximate ψ_T and beats VMC at the same ψ_T; exact-ψ_T identity (mixed = pure = E₀ machine-precision); population-control bias shrinks with N_w at equal budgets; bit-identical replay through lineage-serializing checkpoints |
| **L4** | Reptation QMC (single-polymer path, no population bias) | Free-particle path sampling exact; barrier crossing time consistent with action |
| **L5** | NQS ansätze + t-VMC (autodiff decision point) | RBM on 1D TFIM vs exact ED (Carleo–Troyer benchmark); t-VMC conserves norm/energy to tolerance |

L0–L3 use **hand-derived deterministic adjoints only** — no autodiff
dependency. The trait contract keeps the seam clean so an AD-backed
implementation can slot in at L5 without touching anything below.

## 3. Core traits

### `WaveFunction` — the anti-debt hinge

```rust
pub trait WaveFunction {
    type Config;                       // SoA positions, e.g. flat [f64; 3N]
    fn log_psi(&self, cfg: &Self::Config) -> f64;                 // ln|ψ_T|
    fn log_grad(&self, cfg: &Self::Config, out: &mut GradBuffer); // ∇_i ln|ψ| per particle
    fn log_laplacian(&self, cfg: &Self::Config) -> f64;           // Σ_i ∇²_i ln|ψ|
    fn n_params(&self) -> usize;
    fn log_grad_params(&self, cfg: &Self::Config, out: &mut ParamGradBuffer); // ∂ ln|ψ|/∂p
    fn update_params(&mut self, delta: &[f64]);
    // THE performance hinge — incremental single-particle updates:
    fn delta_log(&self, cfg: &Self::Config, particle: usize, new_pos: &Point)
        -> DeltaLog;                   // O(N) Jastrow / O(N²) rank-1 determinant
    fn commit_move(&mut self, cfg: &mut Self::Config, particle: usize, new_pos: &Point);
    fn rebuild(&mut self, cfg: &Self::Config);  // periodic O(N³) determinant refresh
}
```

- `delta_log` vs full recompute equivalence is a **day-1 machine-precision
  test** — this invariant is what keeps fast paths honest.
- L1 determinants: Sherman-Morrison rank-1 updates with a rebuild every K
  accepted moves (standard practice, guarded by the equivalence test).

### Kernels — populations as solver-internal state

```rust
pub struct VmcKernel<W: WaveFunction> {
    walkers: Vec<Walker>,              // Walker { cfg, log_psi, rng salt }
    wf: W,
    scratch: ScratchBuffers,           // reused, zero-alloc hot path
}
```

- `sweep()` = one epoch of single-particle Metropolis over all walkers in
  deterministic order; `measure()` accumulates E_L, E_L², |∇lnψ|² per walker.
- **Weighted-observable convention (documented once, no ad-hoc sites):**
  weighted quantities are accumulated as `(x·w, w)` pairs in separate
  accumulators; the ratio is formed in postprocessing. Nothing touches
  Carlo.rs.
- RNG: one stream per walker derived via `RngStreamKey` (walker index in the
  replica field) — thread-count independence by construction, future rayon
  fan-out inside `sweep` changes no results.

### `DmcKernel` (L3)

- `sweep()` = one τ-step: drift-diffusion move of every walker (Metropolis
  against the Green-function ratio), branching birth/death from
  `exp(-τ(E_L + E_T))`, population control rescale, E_T feedback.
- Pure estimators via an internal **descendant ledger**: measurements are
  flushed `N_delay` sweeps later when descendants are known — delayed
  `measure` accumulation inside the kernel, checkpoints serialize the
  ledger (versioned `vmc-rs-dmc-v1` tag, loud rejects — lesson from the CMC
  worm v1/v2 history).

### `Optimizer` — an outer loop, never inside sweep

```rust
pub trait Optimizer { fn step(&mut self, stats: &BlockStats) -> ParamUpdate; }
```

- **Mature-crate policy (user directive 2026-08-19): ML machinery is
  adopted from established Rust libraries, never hand-rolled.** Custom code
  is restricted to the physics quantities — E_L statistics, SR force
  vectors `S`, covariance matrices `G` — while the numerical machinery
  comes from: `nalgebra` (already a workspace dep) for the SR solve
  `G Δp = -S` (Cholesky with λ-regularization schedule) and the
  linear-method overlap-matrix eigenproblem; `argmin` for generic
  optimization scaffolding (line search, trust region, CG/L-BFGS) where a
  step is posed as plain optimization (e.g. variance minimization).
- L5 autodiff: adopt a mature AD framework (`candle` / `burn` / `dfdx` —
  decide at L5 against: reverse-mode with double-backward for ∇²
  ln|ψ|, serialization of ansatz weights, maintenance health); never
  hand-build AD. Plain optimizers (Adam etc.) ship with those frameworks.
- Integration via phase transitions (measure block → optimize → next block),
  following the `AdaptiveRunControl` precedent; optimizer state is
  checkpointed alongside walkers. Optimizers run between blocks only —
  zero cost on the hot path.

## 4. Module layout

```
QMC.rs/src/variational/
    wavefunction/   mod.rs, jastrow.rs (L0: McMillan/Pade), determinant.rs (L1),
                    backflow.rs (L1), nqs.rs (L5, later)
    kernel/         vmc.rs (L0), dmc.rs (L3), reptation.rs (L4)
    optimizer/      mod.rs, sr.rs, linear_method.rs, variance.rs (L2)
    estimators/     local_energy.rs, pure.rs (descendant weighting)
    adapters/       carlo.rs (MonteCarlo impls, FromParams)
QMC.rs/tests/variational/   house-style physics tests, input_validation.rs
```

No new workspace crate. The family lives in `qmc-rs` and reuses its deps:
`carlo-rs`, `rand`, `thiserror`, `serde`; `nalgebra` (already in
`[workspace.dependencies]`) is added to `qmc-rs` when L2 optimizers land.
A `criterion` bench target joins QMC.rs when the hot-path budget needs it.
Workspace lints already inherit (`clippy::all = deny`, `unsafe_code = deny`).
Checkpoint tags follow the crate convention: `qmc-rs-vmc-v1`, later
`qmc-rs-dmc-v1`.

## 5. Debt-avoidance decisions (explicit)

1. No autodiff before L5 — and at L5 it is **adopted from a mature AD
   framework**, not built. L0–L4 hand adjoints are physics (∇ ln|ψ| and
   ∇² ln|ψ| are required for E_L regardless), verified against finite
   differences in CI. The optimizer/AD crate boundary sits at the
   `Optimizer` trait and the `WaveFunction` trait respectively — swapping
   crates later touches nothing below.
2. Static dispatch everywhere on the hot path; no `dyn WaveFunction` in the
   inner loop.
3. Weighted-observable convention written once in README + enforced by a
   lint-style test naming convention.
4. Checkpoint version tags from day 1; unknown/corrupt snapshots reject
   loudly.
5. Input validation (criterion G) from L0: invalid walker counts, τ ≤ 0,
   non-finite params, empty populations → `Err`, never panic.
6. Benchmarks before optimization claims: criterion harness guards the
   zero-alloc/O(N) hot-path budgets.

## 6. Canonical references

McMillan (1965) Jastrow He; Ceperley–Chester–Kalos (1977); Anderson (1975)
DMC; Ceperley–Alder (1980); Reynolds et al (1982) fixed-node; Umrigar–
Nightingale–Runge (1993) pure estimators & population control; Umrigar et al
(1988) variance minimization; Sorella (1998) stochastic reconfiguration;
Umrigar–Nightingale (2001) linear method; Baroni–Petrosyan (2002) /
Boninsegni et al (2005) reptation; Foulkes et al RMP 73, 33 (2001) review;
Becca–Sorella (2017) book; Carleo–Troyer (2017) NQS; Carleo et al (2012–14)
t-VMC. Machine-readable copies under `research/vmc/papers/` (see its
MANIFEST.md for what is available on arXiv vs paywalled classics).
