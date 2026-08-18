# CMC.rs — Physics Validation & Validated Domain

> Updated 2026-08-14. Branch: `dev`.

## Test suite summary

| Layer | Tests | Runtime |
|-------|-------|---------|
| Default (`cargo test`) | 230 | ~30s |
| Long stochastic (`--ignored`) | 15 | ~120s |
| **Suite total** | **245** | (+69 lib unit tests) |

## Per-solver validated domain

### Local Metropolis (`MetropolisCore`)
- **Validated:** Ising S=1/2 on chain/square/hypercubic graphs, PBC and OBC
- **Observables:** ⟨E⟩, ⟨m²⟩ vs exact enumeration (N=2,3,4,8)
- **Detailed balance:** Direct DB verified on N=2 (1e-15 precision)
- **Ergodicity:** 16-seed z-score framework (|z|<4, no systematic bias)
- **Connectivity:** Explicit Markov chain strongly connected on N=2
- **Cross-solver:** Agrees with Wolff and SW within 3σ on 8×8 Ising; continuous-spin cross-solver XY (O(2)) and Heisenberg (O(3)) vs Wolff on 8×8 (pooled |z|<4 at β=0.7/0.9 and 0.5/0.9)
- **Frustrated:** AFM XY triangle (J<0): full equilibrium vs spectral quadrature reference — ⟨E⟩ and chirality ⟨κ²⟩ at β=1,3 (per-seed |z|<3, ground-state algebra exact)
- **NOT validated:** — (continuous spins and frustrated systems validated 2026-08-14)

### Wolff cluster (`WolffCore`)
- **Validated:** Ising S=1/2, XY (O(2)), Heisenberg (O(3)) on chain/square
- **Observables:** ⟨E⟩ vs exact enumeration, Langevin/Bessel-ratio analytic limits
- **Detailed balance:** Direct DB verified on N=3
- **Ergodicity:** 16-seed z-score, strongly connected on N=2
- **NOT validated:** Potts q>2 cluster updates

### Swendsen-Wang (`SWCore`)
- **Validated:** Ising S=1/2 on chain/square
- **Detailed balance:** Direct DB verified on N=2 (50k samples/state)
- **Ergodicity:** Strongly connected on N=2
- **Cross-solver:** Agrees with Metropolis within 3σ on 8×8
- **NOT validated:** Potts q>2, continuous spins

### Heat bath — discrete (`HeatBathCore`)
- **Validated:** Ising S=1/2
- **Detailed balance:** Direct DB verified on N=2
- **Exact energies:** ⟨E⟩ and ⟨m²⟩ vs exact enumeration on N=4 and N=8 at β=0.3/0.8 (16-seed z-scores, max|z|<2.5)
- **NOT validated:** Potts q>2

### Heat bath — continuous O(N) (`ContinuousHeatBathCore`)
- **Validated:** O(3) uniform-on-sphere at β→0 (⟨s_α⟩≈0, ⟨s_α²⟩≈1/3)
- **Finite-T conditional distribution:** single-spin ⟨cosθ⟩/⟨cos²θ⟩ vs exact O(2) Bessel ratio I₁(x)/I₀(x) and O(3) Langevin coth(x)−1/x at x ∈ [0.2, 5]; analytic anchors (limits at small/large x, Bessel recurrence) to 1e-4–1e-12; two-spin XY pair equilibrium matches the conditional moment (rotational invariance)
- **NOT validated:** — (finite-T distribution and XY validated 2026-08-14; lattice-scale O(N) cross-solver via Metropolis-vs-Wolff above)

### Microcanonical over-relaxation (`MicrocanonicalCore`)
- **Validated:** Energy conservation to 2e-10, unit norm preservation to 3e-12
- **NOT validated:** Ergodicity (over-relaxation alone is not ergodic)

### Hybrid (`HybridCore<A, B>`)
- **Validated:** Hybrid(Metropolis, Wolff) ⟨E⟩ and ⟨m²⟩ vs exact enumeration on 4-site Ising
- **NOT validated:** Other combinations, continuous spins

### Wang-Landau (`WangLandauCore`)
- **Validated:** DOS matches exact 4×4 Ising enumeration (un-ignored, ~11s)
- **Validated:** 2-site exact DOS levels/degeneracies
- **Validated:** Canonical reweighting recovers exact ⟨E⟩
- **BinnedAxis production run:** binned DOS vs exact binned degeneracies (8-level weighted 6-ring, 14-bin axis; per-bin |Δln g| RMS 0.013, gate 0.05); canonical reweighting ⟨E⟩(T) at 3 temperatures vs exact (|z|≤1.2); flat-histogram-only route error floor documented (RMS ≤0.1)
- **NOT validated:** — (continuous-axis production run validated 2026-08-14)

### Multicanonical / umbrella (`EnergyBiasCore`)
- **Validated:** Transactional rejection, bias algebra, axis boundaries
- **Validated:** Discrete DOS reweighting matches enumeration
- **Full MC-vs-exact:** real EnergyBiasCore sampling with exact-DOS bias on the 6-site ring — ⟨E⟩(T), ⟨m²⟩(T) at β ∈ {0.2, 0.5, 1.0} (max|z|=3.4 over 16 seeds, |z̄|≤0.65), full P(E) at β=0.5 within 4σ binomial band (TV distance 0.005 vs gate 0.02), biased histogram flatness ≥0.90 (canonical would give 0.003)
- **NOT validated:** — (full MC distribution validated 2026-08-14; experimental status lifted)

### Worm (Ising HT graph) (`WormKernel`)
- **Validated:** HT graph partition identity matches spin enumeration
- **Validated:** Graph energy estimator matches exact spin energy
- **Validated:** Hastings reciprocity, endpoint correlation
- **NOT validated:** Multi-component worms

### Kawasaki dynamics (`KawasakiCore`)
- **Validated:** Signed magnetization conservation (exact)
- **Validated:** Energy decreases on cooling (directional)
- **Validated:** High-T equilibration
- **Quantitative equilibrium:** sector-restricted exact validation — BFS over the exchange graph from the initial state pins the reachable fixed-M sector (full sector on 4-ring M=0 and 8-ring M=+2; no hidden invariants), exact sector Boltzmann reference, ⟨E⟩ at β=0.3/0.8 with multi-seed |z|<2.4; sector ⟨E⟩ demonstrably ≠ canonical ⟨E⟩ (0.08 vs −3.20 at β=0.8)
- **NOT validated:** — (quantitative equilibrium validated 2026-08-14; superseded a physically-wrong cross-seed zombie test, see Known issues 3)

### BKL / n-fold-way (`BklIsingKernel`)
- **Validated:** Exact-trajectory reproducibility (bit-exact checkpoint)
- **Validated:** Fixed-time sampling matches exact small-Ising energy
- **Long-time equilibrium:** residence-time-weighted ⟨E⟩ and ⟨m²⟩ vs exact enumeration on N=4 and N=8 at β ∈ {0.2, 0.6, 1.0} (8 seeds; per-seed stderr from time-blocked means — per-visit averaging is length-biased and would be wrong, the estimator is pinned explicitly)
- **NOT validated:** — (long-time equilibrium validated 2026-08-14)

### Gillespie (`GillespieKernel`)
- **Validated:** Rate selection + exponential wait-time mean
- **Validated:** Absorbing-state clock advance
- **Multi-state equilibrium:** 3-state asymmetric CTMC — occupancy fractions vs exact stationary π (solved in-code via Cramer, πQ=0 verified to 1e-12); π = (0.1955, 0.4965, 0.3080); plus Ising-via-Gillespie ⟨E⟩ vs exact at β=0.6 (|z|≤2.1)
- **NOT validated:** — (multi-state equilibrium validated 2026-08-14)

### Event chain (`HardSphereEventChain`)
- **Validated:** Collision geometry, lifting at exact collision, PBC wrapping
- **Validated:** Snapshot restore
- **Equation of state:** contact-value identity Z = 1 + 2η(1−1/N)g(σ⁺) — event-chain collision rate (sausage average + Richardson extrapolation over chain lengths 1σ/2σ) reproduces the exact hard-disk virial series through B₃ (B₂ = πσ²/2, B₃ = 1.9295σ⁴) at η = 0.04/0.07 within 4σ; B₃ constant independently cross-checked by Mayer triple-integral MC (4σ)
- **Cross-solver:** g(σ⁺) at η=0.2 agrees with particle Metropolis NVT shell estimator within pooled 4σ (both in literature window 1.15–1.5); long variant adds η=0.12 and η=0.3
- **NOT validated:** — (EOS and pressure validated 2026-08-14; experimental status lifted)

### Particle NVT (`ParticleMetropolisCore`)
- **Validated:** Two-particle analytic pair energy
- **Validated:** Energy distribution vs quadrature
- **Validated:** Cache integrity, fixed-seed reproducibility
- **NOT validated:** Multi-particle EOS

### Particle NPT (`ParticleNptMetropolisCore`)
- **Validated:** V increases when P decreases (directional)
- **Validated:** V(P1)/V(P2) > 1.02 (non-trivial response)
- **Validated:** Finite-N ideal gas exact: ⟨V⟩ = (N+1)kT/P (long test, `npt_ideal_gas_volume_matches_finite_n_exact`)
- **Resolved:** Earlier "equilibrium volume mismatch" was a missing finite-N correction in the reference formula, not a solver bug

### Particle μVT (`ParticleGrandCanonicalCore`)
- **Validated:** N increases with μ (directional)
- **Validated:** N(μ1)/N(μ2) > 1.02 (non-trivial response)
- **Validated:** Ideal gas Poisson ⟨N⟩ exact (long test, `muvt_ideal_gas_particle_number_matches_poisson_most_probable`; plus `ideal_gas_grand_canonical_number_mean_is_poisson`)
- **Resolved:** Same finite-N reference correction as NPT

### Rigid molecule (`MolecularMetropolisCore`)
- **Validated:** Bond-length preservation, geometry preservation
- **Equilibrium distribution:** three analytic cases against in-code quadrature references through the real solver (translation + plane-rotation moves): two-molecule pair ⟨U⟩ and bound fraction (1D Simpson); dumbbell+atom probe nematic ⟨cos 2α⟩ and ⟨U⟩ (2D midpoint); rotor-pair alignment ⟨cos 2Δθ⟩ across a coupling sweep ε=1→3 (linear response → saturation, the Langevin-x analog; max|z|=1.26 default, 7-coupling long variant max|z|=1.7). Thermalization-length pitfall documented (20k+ sweeps needed at strong coupling)
- **NOT validated:** external-field Langevin case — the API has no one-body external term (report-only; see completion log 2026-08-14)

## Statistical validation framework

- **z-score tests:** 16 independent seeds per solver, |z|<4 per seed, |z̄|<1.5 mean, no one-sided bias. Seed counts scale via `SCUTTLE_ZSCORE_SEEDS` (unset → default 16 unchanged; nightly `zscore-monitor` uses 64; `just nightly-zscore` reproduces locally). Σz thresholds are scale-invariant (−2√n).
- **Cross-solver:** Metropolis vs Wolff pooled z-scores agree (|Δz̄|<2)
- **Connectivity:** Explicit Markov chain enumeration on N=2 Ising (4 states), BFS strong connectivity + aperiodicity check

## Known issues

1. ~~NPT/μVT equilibrium values~~ → Resolved 2026-08-14 (fixed in e1a07e4): the finite-N reference formulas (⟨V⟩=(N+1)kT/P, Poisson ⟨N⟩) match exactly; the old "mismatch" was a test-side formula error.
2. **strict-repro feature:** Was defined in Cargo.toml but had zero implementation. Removed.
3. ~~Kawasaki cross-seed zombie test~~ → Removed 2026-08-14: `kawasaki_2d_ising_energy_converges_same_regardless_of_seed` was #[ignore]'d with a reason stating it cannot pass (random starts land in different fixed-M sectors whose ⟨E⟩ genuinely differ), yet `--ignored` runs executed it anyway and it failed deterministically at HEAD. Physically-wrong criterion; superseded by the sector-restricted exact validation in `kawasaki_exact.rs`.

## Completion log

| Date | Task | Result |
|------|------|--------|
| 2026-07-23 | Gap 1: z-score framework | ✅ 5 tests (16 seeds × 3 solvers + cross-solver) |
| 2026-07-23 | Gap 2: Markov chain connectivity | ✅ 4 tests (Metropolis/Wolff/SW strong connectivity) |
| 2026-07-23 | Gap 3: Quantitative EOS | ✅ NPT/μVT directional + non-trivial response; documented equilibrium issue |
| 2026-07-23 | Gap 4: Validated domain docs | ✅ This document |
| 2026-08-14 | Production hardening: 8 solvers | ✅ 31 new tests: multicanonical full-MC (4), WL binned (4), event-chain EOS (4+1 long), molecule equilibrium (4+1 long), Kawasaki sector-exact (2), heat-bath discrete+O(N) (5), BKL/Gillespie equilibrium (4), continuous-spin cross-solver + frustrated triangle (6) |
| 2026-08-14 | Zombie test removal | ✅ Physically-wrong Kawasaki cross-seed test deleted; superseded by sector-exact validation |
