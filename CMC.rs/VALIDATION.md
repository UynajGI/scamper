# CMC.rs — Physics Validation & Validated Domain

> Updated 2026-07-23. Branch: `dev`.

## Test suite summary

| Layer | Tests | Runtime |
|-------|-------|---------|
| Default (`cargo test`) | 185 | ~12s |
| Long stochastic (`--ignored`) | 8 | ~70s |
| **Total** | **193** | |

## Per-solver validated domain

### Local Metropolis (`MetropolisCore`)
- **Validated:** Ising S=1/2 on chain/square/hypercubic graphs, PBC and OBC
- **Observables:** ⟨E⟩, ⟨m²⟩ vs exact enumeration (N=2,3,4,8)
- **Detailed balance:** Direct DB verified on N=2 (1e-15 precision)
- **Ergodicity:** 16-seed z-score framework (|z|<4, no systematic bias)
- **Connectivity:** Explicit Markov chain strongly connected on N=2
- **Cross-solver:** Agrees with Wolff and SW within 3σ on 8×8 Ising
- **NOT validated:** Continuous spins (XY/Heisenberg) via Metropolis, frustrated systems

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
- **NOT validated:** Exact energy comparison, Potts

### Heat bath — continuous O(N) (`ContinuousHeatBathCore`)
- **Validated:** O(3) uniform-on-sphere at β→0 (⟨s_α⟩≈0, ⟨s_α²⟩≈1/3)
- **NOT validated:** Finite-T distribution, XY (O(2))

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
- **NOT validated:** Continuous-axis (BinnedAxis) physical production run

### Multicanonical / umbrella (`EnergyBiasCore`)
- **Validated:** Transactional rejection, bias algebra, axis boundaries
- **Validated:** Discrete DOS reweighting matches enumeration
- **NOT validated:** Full MC-vs-exact distribution through EnergyBiasCore

### Worm (Ising HT graph) (`WormKernel`)
- **Validated:** HT graph partition identity matches spin enumeration
- **Validated:** Graph energy estimator matches exact spin energy
- **Validated:** Hastings reciprocity, endpoint correlation
- **NOT validated:** Multi-component worms

### Kawasaki dynamics (`KawasakiCore`)
- **Validated:** Signed magnetization conservation (exact)
- **Validated:** Energy decreases on cooling (directional)
- **Validated:** High-T equilibration
- **NOT validated:** Quantitative equilibrium distribution

### BKL / n-fold-way (`BklIsingKernel`)
- **Validated:** Exact-trajectory reproducibility (bit-exact checkpoint)
- **Validated:** Fixed-time sampling matches exact small-Ising energy
- **NOT validated:** Long-time equilibrium distribution

### Gillespie (`GillespieKernel`)
- **Validated:** Rate selection + exponential wait-time mean
- **Validated:** Absorbing-state clock advance
- **NOT validated:** Multi-state equilibrium distribution

### Event chain (`HardSphereEventChain`)
- **Validated:** Collision geometry, lifting at exact collision, PBC wrapping
- **Validated:** Snapshot restore
- **NOT validated:** Equation of state, pressure comparison

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
- **NOT validated:** Equilibrium distribution

## Statistical validation framework

- **z-score tests:** 16 independent seeds per solver, |z|<4 per seed, |z̄|<1.5 mean, no one-sided bias
- **Cross-solver:** Metropolis vs Wolff pooled z-scores agree (|Δz̄|<2)
- **Connectivity:** Explicit Markov chain enumeration on N=2 Ising (4 states), BFS strong connectivity + aperiodicity check

## Known issues

1. ~~NPT/μVT equilibrium values~~ → Resolved 2026-08-14 (fixed in e1a07e4): the finite-N reference formulas (⟨V⟩=(N+1)kT/P, Poisson ⟨N⟩) match exactly; the old "mismatch" was a test-side formula error.
2. **strict-repro feature:** Was defined in Cargo.toml but had zero implementation. Removed.

## Completion log

| Date | Task | Result |
|------|------|--------|
| 2026-07-23 | Gap 1: z-score framework | ✅ 5 tests (16 seeds × 3 solvers + cross-solver) |
| 2026-07-23 | Gap 2: Markov chain connectivity | ✅ 4 tests (Metropolis/Wolff/SW strong connectivity) |
| 2026-07-23 | Gap 3: Quantitative EOS | ✅ NPT/μVT directional + non-trivial response; documented equilibrium issue |
| 2026-07-23 | Gap 4: Validated domain docs | ✅ This document |
