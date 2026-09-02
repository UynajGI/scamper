# CMC.rs — Physics Validation & Validated Domain

> Updated 2026-09-02. Branch: `dev`.

## Test suite summary

| Layer | Tests | Runtime |
|-------|-------|---------|
| Default (`cargo test`) | 289 | ~60s |
| Long stochastic (`--ignored`) | 17 | ~40s |
| **Suite total** | **306** | (+90 lib unit tests) |

## Per-solver validated domain

### Local Metropolis (`MetropolisCore`)
- **Validated:** Ising S=1/2 on chain/square/hypercubic graphs, PBC and OBC; q-state Potts (q=3,4) via cross-solver and β=0/strong-coupling limits
- **Observables:** ⟨E⟩, ⟨m²⟩ vs exact enumeration (N=2,3,4,8)
- **Detailed balance:** Direct DB verified on N=2 (1e-15 precision)
- **Ergodicity:** 16-seed z-score framework (|z|<4, no systematic bias)
- **Connectivity:** Explicit Markov chain strongly connected on N=2
- **Cross-solver:** Agrees with Wolff and SW within 3σ on 8×8 Ising; continuous-spin cross-solver XY (O(2)) and Heisenberg (O(3)) vs Wolff on 8×8 (pooled |z|<4 at β=0.7/0.9 and 0.5/0.9); 8×8 q=3 Potts 4-solver agreement (see Potts section)
- **Frustrated:** AFM XY triangle (J<0): full equilibrium vs spectral quadrature reference — ⟨E⟩ and chirality ⟨κ²⟩ at β=1,3 (per-seed |z|<3, ground-state algebra exact)
- **NOT validated:** —

### Wolff cluster (`WolffCore`)
- **Validated:** Ising S=1/2, XY (O(2)), Heisenberg (O(3)) on chain/square; q-state Potts q=3,4 (2026-08-19)
- **Observables:** ⟨E⟩ vs exact enumeration, Langevin/Bessel-ratio analytic limits; Potts ⟨E⟩, ⟨m²⟩, C vs full q^N enumeration (N=4, 8; β=0.3/0.8; per-seed |z|<4, pooled Σz gate)
- **Detailed balance:** Direct DB verified on N=3
- **Ergodicity:** 16-seed z-score, strongly connected on N=2
- **NOT validated:** —

### Swendsen-Wang (`SWCore`)
- **Validated:** Ising S=1/2 on chain/square; q-state Potts q=3,4 (2026-08-19); continuous O(2)/O(3) cluster updates (2026-08-19, see below)
- **Detailed balance:** Direct DB verified on N=2 (50k samples/state)
- **Exact equilibrium:** Potts ⟨E⟩, ⟨m²⟩, C vs full q^N enumeration on N=4/N=8 at β=0.3/0.8 (16-seed z-scores); β=0 exactly-uniform state distribution (⟨E⟩ = −JΣw/q closed form)
- **Continuous spins (2026-08-19):** XY 4-ring vs exact spectral quadrature — ⟨E⟩, ⟨m²⟩, ⟨cos Δθ⟩ at β ∈ {0.6, 1.2}, 8 seeds, |z| < 4 per seed, |z̄| < 2; O(3) 4-ring analytic limits (β→0: ⟨E⟩ = 0 and ⟨m²⟩ = 1/N exactly; β=8 spin-wave ⟨E⟩ = −J·N + (N−1)/β); cross-solver vs Wolff on 8×8 O(2)/O(3) at β=0.9 (pooled |z| < 4 on ⟨E⟩ and ⟨m²⟩) — `sw_continuous.rs`
- **Ergodicity:** Strongly connected on N=2
- **Cross-solver:** Agrees with Metropolis within 3σ on 8×8; 8×8 q=3 near βc = ln(1+√3) 4-solver agreement
- **NOT validated:** —

### Heat bath — discrete (`HeatBathCore`)
- **Validated:** Ising S=1/2; q-state Potts q=3,4 (2026-08-19)
- **Detailed balance:** Direct DB verified on N=2
- **Exact energies:** ⟨E⟩ and ⟨m²⟩ vs exact enumeration on N=4 and N=8 at β=0.3/0.8 (16-seed z-scores, max|z|<2.5); Potts adds ⟨E⟩, ⟨m²⟩, C on N=4/N=8 (2026-08-19)
- **NOT validated:** —

### Heat bath — continuous O(N) (`ContinuousHeatBathCore`)
- **Validated:** O(3) uniform-on-sphere at β→0 (⟨s_α⟩≈0, ⟨s_α²⟩≈1/3)
- **Finite-T conditional distribution:** single-spin ⟨cosθ⟩/⟨cos²θ⟩ vs exact O(2) Bessel ratio I₁(x)/I₀(x) and O(3) Langevin coth(x)−1/x at x ∈ [0.2, 5]; analytic anchors (limits at small/large x, Bessel recurrence) to 1e-4–1e-12; two-spin XY pair equilibrium matches the conditional moment (rotational invariance)
- **NOT validated:** — (finite-T distribution and XY validated 2026-08-14; lattice-scale O(N) cross-solver via Metropolis-vs-Wolff above)

### Microcanonical over-relaxation (`MicrocanonicalCore`)
- **Validated:** Energy conservation to 2e-10, unit norm preservation to 3e-12 (long sweeps)
- **Reflection-map identities (2026-08-19, machine precision):** R(s) = 2(s·ĥ)ĥ − s is an exact involution (R∘R = id, 1e-15), norm-preserving (1e-15), preserves the local-field projection (bond energies unchanged, 1e-15), and is an isometry with |Jacobian| = 1 — O(2): θ' = 2φ − θ exactly (residual < 1e-14); O(3): the map matrix 2ħħᵀ − I is orthogonal with det = +1 (π rotation about the field axis; the field direction is fixed, the orthogonal plane reversed). A deterministic involutive isometry satisfies detailed balance against any rotation-invariant measure: T(s→s') = T(s'→s) = 1
- **Kernel ≡ the reflection map (2026-08-19):** with the sequential visit schedule the kernel sweep is bit-identical to the manual per-site reflection, and every transactional update reports |ΔE| < 1e-12 (per-update machine-precision conservation; `energy_error` after full sweeps ≤ 1e-12)
- **Composition equilibrium (2026-08-19):** Hybrid(Metropolis, Microcanonical) on the XY 4-ring vs exact spectral quadrature (zero mode factored, grid-doubling anchored) — ⟨E⟩, ⟨m²⟩, ⟨cos(θ1−θ3)⟩ at β ∈ {0.6, 1.2}, from hot and cold initializations, 8 seeds, |z| < 4 per seed, |z̄| < 2 per observable
- **Cross-solver (2026-08-19):** the composition vs Wolff on 8×8 O(2) and O(3) at β=0.9 (pooled |z| < 4 on ⟨E⟩ and ⟨m²⟩)
- **Analytic limits in composition (2026-08-19):** β→0 gives exactly ⟨E⟩ = 0 and ⟨m²⟩ = 1/N; β=8 approaches the spin-wave result ⟨E⟩ = −J·N + (N−1)/(2β) (tolerance 0.05; anharmonic remainder O(β⁻²)) with ⟨m²⟩ → 1
- **Physics note:** over-relaxation alone is deterministic and energy-conserving, hence **not ergodic** — this is physics, not a defect; the validated production mode is composition with an ergodic kernel (`over_relaxation.rs`)
- **NOT validated:** —

### Hybrid (`HybridCore<A, B>`)
- **Validated:** Hybrid(Metropolis, Wolff) ⟨E⟩ and ⟨m²⟩ vs exact enumeration on 4-site Ising (p2_remaining)
- **All Ising pairwise compositions (2026-08-19):** every pairing of the four Ising-capable kernels — Metropolis+Wolff, Metropolis+SW, Metropolis+HeatBath, Wolff+SW, Wolff+HeatBath, SW+HeatBath — vs full 2^4 enumeration on the 4-site ring at β=0.5: ⟨E⟩, ⟨m²⟩ and C = β²(⟨E²⟩−⟨E⟩²), 8 seeds each, |z| < 4 per seed, |z̄| < 2 per combo, pooled one-sided Σz gate across the matrix
- **Composition boundary semantics (2026-08-19):** same-seed determinism; Hybrid(A, B) ≡ manual `A; B` sequencing bit-for-bit; `repetitions(2, 3)` ≡ `A; A; B; B; B`; the nested combinator closure Hybrid(A, Hybrid(B, C)) ≡ `A; B; C`
- **Continuous compositions (2026-08-19):** O(2) Wolff+SW, Metropolis+Wolff and Metropolis+ContinuousHeatBath vs a pure-Wolff reference on 8×8 at β=0.9 (pooled |z| < 4 on ⟨E⟩ and ⟨m²⟩); Metropolis+Microcanonical validated in the microcanonical section (quadrature, limits, Wolff cross-solver)
- **NOT validated:** —

### MultiSpinIsing (64-replica bit-parallel, `MultiSpinIsing`)
- **Validated:** 8-site exact enumeration cross-check (`p2_validation.rs`); acceptance LUT anchors, PT weight ratio, per-replica array observables (lib tests)
- **Multi-seed z coverage (2026-08-19):** 8-site PBC chain, β ∈ {0.4, 0.8}, 8 seeds — ⟨E⟩, ⟨m²⟩ and C = β²(⟨E²⟩−⟨E⟩²) vs full 2^8 enumeration, |z| < 4 per seed, |z̄| < 2 per (β, observable), pooled one-sided Σz gate
- **All 64 replicas ensemble-consistent (2026-08-19):** the `Energy_replica` array observable (averaged over replicas and time) matches the same exact ⟨E⟩ (|z̄| < 1) — replica 0 is not special; the replicas are exchangeable valid Metropolis chains
- **Cross-solver (2026-08-19):** vs scalar per-site `MetropolisCore` on identical physics (same lattice, β, J) at both temperatures — pooled cross-solver z on ⟨E⟩, ⟨m²⟩ and ⟨|m|⟩, all |z| < 4 (`multispin_cross_solver.rs`)
- **NOT validated:** —

### Wang-Landau (`WangLandauCore`)
- **Validated:** DOS matches exact 4×4 Ising enumeration (un-ignored, ~11s)
- **Validated:** 2-site exact DOS levels/degeneracies
- **Validated:** Canonical reweighting recovers exact ⟨E⟩
- **BinnedAxis production run:** binned DOS vs exact binned degeneracies (8-level weighted 6-ring, 14-bin axis; per-bin |Δln g| RMS 0.013, gate 0.05); canonical reweighting ⟨E⟩(T) at 3 temperatures vs exact (|z|≤1.2); flat-histogram-only route error floor documented (RMS ≤0.1)
- **Convergence robustness (2026-08-19):** an unattainable `minimum_visited_fraction` (more visited bins demanded than physically reachable) terminates loudly with `WangLandauTermination::UnreachableBins` after the discovery plateau is established (500 stalled flatness checks), instead of silently burning sweeps to the maximum-sweep guard; checkpoint round-trips the failure and its evidence; version-1 checkpoints without the stall fields still load
- **NOT validated:** —

### Multicanonical / umbrella (`EnergyBiasCore`)
- **Validated:** Transactional rejection, bias algebra, axis boundaries
- **Validated:** Discrete DOS reweighting matches enumeration
- **Full MC-vs-exact:** real EnergyBiasCore sampling with exact-DOS bias on the 6-site ring — ⟨E⟩(T), ⟨m²⟩(T) at β ∈ {0.2, 0.5, 1.0} (max|z|=3.4 over 16 seeds, |z̄|≤0.65), full P(E) at β=0.5 within 4σ binomial band (TV distance 0.005 vs gate 0.02), biased histogram flatness ≥0.90 (canonical would give 0.003)
- **NOT validated:** — (full MC distribution validated 2026-08-14; experimental status lifted)

### Worm (Ising HT graph) (`WormKernel`)
- **Validated:** HT graph partition identity matches spin enumeration
- **Validated:** Graph energy estimator matches exact spin energy
- **Validated:** Hastings reciprocity, endpoint correlation
- **Cross-solver (2026-08-19):** 4×4 square ⟨E⟩ at β=0.44 agrees with spin Metropolis within pooled 4σ
- **Multi-component lattices supported (2026-08-19, second pass):** the HT-graph ensemble factorizes over connected components, so `IsingGraphWormMC::from_lattice` / `IsingGraphWormEnsemble` (`src/worm/ensemble.rs`) run one independent two-defect worm per component on domain-separated derived streams (`carlo_rs::RngStreamKey`, component index in the replica field; one salt per component per sweep from the shared context stream, so no RNG state is hidden from a checkpoint). Observables combine additively; total energy is measured under the all-physical conditioning (the product ensemble is preserved); endpoint correlations are per component (cross-component pairs have no worm estimator — they factorize to zero); isolated sites form trivial components sampled exactly by the empty graph. Validated vs full 2^8 spin enumeration on a 4-ring + 3-chain + isolated-site lattice (total and per-component ⟨E⟩, two-point correlations, and ⟨m²⟩ reconstructed from the worm pair correlations); partition identity Z_spin = 2^N Π cosh(βJ_e) Π_c Z_graph,c at 1e-10; cross-solver vs spin Metropolis on two disjoint 4×4 squares (pooled |z| < 4); v2 multi-component snapshots round-trip bit-exact trajectories, v1 single-component snapshots still load (and are rejected loudly for multi-component ensembles)
- **Input rejection:** genuinely invalid input is still rejected loudly at construction — empty lattice, non-finite/negative β, non-finite coupling, `J · weight < 0` on any edge, self-loops. The raw `IsingGraphWormModel` + `WormKernel` pair additionally requires a **connected** lattice: its single defect pair would silently freeze the other components, so `IsingGraphWormModel::new` rejects disconnected input for direct users while the ensemble adapter handles it by decomposition
- **NOT validated:** multi-defect / multi-leg (multi-component) worm algorithms — not implemented; the kernel is a two-defect kernel by design (documented in `worm` module docs)

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
- **External field (2026-08-19):** one-body dipolar term `DipolarExternalField` (per-atom charges, wrap-safe minimum-image dipoles, non-neutral molecules rejected loudly); free-rotor equilibrium vs the analytic Langevin-dipole answers through the real kernel — 2D: ⟨cosθ⟩=I₁(x)/I₀(x), ⟨cos²θ⟩=(1+I₂/I₀)/2; 3D: ⟨cosθ⟩=L(x), ⟨cos²θ⟩=1−2L(x)/x; x=βpE grid 0.5–5, per-seed |z|<4; machine-precision identity `external_field_energy = −E·μ` (1e-12) every sweep
- **NOT validated:** — (external-field Langevin case validated 2026-08-19)

### Percolation, site + bond (`PercolationMC`, 2026-09-02)
- **Validated:** Ordinary site and bond percolation on arbitrary `CsrLattice` graphs (i.i.d. occupancy resampling, union-find cluster analysis). 2×2 open square: full 16-configuration enumeration vs hand-derived closed-form moments for both modes — ⟨MaxCluster⟩ = 30/16 (site) and 45/16 (bond), ⟨sum(s_i²)⟩ = 76/16 and 164/16, ⟨NClusters⟩ = 17/16 and 33/16, P(spanning) = 7/16 and 12/16 at p = 1/2; site spanning probability matches the polynomial 2p²(1−p)² + 4p³(1−p) + p⁴ across p ∈ {0.2, 0.44, 0.5927, 0.8} — `tests/physics/percolation.rs`
- **Independent algorithm cross-check:** `cluster_stats` (union find) vs an in-test flood-fill reference sharing no algorithmic path, configuration-by-configuration, both modes — exhaustive on chain-8, square-3x3, cubic-2x2x2, triangular-2x2, honeycomb-2x2, kagome-2x2 site, random-graph site (≈22k configurations); seeded random configurations beyond (kagome/random bond)
- **1D exact solution:** open chain P(span) = p^L (site) and p^(L−1) (bond) — exact enumeration at p ∈ {0.3, 0.6, 0.9} and through the full scheduler stack at L = 6 (100k sweeps, |z| < 4)
- **Scheduler end-to-end:** 200k i.i.d. sweeps on the 2×2 square reproduce all four enumerated moments within |z| < 4; p = 0 and p = 1 boundary behavior exact in unit tests (no span / single spanning cluster); fixed seed reproduces bitwise, different seed shifts the stream
- **Statistical:** 4×4 site at p = 0.6 (2¹⁶ = 65536 configurations) and 3×3 bond at p = 0.55 near p_c = 1/2 (2¹² bond configurations) fully enumerated as references; 16-seed z-scores on `Spanning` and `MaxCluster` for both modes (|z| < 4, |z̄| < 1.5, no one-sided bias) — `tests/physics/percolation_zscore.rs`
- **Physics sanity:** crossing probability monotone non-decreasing in p (8×8, both modes, p = 0.1…0.9); overlapping spanning-set degeneracy pinned by unit test
- **Critical-point checks (long, `#[ignore]`, nightly):** 32×32 bond at p_c = 1/2 (exact self-duality) → crossing within 0.06 of 1/2; 16³ cubic bond brackets the critical region — P(cross) < 0.05 at p = 0.12 and > 0.95 at p = 0.40 around p_c ≈ 0.2488 (no unproven 3D duality assumed)
- **NOT validated:** invasion/kinetic percolation variants (not implemented); spanning defaults limited to square/chain (other graphs require explicit site sets, rejected loudly otherwise); critical crossing on pbc triangular/honeycomb/kagome (builders are periodic-only; no clean crossing convention)

## Input-validation coverage (criterion G)

`tests/physics/input_validation.rs` (2026-08-19) sweeps the parameter and
constructor surface of every scheduler-ready solver; invalid input must
return `InvalidConfig` (or the kernel error type), never a silent accept and
never an unintended panic:

- Lattice kernels (Metropolis/Wolff/SW/heat baths/continuous heat
  bath/microcanonical/hybrid): negative/NaN β, unknown lattice names, zero
  or malformed dimensions, OBC triangles, Potts q<2, non-finite J, unknown
  `initial_state`.
- MultiSpinIsing, Wang-Landau (fraction/flatness/log_f/interval ranges,
  24-site exact-axis limit), worm (β, coupling sign, close probability,
  fugacity), Kawasaki/BKL (incl. the fixed β/J validation — previously
  `J = NaN` panicked in an assert-backed constructor), BKL event windows,
  event chain (box/particles/diameter/chain lengths), particle NVT/NPT/μVT
  (density, cutoff, displacement, pressure, activity) and the molecule
  kernel (scales, D≥2, topology corruption).

## Statistical validation framework

- **z-score tests:** 16 independent seeds per solver (8 where chains are more expensive), |z|<4 per seed, |z̄|<2 mean, no one-sided bias. Seed counts scale via `SCUTTLE_ZSCORE_SEEDS` (unset → default unchanged; nightly `zscore-monitor` uses 64; `just nightly-zscore` reproduces locally). Σz thresholds are scale-invariant (−2√n); in test files with many configuration cells (Potts, external field) the Σz gate is pooled over all scores of one solver/test to control the multiple-comparisons false-alarm rate.
- **Cross-solver:** Metropolis vs Wolff pooled z-scores agree (|Δz̄|<2); Potts 4-solver pairwise |z|<4 at βc; worm vs Metropolis pooled |z|<4; MultiSpinIsing vs scalar Metropolis pooled |z|<4; SW-continuous and over-relaxation-composition vs Wolff pooled |z|<4
- **Connectivity:** Explicit Markov chain enumeration on N=2 Ising (4 states), BFS strong connectivity + aperiodicity check

## Known issues

1. ~~NPT/μVT equilibrium values~~ → Resolved 2026-08-14 (fixed in e1a07e4): the finite-N reference formulas (⟨V⟩=(N+1)kT/P, Poisson ⟨N⟩) match exactly; the old "mismatch" was a test-side formula error.
2. **strict-repro feature:** Was defined in Cargo.toml but had zero implementation. Removed.
3. ~~Kawasaki cross-seed zombie test~~ → Removed 2026-08-14: `kawasaki_2d_ising_energy_converges_same_regardless_of_seed` was #[ignore]'d with a reason stating it cannot pass (random starts land in different fixed-M sectors whose ⟨E⟩ genuinely differ), yet `--ignored` runs executed it anyway and it failed deterministically at HEAD. Physically-wrong criterion; superseded by the sector-restricted exact validation in `kawasaki_exact.rs`.
4. ~~Kawasaki/BKL β/J validation~~ → Fixed 2026-08-19: `from_params` previously forwarded user `beta`/`J` into assert-backed constructors, so `J = "NaN"` panicked instead of returning `InvalidConfig`. Both adapters now validate through a shared `validate_kinetic_ising_params`.

## Completion log

| Date | Task | Result |
|------|------|--------|
| 2026-07-23 | Gap 1: z-score framework | ✅ 5 tests (16 seeds × 3 solvers + cross-solver) |
| 2026-07-23 | Gap 2: Markov chain connectivity | ✅ 4 tests (Metropolis/Wolff/SW strong connectivity) |
| 2026-07-23 | Gap 3: Quantitative EOS | ✅ NPT/μVT directional + non-trivial response; documented equilibrium issue |
| 2026-07-23 | Gap 4: Validated domain docs | ✅ This document |
| 2026-08-14 | Production hardening: 8 solvers | ✅ 31 new tests: multicanonical full-MC (4), WL binned (4), event-chain EOS (4+1 long), molecule equilibrium (4+1 long), Kawasaki sector-exact (2), heat-bath discrete+O(N) (5), BKL/Gillespie equilibrium (4), continuous-spin cross-solver + frustrated triangle (6) |
| 2026-08-14 | Zombie test removal | ✅ Physically-wrong Kawasaki cross-seed test deleted; superseded by sector-exact validation |
| 2026-08-19 | Production: Potts q>2 | ✅ 8 tests (`potts_exact.rs`): full q^N enumeration on 2×2 square and N=8 chain, q ∈ {3,4}, β ∈ {0.3, 0.8}, observables ⟨E⟩/⟨m²⟩/C for heat bath + SW + Wolff; enumeration anchors (q=1 degenerate, q=2 ≡ Ising with J/2); β=0 uniform and β=8 frozen analytic limits; 8×8 q=3 4-solver cross-solver at βc; q=4 βc = ln 3 directional anchor |
| 2026-08-19 | Production: molecule external field | ✅ 4 tests (`molecule_external_field.rs`): `DipolarExternalField` API (additive, backward-compatible), 2D von Mises + 3D Langevin free-rotor moments, −E·μ machine-precision identity, loud rejection of non-neutral/short/non-finite charge tables |
| 2026-08-19 | Production: WL convergence robustness | ✅ 4 tests (`wang_landau_convergence.rs`): unattainable `minimum_visited_fraction` → loud `UnreachableBins` termination (not Converged/MaximumSweeps), auto-derived reachable set matches enumeration, checkpoint evidence guards, version-1 back-compat |
| 2026-08-19 | Production: multi-component worm | ✅ Honest limitation route: multi-component (disconnected/isolated-site) lattices now rejected loudly at `IsingGraphWormModel::new` (defect pair confined to one component would silently freeze the rest); multi-defect/multi-leg worms documented as not implemented; + worm-vs-Metropolis cross-solver test |
| 2026-08-19 | Production: input-validation audit | ✅ 10 tests (`input_validation.rs`) across all 19 solvers; fixed three source holes (Kawasaki/BKL β+J panic path, MultiSpinIsing lattice validation, BKL event-window validation) |
| 2026-08-19 | Production: multi-component worm implemented | ✅ `src/worm/ensemble.rs` + `IsingGraphWormMC::from_lattice` (per-component two-defect worms, `RngStreamKey` domain-separated streams, additive observables, v1/v2 checkpoints) + `CsrLattice::connected_components`; 4 suite tests (`worm_multi_component.rs`: exact 2^8 enumeration incl. per-component energies/correlations/worm-reconstructed m², partition identity at 1e-10, snapshot round-trip + loud rejections, cross-solver vs spin Metropolis on two disjoint 4×4 squares) + 3 lib tests; multi-component input no longer rejected at the scheduler surface — genuinely invalid input still is |
| 2026-08-19 | Production: over-relaxation in composition | ✅ 6 tests (`over_relaxation.rs`): reflection-map machine-precision identities + deterministic DB (involution/isometry, O(2) angle form, O(3) orthogonal π-rotation); kernel ≡ manual reflection bit-exact with per-update \|ΔE\|<1e-12; Hybrid(Metropolis, Microcanonical) vs exact XY-ring quadrature from 2 inits; analytic limits (β→0 exact, β=8 spin-wave); cross-solver vs Wolff 8×8 O(2)/O(3) |
| 2026-08-19 | Production: all hybrid compositions | ✅ 3 tests (`hybrid_compositions.rs`): all six Ising pairwise combos vs exact enumeration (E, m², C; multi-seed + pooled Σz); boundary semantics bit-exact (sequencing/repetitions/nesting/determinism); continuous O(2) combos vs Wolff 8×8 |
| 2026-08-19 | Production: MultiSpinIsing | ✅ 2 tests (`multispin_cross_solver.rs`): multi-seed z vs exact enumeration (E, m², C at β=0.4/0.8) + all-64-replica array-observable consistency; cross-solver vs scalar Metropolis (E, m², \|m\| pooled z) |
| 2026-08-19 | Production: SW continuous spins | ✅ 4 tests (`sw_continuous.rs`): XY-ring exact quadrature (E, m², cos Δθ), O(3) analytic limits (β→0, β=8 spin-wave with the zero-mode-counted formula), cross-solver vs Wolff 8×8 O(2)/O(3) |
