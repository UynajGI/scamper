# Maturity Assessment — Physics Validation Coverage

> Updated 2026-07-22. Covers all 4 crates: Carlo.rs, CMC.rs, QMC.rs, MCMC.rs.
> Per-solver assessment against 8 production-readiness criteria.

## Maturity levels

| Level | Meaning |
|-------|---------|
| **Experimental** | Algorithm implemented, minimal validation |
| **Research-grade** | Validated vs ED/analytic on small systems |
| **Production-ready** | Physics + statistics + software sustained guarantees |
| **Stable** | Production + backward-compat commitment |

## Criteria

| # | Criterion | Standard |
|---|-----------|----------|
| A | Deterministic identities (machine precision) | `abs_error < 1e-12` |
| B | ED/exact comparison (≥3 observables) | E, ⟨M²⟩, C(τ) etc. |
| C | Analytic solvable limits | g=0, high-T, strong field |
| D | Per-update balance & invariants | Each update independently |
| E | Ergodicity evidence | Multi-init convergence, sector access |
| F | Cross-solver validation | Two independent solvers agree |
| G | Input validation | Invalid → error, never silent |
| H | Documented limitations | User-facing validated domain |

---

## 1. Summary status table

### Carlo.rs (framework)

| Component | Status | Notes |
|-----------|--------|-------|
| MC loop / Scheduler / Context | **stable** | 294 tests, lifecycle, reproducibility |
| Measurements / Accumulators | **stable** | Binning, covariance, jackknife |
| Checkpoint (in-memory) | **stable** | JSON serde, clock round-trip |
| Checkpoint (HDF5) | **experimental** | Zero tests; I/O paths untested |
| Error analysis / merge | **research-grade** | Regular autocorr tested; decorrelated only non-negativity |
| Parallel tempering (non-MPI) | **research-grade** | Wrapper tested; exchange dynamics untested |
| MPI backend | **experimental** | 2 ignored tests; controller/worker untested |
| RNG stream derivation | **stable** | Domain separation, reproducibility |

### CMC.rs (19 solvers)

| Solver | Status | Key gap |
|--------|--------|---------|
| Local Metropolis | **research-grade** | Ergodicity untested |
| Wolff cluster | **research-grade** | Ergodicity untested |
| Swendsen-Wang | **research-grade** | No direct exact comparison; cross-solver only |
| Heat bath (discrete) | **research-grade** | No exact-energy run |
| Heat bath (continuous O(N)) | **experimental** | Only "physical range" smoke |
| Microcanonical over-relaxation | **research-grade** | Energy conservation tested |
| Hybrid (composed) | **experimental** | Smoke only |
| MultiSpinIsing (64-replica) | **experimental** | No exact-energy or DB test |
| Wang-Landau DOS | **research-grade** | 2/4-site exact; 4×4 is #[ignore] |
| Multicanonical / umbrella | **experimental** | No MC-vs-exact distribution |
| Worm (Ising HT graph) | **research-grade** | Exact energy + endpoint correlation |
| Kawasaki dynamics | **research-grade** | M-conservation; no equilibrium check |
| BKL / n-fold-way | **research-grade** | Exact-trajectory reproducibility |
| Gillespie | **experimental** | 2-event toy only |
| Event chain (hard sphere) | **experimental** | Geometry only; no EOS comparison |
| Particle NVT | **research-grade** | Energy distribution vs quadrature |
| Particle NPT | **experimental** | Jacobian tested; no EOS |
| Particle μVT | **experimental** | Ideal-gas Poisson; no interacting |
| Rigid molecule | **experimental** | Geometry preservation only |

### QMC.rs (4 solvers)

| Solver | Status | Key gap |
|--------|--------|---------|
| Lattice directed-loop | **research-grade** | E+⟨m²⟩+NN Sz vs ED; no χ, no Binder vs ED |
| Wormhole (spin-boson) | **research-grade** | Free-limit only; **no interacting MC-vs-ED** |
| Occupation (cavity-QED) | **research-grade** | Strongest: ⟨σz⟩,⟨σx⟩,E,⟨n⟩ vs ED (Rabi+JC) |
| Cluster (longitudinal SB) | **research-grade** | ⟨σz⟩,⟨σx⟩,C(τ) vs ED; single-mode only |

### MCMC.rs (6 kernels + 3 combinators)

| Solver | Status | Key gap |
|--------|--------|---------|
| NUTS | **research-grade** | Distribution recovery; no detailed-balance test |
| Static HMC | **research-grade** | Same as NUTS |
| Random-walk Metropolis | **experimental** | Coarse moment recovery only |
| Component-wise Metropolis | **experimental** | No distribution-recovery test |
| Slice | **experimental** | Mean recovery only |
| Gibbs | **research-grade** | Exact conditional + atomicity |
| Combinators (Then/Repeat/Mixture) | **research-grade** | Determinism at boundaries |

---

## 2. Per-solver 8-criteria breakdown

### CMC.rs — Local Metropolis

| # | Criterion | Status | Evidence |
|---|-----------|--------|----------|
| A | Deterministic identities | ✅ PASS | `metropolis_cache_matches_exact_energy`; energy-cache audit |
| B | ED/exact (≥3 obs) | ✅ PASS | Onsager E, exact finite-Ising N=2–4, 2×2×2 within 4σ |
| C | Analytic limits | ✅ PASS | g=0, high-T ⟨m²⟩→1/4N, strong field polarizes, dimer exact |
| D | Per-update balance | ✅ PASS | `asymmetric_hastings_detailed_balance_n2` (direct DB) |
| E | Ergodicity | ❌ MISSING | No transition-graph connectivity, no multi-init test |
| F | Cross-solver | ✅ PASS | Metropolis vs Wolff vs SW on 8×8 (3σ) |
| G | Input validation | ✅ PASS | Topology mismatch, malformed params, unknown lattice |
| H | Documented limits | ⚠️ PARTIAL | README mentions energy audit; no explicit validated domain |

### CMC.rs — Wolff cluster

| # | Criterion | Status | Evidence |
|---|-----------|--------|----------|
| A | ✅ PASS | `wolff_batch_cache_matches_exact_energy` |
| B | ✅ PASS | Langevin O(3), Bessel-ratio O(2) exact; 8×8 cross-solver |
| C | ✅ PASS | XY/Heisenberg cluster activation, low-T M→1 |
| D | ✅ PASS | `wolff_detailed_balance_n3` (direct DB) |
| E | ❌ MISSING | No ergodicity test |
| F | ✅ PASS | Cross-solver 8×8 |
| G | ✅ PASS | Sign-problem rejection, model manifold |
| H | ⚠️ PARTIAL | |

### CMC.rs — Swendsen-Wang

| # | Criterion | Status | Evidence |
|---|-----------|--------|----------|
| A | ✅ PASS | Batch delta cache |
| B | ⚠️ PARTIAL | beta=0 ergodicity smoke + 8×8 cross-solver; no direct exact |
| C | ⚠️ PARTIAL | SW assign independent states at beta=0 |
| D | ❌ MISSING | No direct DB test for SW (only Wolff has it) |
| E | ❌ MISSING | |
| F | ✅ PASS | 8×8 cross-solver |
| G | ✅ PASS | |
| H | ❌ MISSING | |

### QMC.rs — Lattice directed-loop

| # | Criterion | Status | Evidence |
|---|-----------|--------|----------|
| A | ⚠️ PARTIAL | Scattering table balance (S=1, 1e-12); generic S missing |
| B | ✅ PASS | 3-site Heisenberg: E, ⟨Sz_iSz_j⟩, ⟨m²⟩ vs ED |
| C | ❌ MISSING | No zero-coupling, high-T, strong-field tests for lattice QMC |
| D | ✅ PASS | Scattering table detailed balance (both policies) |
| E | ❌ MISSING | No multi-init convergence |
| F | ❌ MISSING | No QMC cross-solver |
| G | ✅ PASS | Sign-problem rejection, frustration detection |
| H | ✅ PASS | README documents S>1/2 caveat, validated domain |

### QMC.rs — Wormhole (spin-boson)

| # | Criterion | Status | Evidence |
|---|-----------|--------|----------|
| A | ✅ PASS | Scattering DB both policies (1e-12); vertex weights |
| B | ❌ MISSING | **No interacting MC-vs-ED**; only free-limit (Poisson, tanh) |
| C | ✅ PASS | Poisson expansion, 2-state partition, spin-inversion symmetry |
| D | ✅ PASS | Table-level DB; loop-abort rollback; worldline invariants |
| E | ❌ MISSING | |
| F | ❌ MISSING | cross_solver.rs documents conventions but no numerical check |
| G | ✅ PASS | Channel validation, non-stoquastic rejection |
| H | ✅ PASS | README documents basis rotation, convention differences |

### QMC.rs — Occupation (cavity-QED)

| # | Criterion | Status | Evidence |
|---|-----------|--------|----------|
| A | ✅ PASS | Transfer matrix eigensystem; sqrt(n) matrix elements |
| B | ✅ PASS | ⟨σz⟩,⟨σx⟩,E,⟨n⟩ vs ED (Rabi+JC); long crossover #[ignore] |
| C | ✅ PASS | Free spin tanh, uncoupled Bose distribution, exact Z |
| D | ⚠️ PARTIAL | Non-stoquastic rejection; no per-update DB (transfer matrix) |
| E | ⚠️ PARTIAL | `sampler_only_visits_states_within_basis` (bounds, not connectivity) |
| F | ❌ MISSING | |
| G | ✅ PASS | Beta/slices/cutoff validation |
| H | ✅ PASS | README + Rabi QPT tests document domain |

### QMC.rs — Cluster (longitudinal SB)

| # | Criterion | Status | Evidence |
|---|-----------|--------|----------|
| A | ✅ PASS | Retarded kernel KMS symmetry; quadrature normalization |
| B | ✅ PASS | ⟨σz⟩,⟨σx⟩ (via kinks), C(β/2) vs exact cosh; interacting ED |
| C | ✅ PASS | Free 2-level: kink count vs exact mean; ⟨σz⟩ vs tanh |
| D | ✅ PASS | Worldline invariants (even kink); 10k updates validated |
| E | ❌ MISSING | |
| F | ❌ MISSING | |
| G | ✅ PASS | |
| H | ⚠️ PARTIAL | Single-mode only; multi-mode interacting untested |

### MCMC.rs — NUTS

| # | Criterion | Status | Evidence |
|---|-----------|--------|----------|
| A | ✅ PASS | Leapfrog reversibility (1e-12); metric velocity exact |
| B | ⚠️ PARTIAL | N(0,1) moments, correlated Gaussian cov; only Gaussian targets |
| C | ❌ MISSING | No non-Gaussian recovery (mixture, Student-t, multimodal) |
| D | ❌ MISSING | **No detailed-balance test** for accept/reject step |
| E | ⚠️ PARTIAL | Step-size search rescues bad scale; no sector test |
| F | ⚠️ PARTIAL | NUTS vs HMC step-size search only, not posterior agreement |
| G | ✅ PASS | NaN, dim mismatch, non-PD mass matrix, bad configs |
| H | ❌ MISSING | |

### MCMC.rs — Static HMC

Same as NUTS minus U-turn tests.

### MCMC.rs — Random-walk Metropolis

| # | Criterion | Status | Evidence |
|---|-----------|--------|----------|
| A | ⚠️ PARTIAL | Scale adaptation freeze; no proposal-ratio machine-precision test |
| B | ⚠️ PARTIAL | N(0,1) mean only (30k draws); no covariance |
| C | ❌ MISSING | |
| D | ❌ MISSING | |
| E | ❌ MISSING | |
| F | ❌ MISSING | |
| G | ✅ PASS | Zero-dim, bad covariance rejected |
| H | ❌ MISSING | |

### Carlo.rs — Checkpoint (HDF5)

| # | Criterion | Status | Evidence |
|---|-----------|--------|----------|
| A | ❌ MISSING | No HDF5 round-trip test at all |
| B | ❌ MISSING | |
| C | ❌ MISSING | |
| D | ❌ MISSING | read_checkpoint_hdf5_full drops clocks (attempt/accepted/event_time) |
| E | ❌ MISSING | |
| F | ❌ MISSING | |
| G | ❌ MISSING | |
| H | ❌ MISSING | |

---

## 3. Cross-cutting gaps (universal blockers)

### Every physics solver (CMC + QMC + MCMC)

1. **Ergodicity entirely untested.** No transition-graph connectivity, no multi-initial-state convergence test.
2. **No multi-seed nightly statistical monitoring.** No z-score stability tracking.

### QMC.rs specific

3. **Cross-solver validation absent.** No wormhole↔occupation, wormhole↔cluster numerical comparison.
4. **Susceptibility χ, Binder M⁴, full C(τ)** measured everywhere, validated almost nowhere.

### MCMC.rs specific

5. **No detailed-balance test** for any Metropolis-type kernel (only leapfrog integrator reversibility).
6. **No AR(1) reference test for ESS.** ESS validated only against IID heuristics.
7. **No cross-solver moment agreement.** Samplers never compared on same posterior.

### Carlo.rs specific

8. **HDF5 checkpoint I/O completely untested.**
9. **MPI controller/worker + PT exchange untested** (2 ignored tests only).
10. **Estimate.autocorr_time hardcoded to 1.0** in fast path.

---

## 4. Prioritized roadmap

### P0 — Correctness risks (must fix)

| # | Task | Crate | Why |
|---|------|-------|-----|
| P0.1 | Wormhole interacting MC-vs-ED (≥3 obs) | QMC | Weakest solver; only free-limit validated |
| P0.2 | Lattice QMC analytic limits (g=0, high-T, strong field) | QMC | Basic limits missing |
| P0.3 | Lattice QMC susceptibility χ_z vs ED | QMC | χ measured but unvalidated |
| P0.4 | MCMC detailed-balance test for accept/reject | MCMC | Core MCMC correctness unverified |
| P0.5 | MCMC AR(1) reference test for ESS | MCMC | ESS estimator uncalibrated |
| P0.6 | Carlo.rs HDF5 checkpoint round-trip | Carlo | Production restart path untested |

### P1 — Production blockers

| # | Task | Crate |
|---|------|-------|
| P1.1 | Ergodicity: multi-init convergence test (all solvers) | ALL |
| P1.2 | QMC cross-solver: wormhole↔occupation numerical agreement | QMC |
| P1.3 | QMC cross-solver: wormhole↔cluster numerical agreement | QMC |
| P1.4 | SW detailed-balance test (direct, not just cross-solver) | CMC |
| P1.5 | MCMC cross-solver: NUTS vs RWM vs Slice posterior agreement | MCMC |
| P1.6 | Carlo.rs MPI PT exchange protocol test | Carlo |
| P1.7 | CMC Wang-Landau 4×4 un-#[ignore] | CMC |
| P1.8 | MCMC non-Gaussian recovery (mixture, Student-t) | MCMC |

### P2 — Robustness

| # | Task | Crate |
|---|------|-------|
| P2.1 | Acceptance formula machine-precision validation | CMC, QMC |
| P2.2 | Continuous heat-bath analytic distribution match | CMC |
| P2.3 | NPT/μVT interacting equation-of-state | CMC |
| P2.4 | MultiSpinIsing exact-energy test | CMC |
| P2.5 | Thread-count independence test | ALL |
| P2.6 | Carlo.rs Estimate.autocorr_time real estimation | Carlo |
| P2.7 | strict-repro feature test | Carlo |
| P2.8 | Multi-seed nightly z-score monitoring infrastructure | ALL |

---

## 5. Status statement

**Repository: research-grade (candidate). No solver is production-ready.**

- **Carlo.rs**: stable framework core; HDF5 checkpoint and MPI are experimental
- **CMC.rs**: Metropolis/Wolff/Wang-Landau/Worm are research-grade with strong DB+exact coverage; ergodicity is the universal gap
- **QMC.rs**: Occupation and Cluster are research-grade with genuine MC-vs-ED; Wormhole needs interacting ED comparison; Lattice needs analytic limits
- **MCMC.rs**: NUTS/HMC are research-grade for Gaussian targets; detailed balance and ESS calibration are the key gaps
