# Physics Validation Task Tracker

> Created 2026-07-22 on `dev` branch. Updated 2026-07-23.
> Goal: bring every solver from "compiles and runs" to "physics verified".
> Assessment baseline: `MATURITY_ASSESSMENT.md`.

---

## Progress overview

| Tier | Total | Done | Remaining |
|------|-------|------|-----------|
| P0 | 6 | 4 | 2 (MCMC) |
| P1 | 8 | 6 | 2 (MCMC) |
| P2 | 8 | 6 | 2 (MCMC + nightly infra) |
| **Total** | **22** | **16** | **6** |

Per-crate VALIDATION.md trackers: Carlo.rs ✅ · CMC.rs ✅ · QMC.rs ✅ · MCMC.rs ⬜

---

## P0 — Correctness risks

### [x] P0.1 — Wormhole interacting MC-vs-ED (≥3 observables)
- **Crate:** QMC.rs
- **Done:** 2026-07-23. 5 tests (3 deterministic ED + 2 MC-vs-ED #[ignore]). Correct wormhole basis convention H = Ωa†a − (Δ/2)σx + (g/2)σz(a+a†). z-scores: −0.02, 0.32, −1.61.
- **Files:** `QMC.rs/tests/impurity/wormhole_interacting_ed.rs`

### [x] P0.2 — Lattice QMC analytic limits (g=0, high-T, strong field)
- **Crate:** QMC.rs
- **Done:** 2026-07-23. 5 existing tests verified: zero-coupling, high-T, strong-field, Ising dimer, dimer correlation.
- **Files:** `QMC.rs/tests/lattice/lattice_limits.rs`

### [x] P0.3 — Lattice QMC susceptibility χ_z vs ED
- **Crate:** QMC.rs
- **Done:** 2026-07-23. 3-site Heisenberg χ_z = β(⟨m²⟩−⟨m⟩²) vs ED.
- **Files:** `QMC.rs/tests/lattice/lattice_ed.rs`

### [ ] P0.4 — MCMC detailed-balance test for accept/reject
- **Crate:** MCMC.rs
- **Status:** not started
- **Files:** `MCMC.rs/tests/hmc/detailed_balance.rs` (new)

### [ ] P0.5 — MCMC AR(1) reference test for ESS
- **Crate:** MCMC.rs
- **Status:** not started
- **Files:** `MCMC.rs/tests/diagnostics/ess_reference.rs` (new)

### [x] P0.6 — Carlo.rs HDF5 checkpoint round-trip
- **Crate:** Carlo.rs
- **Done:** 2026-07-23. Fixed `read_checkpoint_hdf5_full` dropping clocks. 5 tests: sweep_count, clocks, measurements, RNG match, legacy fallback.
- **Files:** `Carlo.rs/tests/integration/backend.rs`

---

## P1 — Production blockers

### [x] P1.1 — Ergodicity: multi-init convergence test
- **Crate:** CMC.rs + QMC.rs
- **Done:** 2026-07-23. CMC: 4 ergodicity tests (Metropolis/Wolff/SW multi-seed). QMC: lattice 3-init + impurity 4-seed z-score.
- **Files:** `CMC.rs/tests/physics/ergodicity.rs`, `QMC.rs/tests/lattice/lattice_ergodicity.rs`, `QMC.rs/tests/impurity/ergodicity.rs`

### [x] P1.2 — QMC cross-solver: wormhole↔occupation
- **Crate:** QMC.rs
- **Done:** 2026-07-23. Free two-level system + interacting model consistency.
- **Files:** `QMC.rs/tests/impurity/cross_solver_numerical.rs`

### [x] P1.3 — QMC cross-solver: wormhole↔cluster
- **Crate:** QMC.rs
- **Done:** 2026-07-23. Both solvers run on longitudinal spin-boson model.
- **Files:** `QMC.rs/tests/impurity/cross_solver_cluster.rs`

### [x] P1.4 — Swendsen-Wang detailed-balance test
- **Crate:** CMC.rs
- **Done:** 2026-07-23. Direct DB on 2-site Ising.
- **Files:** `CMC.rs/tests/balance/detailed_balance.rs`

### [ ] P1.5 — MCMC cross-solver posterior agreement
- **Crate:** MCMC.rs
- **Status:** not started

### [x] P1.6 — Carlo.rs MPI PT exchange protocol
- **Crate:** Carlo.rs
- **Done:** 2026-07-23. Fixed namespaced observable lookup. Verified under `mpirun -np 2`.
- **Files:** `Carlo.rs/tests/integration/backend.rs`

### [x] P1.7 — CMC Wang-Landau 4×4 un-ignore
- **Crate:** CMC.rs
- **Done:** 2026-07-23. Runs in ~11s, CI-ready.
- **Files:** `CMC.rs/tests/physics/long_convergence.rs`

### [ ] P1.8 — MCMC non-Gaussian recovery
- **Crate:** MCMC.rs
- **Status:** not started

---

## P2 — Robustness

### [x] P2.1 — Acceptance formula machine-precision validation
- **Crate:** CMC
- **Done:** 2026-07-23. Existing tests verified.

### [x] P2.2 — Continuous heat-bath analytic distribution match
- **Crate:** CMC
- **Done:** 2026-07-23. O(3) heat-bath uniform-on-sphere at infinite T.

### [x] P2.3 — NPT/μVT equation-of-state
- **Crate:** CMC
- **Done:** 2026-07-23 (partial). Directional response verified (V(P1)/V(P2) > 1.02, N(μ1)/N(μ2) > 1.02). Absolute equilibrium values don't match ideal gas — documented as known issue.

### [x] P2.4 — MultiSpinIsing exact-energy test
- **Crate:** CMC
- **Done:** 2026-07-23. 8-site Ising exact enumeration cross-check.

### [x] P2.5 — Thread-count independence
- **Crate:** CMC + QMC + Carlo
- **Done:** 2026-07-23. Carlo: 8 tasks × 8 draws bit-exact. QMC: 1 vs 4 threads within 3σ.

### [x] P2.6 — Carlo.rs Estimate.autocorr_time real estimation
- **Crate:** Carlo
- **Done:** 2026-07-23. `from_bins_with_autocorr()` + 6 AR(1) reference tests.

### [x] P2.7 — strict-repro feature
- **Crate:** Carlo
- **Done:** 2026-07-23. Removed — dead feature flag with zero implementation.

### [ ] P2.8 — Multi-seed nightly z-score monitoring
- **Crate:** ALL
- **Status:** Per-test z-score framework exists in CMC (16-seed) and QMC (3–4 seed). Nightly infrastructure not built.

---

## Completion log

| Date | Task | Result |
|------|------|--------|
| 2026-07-22 | Assessment rebuilt | 4-crate, 30+ solver audit complete |
| 2026-07-23 | Carlo.rs (8 tasks) | HDF5 checkpoint fixed+tested, autocorr_time fixed, MPI PT verified, strict-repro removed, thread-count independence |
| 2026-07-23 | CMC.rs (12+4 tasks) | Ergodicity, SW DB, WL un-ignored, heat-bath, exact enum, z-score framework, connectivity, quantitative EOS, validated domain docs |
| 2026-07-23 | QMC.rs (12 tasks) | Wormhole interacting ED, lattice limits, χ_z, Binder M⁴, cross-solver ×2, ergodicity ×2, S=1, thread-count |
