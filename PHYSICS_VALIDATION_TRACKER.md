# Physics Validation Task Tracker

> Created 2026-07-22 on `dev` branch.
> Goal: bring every solver from "compiles and runs" to "physics verified".
> Assessment baseline: `MATURITY_ASSESSMENT.md`.

---

## Progress overview

| Tier | Total | Done | In Progress | Blocked |
|------|-------|------|-------------|---------|
| P0 | 6 | 0 | 0 | 0 |
| P1 | 8 | 0 | 0 | 0 |
| P2 | 8 | 0 | 0 | 0 |

---

## P0 — Correctness risks

These could mask physics bugs right now.

### [ ] P0.1 — Wormhole interacting MC-vs-ED (≥3 observables)
- **Crate:** QMC.rs
- **Problem:** Wormhole only validated against free/analytic limits (Poisson expansion, 2-state tanh). No interacting test where MC-sampled ⟨σx⟩, ⟨σz⟩, C(τ) are compared to dense diagonalization of the same effective model.
- **Plan:** Build a small interacting single-mode Rabi model. ED with nalgebra Jacobi. Run WormholeEngine at β=20–50, 50k+100k sweeps. Compare ⟨σz⟩ (sampled = physical σx), expansion order (= −βE), and C(β/2). Tolerance: |z| < 4 per seed, 10 seeds.
- **Files:** `QMC.rs/tests/impurity/wormhole_interacting_ed.rs` (new)
- **Status:** not started

### [ ] P0.2 — Lattice QMC analytic limits (g=0, high-T, strong field)
- **Crate:** QMC.rs
- **Problem:** Lattice solver has ED cross-checks but no analytic-limit tests. CMC.rs has these (zero-coupling, high-T, strong field); QMC.rs lattice does not.
- **Plan:** Add 3 tests: (1) zero coupling → zero energy, zero vertices; (2) β=0.1 high-T → ⟨m²⟩ ≈ S(S+1)/3 per site; (3) strong longitudinal field → all spins polarize.
- **Files:** `QMC.rs/tests/lattice/lattice_limits.rs` (new)
- **Status:** not started

### [ ] P0.3 — Lattice QMC susceptibility χ_z vs ED
- **Crate:** QMC.rs
- **Problem:** Susceptibility χ_z = β(⟨m²⟩−⟨m⟩²) is measured but never compared to ED. The connected piece (⟨m²⟩−⟨m⟩²) involves a subtraction that could be wrong.
- **Plan:** Extend `lattice_ed.rs` 3-site Heisenberg test to also compare χ_z vs ED. Requires computing ⟨m²⟩_connected from MC and from density matrix.
- **Files:** `QMC.rs/tests/lattice/lattice_ed.rs` (extend)
- **Status:** not started

### [ ] P0.4 — MCMC detailed-balance test for accept/reject
- **Crate:** MCMC.rs
- **Problem:** No test verifies that the Metropolis accept/reject step satisfies detailed balance π(x)P(x→y) = π(y)P(y→x). Only leapfrog integrator reversibility is tested.
- **Plan:** Small 1D target (e.g., 2-state discrete or small Gaussian grid). Enumerate all transitions for RWM and HMC. Verify DB ratio at machine precision. For stochastic accept, verify empirical transition frequencies match π(x)P(x→y) within statistical error over 10k+ samples.
- **Files:** `MCMC.rs/tests/hmc/detailed_balance.rs` (new)
- **Status:** not started

### [ ] P0.5 — MCMC AR(1) reference test for ESS
- **Crate:** MCMC.rs
- **Problem:** ESS estimator validated only against IID heuristics (rhat<1.02, ess>3000). No closed-form reference. AR(1) process with ρ=0.5 has τ_int = (1+ρ)/(1−ρ) = 3, so ESS = N/3.
- **Plan:** Generate AR(1) series with known ρ, feed to `effective_sample_size`, assert ESS ≈ N(1−ρ)/(1+ρ) within 10%.
- **Files:** `MCMC.rs/tests/diagnostics/ess_reference.rs` (new)
- **Status:** not started

### [ ] P0.6 — Carlo.rs HDF5 checkpoint round-trip
- **Crate:** Carlo.rs
- **Problem:** Zero tests for HDF5 file-level checkpoint I/O. `read_checkpoint_hdf5_full` drops attempt/accepted/event_time clocks. RNG state deserialization untested.
- **Plan:** Create a Run, take HDF5 checkpoint, restore, verify sweep_count + RNG state + measurements match. Test behind `#[cfg(feature = "hdf5")]`.
- **Files:** `Carlo.rs/tests/io/checkpoint_hdf5.rs` (new)
- **Status:** not started

---

## P1 — Production blockers

### [ ] P1.1 — Ergodicity: multi-init convergence test
- **Crate:** ALL
- **Problem:** No test checks that different initial states converge to the same distribution. Could hide sector-trapping bugs.
- **Plan:** For each solver: start from ordered, disordered, random initial states. Run moderate sweeps. Compare ⟨E⟩ and ⟨M²⟩ across initial conditions within statistical error.
- **Files:** per-crate `tests/physics/ergodicity.rs`
- **Status:** not started

### [ ] P1.2 — QMC cross-solver: wormhole↔occupation
- **Crate:** QMC.rs
- **Problem:** `cross_solver.rs` documents 4 convention differences but has no numerical comparison. Need a shared reference both solvers can compare against independently.
- **Plan:** Use the free two-level system (g→0 limit) where both solvers reduce to the same physics. Compare ⟨σz⟩, ⟨n⟩=0, energy.
- **Files:** `QMC.rs/tests/impurity/cross_solver.rs` (extend)
- **Status:** not started

### [ ] P1.3 — QMC cross-solver: wormhole↔cluster
- **Crate:** QMC.rs
- **Problem:** Both solvers handle longitudinal spin-boson but via different algorithms. Never compared numerically.
- **Plan:** Single-mode bath, longitudinal coupling only (no transverse → both sign-free in same basis). Compare ⟨σz⟩, expansion order/kink count, C(β/2).
- **Files:** `QMC.rs/tests/impurity/cross_solver.rs` (extend)
- **Status:** not started

### [ ] P1.4 — Swendsen-Wang detailed-balance test
- **Crate:** CMC.rs
- **Problem:** SW has no direct DB test. Only Wolff has one.
- **Plan:** Enumerate SW transitions on 2–3 site Ising. Verify FK bond activation probabilities satisfy DB.
- **Files:** `CMC.rs/tests/balance/detailed_balance.rs` (extend)
- **Status:** not started

### [ ] P1.5 — MCMC cross-solver posterior agreement
- **Crate:** MCMC.rs
- **Problem:** Each sampler validated in isolation. Never compared on same posterior.
- **Plan:** Run NUTS, RWM, Slice on same correlated Gaussian. Assert recovered mean/covariance agree within statistical error (3–4σ).
- **Files:** `MCMC.rs/tests/kernels/cross_solver.rs` (new)
- **Status:** not started

### [ ] P1.6 — Carlo.rs MPI PT exchange protocol
- **Crate:** Carlo.rs
- **Problem:** MPI parallel-tempering exchange protocol untested (2 ignored smoke tests only).
- **Plan:** Test behind `#[cfg(feature = "mpi")] #[ignore]`. Verify exchange acceptance ratio, chain permutation, measurement synchronization.
- **Files:** `Carlo.rs/tests/mpi/pt_exchange.rs` (new)
- **Status:** not started

### [ ] P1.7 — CMC Wang-Landau 4×4 un-ignore
- **Crate:** CMC.rs
- **Problem:** WL 4×4 Ising DOS comparison is `#[ignore]`. Should be fast enough to run in CI (<10s).
- **Plan:** Evaluate runtime. If <10s, remove `#[ignore]`. If longer, optimize or reduce sweeps.
- **Files:** `CMC.rs/tests/physics/long_convergence.rs`
- **Status:** not started

### [ ] P1.8 — MCMC non-Gaussian recovery
- **Crate:** MCMC.rs
- **Problem:** All recovery targets are Gaussian. Tail behavior on heavy-tailed/multimodal targets untested.
- **Plan:** Add mixture-of-Gaussians (bimodal) and Student-t (heavy tail) recovery tests. Check mean recovery and tail-ESS.
- **Files:** `MCMC.rs/tests/kernels/non_gaussian.rs` (new)
- **Status:** not started

---

## P2 — Robustness

### [ ] P2.1 — Acceptance formula machine-precision validation
- **Crate:** CMC, QMC
- **Plan:** Test `accept_log_probability` and equivalent QMC acceptance ratios at machine precision against hand-computed values.

### [ ] P2.2 — Continuous heat-bath analytic distribution match
- **Crate:** CMC
- **Plan:** O(3) heat bath should produce uniform distribution on S² at infinite T. Test ⟨m⟩ = 0, ⟨m²⟩ = 1/3.

### [ ] P2.3 — NPT/μVT interacting equation-of-state
- **Crate:** CMC
- **Plan:** LJ fluid at known T,P → compare ⟨ρ⟩ to literature values.

### [ ] P2.4 — MultiSpinIsing exact-energy test
- **Crate:** CMC
- **Plan:** 8-spin Ising, compare ⟨E⟩ to exact enumeration.

### [ ] P2.5 — Thread-count independence
- **Crate:** ALL
- **Plan:** Run with 1 vs 4 threads. Assert statistical equivalence (not bitwise).

### [ ] P2.6 — Carlo.rs Estimate.autocorr_time real estimation
- **Crate:** Carlo
- **Plan:** Replace hardcoded τ=1.0 with actual autocorrelation estimate. Test against AR(1).

### [ ] P2.7 — strict-repro feature test
- **Crate:** Carlo
- **Plan:** Exercise the `strict-repro` RNG jump-sequence feature. Verify reproducibility across task counts.

### [ ] P2.8 — Multi-seed nightly z-score monitoring
- **Crate:** ALL
- **Plan:** Infrastructure for nightly runs: 20 seeds × key parameters → z-score report. Flag |z| > 4 trends.

---

## Completion log

| Date | Task | Result | Files |
|------|------|--------|-------|
| 2026-07-22 | Assessment rebuilt | 4-crate, 30+ solver audit complete | `MATURITY_ASSESSMENT.md` |
