# CMC.rs — Physics Validation Task Tracker

> Created 2026-07-22. Branch: `dev`.
> Baseline: `../MATURITY_ASSESSMENT.md`.

## Current status

**Metropolis/Wolff/Wang-Landau/Worm: research-grade with strong DB+exact coverage.**
**Heat-bath(continuous)/Hybrid/MultiSpin/Umbrella/Gillespie/EventChain: experimental.**

~110 test functions. Detailed balance is directly tested for Metropolis, Wolff, batch moves, kinetic rates. Exact comparisons include Onsager, finite-Ising enumeration, Langevin, Bessel-ratio, exact DOS. Cross-solver agreement (Metropolis vs Wolff vs SW) on 8×8 Ising.

## Tasks

### [ ] CMC-P0.1 — Ergodicity: multi-init convergence for all lattice solvers
- **Problem:** No test checks convergence from different initial states. Could hide sector-trapping or non-ergodic bugs.
- **Plan:** For Metropolis, Wolff, SW, heat-bath: start from (a) all-up ordered, (b) all-down, (c) random. Run 10k sweeps on 4×4 Ising at Tc. Compare ⟨E⟩ and ⟨m²⟩ across init conditions within 3σ.
- **File:** `CMC.rs/tests/physics/ergodicity.rs` (new)
- **Status:** not started

### [ ] CMC-P0.2 — Continuous heat-bath: infinite-T uniform distribution
- **Problem:** O(N) continuous heat bath only has "energy in physical range" smoke test. No analytic distribution match.
- **Plan:** β→0 (infinite T). O(3) spins. Assert ⟨m_x⟩≈0, ⟨m_y⟩≈0, ⟨m_z⟩≈0, ⟨m²⟩≈1/3 per site. Compare to analytic uniform-on-sphere distribution.
- **File:** `CMC.rs/tests/exact/classic_models.rs` (extend)
- **Status:** not started

### [ ] CMC-P1.1 — Swendsen-Wang detailed-balance test
- **Problem:** SW has no direct DB test. Only cross-solver agreement on 8×8. An error in SW's FK bond construction could be masked by Wolff's correctness.
- **Plan:** Enumerate SW transitions on 2-site and 3-site Ising. Verify FK bond activation probability P(bond active | s_i, s_j) = 1 − exp(−2βJ) for aligned spins, and that cluster flip preserves DB.
- **File:** `CMC.rs/tests/balance/detailed_balance.rs` (extend)
- **Status:** not started

### [ ] CMC-P1.2 — Wang-Landau 4×4 DOS: evaluate for un-#[ignore]
- **Problem:** WL 4×4 Ising DOS comparison is `#[ignore]`. If fast enough, should be in CI.
- **Plan:** Time the test. If <10s, remove `#[ignore]`. If longer, optimize (reduce refinement sweeps) or keep ignored but document runtime.
- **File:** `CMC.rs/tests/physics/long_convergence.rs`
- **Status:** not started

### [ ] CMC-P1.3 — Multicanonical/umbrella MC-vs-exact distribution
- **Problem:** `EnergyBiasCore` only has transactional-rejection and bias-algebra tests. No test verifies sampling actually reproduces the target distribution.
- **Plan:** 4-site Ising. Apply known bias → sample → reweight → compare reweighted ⟨E⟩ to exact canonical ⟨E⟩ at several temperatures.
- **File:** `CMC.rs/tests/generalized/umbrella_exact.rs` (new)
- **Status:** not started

### [ ] CMC-P1.4 — Wang-Landau continuous-axis production test
- **Problem:** `BinnedAxis` WL has no physics run test (only state-machine unit test).
- **Plan:** Small system with continuous energy axis. Verify WL converges to flat histogram, then 1/t refinement reduces DOS error.
- **File:** `CMC.rs/tests/generalized/wl_continuous.rs` (new)
- **Status:** not started

### [ ] CMC-P2.1 — MultiSpinIsing exact-energy test
- **Problem:** Bit-packed 64-replica solver has no exact-energy comparison.
- **Plan:** 8-spin Ising. Run MultiSpinIsing. Compare ⟨E⟩ to exact enumeration (256 states).
- **File:** `CMC.rs/tests/exact/multispin_exact.rs` (new)
- **Status:** not started

### [ ] CMC-P2.2 — HybridCore correctness
- **Problem:** Only smoke-tested (alternation/repetitions). No DB or exact comparison.
- **Plan:** Hybrid(Metropolis, Wolff) on 4×4 Ising. Verify ⟨E⟩ matches exact within statistical error. Verify both updates contribute.
- **File:** `CMC.rs/tests/balance/hybrid_correctness.rs` (new)
- **Status:** not started

### [ ] CMC-P2.3 — NPT equation-of-state
- **Problem:** LJ NPT has Jacobian/pressure terms tested algebraically but no equilibrium comparison.
- **Plan:** LJ at T=1.0, P=0.1, N=32. Compare ⟨ρ⟩ to known values (e.g., Johnson et al. tables).
- **File:** `CMC.rs/tests/physics/long_convergence.rs` (extend, #[ignore])
- **Status:** not started

### [ ] CMC-P2.4 — μVT interacting: particle number vs grand-canonical
- **Problem:** μVT only tested for ideal gas (Poisson). No interacting test.
- **Plan:** LJ μVT at small N (≤10). Compare ⟨N⟩ to NVT at same chemical potential.
- **File:** `CMC.rs/tests/physics/long_convergence.rs` (extend, #[ignore])
- **Status:** not started

### [ ] CMC-P2.5 — Gillespie equilibrium distribution
- **Problem:** Only 2-event toy model tested.
- **Plan:** Simple 3-state model with known rates. Run Gillespie. Compare occupation probabilities to stationary distribution π_i = r_i / Σr.
- **File:** `CMC.rs/tests/dynamics/gillespie_exact.rs` (new)
- **Status:** not started

### [ ] CMC-P2.6 — Event-chain pressure comparison
- **Problem:** Event chain only tested for collision geometry. No thermodynamic comparison.
- **Plan:** Hard-sphere NVT at η=0.3. Compare P to Carnahan-Starling or Metropolis result.
- **File:** `CMC.rs/tests/physics/long_convergence.rs` (extend, #[ignore])
- **Status:** not started

## Completion log

| Date | Task | Result |
|------|------|--------|
