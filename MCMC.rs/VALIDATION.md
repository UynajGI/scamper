# MCMC.rs — Physics Validation Task Tracker

> Created 2026-07-22. Branch: `dev`.
> Baseline: `../MATURITY_ASSESSMENT.md`.

## Current status

**NUTS/HMC: research-grade for Gaussian targets. RWM/ComponentWise/Slice: experimental.**

44 test functions (3 inline + 41 integration). Leapfrog reversibility, energy conservation, U-turn detection, gradient correctness, distribution recovery (Gaussian only), checkpoint reproducibility, input validation — all solid. But detailed balance, ESS calibration, cross-solver agreement, and non-Gaussian recovery are missing.

## Tasks

### [ ] MC-P0.1 — Detailed-balance test for Metropolis accept/reject
- **Problem:** No test verifies π(x)P(x→y) = π(y)P(y→x) for any kernel. Leapfrog reversibility ≠ detailed balance (the accept/reject step introduces asymmetry that could be wrong).
- **Plan:** Two approaches:
  - **Analytic:** Small 1D target on a grid (e.g., 5-point {0.1, 0.2, 0.4, 0.2, 0.1}). For RWM with known proposal width: enumerate transition probabilities analytically, verify P(x→y)π(x) = P(y→x)π(y) at machine precision.
  - **Empirical:** Run RWM 100k steps from each state. Measure empirical transition frequency f(x→y). Verify f(x→y)/f(y→x) ≈ π(y)/π(x) within statistical error.
  - Repeat for HMC accept step.
- **File:** `MCMC.rs/tests/hmc/detailed_balance.rs` (new)
- **Status:** not started

### [ ] MC-P0.2 — AR(1) reference test for ESS estimator
- **Problem:** ESS (`effective_sample_size`) only validated against IID data (ess≈N for uncorrelated). The estimator's behavior on correlated data is untested. AR(1) process with parameter ρ has exact integrated autocorrelation time τ_int = (1+ρ)/(1−ρ), so ESS = N/τ_int.
- **Plan:** Generate AR(1) series: x_{t+1} = ρ·x_t + √(1−ρ²)·ε_t, ε~N(0,1). For ρ ∈ {0.0, 0.3, 0.5, 0.7, 0.9}: feed N=10000 samples to `effective_sample_size`. Assert ESS ≈ N(1−ρ)/(1+ρ) within 10%.
- **File:** `MCMC.rs/tests/diagnostics/ess_reference.rs` (new)
- **Status:** not started

### [ ] MC-P0.3 — R-hat reference test on known non-converged chains
- **Problem:** R-hat only tested qualitatively (IID → rhat<1.02, shifted → rhat>1.05). No closed-form reference for rank-normalized R-hat.
- **Plan:** Two chains from N(0,1) and N(μ,1) with known separation. Compute theoretical R-hat (depends on between/within variance ratio). Compare to `rhat()` output.
- **File:** `MCMC.rs/tests/diagnostics/convergence.rs` (extend)
- **Status:** not started

### [ ] MC-P1.1 — Cross-solver posterior agreement
- **Problem:** NUTS, HMC, RWM, Slice each recover N(0,1) independently. Never compared to each other on the same target.
- **Plan:** Correlated Gaussian (ρ=0.8, dim=2). Run all 4 samplers: 10k draws each. Assert recovered mean and covariance agree within 3σ across all pairs. This catches systematic bias in any single sampler.
- **File:** `MCMC.rs/tests/kernels/cross_solver.rs` (new)
- **Status:** not started

### [ ] MC-P1.2 — Non-Gaussian recovery: bimodal mixture
- **Problem:** All targets are Gaussian. Bimodal targets expose tunneling failures.
- **Plan:** 0.5·N(−3,1) + 0.5·N(3,1). Run NUTS with default warmup. Assert recovered mean ≈ 0 (both modes visited). Check tail-ESS.
- **File:** `MCMC.rs/tests/kernels/non_gaussian.rs` (new)
- **Status:** not started

### [ ] MC-P1.3 — Non-Gaussian recovery: Student-t heavy tails
- **Problem:** Gaussian targets don't stress tail-ESS.
- **Plan:** Student-t with ν=3 (finite variance, heavy tail). Run NUTS. Assert ⟨x²⟩ ≈ ν/(ν−2) = 3. Check tail-ESS > 100 for 10k draws.
- **File:** `MCMC.rs/tests/kernels/non_gaussian.rs` (extend)
- **Status:** not started

### [ ] MC-P1.4 — ComponentWiseMetropolis distribution recovery
- **Problem:** No test verifies ComponentWiseMetropolis recovers the correct distribution. Only iteration-advancement and NaN-rejection tested.
- **Plan:** Correlated 2D Gaussian. Run CWM with adaptation. Assert ⟨x⟩≈0, ⟨y⟩≈0, Var[x]≈1, Var[y]≈1, Cov[x,y]≈0.8 within 0.15.
- **File:** `MCMC.rs/tests/diagnostics/gaussian_moments.rs` (extend)
- **Status:** not started

### [ ] MC-P1.5 — Replica exchange: exchange acceptance ratio vs theory
- **Problem:** PT exchange validated only for trace recording and reproducibility. Exchange acceptance ratio never compared to theoretical value.
- **Plan:** 3-replica ladder, geometric β spacing. Run. Compare empirical exchange acceptance rate to theoretical: P_acc = min(1, exp(Δβ·ΔE)).
- **File:** `MCMC.rs/tests/tempering/replica_exchange.rs` (extend)
- **Status:** not started

### [ ] MC-P2.1 — Gibbs sampler: multivariate recovery
- **Problem:** Gibbs tested for exact conditional draw + atomicity, but not for full distribution recovery on a multivariate target.
- **Plan:** 3D correlated Gaussian. Gibbs update (sample x|y,z; y|x,z; z|x,y). Assert covariance recovery.
- **File:** `MCMC.rs/tests/kernels/composition_gibbs.rs` (extend)
- **Status:** not started

### [ ] MC-P2.2 — Autocorrelation time: public API exposure
- **Problem:** Autocorrelation is internal to `ess.rs`, not publicly exposed. Users can't compute ACF directly.
- **Plan:** Expose `autocovariance(series, max_lag)` as public function. Test against AR(1): ρ(τ) = ρ^τ for AR(1).
- **File:** `MCMC.rs/src/diagnostics/mod.rs` (extend)
- **Status:** not started

### [ ] MC-P2.3 — Proposal ratio machine-precision validation
- **Problem:** Acceptance probability computed in log space but never tested at machine precision against hand-computed values.
- **Plan:** Known log-target values: construct artificial log-density ratio. Feed to `accept_log_probability`. Assert exact acceptance probability.
- **File:** `MCMC.rs/tests/hmc/acceptance_formula.rs` (new)
- **Status:** not started

### [ ] MC-P2.4 — NUTS multinomial sampling validity
- **Problem:** NUTS uses multinomial sampling for trajectory selection. Never tested that the multinomial weights are proportional to exp(−H).
- **Plan:** Small trajectory (4 leaf nodes). Compute weights. Verify multinomial selection probabilities match exp(−H_i)/Σexp(−H_j).
- **File:** `MCMC.rs/tests/hmc/nuts_weights.rs` (new)
- **Status:** not started

## Completion log

| Date | Task | Result |
|------|------|--------|
