# MCMC.rs — Physics Validation Task Tracker

> Created 2026-07-22. Branch: `dev`.
> Baseline: `../MATURITY_ASSESSMENT.md`.

## Current status

**All 6 kernels research-grade (2026-08-14). Remaining fine-grained items: MC-P0.3, MC-P1.5, MC-P2.2, MC-P2.4.**

69 test functions. Leapfrog reversibility, energy conservation, U-turn detection, gradient correctness, checkpoint reproducibility, input validation; detailed balance (machine-precision + statistical), ESS calibration on AR(1), 6-solver posterior agreement, and non-Gaussian recovery (bimodal + Student-t) all landed 2026-08-14.

## Tasks

### [x] MC-P0.1 — Detailed-balance test for Metropolis accept/reject
- **Done:** 2026-08-14. Machine precision: the kernel's reported acceptance statistic is bit-for-bit `exp(min(0, Δlog p))` (tied to the real code path via an instrumented target), and the log-Metropolis rule satisfies p(x)A(x,y) = p(y)A(y,x) to round-off incl. ulp-scale Δ. Empirical: 8-bucket transition counting over 300k draws, flow symmetry n_ij ≈ n_ji within a Poisson envelope + occupancy vs Simpson quadrature, for RandomWalk and ComponentWise kernels.
- **File:** `MCMC.rs/tests/diagnostics/detailed_balance.rs`

### [x] MC-P0.2 — AR(1) reference test for ESS estimator
- **Done:** 2026-08-14. AR(1) chains ρ ∈ {0.0, 0.5, 0.9, 0.99} (4 chains, burn-in 1000) through `diagnose()`; ess_bulk vs closed-form ESS/N = (1−ρ)/(1+ρ), measured tolerances 5/5/15/30%. Note: rank-normalization is an identity in distribution for the exactly-N(0,1) AR(1) marginal, so the closed form carries over.
- **File:** `MCMC.rs/tests/diagnostics/ess_ar1.rs`

### [ ] MC-P0.3 — R-hat reference test on known non-converged chains
- **Problem:** R-hat only tested qualitatively (IID → rhat<1.02, shifted → rhat>1.05). No closed-form reference for rank-normalized R-hat.
- **Plan:** Two chains from N(0,1) and N(μ,1) with known separation. Compute theoretical R-hat (depends on between/within variance ratio). Compare to `rhat()` output.
- **File:** `MCMC.rs/tests/diagnostics/convergence.rs` (extend)
- **Status:** not started

### [x] MC-P1.1 — Cross-solver posterior agreement
- **Done:** 2026-08-14. Six solvers (RW-Metropolis, ComponentWise, Slice, StaticHMC, NUTS, Gibbs) on the same correlated 2D Gaussian (Σ=[[2.0,0.6],[0.6,1.0]], μ=(1,−2)): pairwise agreement 15 pairs × 5 moments, |z| < 4, MC errors from 8 independent seeds; each solver also vs analytic moments. `#[ignore]` long-run variant (12 seeds, 3× chains) for nightly.
- **File:** `MCMC.rs/tests/integration/cross_solver.rs`

### [x] MC-P1.2 — Non-Gaussian recovery: bimodal mixture
- **Done:** 2026-08-14. 0.5·N(−1.5, 0.6²) + 0.5·N(1.5, 0.6²) (mean 0, var 2.61): 4 solvers recover moments within |z| < 4 + mode-occupancy symmetry; 16 seeds alternating mode inits; `#[ignore]` long variant. Note: the originally planned deeper valley (e.g. ±3/σ=1 → valley ≈ e⁻⁴·⁵ below the peaks) sits in the rare-tunneling regime where finite-chain variance carries a −Var(chain mean)/2 bias — ±1.5/σ=0.6 keeps a ≈23× deep valley while all solvers genuinely tunnel.
- **File:** `MCMC.rs/tests/covariance/non_gaussian.rs`

### [x] MC-P1.3 — Non-Gaussian recovery: Student-t heavy tails
- **Done:** 2026-08-14. 2D Student-t ν=5, Σ₀=[[1,0.3],[0.3,0.5]], μ=(0.5,−1): 4 solvers recover mean + 3 cov elements (analytic ν/(ν−2)·Σ₀) within |z| < 4; 8 seeds; `#[ignore]` long variant.
- **File:** `MCMC.rs/tests/covariance/non_gaussian.rs`

### [x] MC-P1.4 — ComponentWiseMetropolis distribution recovery
- **Done:** 2026-08-14 (via cross-solver + non-Gaussian suites). ComponentWise is one of the six solvers in `cross_solver.rs` (correlated Gaussian moments vs analytic, |z| < 4) and one of the four in `non_gaussian.rs`.
- **File:** `MCMC.rs/tests/integration/cross_solver.rs`

### [ ] MC-P1.5 — Replica exchange: exchange acceptance ratio vs theory
- **Problem:** PT exchange validated only for trace recording and reproducibility. Exchange acceptance ratio never compared to theoretical value.
- **Plan:** 3-replica ladder, geometric β spacing. Run. Compare empirical exchange acceptance rate to theoretical: P_acc = min(1, exp(Δβ·ΔE)).
- **File:** `MCMC.rs/tests/tempering/replica_exchange.rs` (extend)
- **Status:** not started

### [x] MC-P2.1 — Gibbs sampler: multivariate recovery
- **Done:** 2026-08-14 (via cross-solver). Gibbs (two exact-Gaussian-conditional `GibbsKernel`s composed with `Then`) is one of the six solvers agreeing on the correlated 2D Gaussian — moment recovery vs analytic and vs 5 other solvers, |z| < 4. (dim 2, not the planned dim 3.)
- **File:** `MCMC.rs/tests/integration/cross_solver.rs`

### [ ] MC-P2.2 — Autocorrelation time: public API exposure
- **Problem:** Autocorrelation is internal to `ess.rs`, not publicly exposed. Users can't compute ACF directly.
- **Plan:** Expose `autocovariance(series, max_lag)` as public function. Test against AR(1): ρ(τ) = ρ^τ for AR(1).
- **File:** `MCMC.rs/src/diagnostics/mod.rs` (extend)
- **Status:** not started

### [x] MC-P2.3 — Proposal ratio machine-precision validation
- **Done:** 2026-08-14 (as part of MC-P0.1). `acceptance_statistic_follows_log_metropolis_formula` asserts the reported log-acceptance equals `min(0, Δlog p)` bit-for-bit over 4000 live proposals (uphill accepted deterministically without a uniform draw; rejections leave the state untouched).
- **File:** `MCMC.rs/tests/diagnostics/detailed_balance.rs`

### [ ] MC-P2.4 — NUTS multinomial sampling validity
- **Problem:** NUTS uses multinomial sampling for trajectory selection. Never tested that the multinomial weights are proportional to exp(−H).
- **Plan:** Small trajectory (4 leaf nodes). Compute weights. Verify multinomial selection probabilities match exp(−H_i)/Σexp(−H_j).
- **File:** `MCMC.rs/tests/hmc/nuts_weights.rs` (new)
- **Status:** not started

## Completion log

| Date | Task | Result |
|------|------|--------|
| 2026-08-14 | MC-P0.1 + MC-P2.3 | ✅ Machine-precision log-Metropolis rule (bit-for-bit, tied to live code path) + binned empirical flow balance (RW + ComponentWise) |
| 2026-08-14 | MC-P0.2 | ✅ AR(1) ESS calibration, ρ ∈ {0, 0.5, 0.9, 0.99}, tiered tolerances |
| 2026-08-14 | MC-P1.1 | ✅ 6-solver posterior agreement (15 pairs × 5 moments, \|z\| < 4) + long variant |
| 2026-08-14 | MC-P1.2/P1.3/P1.4 | ✅ Bimodal mixture + Student-t ν=5 recovery, 4 solvers, \|z\| < 4 + long variants |
| 2026-08-14 | MC-P2.1 | ✅ Gibbs multivariate recovery via cross-solver suite |
