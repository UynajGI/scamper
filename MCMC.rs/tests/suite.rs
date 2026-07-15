//! MCMC.rs integration test suite.
//!
//! Each module lives in a subdirectory and is wired via explicit `#[path]`
//! attributes so the directory tree mirrors the test taxonomy without
//! requiring `mod.rs` files.

// ── Kernels — transition kernel correctness ──────────────────────────────

#[path = "kernels/composition_gibbs.rs"]
mod composition_gibbs;
#[path = "kernels/kernel_state_invariants.rs"]
mod kernel_state_invariants;
#[path = "kernels/slice_sampling.rs"]
mod slice_sampling;
#[path = "kernels/transforms.rs"]
mod transforms;

// ── HMC — Hamiltonian Monte Carlo (static + dynamic) ─────────────────────

#[path = "hmc/dynamic_hmc.rs"]
mod dynamic_hmc;
#[path = "hmc/energy_diagnostics.rs"]
mod energy_diagnostics;
#[path = "hmc/gradient_check.rs"]
mod gradient_check;
#[path = "hmc/metrics_integrator.rs"]
mod metrics_integrator;
#[path = "hmc/nuts.rs"]
mod nuts;
#[path = "hmc/static_hmc.rs"]
mod static_hmc;
#[path = "hmc/warmup.rs"]
mod warmup;

// ── Transforms — constrained-to-unconstrained bijections ─────────────────

#[path = "transforms/transform_gradients.rs"]
mod transform_gradients;

// ── Covariance — proposal adaptation ─────────────────────────────────────

#[path = "covariance/dense_covariance.rs"]
mod dense_covariance;

// ── Tempering — replica exchange ─────────────────────────────────────────

#[path = "tempering/replica_exchange.rs"]
mod replica_exchange;

// ── Adaptation ───────────────────────────────────────────────────────────

#[path = "adaptation/adaptation_freeze.rs"]
mod adaptation_freeze;

// ── Diagnostics — convergence, moments ───────────────────────────────────

#[path = "diagnostics/convergence.rs"]
mod convergence;
#[path = "diagnostics/gaussian_moments.rs"]
mod gaussian_moments;

// ── Integration — Carlo.rs adapter, checkpoint, multi-chain ──────────────

#[path = "integration/carlo_adapter.rs"]
mod carlo_adapter;
#[path = "integration/checkpoint_equivalence.rs"]
mod checkpoint_equivalence;
#[path = "integration/multichain_reproducibility.rs"]
mod multichain_reproducibility;
