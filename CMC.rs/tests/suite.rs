//! CMC.rs integration test suite.
//!
//! Each module lives in a subdirectory and is wired via explicit `#[path]`
//! attributes so the directory tree mirrors the test taxonomy without
//! requiring `mod.rs` files.

// ── Algebra — deterministic math correctness ─────────────────────────────

#[path = "algebra/model_smoke.rs"]
mod model_smoke;
#[path = "algebra/sampling_core.rs"]
mod sampling_core;
#[path = "algebra/statistical_correctness.rs"]
mod statistical_correctness;

// ── Balance — Markov chain detailed-balance checks ───────────────────────

#[path = "balance/detailed_balance.rs"]
mod detailed_balance;

// ── Exact — exact enumeration / analytical reference solutions ───────────

#[path = "exact/classic_models.rs"]
mod classic_models;
#[path = "exact/onsager.rs"]
mod onsager;

// ── Integration — algorithm-level Monte Carlo tests ──────────────────────

#[path = "integration/checkpoint.rs"]
mod checkpoint;
#[path = "integration/dynamics_stage6.rs"]
mod dynamics_stage6;
#[path = "integration/generalized_stage4.rs"]
mod generalized_stage4;
#[path = "integration/particle_core.rs"]
mod particle_core;
#[path = "integration/particle_metropolis.rs"]
mod particle_metropolis;
#[path = "integration/particle_stage3.rs"]
mod particle_stage3;
#[path = "integration/usage.rs"]
mod usage;
#[path = "integration/worm_stage5.rs"]
mod worm_stage5;

// ── Physics — deterministic exact validation (was physics_validation.rs) ─

#[path = "physics/bkl_gillespie_equilibrium.rs"]
mod bkl_gillespie_equilibrium;
#[path = "physics/common.rs"]
mod common;
#[path = "physics/connectivity.rs"]
mod connectivity;
#[path = "physics/continuous_cross_solver.rs"]
mod continuous_cross_solver;
#[path = "physics/continuous_spins.rs"]
mod continuous_spins;
#[path = "physics/dynamics_exact.rs"]
mod dynamics_exact;
#[path = "physics/ergodicity.rs"]
mod ergodicity;
#[path = "physics/ergodicity_extended.rs"]
mod ergodicity_extended;
#[path = "physics/event_chain_eos.rs"]
mod event_chain_eos;
#[path = "physics/generalized_exact.rs"]
mod generalized_exact;
#[path = "physics/heat_bath_exact.rs"]
mod heat_bath_exact;
#[path = "physics/kawasaki_exact.rs"]
mod kawasaki_exact;
#[path = "physics/lattice_exact.rs"]
mod lattice_exact;
#[path = "physics/long_convergence.rs"]
mod long_convergence;
#[path = "physics/molecule_equilibrium.rs"]
mod molecule_equilibrium;
#[path = "physics/molecule_external_field.rs"]
mod molecule_external_field;
#[path = "physics/multicanonical_exact.rs"]
mod multicanonical_exact;
#[path = "physics/observables_exact.rs"]
mod observables_exact;
#[path = "physics/p2_remaining.rs"]
mod p2_remaining;
#[path = "physics/p2_validation.rs"]
mod p2_validation;
#[path = "physics/particle_zscore.rs"]
mod particle_zscore;
#[path = "physics/particles_exact.rs"]
mod particles_exact;
#[path = "physics/statistical_regression.rs"]
mod statistical_regression;
#[path = "physics/transition_balance.rs"]
mod transition_balance;
#[path = "physics/usage_exact.rs"]
mod usage_exact;
#[path = "physics/wang_landau_binned.rs"]
mod wang_landau_binned;
#[path = "physics/worm_exact.rs"]
mod worm_exact;
#[path = "physics/zscore_extended.rs"]
mod zscore_extended;
#[path = "physics/zscore_validation.rs"]
mod zscore_validation;
