//! Deterministic and exact physical validation for CMC.rs.
//!
//! Default tests use finite-state enumeration, analytical identities or
//! transaction invariants. Long stochastic convergence checks live in
//! `physics/statistical_regression.rs` and are ignored by default.

#[path = "physics/common.rs"]
mod common;
#[path = "physics/continuous_spins.rs"]
mod continuous_spins;
#[path = "physics/dynamics_exact.rs"]
mod dynamics_exact;
#[path = "physics/generalized_exact.rs"]
mod generalized_exact;
#[path = "physics/lattice_exact.rs"]
mod lattice_exact;
#[path = "physics/observables_exact.rs"]
mod observables_exact;
#[path = "physics/particles_exact.rs"]
mod particles_exact;
#[path = "physics/statistical_regression.rs"]
mod statistical_regression;
#[path = "physics/transition_balance.rs"]
mod transition_balance;
#[path = "physics/worm_exact.rs"]
mod worm_exact;
