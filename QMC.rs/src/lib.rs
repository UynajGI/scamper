//! # QMC.rs — Quantum Monte Carlo algorithm toolbox
//!
//! Built on [Carlo.rs](../carlo_rs/index.html) for scheduling, measurement accumulation,
//! and result analysis. QMC.rs provides the QMC-specific data structures and algorithms.
//!
//! ## Modules
//!
//! - [`hamiltonian`] — [`QuantumHamiltonian`] trait + concrete models (Heisenberg, …)
//! - [`lattice`] — minimal lattice topology ([`ChainLattice`])
//! - [`worldline`] — single-site worldline objects (path-integral / worm foundation)
//! - [`discrete`] — discrete-time Suzuki-Trotter worm algorithm

pub mod discrete;
pub mod hamiltonian;
pub mod heisenberg_chain;
pub mod lattice;
pub mod worldline;

pub use discrete::{local_metropolis_sweep, worm_sweep, SpaceTimeConfig, Spin};
pub use hamiltonian::{
    heisenberg_chain_ground_energy_per_site, HeisenbergChain, QuantumHamiltonian,
};
pub use heisenberg_chain::HeisenbergChainMC;
pub use lattice::ChainLattice;
