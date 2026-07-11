//! # QMC.rs — Quantum Monte Carlo algorithm toolbox
//!
//! QMC.rs is the physics layer built on top of Carlo.rs.  Carlo.rs owns
//! scheduling, random-number seeding, accumulation, error analysis,
//! checkpoint orchestration, and parallel execution.  QMC.rs owns quantum
//! representations, model catalogs, update kernels, and estimators.
//!
//! ## Architecture
//!
//! - [`algorithm`] — reusable update-kernel contracts and fixed sweep schedules
//! - [`spin_boson`] — continuous-time retarded-interaction wormhole QMC
//! - [`worldline`] — reusable continuous and discrete worldline containers
//! - [`discrete`] — discrete-time path-integral update prototypes
//! - [`hamiltonian`] / [`lattice`] — lattice-model foundations
//! - [`heisenberg_chain`] — existing Carlo.rs adapter for the Heisenberg chain

pub mod algorithm;
pub mod discrete;
pub mod hamiltonian;
pub mod heisenberg_chain;
pub mod lattice;
pub mod spin_boson;
pub mod worldline;

pub use algorithm::{QmcKernel, UpdateSchedule};
pub use discrete::{local_metropolis_sweep, worm_sweep, SpaceTimeConfig, Spin as DiscreteSpin};
pub use hamiltonian::{
    heisenberg_chain_ground_energy_per_site, HeisenbergChain, QuantumHamiltonian,
};
pub use heisenberg_chain::HeisenbergChainMC;
pub use lattice::ChainLattice;
pub use spin_boson::{
    Bath, PowerLawBath, SingleModeBath, SpinBosonModel, SpinBosonModelKind, SpinBosonQmc,
    TabulatedBath, WormholeConfiguration, WormholeEngine,
};
