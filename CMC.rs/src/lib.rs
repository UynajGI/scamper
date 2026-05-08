//! CMC.rs - Classical Monte Carlo algorithm toolbox
//!
//! This crate provides classical Monte Carlo algorithms
//! for lattice simulations.

pub mod algorithms;
pub mod lattice;
pub mod models;

// Re-export Carlo.rs types for convenience
pub use carlo_rs::{
    CarloError, Context, Estimate, FromParams, MonteCarlo, Params, RayonBackend, Results,
    RunConfig, Scheduler,
};

// Re-export key types
pub use algorithms::{MetropolisCore, OPSSStrategy, ProposalStrategy, StandardStrategy, SWCore, WolffCore};
pub use lattice::{BondType, Lattice, LatticeMC, Neighbor};
pub use models::{HeisenbergModel, IsingModel, IsingModel2D, ModelMC, PottsModel, XYModel};
