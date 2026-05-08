//! QMC.rs - Quantum Monte Carlo algorithm toolbox
//!
//! This crate provides SSE (Stochastic Series Expansion) algorithms
//! for lattice quantum Monte Carlo simulations.

pub mod ed;
pub mod hilbert;
pub mod lattice;
pub mod models;
pub mod sse;

// Re-export Carlo.rs types for convenience
pub use carlo_rs::{
    CarloError, Context, Estimate, FromParams, MonteCarlo, Params, RayonBackend, Results,
    RunConfig, Scheduler,
};

// Re-export key types
pub use hilbert::{HilbertSpace, LocalState, OpType, SpinHalfHS};
pub use lattice::{BondType, Lattice, Neighbor};
pub use models::heisenberg::HeisenbergModel;
pub use models::xxz::XxzModel;
pub use sse::{LatticeQMC, OperatorSequence, SSECore, SSEEngine, SSEMonteCarlo, Vertex};