//! Persistent classical worm algorithms in an extended configuration space.
//!
//! The generic driver samples explicit [`WormSector::Physical`] and
//! [`WormSector::Worm`] sectors with Metropolis-Hastings open, close and local
//! head moves. The first model is the ferromagnetic Ising high-temperature
//! graph representation on an arbitrary loop-free [`crate::CsrLattice`].

mod error;
mod ising;
mod kernel;
mod mc;
mod model;
mod state;

pub use error::WormError;
pub use ising::{
    enumerate_ising_graph_expansion, ExactIsingGraphExpansion, IsingGraphConfiguration,
    IsingGraphPatch, IsingGraphWormModel, IsingWormStep,
};
pub use kernel::{
    EndpointPairHistogram, WormConfig, WormKernel, WormTransition, WormTransitionStatistics,
};
pub use mc::IsingGraphWormMC;
pub use model::{WormModel, WormStepDelta, WormStepProposal};
pub use state::{WormSector, WormState};
