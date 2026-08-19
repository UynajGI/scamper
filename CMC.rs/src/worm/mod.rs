//! Persistent classical worm algorithms in an extended configuration space.
//!
//! The generic driver samples explicit [`WormSector::Physical`] and
//! [`WormSector::Worm`] sectors with Metropolis-Hastings open, close and local
//! head moves. The first model is the ferromagnetic Ising high-temperature
//! graph representation on an arbitrary loop-free [`crate::CsrLattice`].
//!
//! # Multi-component lattices
//!
//! The kernel is a **two-defect (single-worm) kernel**: exactly one head and
//! one tail. Multi-defect/multi-leg worms are not implemented. Because the
//! Ising high-temperature graph ensemble factorizes over connected components,
//! multi-component lattices are still sampled exactly by
//! [`IsingGraphWormEnsemble`] (and the scheduler-ready
//! [`IsingGraphWormMC::from_lattice`]): one independent two-defect worm per
//! component, domain-separated RNG streams, observables combined additively.
//! A raw [`IsingGraphWormModel`] + [`WormKernel`] pair remains restricted to
//! connected lattices — its single defect pair would otherwise silently freeze
//! the other components — so [`IsingGraphWormModel::new`] rejects disconnected
//! input loudly for direct users.

mod ensemble;
mod error;
mod ising;
mod kernel;
mod mc;
mod model;
mod state;

pub use ensemble::{IsingComponentWorm, IsingGraphWormEnsemble};
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
