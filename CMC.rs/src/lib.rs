//! CMC.rs — Classical Monte Carlo algorithm toolbox.
//!
//! Built on [Carlo.rs] for scheduling, measurement, and result analysis.
//!
//! # Architecture
//!
//! ```text
//! ClassicalMC<M, A>  ← impl MonteCarlo + FromParams (pre-built)
//!   ├── System       ← mutable state: spins, energy
//!   ├── Model (M)    ← stateless physics (Ising, Potts, XY, Heisenberg)
//!   └── Algorithm (A)← update strategy (Metropolis, Wolff, Swendsen-Wang)
//! ```
//!
//! # Usage
//!
//! ```ignore
//! use cmc_rs::{ClassicalMC, IsingModel, MetropolisCore};
//! use carlo_rs::{Scheduler, RayonBackend, RunConfig, Params};
//!
//! type IsingMetro = ClassicalMC<IsingModel, MetropolisCore>;
//!
//! let mut params = Params::new();
//! params.set("L", 16);
//! params.set("beta", 0.5);
//!
//! let config = RunConfig {
//!     thermalization_sweeps: 1000,
//!     measurement_sweeps: 10_000,
//!     binsize: 100,
//!     base_seed: 42,
//!     ..Default::default()
//! };
//!
//! let scheduler = Scheduler::new(RayonBackend::new(1), config);
//! let results = scheduler.run_one::<IsingMetro>(&params);
//! ```

pub mod algorithm;
pub mod classical_mc;
pub mod lattice;
pub mod model;
pub mod proposal;
pub mod system;

// Re-export key types
pub use algorithm::{Algorithm, MetropolisCore, SWCore, WolffCore};
pub use classical_mc::{ClassicalMC, FromModelParams};
pub use lattice::{build_chain, build_hypercubic, build_square, BondType, Lattice, Neighbor};
pub use model::{HeisenbergModel, IsingModel, Model, PottsModel, XYModel};
pub use proposal::{OPSSStrategy, ProposalStrategy, StandardStrategy};
pub use system::System;
