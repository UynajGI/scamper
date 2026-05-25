//! CMC.rs — Classical Monte Carlo algorithm toolbox.
//!
//! Built on [Carlo.rs] for scheduling, measurement, and result analysis.
//!
//! # Architecture
//!
//! ```text
//! ClassicalMC<H, A>  ← impl MonteCarlo + FromParams (pre-built)
//!   ├── System       ← mutable state: spins, energy, beta
//!   ├── H: Hamiltonian ← stateless physics (Ising, Potts, XY, Heisenberg)
//!   │   + ClusterModel  ← cluster algorithm support
//!   │   + Proposable    ← spin proposal
//!   │   + Measurable    ← magnetization
//!   ├── A: Algorithm ← update strategy (Metropolis, Wolff, Swendsen-Wang)
//!   └── observables  ← pluggable measurement system
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
pub mod hamiltonian;
pub mod lattice;
pub mod models;
pub mod multi_spin;
pub mod observables;
pub mod postprocess;
pub mod proposal;
pub mod system;

// Re-export key types from hamiltonian (traits)
pub use hamiltonian::{ClusterModel, Hamiltonian, HeatBathable, Measurable, Proposable};

// Re-export models
pub use models::{HeisenbergModel, IsingModel, PottsModel, XYModel};

// Re-export algorithms
pub use algorithm::{
    Algorithm, HeatBathCore, MetropolisCore, MicrocanonicalCore, SWCore, WolffCore,
};

// Re-export observables
pub use observables::{
    DefaultObservableSet, EnergyPerSite, Magnetization, Observable, TotalEnergy,
};

// Re-export classical_mc
pub use classical_mc::{ClassicalMC, FromHamiltonianParams};

// Re-export lattice
pub use lattice::{
    build_chain, build_honeycomb, build_hypercubic, build_kagome, build_square, build_triangular,
    BondType, CsrLattice,
};

// Re-export multi-spin
pub use multi_spin::{MultiSpinIsing, N_REPLICAS};

// Re-export proposal
pub use proposal::{OPSSStrategy, ProposalStrategy, StandardStrategy};

// Re-export system
pub use system::System;

// Re-export postprocess
pub use postprocess::{binder_cumulant, specific_heat, susceptibility};
