//! CMC.rs — reusable classical lattice Monte Carlo kernels on top of Carlo.rs.
//!
//! Carlo.rs owns execution concerns (RNG contexts, thermalization/measurement
//! scheduling, parallel backends, accumulation and parallel tempering). CMC.rs
//! owns graph-based physical models, state transitions and observables.
//!
//! # Quick start
//!
//! ```ignore
//! use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
//! use cmc_rs::{ClassicalMC, IsingModel, MetropolisCore};
//!
//! type Simulation = ClassicalMC<IsingModel, MetropolisCore>;
//! let mut params = Params::new();
//! params.set("Lx", 16);
//! params.set("Ly", 16);
//! params.set("beta", 0.44);
//! let results = Scheduler::new(RayonBackend::new(1), RunConfig::default())
//!     .run_one::<Simulation>(&params);
//! ```
//!
//! # Extensibility
//!
//! * [`CsrLattice`] is an arbitrary weighted undirected multigraph.
//! * [`Hamiltonian`] defines physical onsite/bond energy.
//! * Capability traits opt models into Metropolis, cluster, heat-bath or
//!   over-relaxation kernels without a monolithic model interface.
//! * [`ClassicalMC`] composes model + kernel + observable set into Carlo.rs's
//!   [`carlo_rs::MonteCarlo`] trait.

pub mod algorithm;
pub mod classical_mc;
pub mod ensemble;
pub mod hamiltonian;
pub mod lattice;
pub mod models;
pub mod moves;
pub mod multi_spin;
pub mod observables;
pub mod postprocess;
pub mod proposal;
pub mod system;
pub mod trial;
pub mod visit;

pub use algorithm::{
    Algorithm, ContinuousHeatBathCore, HeatBathCore, HybridCore, MetropolisCore,
    MicrocanonicalCore, SWCore, SimulationPhase, WolffCore,
};
pub use classical_mc::{build_lattice_from_params, ClassicalMC, FromHamiltonianParams};
pub use ensemble::{CanonicalEnsemble, Ensemble, ThermodynamicDelta};
pub use hamiltonian::{
    ClusterAuxiliary, ClusterModel, ContinuousHeatBathable, Hamiltonian, HeatBathable,
    Initializable, LocalFieldModel, Measurable, PairInteraction, Proposable,
};
pub use lattice::{
    build_chain, build_honeycomb, build_hypercubic, build_kagome, build_square, build_triangular,
    Bond, BondType, CsrLattice,
};
pub use models::{HeisenbergModel, IsingModel, ONModel, PottsModel, XYModel};
pub use moves::{
    BatchEnergyPatch, BatchEnergyWorkspace, BatchSpinMove, EnergyPatch, SiteSpinMove, Spin,
};
pub use multi_spin::{MultiSpinIsing, N_REPLICAS};
pub use observables::{
    DefaultObservableSet, EmptyObservableSet, EnergyPerSite, Magnetization, MomentSpec, Observable,
    ObservableSet, TotalEnergy,
};
pub use postprocess::{binder_cumulant, compute_correlation_1d, specific_heat, susceptibility};
pub use proposal::{OPSSStrategy, ProposalStrategy, ProposedSpin, StandardStrategy};
pub use system::{SiteChange, System};
pub use trial::{metropolis_hastings_step, ProposedMove, TrialEvaluator, TrialOutcome};
pub use visit::{SiteOrder, VisitSchedule};
