//! Generalized ensembles, density-of-states estimation and reweighting.

mod axis;
mod bias;
mod error;
mod exact;
mod histogram;
mod macrostate;
mod multicanonical;
mod reweight;
mod wang_landau;

pub use axis::{BinnedAxis, DiscreteAxis, MacrostateAxis};
pub use bias::{FixedBias, HarmonicUmbrellaBias, LogBias, MulticanonicalBias};
pub use error::GeneralizedError;
pub use exact::{enumerate_ising_density_of_states, ExactIsingDensityOfStates};
pub use histogram::{Histogram, LogDensityOfStates};
pub use macrostate::{
    EnergyMacrostate, Macrostate, MagnetizationMacrostate, ParticleNumberMacrostate,
};
pub use multicanonical::EnergyBiasCore;
pub use reweight::{canonical_reweight, CanonicalReweighting};
pub use wang_landau::{
    IsingWangLandau, WangLandauConfig, WangLandauCore, WangLandauPhase, WangLandauRefinement,
    WangLandauRunControl, WangLandauState, WangLandauTermination,
};
