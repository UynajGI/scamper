//! Algorithm engines for classical Monte Carlo.

mod metropolis;
mod opss_strategy;
mod proposal_strategy;
mod standard_strategy;
mod swendsen_wang;
mod wolff;

pub use metropolis::MetropolisCore;
pub use opss_strategy::OPSSStrategy;
pub use proposal_strategy::ProposalStrategy;
pub use standard_strategy::StandardStrategy;
pub use swendsen_wang::SWCore;
pub use wolff::WolffCore;
