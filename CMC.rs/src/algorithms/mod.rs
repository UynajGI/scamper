pub mod heat_bath;
pub mod hybrid;
pub mod metropolis;
pub mod microcanonical;
pub mod swendsen_wang;
pub mod wolff;

mod algorithm_tests;
mod common;

pub use common::{checked_probability, Algorithm, CanonicalLatticeKernel, SimulationPhase};
pub use heat_bath::{ContinuousHeatBathCore, HeatBathCore};
pub use hybrid::HybridCore;
pub use metropolis::MetropolisCore;
pub use microcanonical::MicrocanonicalCore;
pub use swendsen_wang::SWCore;
pub use wolff::WolffCore;

// These built-in kernels all target the ordinary canonical lattice weight.
impl CanonicalLatticeKernel for MetropolisCore {}
impl CanonicalLatticeKernel for HeatBathCore {}
impl CanonicalLatticeKernel for ContinuousHeatBathCore {}
impl CanonicalLatticeKernel for WolffCore {}
impl CanonicalLatticeKernel for SWCore {}
impl CanonicalLatticeKernel for MicrocanonicalCore {}
impl<A, B> CanonicalLatticeKernel for HybridCore<A, B>
where
    A: CanonicalLatticeKernel,
    B: CanonicalLatticeKernel,
{
}

impl CanonicalLatticeKernel for crate::dynamics::KawasakiCore {}
