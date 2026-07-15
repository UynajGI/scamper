//! Explicit finite-occupation spin-boson solver for one or a few cavity modes.
//!
//! This backend keeps boson occupations instead of integrating out the bath.
//! It uses a closed-worldline bridge heat-bath sampler built from the exact finite-basis
//! propagator, so there is no Trotter error. The only controlled physical
//! approximation is the per-mode occupation cutoff.

pub mod basis;
pub mod mc;
pub mod model;
pub mod transfer;
pub mod worldline;

pub use basis::{OccupationBasis, SpinState};
pub use mc::OccupationWorldlineQmc;
pub use model::{CavityMode, OccupationModelKind, OccupationSpinBosonModel};
pub use worldline::{OccupationObservables, OccupationWorldlineSampler};
