//! Discrete-time Suzuki-Trotter worm / path-integral QMC.
//!
//! Layered:
//! - [`config`] — `SpaceTimeConfig` (N sites × M slices) + observables
//! - [`worm`]   — update kernels: [`worm_sweep`](worm::worm_sweep),
//!   [`local_metropolis_sweep`](worm::local_metropolis_sweep)
//!
//! The future continuous-time worm will share [`hamiltonian`](crate::hamiltonian)
//! and the observable layer, differing only in the time representation and
//! its move set.

pub mod config;
pub mod worm;

pub use config::{SpaceTimeConfig, Spin};
pub use worm::{local_metropolis_sweep, worm_sweep};
