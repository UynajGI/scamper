//! Scamper — a Monte Carlo simulation workspace.
//!
//! Umbrella crate for the scamper workspace: it re-exports every
//! registry-published member so `cargo add scamper` pulls in the whole
//! stack. Members stay independently usable as individual crates.
//!
//! - [`carlo_rs`] — core framework: run scheduling, deterministic RNG
//!   setup, thermalization/measurement phases, accumulation and error
//!   analysis, checkpointing, rayon and MPI parallel backends.
//! - [`cmc_rs`] — classical kernels: lattice (Metropolis/Wolff/SW,
//!   multi-spin coding), particle (NVT/NPT/μVT), generalized ensembles
//!   (Wang-Landau, multicanonical, umbrella), worm and dynamics.
//! - [`qmc_rs`] — quantum kernels: continuous-time interaction-expansion
//!   lattice QMC, wormhole impurity QMC, variational Monte Carlo.
//!
//! The statistical MCMC member is not re-exported: its crate name is
//! unavailable on crates.io, so it stays workspace-only under a path
//! dependency.
//!
//! The `hdf5` and `mpi` features forward to every member crate.

pub use carlo_rs;
pub use cmc_rs;
pub use qmc_rs;
