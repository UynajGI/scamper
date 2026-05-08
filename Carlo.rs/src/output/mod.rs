//! Output serialization for simulation results.
//!
//! Provides functions for saving results to different formats:
//!
//! - [`save_hdf5()`]: Save to HDF5 format (requires `hdf5` feature)
//! - [`save_json()`]: Save to JSON format (always available)
//!
//! HDF5 is recommended for large simulations with many observables,
//! while JSON is suitable for smaller simulations or human-readable output.

pub mod hdf5;
pub mod json;

pub use hdf5::save_hdf5;
pub use json::save_json;
