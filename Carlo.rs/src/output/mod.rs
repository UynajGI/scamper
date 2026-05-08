//! Output serialization for simulation results.
//!
//! Provides functions for saving and reading results:
//!
//! - [`save_hdf5()`]: Save to HDF5 format (requires `hdf5` feature)
//! - [`save_json()`]: Save to JSON format (always available)
//! - [`dataframe()`]: Load `*.results.json` for analysis
//! - [`measurement_from_obs()`]: Parse a single observable dict
//!
//! HDF5 is recommended for large simulations with many observables,
//! while JSON is suitable for smaller simulations or human-readable output.

pub mod hdf5;
pub mod json;
pub mod resulttools;

pub use hdf5::save_hdf5;
pub use json::save_json;
pub use resulttools::{
    dataframe, make_scalar, make_scalar_owned, measurement_from_obs, recursive_stack, ResultRow,
};
