//! HDF5 serialization for simulation results.
//!
//! Requires the `hdf5` feature. Without it, [`save_hdf5`] returns an error.

use std::path::Path;

use crate::{CarloError, Results};

/// Write `results` to an HDF5 file at `path`.
///
/// Requires the `hdf5` feature; returns [`CarloError::InvalidConfig`] otherwise.
#[cfg(feature = "hdf5")]
pub fn save_hdf5(results: &Results, path: &Path) -> Result<(), CarloError> {
    use hdf5::{File, H5};

    let file = File::create(path)?;

    // Write observables
    let obs_group = file.create_group("observables")?;
    for (name, est) in results.estimates() {
        let obs = obs_group.create_group(name)?;
        obs.write_scalar("mean", &est.mean)?;
        obs.write_scalar("stderr", &est.stderr)?;
        obs.write_scalar("autocorr_time", &est.autocorr_time)?;
        obs.write_scalar("n_bins", &(est.n_bins as i64))?;
    }

    // Write metadata
    let meta = results.metadata();
    let meta_group = file.create_group("metadata")?;
    meta_group.write_scalar("version", &meta.version)?;
    meta_group.write_scalar("timestamp", &meta.timestamp.to_rfc3339())?;
    meta_group.write_scalar("base_seed", &(meta.base_seed as i64))?;
    meta_group.write_scalar(
        "thermalization_sweeps",
        &(meta.thermalization_sweeps as i64),
    )?;
    meta_group.write_scalar("measurement_sweeps", &(meta.measurement_sweeps as i64))?;

    Ok(())
}

/// Stub returned when the `hdf5` feature is not enabled.
#[cfg(not(feature = "hdf5"))]
pub fn save_hdf5(_results: &Results, _path: &Path) -> Result<(), CarloError> {
    Err(CarloError::InvalidConfig {
        field: "hdf5".into(),
        reason: "hdf5 feature not enabled".into(),
    })
}
