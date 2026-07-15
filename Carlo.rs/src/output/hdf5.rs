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
    use hdf5::File;

    let file = File::create(path)?;

    // Write observables
    let obs_group = file.create_group("observables")?;
    for (name, est) in results.estimates() {
        let obs = obs_group.create_group(name)?;

        let ds = obs.new_dataset::<f64>().create("mean")?;
        ds.write_scalar(&est.mean)?;
        let ds = obs.new_dataset::<f64>().create("stderr")?;
        ds.write_scalar(&est.stderr)?;
        let ds = obs.new_dataset::<f64>().create("autocorr_time")?;
        ds.write_scalar(&est.autocorr_time)?;
        let ds = obs.new_dataset::<i64>().create("n_bins")?;
        ds.write_scalar(&(est.n_bins as i64))?;
    }

    // Write metadata
    let meta = results.metadata();
    let meta_group = file.create_group("metadata")?;

    // Store strings as byte arrays
    let version_bytes = meta.version.as_bytes();
    meta_group
        .new_dataset_builder()
        .with_data(version_bytes)
        .create("version")?;

    let timestamp_bytes = meta.timestamp.to_rfc3339();
    let timestamp_bytes = timestamp_bytes.as_bytes();
    meta_group
        .new_dataset_builder()
        .with_data(timestamp_bytes)
        .create("timestamp")?;

    let ds = meta_group.new_dataset::<i64>().create("base_seed")?;
    ds.write_scalar(&(meta.base_seed as i64))?;
    let ds = meta_group
        .new_dataset::<i64>()
        .create("thermalization_sweeps")?;
    ds.write_scalar(&(meta.thermalization_sweeps as i64))?;
    let ds = meta_group
        .new_dataset::<i64>()
        .create("measurement_sweeps")?;
    ds.write_scalar(&(meta.measurement_sweeps as i64))?;

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
