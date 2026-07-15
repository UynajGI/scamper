//! Version information for checkpoint files.

use serde::{Deserialize, Serialize};

/// Version information for HDF5 checkpoints.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Version {
    /// Carlo.rs version.
    pub carlo_version: String,

    /// Monte Carlo implementation version (optional).
    pub mc_version: Option<String>,

    /// RNG version (for Xoshiro compatibility).
    pub rng_version: u64,
}

impl Version {
    /// Create new version info.
    pub fn new(mc_version: Option<&str>) -> Self {
        Self {
            carlo_version: env!("CARGO_PKG_VERSION").to_string(),
            mc_version: mc_version.map(|s| s.to_string()),
            rng_version: crate::RNG_VERSION,
        }
    }

    /// Create version with default values.
    pub fn current() -> Self {
        Self::new(None)
    }
}

#[cfg(feature = "hdf5")]
impl Version {
    /// Write version to HDF5 group.
    pub fn write_hdf5(&self, group: &mut hdf5::Group) -> Result<(), crate::CarloError> {
        group
            .new_dataset_builder()
            .with_data(self.carlo_version.as_bytes())
            .create("carlo_version")
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "version".into(),
                reason: format!("Cannot write carlo_version: {}", e),
            })?;

        if let Some(ref mc_ver) = self.mc_version {
            group
                .new_dataset_builder()
                .with_data(mc_ver.as_bytes())
                .create("mc_version")
                .map_err(|e| crate::CarloError::InvalidConfig {
                    field: "version".into(),
                    reason: format!("Cannot write mc_version: {}", e),
                })?;
        }

        let rng_bytes = self.rng_version.to_ne_bytes();
        group
            .new_dataset_builder()
            .with_data(&rng_bytes)
            .create("rng_version")
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "version".into(),
                reason: format!("Cannot write rng_version: {}", e),
            })?;

        Ok(())
    }

    /// Read version from HDF5 group.
    pub fn read_hdf5(group: &hdf5::Group) -> Result<Self, crate::CarloError> {
        let carlo_bytes: Vec<u8> = group
            .dataset("carlo_version")
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "version".into(),
                reason: format!("Cannot read carlo_version: {}", e),
            })?
            .read_1d::<u8>()
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "version".into(),
                reason: format!("Cannot parse carlo_version: {}", e),
            })?
            .to_vec();
        let carlo_version = String::from_utf8_lossy(&carlo_bytes).to_string();

        let mc_version = group
            .dataset("mc_version")
            .ok()
            .and_then(|ds| ds.read_1d::<u8>().ok())
            .map(|bytes| String::from_utf8_lossy(&bytes.to_vec()).to_string());

        let rng_bytes: Vec<u8> = group
            .dataset("rng_version")
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "version".into(),
                reason: format!("Cannot read rng_version: {}", e),
            })?
            .read_1d::<u8>()
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "version".into(),
                reason: format!("Cannot parse rng_version: {}", e),
            })?
            .to_vec();

        // Reconstruct u64 from NE bytes
        let mut arr = [0u8; 8];
        let len = rng_bytes.len().min(8);
        arr[..len].copy_from_slice(&rng_bytes[..len]);
        let rng_version = u64::from_ne_bytes(arr);

        Ok(Self {
            carlo_version,
            mc_version,
            rng_version,
        })
    }
}
