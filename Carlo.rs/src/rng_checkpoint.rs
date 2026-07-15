//! RNG checkpoint serialization for Xoshiro family.
//!
//! Matches Carlo.jl's random_wrap.jl format for interoperability.

#[cfg(feature = "hdf5")]
use hdf5::Group;

/// Xoshiro RNG type identifier for HDF5.
pub const RNG_TYPE: &str = "xoroshiro256++";

/// RNG version: 1 = 4-field Xoshiro (standard), 2 = 5-field.
pub const RNG_VERSION: u64 = 1;

/// Trait for RNG checkpoint serialization via HDF5.
pub trait RngCheckpointHdf5: rand_core::Rng + rand_core::SeedableRng {
    /// Write RNG state to HDF5 group.
    #[cfg(feature = "hdf5")]
    fn write_checkpoint(&self, group: &mut Group) -> Result<(), crate::CarloError>;

    /// Read RNG from HDF5 group.
    #[cfg(feature = "hdf5")]
    fn read_checkpoint(group: &Group) -> Result<Self, crate::CarloError>;
}

#[cfg(feature = "hdf5")]
impl RngCheckpointHdf5 for rand_xoshiro::Xoshiro256PlusPlus {
    fn write_checkpoint(&self, group: &mut Group) -> Result<(), crate::CarloError> {
        // Write type identifier (stored as byte array)
        group
            .new_dataset_builder()
            .with_data(RNG_TYPE.as_bytes())
            .create("rng_type")
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot write rng_type: {}", e),
            })?;

        // Write version (stored as NE bytes)
        let ver_bytes = RNG_VERSION.to_ne_bytes();
        group
            .new_dataset_builder()
            .with_data(&ver_bytes)
            .create("rng_version")
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot write rng_version: {}", e),
            })?;

        // Serialize RNG state via serde
        let state_json =
            serde_json::to_string(self).map_err(|e| crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot serialize RNG state: {}", e),
            })?;

        group
            .new_dataset_builder()
            .with_data(state_json.as_bytes())
            .create("rng_state_json")
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot write rng_state_json: {}", e),
            })?;

        Ok(())
    }

    fn read_checkpoint(group: &Group) -> Result<Self, crate::CarloError> {
        // Verify type
        let type_bytes: Vec<u8> = group
            .dataset("rng_type")
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot read rng_type: {}", e),
            })?
            .read_1d::<u8>()
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot parse rng_type: {}", e),
            })?
            .to_vec();

        let rng_type = String::from_utf8_lossy(&type_bytes);
        if rng_type != RNG_TYPE {
            return Err(crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("RNG type mismatch: expected {}, got {}", RNG_TYPE, rng_type),
            });
        }

        // Check version
        let ver_bytes: Vec<u8> = group
            .dataset("rng_version")
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot read rng_version: {}", e),
            })?
            .read_1d::<u8>()
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot parse rng_version: {}", e),
            })?
            .to_vec();
        let mut arr = [0u8; 8];
        let len = ver_bytes.len().min(8);
        arr[..len].copy_from_slice(&ver_bytes[..len]);
        let rng_version = u64::from_ne_bytes(arr);

        if rng_version != RNG_VERSION {
            return Err(crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!(
                    "RNG version mismatch: checkpoint was done with version {}",
                    rng_version
                ),
            });
        }

        // Read and deserialize state
        let state_bytes: Vec<u8> = group
            .dataset("rng_state_json")
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot read rng_state_json: {}", e),
            })?
            .read_1d::<u8>()
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot parse rng_state_json: {}", e),
            })?
            .to_vec();

        let state_json = String::from_utf8_lossy(&state_bytes);
        serde_json::from_str(&state_json).map_err(|e| crate::CarloError::InvalidConfig {
            field: "checkpoint".into(),
            reason: format!("Cannot deserialize RNG state: {}", e),
        })
    }
}
