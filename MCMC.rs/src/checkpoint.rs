use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;

use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};

use crate::{EuclideanState, McmcError, SamplingPhase, TransitionReport};

pub const CHECKPOINT_FORMAT: &str = "mcmc-rs-chain-v1";

/// Stable metadata used to reject checkpoints restored against another model.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct TargetFingerprint {
    pub name: String,
    pub version: String,
    pub dimension: usize,
    pub parameter_names: Vec<String>,
}

impl TargetFingerprint {
    pub fn validate(&self) -> Result<(), McmcError> {
        if self.name.is_empty() || self.version.is_empty() || self.dimension == 0 {
            return Err(McmcError::InvalidConfig(
                "target fingerprint contains empty required fields".to_string(),
            ));
        }
        if !self.parameter_names.is_empty() && self.parameter_names.len() != self.dimension {
            return Err(McmcError::DimensionMismatch {
                expected: self.dimension,
                actual: self.parameter_names.len(),
            });
        }
        Ok(())
    }
}

/// Serializable chain state including RNG, kernel adaptation and retained trace.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChainCheckpoint<K, R, Tr> {
    pub format: String,
    pub chain_id: usize,
    pub phase: SamplingPhase,
    pub target: TargetFingerprint,
    pub state: EuclideanState,
    pub kernel: K,
    pub rng: R,
    pub trace: Tr,
    pub last_report: TransitionReport,
}

impl<K, R, Tr> ChainCheckpoint<K, R, Tr>
where
    K: Serialize + DeserializeOwned,
    R: Serialize + DeserializeOwned,
    Tr: Serialize + DeserializeOwned,
{
    pub fn save_json(&self, path: impl AsRef<Path>) -> Result<(), McmcError> {
        self.validate_format()?;
        let writer = BufWriter::new(File::create(path)?);
        serde_json::to_writer(writer, self)?;
        Ok(())
    }

    pub fn load_json(path: impl AsRef<Path>) -> Result<Self, McmcError> {
        let reader = BufReader::new(File::open(path)?);
        let checkpoint: Self = serde_json::from_reader(reader)?;
        checkpoint.validate_format()?;
        Ok(checkpoint)
    }

    pub fn validate_target(&self, expected: &TargetFingerprint) -> Result<(), McmcError> {
        if &self.target == expected {
            Ok(())
        } else {
            Err(McmcError::TargetMismatch)
        }
    }

    fn validate_format(&self) -> Result<(), McmcError> {
        if self.format != CHECKPOINT_FORMAT {
            return Err(McmcError::CheckpointFormat {
                expected: CHECKPOINT_FORMAT.to_string(),
                found: self.format.clone(),
            });
        }
        self.target.validate()?;
        self.last_report.validate()?;
        self.state.validate()?;
        if self.state.dimension() != self.target.dimension {
            return Err(McmcError::DimensionMismatch {
                expected: self.target.dimension,
                actual: self.state.dimension(),
            });
        }
        Ok(())
    }
}
