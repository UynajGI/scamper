use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;

use serde::{Deserialize, Serialize};

use super::{TraceStore, TraceView};
use crate::{EuclideanState, McmcError, TransitionReport};

/// Contiguous row-major in-memory posterior trace.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct MemoryTrace {
    dimension: usize,
    thinning: usize,
    positions: Vec<f64>,
    log_density: Vec<f64>,
    accepted: Vec<i8>,
    acceptance_rate: Vec<Option<f64>>,
    divergent: Vec<u8>,
    energy_error: Vec<Option<f64>>,
    #[serde(default)]
    energy: Vec<Option<f64>>,
    #[serde(default)]
    tree_depth: Vec<Option<u16>>,
    #[serde(default)]
    max_tree_depth_reached: Vec<u8>,
    seen_iterations: usize,
    chain_id: Option<usize>,
}

impl MemoryTrace {
    pub fn new(dimension: usize, thinning: usize) -> Result<Self, McmcError> {
        if dimension == 0 {
            return Err(McmcError::InvalidConfig(
                "trace dimension must be positive".to_string(),
            ));
        }
        if thinning == 0 {
            return Err(McmcError::InvalidConfig(
                "trace thinning must be at least one".to_string(),
            ));
        }
        Ok(Self {
            dimension,
            thinning,
            positions: Vec::new(),
            log_density: Vec::new(),
            accepted: Vec::new(),
            acceptance_rate: Vec::new(),
            divergent: Vec::new(),
            energy_error: Vec::new(),
            energy: Vec::new(),
            tree_depth: Vec::new(),
            max_tree_depth_reached: Vec::new(),
            seen_iterations: 0,
            chain_id: None,
        })
    }

    pub fn reserve_draws(&mut self, draws: usize) {
        self.positions.reserve(draws.saturating_mul(self.dimension));
        self.log_density.reserve(draws);
        self.accepted.reserve(draws);
        self.acceptance_rate.reserve(draws);
        self.divergent.reserve(draws);
        self.energy_error.reserve(draws);
        self.energy.reserve(draws);
        self.tree_depth.reserve(draws);
        self.max_tree_depth_reached.reserve(draws);
    }

    pub const fn dimension(&self) -> usize {
        self.dimension
    }

    pub const fn thinning(&self) -> usize {
        self.thinning
    }

    pub const fn seen_iterations(&self) -> usize {
        self.seen_iterations
    }

    pub const fn chain_id(&self) -> Option<usize> {
        self.chain_id
    }

    pub fn positions(&self) -> &[f64] {
        &self.positions
    }

    pub fn log_densities(&self) -> &[f64] {
        &self.log_density
    }

    pub fn accepted(&self) -> &[i8] {
        &self.accepted
    }

    pub fn acceptance_rates(&self) -> &[Option<f64>] {
        &self.acceptance_rate
    }

    pub fn divergences(&self) -> &[u8] {
        &self.divergent
    }

    pub fn energy_errors(&self) -> &[Option<f64>] {
        &self.energy_error
    }

    pub fn energies(&self) -> &[Option<f64>] {
        &self.energy
    }

    pub fn tree_depths(&self) -> &[Option<u16>] {
        &self.tree_depth
    }

    pub fn max_tree_depth_reached(&self) -> &[u8] {
        &self.max_tree_depth_reached
    }

    pub fn view(&self) -> TraceView<'_> {
        TraceView::new(&self.positions, self.len(), self.dimension)
    }

    pub fn draw(&self, index: usize) -> Option<&[f64]> {
        self.view().draw(index)
    }

    pub fn parameter(&self, index: usize) -> Option<super::ParameterIter<'_>> {
        self.view().parameter(index)
    }

    pub fn save_json(&self, path: impl AsRef<Path>) -> Result<(), McmcError> {
        self.validate()?;
        let writer = BufWriter::new(File::create(path)?);
        serde_json::to_writer(writer, self)?;
        Ok(())
    }

    pub fn load_json(path: impl AsRef<Path>) -> Result<Self, McmcError> {
        let reader = BufReader::new(File::open(path)?);
        let trace: Self = serde_json::from_reader(reader)?;
        trace.validate()?;
        Ok(trace)
    }

    /// Export the complete retained trace to a flat HDF5 schema.
    #[cfg(feature = "hdf5")]
    pub fn save_hdf5(&self, path: impl AsRef<Path>) -> Result<(), McmcError> {
        self.validate()?;
        let file = hdf5::File::create(path)?;
        file.new_dataset::<u64>()
            .shape([1])
            .create("dimension")?
            .write(&[self.dimension as u64])?;
        file.new_dataset::<u64>()
            .shape([1])
            .create("thinning")?
            .write(&[self.thinning as u64])?;
        file.new_dataset::<u64>()
            .shape([1])
            .create("seen_iterations")?
            .write(&[self.seen_iterations as u64])?;
        file.new_dataset::<u64>()
            .shape([1])
            .create("chain_id")?
            .write(&[self.chain_id.map_or(u64::MAX, |value| value as u64)])?;
        file.new_dataset::<f64>()
            .shape([self.positions.len()])
            .create("positions")?
            .write(self.positions.as_slice())?;
        file.new_dataset::<f64>()
            .shape([self.log_density.len()])
            .create("log_density")?
            .write(self.log_density.as_slice())?;
        file.new_dataset::<i8>()
            .shape([self.accepted.len()])
            .create("accepted")?
            .write(self.accepted.as_slice())?;
        let acceptance_rate = self
            .acceptance_rate
            .iter()
            .map(|value| value.unwrap_or(f64::NAN))
            .collect::<Vec<_>>();
        file.new_dataset::<f64>()
            .shape([acceptance_rate.len()])
            .create("acceptance_rate")?
            .write(acceptance_rate.as_slice())?;
        file.new_dataset::<u8>()
            .shape([self.divergent.len()])
            .create("divergent")?
            .write(self.divergent.as_slice())?;
        let energy_error = self
            .energy_error
            .iter()
            .map(|value| value.unwrap_or(f64::NAN))
            .collect::<Vec<_>>();
        file.new_dataset::<f64>()
            .shape([energy_error.len()])
            .create("energy_error")?
            .write(energy_error.as_slice())?;
        let energy = normalized_optional_f64(&self.energy, self.len());
        file.new_dataset::<f64>()
            .shape([energy.len()])
            .create("energy")?
            .write(energy.as_slice())?;
        let tree_depth = normalized_optional_u16(&self.tree_depth, self.len());
        file.new_dataset::<u16>()
            .shape([tree_depth.len()])
            .create("tree_depth")?
            .write(tree_depth.as_slice())?;
        let max_tree_depth_reached = normalized_u8(&self.max_tree_depth_reached, self.len());
        file.new_dataset::<u8>()
            .shape([max_tree_depth_reached.len()])
            .create("max_tree_depth_reached")?
            .write(max_tree_depth_reached.as_slice())?;
        Ok(())
    }

    /// Load a trace previously written by [`Self::save_hdf5`].
    #[cfg(feature = "hdf5")]
    pub fn load_hdf5(path: impl AsRef<Path>) -> Result<Self, McmcError> {
        let file = hdf5::File::open(path)?;
        let dimension = file.dataset("dimension")?.read_1d::<u64>()?[0] as usize;
        let thinning = file.dataset("thinning")?.read_1d::<u64>()?[0] as usize;
        let seen_iterations = file.dataset("seen_iterations")?.read_1d::<u64>()?[0] as usize;
        let chain_id_value = file.dataset("chain_id")?.read_1d::<u64>()?[0];
        let chain_id = (chain_id_value != u64::MAX).then_some(chain_id_value as usize);
        let trace = Self {
            dimension,
            thinning,
            positions: file.dataset("positions")?.read_1d::<f64>()?.to_vec(),
            log_density: file.dataset("log_density")?.read_1d::<f64>()?.to_vec(),
            accepted: file.dataset("accepted")?.read_1d::<i8>()?.to_vec(),
            acceptance_rate: file
                .dataset("acceptance_rate")?
                .read_1d::<f64>()?
                .iter()
                .map(|value| (!value.is_nan()).then_some(*value))
                .collect(),
            divergent: file.dataset("divergent")?.read_1d::<u8>()?.to_vec(),
            energy_error: file
                .dataset("energy_error")?
                .read_1d::<f64>()?
                .iter()
                .map(|value| (!value.is_nan()).then_some(*value))
                .collect(),
            energy: file
                .dataset("energy")
                .ok()
                .map_or_else(Vec::new, |dataset| {
                    dataset
                        .read_1d::<f64>()
                        .map(|values| {
                            values
                                .iter()
                                .map(|value| (!value.is_nan()).then_some(*value))
                                .collect()
                        })
                        .unwrap_or_default()
                }),
            tree_depth: file
                .dataset("tree_depth")
                .ok()
                .map_or_else(Vec::new, |dataset| {
                    dataset
                        .read_1d::<u16>()
                        .map(|values| {
                            values
                                .iter()
                                .map(|value| (*value != u16::MAX).then_some(*value))
                                .collect()
                        })
                        .unwrap_or_default()
                }),
            max_tree_depth_reached: file
                .dataset("max_tree_depth_reached")
                .ok()
                .and_then(|dataset| dataset.read_1d::<u8>().ok())
                .map_or_else(Vec::new, |values| values.to_vec()),
            seen_iterations,
            chain_id,
        };
        trace.validate()?;
        Ok(trace)
    }

    pub fn validate(&self) -> Result<(), McmcError> {
        let draws = self.log_density.len();
        if self.positions.len() != draws.saturating_mul(self.dimension)
            || self.accepted.len() != draws
            || self.acceptance_rate.len() != draws
            || self.divergent.len() != draws
            || self.energy_error.len() != draws
            || !optional_column_length_is_valid(self.energy.len(), draws)
            || !optional_column_length_is_valid(self.tree_depth.len(), draws)
            || !optional_column_length_is_valid(self.max_tree_depth_reached.len(), draws)
        {
            return Err(McmcError::InvalidConfig(
                "trace columns have inconsistent lengths".to_string(),
            ));
        }
        if self.dimension == 0 || self.thinning == 0 {
            return Err(McmcError::InvalidConfig(
                "trace metadata is invalid".to_string(),
            ));
        }
        let expected_draws = self.seen_iterations.div_ceil(self.thinning);
        if draws != expected_draws
            || self.chain_id.is_some() != (self.seen_iterations > 0)
            || self.accepted.iter().any(|value| !(-1..=1).contains(value))
            || self.divergent.iter().any(|value| *value > 1)
            || self.max_tree_depth_reached.iter().any(|value| *value > 1)
        {
            return Err(McmcError::InvalidConfig(
                "trace metadata or discrete columns are inconsistent".to_string(),
            ));
        }
        if self.positions.iter().any(|value| !value.is_finite())
            || self.log_density.iter().any(|value| !value.is_finite())
            || self
                .acceptance_rate
                .iter()
                .flatten()
                .any(|value| !value.is_finite() || !(0.0..=1.0).contains(value))
            || self
                .energy_error
                .iter()
                .flatten()
                .any(|value| !value.is_finite())
            || self.energy.iter().flatten().any(|value| !value.is_finite())
        {
            return Err(McmcError::InvalidConfig(
                "trace contains invalid floating-point values".to_string(),
            ));
        }
        Ok(())
    }
}

fn optional_column_length_is_valid(length: usize, draws: usize) -> bool {
    length == 0 || length == draws
}

fn normalize_optional_column<T: Clone>(column: &mut Vec<T>, draws: usize, missing: T) {
    if column.is_empty() && draws > 0 {
        column.resize(draws, missing);
    }
}

#[cfg(feature = "hdf5")]
fn normalized_optional_f64(column: &[Option<f64>], draws: usize) -> Vec<f64> {
    if column.is_empty() {
        vec![f64::NAN; draws]
    } else {
        column
            .iter()
            .map(|value| value.unwrap_or(f64::NAN))
            .collect()
    }
}

#[cfg(feature = "hdf5")]
fn normalized_optional_u16(column: &[Option<u16>], draws: usize) -> Vec<u16> {
    if column.is_empty() {
        vec![u16::MAX; draws]
    } else {
        column
            .iter()
            .map(|value| value.unwrap_or(u16::MAX))
            .collect()
    }
}

#[cfg(feature = "hdf5")]
fn normalized_u8(column: &[u8], draws: usize) -> Vec<u8> {
    if column.is_empty() {
        vec![0; draws]
    } else {
        column.to_vec()
    }
}

impl TraceStore for MemoryTrace {
    fn record(
        &mut self,
        chain_id: usize,
        state: &EuclideanState,
        report: &TransitionReport,
    ) -> Result<bool, McmcError> {
        state.validate()?;
        report.validate()?;
        if state.dimension() != self.dimension {
            return Err(McmcError::DimensionMismatch {
                expected: self.dimension,
                actual: state.dimension(),
            });
        }
        if self.chain_id.is_some_and(|stored| stored != chain_id) {
            return Err(McmcError::InvalidConfig(
                "one MemoryTrace cannot mix multiple chain IDs".to_string(),
            ));
        }
        self.chain_id = Some(chain_id);
        self.seen_iterations = self.seen_iterations.saturating_add(1);
        if !(self.seen_iterations - 1).is_multiple_of(self.thinning) {
            return Ok(false);
        }
        normalize_optional_column(&mut self.energy, self.log_density.len(), None);
        normalize_optional_column(&mut self.tree_depth, self.log_density.len(), None);
        normalize_optional_column(&mut self.max_tree_depth_reached, self.log_density.len(), 0);
        self.positions.extend_from_slice(state.position());
        self.log_density.push(state.log_density());
        self.accepted.push(match report.accepted {
            Some(true) => 1,
            Some(false) => 0,
            None => -1,
        });
        self.acceptance_rate.push(
            report
                .acceptance_statistic
                .or_else(|| report.acceptance_rate()),
        );
        self.divergent.push(u8::from(report.divergent));
        self.energy_error.push(report.energy_error);
        self.energy.push(report.energy);
        self.tree_depth.push(report.tree_depth);
        self.max_tree_depth_reached
            .push(u8::from(report.max_tree_depth_reached));
        Ok(true)
    }

    fn clear(&mut self) {
        self.positions.clear();
        self.log_density.clear();
        self.accepted.clear();
        self.acceptance_rate.clear();
        self.divergent.clear();
        self.energy_error.clear();
        self.energy.clear();
        self.tree_depth.clear();
        self.max_tree_depth_reached.clear();
        self.seen_iterations = 0;
        self.chain_id = None;
    }

    fn len(&self) -> usize {
        self.log_density.len()
    }
}
