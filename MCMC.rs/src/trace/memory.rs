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
        file.create_dataset_simple("dimension", &[1], &(self.dimension as u64))?;
        file.create_dataset_simple("thinning", &[1], &(self.thinning as u64))?;
        file.create_dataset_simple("seen_iterations", &[1], &(self.seen_iterations as u64))?;
        file.create_dataset_simple(
            "chain_id",
            &[1],
            &(self.chain_id.map_or(u64::MAX, |value| value as u64)),
        )?;
        file.create_dataset_simple("positions", &[self.positions.len() as u64], &self.positions)?;
        file.create_dataset_simple(
            "log_density",
            &[self.log_density.len() as u64],
            &self.log_density,
        )?;
        file.create_dataset_simple("accepted", &[self.accepted.len() as u64], &self.accepted)?;
        let acceptance_rate = self
            .acceptance_rate
            .iter()
            .map(|value| value.unwrap_or(f64::NAN))
            .collect::<Vec<_>>();
        file.create_dataset_simple(
            "acceptance_rate",
            &[acceptance_rate.len() as u64],
            &acceptance_rate,
        )?;
        file.create_dataset_simple("divergent", &[self.divergent.len() as u64], &self.divergent)?;
        let energy_error = self
            .energy_error
            .iter()
            .map(|value| value.unwrap_or(f64::NAN))
            .collect::<Vec<_>>();
        file.create_dataset_simple("energy_error", &[energy_error.len() as u64], &energy_error)?;
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
        {
            return Err(McmcError::InvalidConfig(
                "trace contains invalid floating-point values".to_string(),
            ));
        }
        Ok(())
    }
}

impl TraceStore for MemoryTrace {
    fn record(
        &mut self,
        chain_id: usize,
        state: &EuclideanState,
        report: &TransitionReport,
    ) -> Result<bool, McmcError> {
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
        self.positions.extend_from_slice(state.position());
        self.log_density.push(state.log_density());
        self.accepted.push(match report.accepted {
            Some(true) => 1,
            Some(false) => 0,
            None => -1,
        });
        self.acceptance_rate.push(report.acceptance_rate());
        self.divergent.push(u8::from(report.divergent));
        self.energy_error.push(report.energy_error);
        Ok(true)
    }

    fn clear(&mut self) {
        self.positions.clear();
        self.log_density.clear();
        self.accepted.clear();
        self.acceptance_rate.clear();
        self.divergent.clear();
        self.energy_error.clear();
        self.seen_iterations = 0;
        self.chain_id = None;
    }

    fn len(&self) -> usize {
        self.log_density.len()
    }
}
