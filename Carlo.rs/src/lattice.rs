//! Basic 2-D lattice parameter parsing.

use crate::{CarloError, Params};

/// Simple 2D lattice parameters.
///
/// `ly` defaults to `lx` when only `lx` is provided.
#[derive(Debug, Clone)]
pub struct LatticeParams {
    /// Number of sites along x.
    pub lx: usize,
    /// Number of sites along y.
    pub ly: usize,
}

impl LatticeParams {
    /// Parse `lx` and `ly` from `params`.
    ///
    /// `lx` is required; `ly` defaults to `lx`.  Returns
    /// [`CarloError::InvalidConfig`] if `lx` is missing or either dimension is zero.
    pub fn from_params(params: &Params) -> Result<Self, CarloError> {
        let lx = params.get("lx").ok_or_else(|| CarloError::InvalidConfig {
            field: "lx".into(),
            reason: "required lattice dimension".into(),
        })?;
        let ly = params.get("ly").unwrap_or(lx);

        if lx == 0 || ly == 0 {
            return Err(CarloError::InvalidConfig {
                field: "lattice".into(),
                reason: "dimensions must be positive".into(),
            });
        }

        Ok(Self { lx, ly })
    }

    /// Total number of lattice sites (`lx * ly`).
    pub fn n_sites(&self) -> usize {
        self.lx * self.ly
    }
}
