use crate::{CarloError, Params};

/// Simple 2D lattice parameters.
#[derive(Debug, Clone)]
pub struct LatticeParams {
    pub lx: usize,
    pub ly: usize,
}

impl LatticeParams {
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

    pub fn n_sites(&self) -> usize {
        self.lx * self.ly
    }
}
