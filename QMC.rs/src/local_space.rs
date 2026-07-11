//! Local Hilbert-space abstractions.
//!
//! Continuous-time lattice QMC only needs a finite local basis and sparse
//! operator matrix elements. The update engine is therefore independent of
//! spin algebra. [`SpinSpace`] is the first production implementation; a
//! future fermionic local space can implement [`LocalHilbertSpace`] and reuse
//! the operator-catalog layer while supplying the required sign backend.

use thiserror::Error;

/// Compact local basis index.
pub type BasisState = u16;

/// Exchange statistics associated with a local space.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ParticleStatistics {
    /// Spin or other distinguishable finite-state degree of freedom.
    Spin,
    /// Bosonic occupation basis.
    Boson,
    /// Fermionic occupation basis. The current lattice engine rejects signed
    /// configurations but reserves this category for determinant/worldline
    /// backends.
    Fermion,
}

/// Error raised by a local Hilbert-space implementation.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum LocalSpaceError {
    /// Site index is invalid.
    #[error("site {site} is outside 0..{site_count}")]
    SiteOutOfRange {
        /// Invalid site.
        site: usize,
        /// Number of sites.
        site_count: usize,
    },
    /// Basis state is invalid for a site.
    #[error("state {state} is outside local dimension {dimension} at site {site}")]
    StateOutOfRange {
        /// Site index.
        site: usize,
        /// Invalid basis state.
        state: BasisState,
        /// Local dimension.
        dimension: usize,
    },
    /// A spin magnitude is invalid.
    #[error("2S must be positive")]
    InvalidSpin,
    /// Site-dependent spin data have the wrong length.
    #[error("expected {expected} site spin values, got {actual}")]
    SiteCountMismatch {
        /// Required count.
        expected: usize,
        /// Supplied count.
        actual: usize,
    },
}

/// Finite local basis required by the sparse-operator engine.
pub trait LocalHilbertSpace: Send + Sync {
    /// Number of sites represented by this object.
    fn site_count(&self) -> usize;
    /// Local dimension at a site.
    fn dimension(&self, site: usize) -> usize;
    /// Exchange statistics.
    fn statistics(&self) -> ParticleStatistics;
    /// Validate one local basis state.
    fn validate_state(&self, site: usize, state: BasisState) -> Result<(), LocalSpaceError> {
        if site >= self.site_count() {
            return Err(LocalSpaceError::SiteOutOfRange {
                site,
                site_count: self.site_count(),
            });
        }
        let dimension = self.dimension(site);
        if usize::from(state) >= dimension {
            return Err(LocalSpaceError::StateOutOfRange {
                site,
                state,
                dimension,
            });
        }
        Ok(())
    }
}

/// Site-resolved quantum-spin Hilbert space.
///
/// A local basis state `q=0..2S` represents `m=-S+q`. Storing `2S` as an
/// integer supports integer and half-integer spins without floating-point
/// ambiguity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SpinSpace {
    two_s: Vec<u16>,
}

impl SpinSpace {
    /// Uniform spin magnitude on all sites.
    pub fn uniform(site_count: usize, two_s: u16) -> Result<Self, LocalSpaceError> {
        if two_s == 0 {
            return Err(LocalSpaceError::InvalidSpin);
        }
        Ok(Self {
            two_s: vec![two_s; site_count],
        })
    }

    /// Site-dependent spin magnitudes.
    pub fn site_resolved(two_s: Vec<u16>) -> Result<Self, LocalSpaceError> {
        if two_s.is_empty() || two_s.contains(&0) {
            return Err(LocalSpaceError::InvalidSpin);
        }
        Ok(Self { two_s })
    }

    /// Validate that this space matches a graph.
    pub fn require_site_count(&self, expected: usize) -> Result<(), LocalSpaceError> {
        if self.two_s.len() != expected {
            return Err(LocalSpaceError::SiteCountMismatch {
                expected,
                actual: self.two_s.len(),
            });
        }
        Ok(())
    }

    /// Integer `2S` at one site.
    pub fn two_s(&self, site: usize) -> u16 {
        self.two_s[site]
    }

    /// Spin magnitude `S`.
    pub fn spin(&self, site: usize) -> f64 {
        0.5 * f64::from(self.two_s(site))
    }

    /// Integer `2m` for a basis state.
    pub fn m_twice(&self, site: usize, state: BasisState) -> i32 {
        2 * i32::from(state) - i32::from(self.two_s(site))
    }

    /// Magnetic quantum number `m`.
    pub fn m(&self, site: usize, state: BasisState) -> f64 {
        0.5 * f64::from(self.m_twice(site, state))
    }

    /// Matrix element of `S+` from `state` to `state+1`.
    pub fn raising_amplitude(&self, site: usize, state: BasisState) -> Option<f64> {
        if state >= self.two_s(site) {
            return None;
        }
        let s = self.spin(site);
        let m = self.m(site, state);
        Some((s * (s + 1.0) - m * (m + 1.0)).sqrt())
    }

    /// Matrix element of `S-` from `state` to `state-1`.
    pub fn lowering_amplitude(&self, site: usize, state: BasisState) -> Option<f64> {
        if state == 0 {
            return None;
        }
        let s = self.spin(site);
        let m = self.m(site, state);
        Some((s * (s + 1.0) - m * (m - 1.0)).sqrt())
    }

    /// Apply a worm step `delta=±1` to a basis state.
    pub fn shifted_state(&self, site: usize, state: BasisState, delta: i8) -> Option<BasisState> {
        match delta {
            1 if state < self.two_s(site) => Some(state + 1),
            -1 if state > 0 => Some(state - 1),
            _ => None,
        }
    }
}

impl LocalHilbertSpace for SpinSpace {
    fn site_count(&self) -> usize {
        self.two_s.len()
    }

    fn dimension(&self, site: usize) -> usize {
        usize::from(self.two_s[site]) + 1
    }

    fn statistics(&self) -> ParticleStatistics {
        ParticleStatistics::Spin
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arbitrary_spin_ladders_are_exact() {
        let space = SpinSpace::uniform(1, 4).expect("S=2");
        assert_eq!(space.dimension(0), 5);
        assert_eq!(space.m(0, 0), -2.0);
        assert_eq!(space.m(0, 4), 2.0);
        let amplitude = space.raising_amplitude(0, 2).expect("raise m=0");
        assert!((amplitude - 6.0_f64.sqrt()).abs() < 1.0e-12);
    }

    #[test]
    fn half_integer_spin_has_no_rounding_convention() {
        let space = SpinSpace::uniform(2, 3).expect("S=3/2");
        assert_eq!(space.m_twice(0, 0), -3);
        assert_eq!(space.m_twice(0, 3), 3);
    }
}
