//! Errors reported while constructing or validating particle systems.

use std::fmt::{Display, Formatter};

/// Configuration and cache errors for continuous particle systems.
#[derive(Debug, Clone, PartialEq)]
pub enum ParticleError {
    /// Particle simulations require at least one spatial dimension.
    ZeroDimension,
    /// A simulation-cell length was not finite and strictly positive.
    InvalidCellLength { axis: usize, length: f64 },
    /// The product of valid side lengths was not a finite cell volume.
    NonFiniteCellVolume,
    /// Position and species buffers have different lengths.
    BufferLengthMismatch { positions: usize, species: usize },
    /// A particle coordinate was not finite.
    NonFinitePosition { particle: usize, axis: usize },
    /// A pair-potential parameter was invalid.
    InvalidPotential(String),
    /// A particle species is unsupported by the selected potential.
    UnsupportedSpecies { particle: usize, species: u16 },
    /// The cutoff exceeds the minimum-image limit of the simulation cell.
    CutoffTooLarge { cutoff: f64, maximum: f64 },
    /// The accepted configuration has a non-finite physical energy.
    NonFiniteAcceptedEnergy,
    /// Cached and exact physical energies differ beyond roundoff.
    EnergyCacheMismatch { cached: f64, exact: f64 },
    /// The packed cell-list cache is inconsistent with the configuration.
    InvalidCellList(String),
    /// A particle or molecule move was malformed.
    InvalidMove(String),
    /// Molecular topology data was inconsistent.
    InvalidTopology(String),
    /// A weighted move mixture was invalid.
    InvalidMoveMixture(String),
}

impl Display for ParticleError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ZeroDimension => write!(formatter, "particle dimension must be positive"),
            Self::InvalidCellLength { axis, length } => write!(
                formatter,
                "cell length on axis {axis} must be positive and yield finite periodic geometry, got {length}"
            ),
            Self::NonFiniteCellVolume => write!(formatter, "simulation-cell volume is non-finite"),
            Self::BufferLengthMismatch { positions, species } => write!(
                formatter,
                "position/species length mismatch: {positions} positions, {species} species"
            ),
            Self::NonFinitePosition { particle, axis } => write!(
                formatter,
                "particle {particle} coordinate {axis} is non-finite"
            ),
            Self::InvalidPotential(reason) => write!(formatter, "invalid pair potential: {reason}"),
            Self::UnsupportedSpecies { particle, species } => write!(
                formatter,
                "particle {particle} uses unsupported species {species}"
            ),
            Self::CutoffTooLarge { cutoff, maximum } => write!(
                formatter,
                "cutoff {cutoff} exceeds minimum-image maximum {maximum}"
            ),
            Self::NonFiniteAcceptedEnergy => {
                write!(formatter, "accepted particle configuration has non-finite energy")
            }
            Self::EnergyCacheMismatch { cached, exact } => write!(
                formatter,
                "particle energy cache mismatch: cached {cached}, exact {exact}"
            ),
            Self::InvalidCellList(reason) => write!(formatter, "invalid cell list: {reason}"),
            Self::InvalidMove(reason) => write!(formatter, "invalid particle move: {reason}"),
            Self::InvalidTopology(reason) => write!(formatter, "invalid molecular topology: {reason}"),
            Self::InvalidMoveMixture(reason) => write!(formatter, "invalid move mixture: {reason}"),
        }
    }
}

impl std::error::Error for ParticleError {}
