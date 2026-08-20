//! Continuum Hamiltonians sampled by the variational kernels.
//!
//! Units: ħ = m = 1 throughout the L0 family. Open box: no periodic
//! boundaries and no minimum-image convention (documented L0 debt; the
//! homogeneous-bulk layer will add them behind this same type).

use super::error::VariationalError;
use super::wavefunction::{read_particle, Point, DIM};

/// One-body harmonic trap `V_trap = Σ_i ½ ω² |r_i − r₀|²`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HarmonicTrap {
    omega: f64,
    center: Point,
}

impl HarmonicTrap {
    /// Construct with trap frequency `omega` (finite, > 0) and center `r₀`.
    pub fn new(omega: f64, center: Point) -> Result<Self, VariationalError> {
        VariationalError::require_positive("omega", omega)?;
        if !center.iter().all(|x| x.is_finite()) {
            return Err(VariationalError::invalid(
                "center",
                "trap center must be finite in every coordinate",
            ));
        }
        Ok(Self { omega, center })
    }

    /// Trap frequency `ω`.
    #[inline]
    pub const fn omega(&self) -> f64 {
        self.omega
    }

    /// Trap center `r₀`.
    #[inline]
    pub const fn center(&self) -> Point {
        self.center
    }

    /// Potential energy of one particle.
    #[inline]
    fn single_particle(&self, r: &Point) -> f64 {
        let d = [
            r[0] - self.center[0],
            r[1] - self.center[1],
            r[2] - self.center[2],
        ];
        0.5 * self.omega * self.omega * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2])
    }
}

/// Pair interaction `V_pair = Σ_{i<j} v(r_ij)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PairPotential {
    /// Lennard-Jones `4ε[(σ/r)¹² − (σ/r)⁶]` — the standard He-4 pair model.
    LennardJones {
        /// Well depth `ε` (finite, > 0).
        epsilon: f64,
        /// Core length `σ` (finite, > 0).
        sigma: f64,
    },
    /// Harmonic confinement `½ k r²` — the "right trap" that makes
    /// [`HarmonicJastrow`](super::wavefunction::HarmonicJastrow) an exact
    /// ground state for `k = 4a²N`.
    Harmonic {
        /// Spring constant `k` (finite, > 0).
        spring_constant: f64,
    },
}

impl PairPotential {
    /// Validate the variant fields (criterion G).
    pub fn validate(&self) -> Result<(), VariationalError> {
        match *self {
            Self::LennardJones { epsilon, sigma } => {
                VariationalError::require_positive("epsilon", epsilon)?;
                VariationalError::require_positive("sigma", sigma)?;
            }
            Self::Harmonic { spring_constant } => {
                VariationalError::require_positive("spring_constant", spring_constant)?;
            }
        }
        Ok(())
    }

    /// Pair energy `v(r)`.
    #[inline]
    pub fn energy(&self, r: f64) -> f64 {
        match *self {
            Self::LennardJones { epsilon, sigma } => {
                let x = (sigma / r).powi(6);
                4.0 * epsilon * (x * x - x)
            }
            Self::Harmonic { spring_constant } => 0.5 * spring_constant * r * r,
        }
    }

    /// Stable label for checkpoint snapshots.
    #[inline]
    pub const fn kind(&self) -> &'static str {
        match self {
            Self::LennardJones { .. } => "lennard_jones",
            Self::Harmonic { .. } => "harmonic",
        }
    }

    /// Variant parameters in a stable order (checkpoint snapshots, logging).
    /// Off the hot path; allocating here is fine.
    pub fn parameters(&self) -> Vec<f64> {
        match *self {
            Self::LennardJones { epsilon, sigma } => vec![epsilon, sigma],
            Self::Harmonic { spring_constant } => vec![spring_constant],
        }
    }
}

/// Continuum many-body Hamiltonian: optional trap plus optional pair terms.
///
/// `H = −½ Σ_i ∇_i² + Σ_i ½ ω² |r_i − r₀|² + Σ_{i<j} v(r_ij)`.
/// At least one term must be present: a free-particle-only Hamiltonian has
/// no normalizable sampling problem for these open-box ansätze.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ContinuumHamiltonian {
    trap: Option<HarmonicTrap>,
    pair: Option<PairPotential>,
}

impl ContinuumHamiltonian {
    /// Assemble from validated parts; rejects the empty Hamiltonian.
    pub fn new(
        trap: Option<HarmonicTrap>,
        pair: Option<PairPotential>,
    ) -> Result<Self, VariationalError> {
        if let Some(trap) = trap {
            HarmonicTrap::new(trap.omega, trap.center)?;
        }
        if let Some(pair) = pair {
            pair.validate()?;
        }
        if trap.is_none() && pair.is_none() {
            return Err(VariationalError::invalid(
                "hamiltonian",
                "at least one of trap or pair potential is required",
            ));
        }
        Ok(Self { trap, pair })
    }

    /// Trap-free pair Hamiltonian.
    pub fn pair_only(pair: PairPotential) -> Result<Self, VariationalError> {
        Self::new(None, Some(pair))
    }

    /// Pair-free trap Hamiltonian.
    pub fn trap_only(trap: HarmonicTrap) -> Result<Self, VariationalError> {
        Self::new(Some(trap), None)
    }

    /// The one-body trap, if present.
    #[inline]
    pub const fn trap(&self) -> Option<HarmonicTrap> {
        self.trap
    }

    /// The pair interaction, if present.
    #[inline]
    pub fn pair(&self) -> Option<&PairPotential> {
        self.pair.as_ref()
    }

    /// Total potential energy `V(cfg)` of any flat configuration. O(N) with
    /// a trap, O(N²) with a pair term; stack-only, no allocation.
    pub fn potential_energy(&self, cfg: &impl AsRef<[f64]>) -> f64 {
        let coords = cfg.as_ref();
        let n = coords.len() / DIM;
        let mut energy = 0.0;
        if let Some(trap) = self.trap {
            for index in 0..n {
                energy += trap.single_particle(&read_particle(cfg, index));
            }
        }
        if let Some(pair) = self.pair {
            for i in 0..n {
                for j in (i + 1)..n {
                    let a = read_particle(cfg, i);
                    let b = read_particle(cfg, j);
                    let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
                    let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                    energy += pair.energy(r);
                }
            }
        }
        energy
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variational::Positions;

    #[test]
    fn rejects_empty_and_invalid_hamiltonians() {
        assert!(ContinuumHamiltonian::new(None, None).is_err());
        assert!(HarmonicTrap::new(-1.0, [0.0; 3]).is_err());
        assert!(HarmonicTrap::new(1.0, [f64::NAN, 0.0, 0.0]).is_err());
        assert!(PairPotential::LennardJones {
            epsilon: 0.0,
            sigma: 1.0
        }
        .validate()
        .is_err());
        assert!(PairPotential::LennardJones {
            epsilon: 1.0,
            sigma: f64::INFINITY
        }
        .validate()
        .is_err());
        assert!(PairPotential::Harmonic {
            spring_constant: -1.0
        }
        .validate()
        .is_err());
    }

    #[test]
    fn lj_minimum_is_minus_epsilon() {
        // d/dr 4eps[(sigma/r)^12 - (sigma/r)^6] = 0 at r = 2^(1/6) sigma
        // with value -eps: every LJ pair is bounded below by -eps.
        let lj = PairPotential::LennardJones {
            epsilon: 1.5,
            sigma: 2.0,
        };
        let r_min = 2f64.powf(1.0 / 6.0) * 2.0;
        assert!((lj.energy(r_min) + 1.5).abs() < 1e-12);
        for r in [0.7 * r_min, r_min, 1.4 * r_min, 4.0] {
            assert!(lj.energy(r) >= -1.5 - 1e-12);
        }
    }

    #[test]
    fn trap_potential_energy_matches_closed_form() {
        let trap = HarmonicTrap::new(2.0, [0.5, -0.5, 0.0]).unwrap();
        let hamiltonian = ContinuumHamiltonian::trap_only(trap).unwrap();
        let cfg = Positions::from_flat(vec![1.5, 0.5, 0.0, -0.5, 0.5, 1.0]).unwrap();
        // Particle 1 sits at distance sqrt(2); particle 2 at sqrt(3).
        let expected = 0.5 * 4.0 * (2.0 + 3.0);
        assert!((hamiltonian.potential_energy(&cfg) - expected).abs() < 1e-12);
    }
}
