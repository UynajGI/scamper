//! Cross-kernel cache auditing.
//!
//! Debug and test builds audit automatically at a conservative cadence. The
//! `cache-audit` feature enables the same checks in optimized builds. A
//! non-zero per-kernel interval always requests explicit auditing.

use crate::generalized::MacrostateAxis;
use crate::lattice::interaction::Hamiltonian;
use crate::lattice::state::System;
use crate::particle::{PairPotential, ParticleSystem};

/// Default number of completed sweeps between automatic audits.
pub const DEFAULT_CACHE_AUDIT_INTERVAL: u64 = 1_024;

/// Whether automatic audits are enabled for this build.
#[inline]
pub const fn automatic_cache_audit_enabled() -> bool {
    cfg!(any(debug_assertions, test, feature = "cache-audit"))
}

/// Resolve an explicit interval or the build-mode default.
#[inline]
pub const fn effective_cache_audit_interval(configured: u64) -> u64 {
    if configured > 0 {
        configured
    } else if automatic_cache_audit_enabled() {
        DEFAULT_CACHE_AUDIT_INTERVAL
    } else {
        0
    }
}

/// Return true after the requested number of completed sweeps.
#[inline]
pub fn should_audit_cache(completed_sweeps: u64, configured: u64) -> bool {
    let interval = effective_cache_audit_interval(configured);
    interval > 0 && completed_sweeps > 0 && completed_sweeps.is_multiple_of(interval)
}

/// Validate graph/configuration shape and cached physical energy.
pub fn audit_lattice_cache<H: Hamiltonian>(system: &System, model: &H) -> Result<(), String> {
    system.validate(model)?;
    let exact = model.compute_total_energy(&system.spins, &system.lattice, system.beta);
    let tolerance = 1e-10 * (1.0 + exact.abs());
    if (system.energy - exact).abs() > tolerance {
        return Err(format!(
            "lattice energy cache mismatch: cached={}, exact={exact}",
            system.energy
        ));
    }
    Ok(())
}

/// Validate accepted coordinates, species, cell membership and cached energy.
pub fn audit_particle_cache<const D: usize, P: PairPotential>(
    system: &ParticleSystem<D>,
    potential: &P,
) -> Result<(), String> {
    system
        .validate(potential)
        .map_err(|error| error.to_string())
}

/// Validate a cached scalar macrostate bin against the accepted value.
pub fn audit_macrostate_bin<A: MacrostateAxis>(
    axis: &A,
    accepted_value: f64,
    cached_bin: usize,
) -> Result<(), String> {
    if cached_bin >= axis.bins() {
        return Err(format!(
            "macrostate cache bin {cached_bin} is outside {} bins",
            axis.bins()
        ));
    }
    let exact_bin = axis.bin(accepted_value).ok_or_else(|| {
        format!("accepted macrostate {accepted_value} lies outside its configured axis")
    })?;
    if exact_bin != cached_bin {
        return Err(format!(
            "macrostate cache mismatch: cached bin={cached_bin}, exact bin={exact_bin}"
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{build_chain, DiscreteAxis, IsingModel, System};

    #[test]
    fn automatic_policy_is_enabled_in_tests() {
        assert!(automatic_cache_audit_enabled());
        assert_eq!(
            effective_cache_audit_interval(0),
            DEFAULT_CACHE_AUDIT_INTERVAL
        );
        assert_eq!(effective_cache_audit_interval(7), 7);
    }

    #[test]
    fn lattice_audit_detects_energy_drift_without_repairing_it() {
        let model = IsingModel::new(1.0);
        let mut system = System::new(build_chain(4, true), 1, 1.0, 1.0);
        system.recompute_energy(&model);
        assert!(audit_lattice_cache(&system, &model).is_ok());
        let corrupted = system.energy + 1.0;
        system.energy = corrupted;
        assert!(audit_lattice_cache(&system, &model).is_err());
        assert_eq!(system.energy, corrupted);
    }

    #[test]
    fn macrostate_audit_detects_stale_bin() {
        let axis = DiscreteAxis::new(vec![-4.0, 0.0, 4.0]).unwrap();
        assert!(audit_macrostate_bin(&axis, 0.0, 1).is_ok());
        assert!(audit_macrostate_bin(&axis, 0.0, 0).is_err());
    }
}
