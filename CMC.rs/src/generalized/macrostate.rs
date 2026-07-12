//! Reusable scalar macrostates for lattice and particle configurations.

use crate::lattice::interaction::Measurable;
use crate::lattice::state::System;
use crate::particle::ParticleSystem;

/// Extract a scalar macrostate from an accepted configuration.
pub trait Macrostate<State, Model>: Send + Sync {
    fn value(&self, state: &State, model: &Model) -> f64;
}

/// Cached physical energy.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct EnergyMacrostate;

impl<H> Macrostate<System, H> for EnergyMacrostate {
    #[inline]
    fn value(&self, state: &System, _model: &H) -> f64 {
        state.energy
    }
}

impl<const D: usize, P> Macrostate<ParticleSystem<D>, P> for EnergyMacrostate {
    #[inline]
    fn value(&self, state: &ParticleSystem<D>, _model: &P) -> f64 {
        state.energy
    }
}

/// Model-native magnetization magnitude.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct MagnetizationMacrostate;

impl<H: Measurable> Macrostate<System, H> for MagnetizationMacrostate {
    #[inline]
    fn value(&self, state: &System, model: &H) -> f64 {
        model.magnetization(&state.spins)
    }
}

/// Accepted particle count.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ParticleNumberMacrostate;

impl<const D: usize, P> Macrostate<ParticleSystem<D>, P> for ParticleNumberMacrostate {
    #[inline]
    fn value(&self, state: &ParticleSystem<D>, _model: &P) -> f64 {
        state.len() as f64
    }
}
