//! Variational quantum Monte Carlo for continuum systems (L0 + L1).
//!
//! Layer plan and contracts live in `research/vmc/DESIGN.md` (authoritative).
//! This module family is the L0 layer: three nodeless trial states
//! ([`GaussianTrap`], [`McMillanJastrow`], [`HarmonicJastrow`]) with
//! hand-derived analytic `∇ ln|ψ|`, `∇² ln|ψ|` and parameter gradients, a
//! single-particle Metropolis [`VmcKernel`] hosting a walker population as
//! solver-internal state inside the Carlo.rs `sweep`/`measure` contract,
//! and versioned JSON checkpoints (`qmc-rs-vmc-v1`). L1 adds the fermionic
//! family behind the same [`WaveFunction`] trait: two-spin-block
//! [`SlaterDeterminant`]s of contracted cartesian Gaussians with
//! Sherman–Morrison single-particle updates and the
//! Kwon–Ceperley–Martin [`Backflow`] quasiparticle transformation.
//!
//! Layer boundaries (the anti-debt firewall): all wave-function physics in
//! `wavefunction/`, all sampling statistics in `kernel/`, all potential
//! physics in `hamiltonian.rs`, all estimator combination logic in
//! `estimators/`, all Carlo.rs glue in `adapters/`. Parameter optimizers
//! (L2), DMC (L3) and NQS (L5) plug in later behind the same
//! [`WaveFunction`] trait without touching this layer.

pub mod adapters;
pub mod error;
pub mod estimators;
pub mod hamiltonian;
pub mod kernel;
pub mod optimizer;
pub mod wavefunction;

pub use error::VariationalError;
pub use estimators::{local_energy, LocalEnergy};
pub use hamiltonian::{ContinuumHamiltonian, HarmonicTrap, PairPotential};
pub use kernel::{
    DmcKernel, DmcStats, VmcKernel, VmcStats, Walker, DMC_CHECKPOINT_FORMAT, VMC_CHECKPOINT_FORMAT,
};
pub use optimizer::{
    BlockStats, LinearMethod, Optimizer, ReferenceSample, StochasticReconfiguration,
    VarianceMinimization, VarianceMinimizationResult, VarianceObjective,
};
pub use wavefunction::{
    harmonic_closed_shell_electrons, harmonic_closed_shell_energy, harmonic_trap_orbitals,
    Backflow, DeltaLog, GaussianTrap, GradBuffer, GtoOrbital, HarmonicJastrow, McMillanJastrow,
    ParamGradBuffer, Point, Positions, Product, SlaterDeterminant, WaveFunction,
    WaveFunctionParams, DIM,
};
