//! X4: dynamic percolation on the rejection-free event substrate.
//!
//! The absorbing-state contact process is the cleanest directed dynamic
//! percolation model: active vertices recover at rate `mu`, and each
//! active-inactive edge transmits at rate `lambda`. The all-inactive state is
//! absorbing, and the transition belongs to the directed-percolation
//! universality class (anchor: Hinrichsen, Adv. Phys. 49, 815 (2000),
//! archived). This module is feature-gated (`temporal`) and CPU-first.
//!
//! Scope of X4a: static topology with a contact process driven through the
//! existing [`crate::dynamics::RejectionFreeModel`] /
//! [`crate::dynamics::GillespieKernel`] substrate. Every rate is recomputed
//! per event, so the recorded event order is the direct-method reference
//! order. Temporal-network input (time-varying edge activity) and a
//! scheduler adapter remain X4b.
//!
//! Non-equilibrium caveats: the reported quantities are finite-time,
//! finite-size survival fractions and densities, never a
//! thermodynamic-limit order parameter. An absorbing run is a true zero,
//! but a surviving run at finite time is only a lower bracket on
//! quasi-stationary activity.

mod contact;

pub use contact::{
    contact_ensemble, contact_process_run, ContactEnsemble, ContactError, ContactPatch,
    ContactProcessModel, ContactRun, ContactState,
};

/// Alias used by the percolation module re-exports.
pub use contact::ContactError as TemporalError;
