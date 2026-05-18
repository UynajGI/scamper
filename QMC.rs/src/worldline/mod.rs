//! Single-site worldline in imaginary time [0, β).
//!
//! A worldline is a piecewise-constant function `state(τ)` with kinks
//! at points where the state changes. Periodic boundary: state(0) = state(β).

pub mod continuous;
pub mod discrete;

pub use continuous::ContinuousWorldline;
pub use discrete::DiscreteWorldline;

/// State index: 0..dim, where dim is the site's Hilbert space dimension.
/// Physical mapping (e.g. Sz = -1/0/+1 ↔ idx 0/1/2) is handled by the model layer.
type State = u8;

/// A worldline for a single site in imaginary time [0, β).
pub trait Worldline {
    /// Inverse temperature β.
    fn beta(&self) -> f64;

    /// Hilbert space dimension at this site.
    fn dim(&self) -> u8;

    /// Number of kinks (state-change points) in the worldline.
    fn num_kinks(&self) -> usize;

    /// State at imaginary time τ ∈ [0, β).
    fn state_at(&self, tau: f64) -> State;

    /// Iterate kinks with zero-copy callback: (tau, from_state, to_state).
    fn for_each_kink(&self, f: impl FnMut(f64, State, State));

    /// Insert a kink at τ, flipping to `to`. If state_at(τ) == to, this is a no-op.
    fn insert_kink(&mut self, tau: f64, to: State);

    /// Remove the kink at index `idx` (panic if out of bounds).
    fn remove_kink(&mut self, idx: usize);

    /// Time-averaged state: ∫₀ᵝ state(τ) dτ / β.
    fn diagonal(&self) -> f64;

    /// Build a Vec of all kinks. Allocates — prefer `for_each_kink` in hot paths.
    fn kinks_vec(&self) -> Vec<(f64, State, State)> {
        let mut v = Vec::with_capacity(self.num_kinks());
        self.for_each_kink(|tau, from, to| v.push((tau, from, to)));
        v
    }
}
