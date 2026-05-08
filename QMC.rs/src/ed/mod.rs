//! Exact diagonalization for QMC validation.
//!
//! Provides a sparse Hamiltonian builder and Lanczos eigensolver
//! to compute exact ground state energies for small systems.
//! Used to validate SSE results on N <= 16 sites.

mod hamiltonian;
mod lanczos;

pub use hamiltonian::SparseHamiltonian;
pub use lanczos::lanczos_ground_state;
