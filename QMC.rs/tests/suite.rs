//! QMC.rs integration test suite.
//!
//! Each module lives in a subdirectory and is wired via explicit `#[path]`
//! attributes so the directory tree mirrors the test taxonomy without
//! requiring `mod.rs` files.

// ── Lattice — continuous-time lattice QMC ────────────────────────────────

#[path = "lattice/lattice_continuous.rs"]
mod lattice_continuous;
#[path = "lattice/lattice_ed.rs"]
mod lattice_ed;
#[path = "lattice/lattice_limits.rs"]
mod lattice_limits;

// ── Impurity — spin-boson wormhole QMC ───────────────────────────────────

#[path = "impurity/occupation.rs"]
mod occupation;
#[path = "impurity/physics.rs"]
mod physics;
#[path = "impurity/wormhole.rs"]
mod wormhole;
