//! QMC.rs integration test suite.
//!
//! Each module lives in a subdirectory and is wired via explicit `#[path]`
//! attributes so the directory tree mirrors the test taxonomy without
//! requiring `mod.rs` files.

// ── Shared seed management for multi-seed z-score tests ──────────────────

#[path = "zscore_seeds.rs"]
mod zscore_seeds;

// ── Lattice — continuous-time lattice QMC ────────────────────────────────

#[path = "lattice/lattice_continuous.rs"]
mod lattice_continuous;
#[path = "lattice/lattice_ed.rs"]
mod lattice_ed;
#[path = "lattice/lattice_ergodicity.rs"]
mod lattice_ergodicity;
#[path = "lattice/lattice_limits.rs"]
mod lattice_limits;
#[path = "lattice/lattice_scattering_generic_s.rs"]
mod lattice_scattering_generic_s;
#[path = "lattice/lattice_spin1.rs"]
mod lattice_spin1;
#[path = "lattice/lattice_zscore.rs"]
mod lattice_zscore;

// ── Impurity — spin-boson wormhole QMC ───────────────────────────────────

#[path = "impurity/cluster_ergodicity.rs"]
mod cluster_ergodicity;
#[path = "impurity/cluster_ergodicity_ed.rs"]
mod cluster_ergodicity_ed;
#[path = "impurity/cluster_multimode.rs"]
mod cluster_multimode;
#[path = "impurity/cross_solver.rs"]
mod cross_solver;
#[path = "impurity/cross_solver_cluster.rs"]
mod cross_solver_cluster;
#[path = "impurity/cross_solver_numerical.rs"]
mod cross_solver_numerical;
#[path = "impurity/ergodicity.rs"]
mod impurity_ergodicity;
#[path = "impurity/occupation.rs"]
mod occupation;
#[path = "impurity/occupation_detailed_balance.rs"]
mod occupation_detailed_balance;
#[path = "impurity/occupation_ergodicity.rs"]
mod occupation_ergodicity;
#[path = "impurity/physics.rs"]
mod physics;
#[path = "impurity/rabi_long.rs"]
mod rabi_long;
#[path = "impurity/rabi_qpt.rs"]
mod rabi_qpt;
#[path = "impurity/thread_count.rs"]
mod thread_count;
#[path = "impurity/wormhole.rs"]
mod wormhole;
#[path = "impurity/wormhole_interacting_ed.rs"]
mod wormhole_interacting_ed;

// ── Variational — continuum VMC (L0) ─────────────────────────────────────

#[path = "variational/dmc.rs"]
mod variational_dmc;
#[path = "variational/machine_precision.rs"]
mod variational_machine_precision;
#[path = "variational/optimizer.rs"]
mod variational_optimizer;
#[path = "variational/software.rs"]
mod variational_software;
#[path = "variational/statistical.rs"]
mod variational_statistical;
