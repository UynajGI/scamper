//! Variational Metropolis kernel: a walker population as solver-internal
//! state.
//!
//! One [`sweep`](VmcKernel::sweep_with_phase) is one epoch of single-particle
//! Metropolis updates over every walker in deterministic order (walker 0..W,
//! particle 0..N inside each walker), proposing symmetric uniform-cube
//! displacements of half-width `proposal_width` and accepting against
//! `|ψ_T|²`. All derivative physics (`delta_log`) and all potential physics
//! live in the [`WaveFunction`](crate::variational::WaveFunction) and
//! [`ContinuumHamiltonian`](crate::variational::ContinuumHamiltonian) the
//! kernel is generic over — static dispatch, no `dyn` on the hot path.
//!
//! # RNG discipline (CMC worm-ensemble pattern)
//!
//! Each sweep draws one `u64` salt per walker from the caller's stream and
//! derives that walker's stream through [`carlo_rs::RngStreamKey`] with the
//! walker index in the `replica` field and the current lifecycle phase
//! folded in. Streams are therefore domain-separated by walker identity,
//! not by execution order: a future intra-sweep rayon fan-out changes no
//! results, and nothing is hidden from a checkpointed context — a restored
//! run replays the exact same per-walker streams.
//!
//! # Hot-path budget
//!
//! Zero heap allocation per single-particle move: proposals and deltas live
//! on the stack, and the only owned scratch ([`GradBuffer`]) is allocated
//! once at construction and reused by every measurement.

use rand::{Rng, RngExt};
use rand_xoshiro::Xoshiro256PlusPlus;
use serde_json::{json, Value as Json};

use super::super::error::VariationalError;
use super::super::estimators::{local_energy, LocalEnergy};
use super::super::hamiltonian::ContinuumHamiltonian;
use super::super::wavefunction::{GradBuffer, Positions, WaveFunctionParams, DIM};

/// Checkpoint format tag for VMC population snapshots.
///
/// Unknown tags (future formats) are rejected loudly on load — the lesson of
/// the CMC worm v1/v2 history: silent format guessing corrupts runs.
pub const VMC_CHECKPOINT_FORMAT: &str = "qmc-rs-vmc-v1";

/// One walker: a configuration plus its cached `ln|ψ_T|`.
#[derive(Debug, Clone, PartialEq)]
pub struct Walker {
    cfg: Positions,
    log_psi: f64,
}

impl Walker {
    /// The walker's configuration.
    #[inline]
    pub fn configuration(&self) -> &Positions {
        &self.cfg
    }

    /// Cached `ln|ψ_T|(cfg)`, maintained incrementally by the kernel.
    #[inline]
    pub const fn log_psi(&self) -> f64 {
        self.log_psi
    }
}

/// Kernel run-counters snapshot.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct VmcStats {
    /// Completed sweeps (epochs over the full population).
    pub sweeps: u64,
    /// Proposed single-particle moves since construction.
    pub attempted_moves: u64,
    /// Accepted single-particle moves since construction.
    pub accepted_moves: u64,
}

impl VmcStats {
    /// Fraction of accepted single-particle moves (`0.0` before the first
    /// proposal).
    #[inline]
    pub fn acceptance_ratio(&self) -> f64 {
        if self.attempted_moves == 0 {
            0.0
        } else {
            self.accepted_moves as f64 / self.attempted_moves as f64
        }
    }
}

/// Metropolis-in-`|ψ_T|²` population kernel over a trial wave function.
///
/// Generic over any stateless-or-cached ansatz whose configurations are the
/// canonical flat [`Positions`] layout — the layout the walker population
/// stores. (An ansatz with an exotic `Config` type would carry its own
/// kernel; none is planned.)
pub struct VmcKernel<W: WaveFunctionParams<Config = Positions>> {
    wave_function: W,
    hamiltonian: ContinuumHamiltonian,
    walkers: Vec<Walker>,
    proposal_width: f64,
    grad_scratch: GradBuffer,
    stats: VmcStats,
    reported_attempts: u64,
    reported_accepts: u64,
}

impl<W: WaveFunctionParams<Config = Positions>> VmcKernel<W> {
    /// Build a population of `n_walkers` walkers with `n_particles` particles
    /// each, initial positions drawn uniformly from a cube of half-width
    /// `initial_spread` around the origin (thermalization is expected, as
    /// with any Markov sampler), proposing moves of uniform half-width
    /// `proposal_width`.
    ///
    /// Validation (criterion G — reject, never panic): zero walker or
    /// particle counts, non-finite or non-positive widths, non-finite
    /// wave-function parameters.
    pub fn new(
        wave_function: W,
        hamiltonian: ContinuumHamiltonian,
        n_walkers: usize,
        n_particles: usize,
        initial_spread: f64,
        proposal_width: f64,
        rng: &mut impl Rng,
    ) -> Result<Self, VariationalError> {
        if n_walkers == 0 {
            return Err(VariationalError::invalid("n_walkers", "must be at least 1"));
        }
        if n_particles == 0 {
            return Err(VariationalError::invalid(
                "n_particles",
                "must be at least 1",
            ));
        }
        VariationalError::require_positive("initial_spread", initial_spread)?;
        VariationalError::require_positive("proposal_width", proposal_width)?;
        let params = wave_function.param_values();
        if let Some(bad) = params.iter().position(|p| !p.is_finite()) {
            return Err(VariationalError::invalid(
                "params",
                format!("parameter {bad} is non-finite"),
            ));
        }

        let mut walkers = Vec::with_capacity(n_walkers);
        for _ in 0..n_walkers {
            let coords: Vec<f64> = (0..n_particles * DIM)
                .map(|_| rng.random_range(-initial_spread..initial_spread))
                .collect();
            let cfg = Positions::from_flat(coords)?;
            walkers.push(Walker {
                log_psi: wave_function.log_psi(&cfg),
                cfg,
            });
        }
        Ok(Self {
            wave_function,
            hamiltonian,
            walkers,
            proposal_width,
            grad_scratch: GradBuffer::new(n_particles),
            stats: VmcStats::default(),
            reported_attempts: 0,
            reported_accepts: 0,
        })
    }

    /// Number of walkers in the population.
    #[inline]
    pub fn n_walkers(&self) -> usize {
        self.walkers.len()
    }

    /// Number of particles per walker.
    #[inline]
    pub fn n_particles(&self) -> usize {
        self.grad_scratch.len() / DIM
    }

    /// The sampled walker population.
    #[inline]
    pub fn walkers(&self) -> &[Walker] {
        &self.walkers
    }

    /// The trial wave function.
    #[inline]
    pub fn wave_function(&self) -> &W {
        &self.wave_function
    }

    /// The Hamiltonian defining the local energy.
    #[inline]
    pub const fn hamiltonian(&self) -> &ContinuumHamiltonian {
        &self.hamiltonian
    }

    /// Proposal half-width.
    #[inline]
    pub const fn proposal_width(&self) -> f64 {
        self.proposal_width
    }

    /// Change the proposal half-width (validated).
    ///
    /// House discipline (see `UpdateSchedule`): adapt only during
    /// thermalization and freeze before measurement — changing the Markov
    /// kernel mid-measurement reweights sweep-sampled observables.
    pub fn set_proposal_width(&mut self, width: f64) -> Result<(), VariationalError> {
        VariationalError::require_positive("proposal_width", width)?;
        self.proposal_width = width;
        Ok(())
    }

    /// Run counters.
    #[inline]
    pub const fn stats(&self) -> VmcStats {
        self.stats
    }

    /// Attempted/accepted move counts accumulated since the last call
    /// (incremental feed for Carlo.rs's attempt clocks).
    pub fn drain_move_counters(&mut self) -> (u64, u64) {
        let attempts = self.stats.attempted_moves - self.reported_attempts;
        let accepts = self.stats.accepted_moves - self.reported_accepts;
        self.reported_attempts = self.stats.attempted_moves;
        self.reported_accepts = self.stats.accepted_moves;
        (attempts, accepts)
    }

    /// One sweep: one epoch of single-particle Metropolis per walker, on
    /// per-walker streams derived from `rng` (see the module docs).
    pub fn sweep_with_phase(&mut self, rng: &mut impl Rng, phase: carlo_rs::RngPhase) {
        let n_particles = self.n_particles();
        let width = self.proposal_width;
        let Self {
            wave_function,
            walkers,
            stats,
            ..
        } = self;
        for (index, walker) in walkers.iter_mut().enumerate() {
            let salt = rng.random::<u64>();
            let mut stream: Xoshiro256PlusPlus = carlo_rs::RngStreamKey::new(salt)
                .with_replica(index as u64)
                .with_phase(phase)
                .seeded();
            for particle in 0..n_particles {
                let old = walker.cfg.particle(particle);
                let mut new_pos = old;
                for coordinate in new_pos.iter_mut() {
                    *coordinate += stream.random_range(-width..width);
                }
                let delta = wave_function.delta_log(&walker.cfg, particle, &new_pos);
                stats.attempted_moves += 1;
                // Metropolis for |psi|^2 with a symmetric proposal:
                // accept iff ln(u) < 2 * (ln|psi'| - ln|psi|). A -inf
                // delta (proposal onto an occupied site) is rejected; a
                // hypothetical NaN delta also fails the comparison.
                if stream.random::<f64>().ln() < 2.0 * delta.log_ratio {
                    wave_function.commit_move(&mut walker.cfg, particle, &new_pos);
                    walker.log_psi += delta.log_ratio;
                    stats.accepted_moves += 1;
                }
            }
        }
        stats.sweeps += 1;
    }

    /// Evaluate the per-walker local energies of the current population,
    /// invoking `sink` once per walker (measurement path; reuses the internal
    /// gradient buffer, allocates nothing).
    pub fn measure_population(&mut self, mut sink: impl FnMut(LocalEnergy)) {
        let Self {
            wave_function,
            hamiltonian,
            walkers,
            grad_scratch,
            ..
        } = self;
        for walker in walkers.iter() {
            sink(local_energy(
                wave_function,
                hamiltonian,
                &walker.cfg,
                grad_scratch,
            ));
        }
    }

    /// Population-mean local energy (and drift norm squared) of the current
    /// walker configurations. Measurement-side convenience; the Carlo.rs
    /// adapter performs the same evaluation per walker into the context.
    pub fn population_mean_local_energy(&mut self) -> LocalEnergy {
        let mut total = LocalEnergy {
            value: 0.0,
            log_grad_squared: 0.0,
        };
        let n_walkers = self.n_walkers() as f64;
        self.measure_population(|sample| {
            total.value += sample.value;
            total.log_grad_squared += sample.log_grad_squared;
        });
        total.value /= n_walkers;
        total.log_grad_squared /= n_walkers;
        total
    }

    /// Versioned population snapshot (walkers + parameters + counters).
    ///
    /// Layout (`format: "qmc-rs-vmc-v1"`): `n_walkers`, `n_particles`,
    /// `proposal_width`, `sweeps`, `attempted_moves`, `accepted_moves`,
    /// `params` (ansatz parameters), `hamiltonian` (trap + pair identity,
    /// bit-checked on load) and `walkers` (positions plus the incrementally
    /// maintained cached `log_psi`). RNG and measurement state stay owned
    /// by the Carlo.rs context checkpoint, per the CMC worm convention.
    pub fn save_snapshot(&self) -> Json {
        let trap = self.hamiltonian.trap().map(|trap| {
            json!({
                "omega": trap.omega(),
                "center": trap.center(),
            })
        });
        let pair = self.hamiltonian.pair().map(|pair| {
            json!({
                "kind": pair.kind(),
                "parameters": pair.parameters(),
            })
        });
        let walkers = self
            .walkers
            .iter()
            .map(|walker| {
                json!({
                    "positions": walker.cfg.as_slice(),
                    "log_psi": walker.log_psi,
                })
            })
            .collect::<Vec<_>>();
        json!({
            "format": VMC_CHECKPOINT_FORMAT,
            "n_walkers": self.n_walkers(),
            "n_particles": self.n_particles(),
            "proposal_width": self.proposal_width,
            "sweeps": self.stats.sweeps,
            "attempted_moves": self.stats.attempted_moves,
            "accepted_moves": self.stats.accepted_moves,
            "params": self.wave_function.param_values(),
            "hamiltonian": {"trap": trap, "pair": pair},
            "walkers": walkers,
        })
    }

    /// Restore population, parameters and counters from a snapshot.
    ///
    /// Rejects loudly (never panics) on: unknown/missing format tag, wrong
    /// field types, walker/particle-count mismatch, proposal-width or
    /// parameter or Hamiltonian mismatch against the constructed kernel,
    /// non-finite stored values, and stored `log_psi` inconsistent with the
    /// recomputation from the stored positions (tamper detection).
    pub fn load_snapshot(&mut self, snapshot: &Json) -> Result<(), VariationalError> {
        if snapshot["format"].as_str() != Some(VMC_CHECKPOINT_FORMAT) {
            return Err(VariationalError::checkpoint(format!(
                "unknown VMC snapshot format {:?}, expected {VMC_CHECKPOINT_FORMAT:?}",
                snapshot["format"]
            )));
        }
        if required_usize(snapshot, "n_walkers")? != self.n_walkers()
            || required_usize(snapshot, "n_particles")? != self.n_particles()
        {
            return Err(VariationalError::checkpoint(
                "snapshot population shape mismatch",
            ));
        }
        require_same_f64(snapshot, "proposal_width", self.proposal_width)?;

        // Parameters must match the constructed ansatz bit-for-bit: silently
        // adopting snapshot parameters would mask loading a snapshot into a
        // differently parameterized run.
        let params = f64_array(&snapshot["params"], "params")?;
        let current = self.wave_function.param_values();
        if params != current {
            return Err(VariationalError::checkpoint(
                "snapshot parameter mismatch: construct the kernel with the \
                 snapshot's ansatz parameters before loading",
            ));
        }
        self.validate_hamiltonian_snapshot(&snapshot["hamiltonian"])?;

        let walker_json = snapshot["walkers"]
            .as_array()
            .ok_or_else(|| VariationalError::checkpoint("walkers must be an array"))?;
        if walker_json.len() != self.n_walkers() {
            return Err(VariationalError::checkpoint(
                "snapshot walker-count mismatch",
            ));
        }
        let expected_len = DIM * self.n_particles();
        let mut walkers = Vec::with_capacity(walker_json.len());
        for (index, entry) in walker_json.iter().enumerate() {
            let coords = f64_array(&entry["positions"], "walkers[].positions")?;
            if coords.len() != expected_len {
                return Err(VariationalError::checkpoint(format!(
                    "walker {index}: expected {expected_len} coordinates, got {}",
                    coords.len()
                )));
            }
            if let Some(bad) = coords.iter().position(|x| !x.is_finite()) {
                return Err(VariationalError::checkpoint(format!(
                    "walker {index}: coordinate {bad} is non-finite"
                )));
            }
            let stored_log_psi = entry["log_psi"].as_f64().ok_or_else(|| {
                VariationalError::checkpoint(format!("walker {index}: log_psi must be a number"))
            })?;
            if !stored_log_psi.is_finite() {
                return Err(VariationalError::checkpoint(format!(
                    "walker {index}: log_psi is non-finite"
                )));
            }
            let cfg = Positions::from_flat(coords)?;
            // Tamper detection: the cached value must agree with a fresh
            // full recompute (pure deterministic function of positions and
            // parameters) to tight relative tolerance.
            let recomputed = self.wave_function.log_psi(&cfg);
            let scale = stored_log_psi.abs().max(recomputed.abs()).max(1.0);
            if (stored_log_psi - recomputed).abs() > 1e-12 * scale {
                return Err(VariationalError::checkpoint(format!(
                    "walker {index}: stored log_psi {stored_log_psi} inconsistent with \
                     recomputed {recomputed}"
                )));
            }
            walkers.push(Walker {
                cfg,
                log_psi: stored_log_psi,
            });
        }

        self.stats = VmcStats {
            sweeps: required_u64(snapshot, "sweeps")?,
            attempted_moves: required_u64(snapshot, "attempted_moves")?,
            accepted_moves: required_u64(snapshot, "accepted_moves")?,
        };
        self.reported_attempts = self.stats.attempted_moves;
        self.reported_accepts = self.stats.accepted_moves;
        self.walkers = walkers;
        Ok(())
    }

    /// Bit-exact Hamiltonian identity check (CMC `validate_snapshot_model`
    /// style): a snapshot from a different trap or pair model is a loud
    /// error, not a silent state transplant.
    fn validate_hamiltonian_snapshot(&self, value: &Json) -> Result<(), VariationalError> {
        match self.hamiltonian.trap() {
            Some(trap) => {
                let trap_json = &value["trap"];
                if trap_json.is_null() {
                    return Err(VariationalError::checkpoint(
                        "snapshot has no trap but this kernel samples one",
                    ));
                }
                require_same_f64(trap_json, "omega", trap.omega())?;
                let center = f64_array(&trap_json["center"], "trap.center")?;
                if center.as_slice() != trap.center() {
                    return Err(VariationalError::checkpoint("trap.center mismatch"));
                }
            }
            None => {
                if !value["trap"].is_null() {
                    return Err(VariationalError::checkpoint(
                        "snapshot has a trap but this kernel samples none",
                    ));
                }
            }
        }
        match self.hamiltonian.pair() {
            Some(pair) => {
                let pair_json = &value["pair"];
                if pair_json["kind"].as_str() != Some(pair.kind()) {
                    return Err(VariationalError::checkpoint(
                        "snapshot pair-potential kind mismatch",
                    ));
                }
                let parameters = f64_array(&pair_json["parameters"], "pair.parameters")?;
                if parameters != pair.parameters() {
                    return Err(VariationalError::checkpoint(
                        "snapshot pair-potential parameters mismatch",
                    ));
                }
            }
            None => {
                if !value["pair"].is_null() {
                    return Err(VariationalError::checkpoint(
                        "snapshot has a pair potential but this kernel samples none",
                    ));
                }
            }
        }
        Ok(())
    }
}

fn required_u64(value: &Json, name: &str) -> Result<u64, VariationalError> {
    value[name]
        .as_u64()
        .ok_or_else(|| VariationalError::checkpoint(format!("{name} must be u64")))
}

fn required_usize(value: &Json, name: &str) -> Result<usize, VariationalError> {
    usize::try_from(required_u64(value, name)?)
        .map_err(|_| VariationalError::checkpoint(format!("{name} does not fit usize")))
}

fn require_same_f64(value: &Json, name: &str, expected: f64) -> Result<(), VariationalError> {
    if value[name].as_f64().map(f64::to_bits) == Some(expected.to_bits()) {
        Ok(())
    } else {
        Err(VariationalError::checkpoint(format!("{name} mismatch")))
    }
}

fn f64_array(value: &Json, name: &str) -> Result<Vec<f64>, VariationalError> {
    value
        .as_array()
        .ok_or_else(|| VariationalError::checkpoint(format!("{name} must be an array")))?
        .iter()
        .map(|entry| {
            entry
                .as_f64()
                .ok_or_else(|| VariationalError::checkpoint(format!("{name} has a non-number")))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variational::hamiltonian::HarmonicTrap;
    use crate::variational::wavefunction::GaussianTrap;
    use rand::SeedableRng;

    fn build_kernel(seed: u64, n_walkers: usize, n_particles: usize) -> VmcKernel<GaussianTrap> {
        let wave_function = GaussianTrap::new(0.5, [0.0; DIM]).unwrap();
        let hamiltonian =
            ContinuumHamiltonian::trap_only(HarmonicTrap::new(1.0, [0.0; DIM]).unwrap()).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
        VmcKernel::new(
            wave_function,
            hamiltonian,
            n_walkers,
            n_particles,
            1.5,
            0.5,
            &mut rng,
        )
        .expect("valid kernel inputs")
    }

    #[test]
    fn constructor_rejects_invalid_populations_and_widths() {
        let build = |walkers: usize, particles: usize, spread: f64, width: f64| {
            let wave_function = GaussianTrap::new(0.5, [0.0; DIM]).unwrap();
            let hamiltonian =
                ContinuumHamiltonian::trap_only(HarmonicTrap::new(1.0, [0.0; DIM]).unwrap())
                    .unwrap();
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(1);
            VmcKernel::new(
                wave_function,
                hamiltonian,
                walkers,
                particles,
                spread,
                width,
                &mut rng,
            )
        };
        assert!(build(0, 4, 1.0, 0.5).is_err());
        assert!(build(8, 0, 1.0, 0.5).is_err());
        for bad_width in [0.0, -0.5, f64::NAN, f64::INFINITY] {
            assert!(build(8, 4, 1.0, bad_width).is_err());
        }
        for bad_spread in [0.0, -1.0, f64::NAN] {
            assert!(build(8, 4, bad_spread, 0.5).is_err());
        }
        assert!(build(8, 4, 1.0, 0.5).is_ok());
        assert!(build_kernel(1, 4, 4).set_proposal_width(0.0).is_err());
        assert!(build_kernel(1, 4, 4).set_proposal_width(0.7).is_ok());
    }

    #[test]
    fn sweep_is_deterministic_and_updates_counters() {
        let mut left = build_kernel(0xA11CE, 4, 4);
        let mut right = build_kernel(0xA11CE, 4, 4);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x5EED);
        let mut replay = rng.clone();
        for _ in 0..25 {
            left.sweep_with_phase(&mut rng, carlo_rs::RngPhase::Thermalization);
            right.sweep_with_phase(&mut replay, carlo_rs::RngPhase::Thermalization);
        }
        for (a, b) in left.walkers().iter().zip(right.walkers().iter()) {
            assert_eq!(a.configuration().as_slice(), b.configuration().as_slice());
            assert_eq!(a.log_psi().to_bits(), b.log_psi().to_bits());
        }
        assert_eq!(left.stats(), right.stats());
        assert_eq!(left.stats().sweeps, 25);
        assert_eq!(left.stats().attempted_moves, 25 * 4 * 4);
        let (attempts, accepts) = left.drain_move_counters();
        assert_eq!(attempts, left.stats().attempted_moves);
        assert_eq!(accepts, left.stats().accepted_moves);
        assert_eq!(left.drain_move_counters(), (0, 0));
        assert!(left.stats().acceptance_ratio() > 0.0);
    }

    #[test]
    fn checkpoint_round_trip_rejects_unknown_format() {
        let mut kernel = build_kernel(0xBEEF, 3, 3);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(7);
        kernel.sweep_with_phase(&mut rng, carlo_rs::RngPhase::Measurement);
        let snapshot = kernel.save_snapshot();
        assert_eq!(snapshot["format"], VMC_CHECKPOINT_FORMAT);

        let mut restored = build_kernel(0xBEEF, 3, 3);
        restored.load_snapshot(&snapshot).unwrap();
        assert_eq!(restored.stats(), kernel.stats());
        for (a, b) in restored.walkers().iter().zip(kernel.walkers()) {
            assert_eq!(a.configuration().as_slice(), b.configuration().as_slice());
            assert_eq!(a.log_psi().to_bits(), b.log_psi().to_bits());
        }

        let mut wrong_tag = snapshot.clone();
        wrong_tag["format"] = json!("qmc-rs-vmc-v0");
        assert!(restored.load_snapshot(&wrong_tag).is_err());
        assert!(restored.load_snapshot(&json!({})).is_err());
        let mut tampered = snapshot;
        tampered["walkers"][0]["log_psi"] = json!(123.456);
        assert!(restored.load_snapshot(&tampered).is_err());
    }
}
