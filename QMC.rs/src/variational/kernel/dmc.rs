//! Diffusion Monte Carlo (L3): drift-diffusion + branching + population
//! control, with forward-walking pure estimators.
//!
//! The walker population is solver-internal state inside the Carlo.rs
//! `sweep`/`measure` contract (DESIGN.md §3): one `sweep` is one
//! imaginary-time step over the whole population. For a time step `tau`:
//!
//! 1. **Drift-diffusion move** per walker, `R' = R + tau * b(R) + sqrt(tau)
//!    * chi` with drift `b = grad ln|psi_T|` (the importance-transformed
//!    Fokker-Planck drift for kinetic `-1/2 lap^2`; the stationary density
//!    is `psi_T^2`, the same target VMC samples) and Gaussian `chi`,
//!    Metropolis-accepted against the Green-function ratio
//!    `ln A = 2(ln|psi(R')| - ln|psi(R)|) + (|D_fwd|^2 - |D_bwd|^2)/(2 tau)`
//!    with `D_fwd = R' - R - tau b(R)` and `D_bwd = R - R' - tau b(R')`
//!    — the standard adjustable form (Umrigar–Nightingale–Runge 1993).
//!    Per-particle drift norms are clamped at `CLAMP_COEFF/sqrt(tau)`
//!    (UMR finite-`tau` stabilization; inactive in the validated regime).
//! 2. **Branching**: each moved walker is replaced by `floor(g + u)`
//!    copies, `g = exp(-tau (E_L - E_T))`, `u` uniform — walkers in
//!    low-energy regions multiply.
//! 3. **Population control**: `E_T = E_ref - ln(N/N_target)/tau` with
//!    `E_ref` an exponential moving average of the population-mean
//!    `E_L` (classic feedback; the mixed estimator is the pre-branching
//!    population mean).
//! 4. **Pure estimators by forward walking**: every walker carries its
//!    lineage as a ring of the last `n_delay` ids; a measurement `(id,
//!    value)` recorded at step `s` is credited to that id when it leaves
//!    every walker's ring at `s + n_delay` — i.e. once per descendant,
//!    which is exactly descendant weighting.
//!
//! Checkpoints (`qmc-rs-dmc-v1`) serialize walkers, lineage rings, the
//! pending-measurement ring and the energy-feedback state so restored
//! runs replay bit-identically.

use std::collections::VecDeque;

use super::super::error::VariationalError;
use super::super::estimators::local_energy;
use super::super::hamiltonian::ContinuumHamiltonian;
use super::super::wavefunction::{GradBuffer, Positions, WaveFunctionParams, DIM};
use carlo_rs::RngStreamKey;
use rand::{Rng, RngExt};
use rand_xoshiro::Xoshiro256PlusPlus;

/// Per-particle drift clamp `|b| <= CLAMP_COEFF / sqrt(tau)` (inactive for
/// the smooth Gaussian states of the validated domain at small `tau`).
const CLAMP_COEFF: f64 = 5.0;

/// Population-safety cap on single-walker branching (`floor(g+u)` copies).
const MAX_BRANCH: u64 = 10;

/// EMA rate of the reference energy feedback.
const EMA_RATE: f64 = 0.05;

/// Checkpoint format tag.
pub const DMC_CHECKPOINT_FORMAT: &str = "qmc-rs-dmc-v1";

/// One DMC walker: configuration, cached amplitude, drift, and the
/// lineage ring for forward-walking estimators.
#[derive(Debug, Clone)]
struct DmcWalker {
    cfg: Positions,
    log_psi: f64,
    /// Flat `3N` drift `grad ln|psi_T|` at `cfg`.
    drift: Vec<f64>,
    /// Local energy of the last moved configuration (the value the
    /// forward-walking ledger records for this walker's fresh id).
    last_energy: f64,
    /// Ancestry: one id per completed step (children copy the parent's
    /// ring and then push their own fresh id, so `front()` is the id of
    /// the ancestor `n_delay` steps back).
    id_history: VecDeque<u64>,
    /// Fresh id assigned at the end of the last step — the key the
    /// pending-measurement ring records values under.
    current_id: u64,
}

/// A measurement recorded for forward walking: `(walker id, value)`.
type PendingMeasurement = (u64, f64);

/// Population statistics of a DMC run.
#[derive(Debug, Clone, Copy, Default)]
pub struct DmcStats {
    pub steps: u64,
    /// Sum of pre-branching population-mean local energies (mixed
    /// estimator accumulators).
    pub energy_sum: f64,
    pub energy_squared_sum: f64,
    pub n_energy_samples: u64,
    /// Forward-walking pure estimator accumulators (credited `n_delay`
    /// steps after measurement, descendant-weighted).
    pub pure_energy_sum: f64,
    pub pure_energy_squared_sum: f64,
    pub n_pure_samples: u64,
    /// Branching bookkeeping.
    pub total_births: u64,
    pub total_deaths: u64,
    /// Drift-diffusion acceptance bookkeeping.
    pub attempted_moves: u64,
    pub accepted_moves: u64,
}

impl DmcStats {
    /// Mixed estimator of the energy and its naive stderr.
    pub fn mixed_energy(&self) -> f64 {
        self.energy_sum / self.n_energy_samples as f64
    }

    /// Pure (forward-walking, descendant-weighted) estimator of the
    /// energy; available once at least `n_delay + 1` steps have run.
    pub fn pure_energy(&self) -> f64 {
        self.pure_energy_sum / self.n_pure_samples as f64
    }
}

/// Diffusion-Monte-Carlo population kernel.
pub struct DmcKernel<W: WaveFunctionParams<Config = Positions>> {
    wave_function: W,
    hamiltonian: ContinuumHamiltonian,
    walkers: Vec<DmcWalker>,
    time_step: f64,
    target_population: usize,
    /// Forward-walking delay (steps between measurement and credit).
    n_delay: usize,
    reference_energy: f64,
    trial_energy: f64,
    /// Ring buffer of per-step measurements awaiting credit.
    pending: VecDeque<Vec<PendingMeasurement>>,
    next_id: u64,
    stats: DmcStats,
    grad_scratch: GradBuffer,
    /// Fresh proposal position scratch (`3N`).
    proposal_scratch: Vec<f64>,
}

impl<W: WaveFunctionParams<Config = Positions>> DmcKernel<W> {
    /// Build a DMC population of `n_walkers` walkers with `n_particles`
    /// particles, time step `tau`, forward-walking delay `n_delay`, and
    /// initial reference/trial energy estimate `initial_energy` (the
    /// population-mean `E_L` of the initial spread; pass a sensible
    /// value to keep the first branching steps tame).
    ///
    /// Validation (criterion G): zero walker/particle counts, non-finite
    /// or non-positive `tau`, `n_delay == 0`, non-finite energies,
    /// non-finite wave-function parameters.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        wave_function: W,
        hamiltonian: ContinuumHamiltonian,
        n_walkers: usize,
        n_particles: usize,
        time_step: f64,
        n_delay: usize,
        initial_energy: f64,
        initial_spread: f64,
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
        VariationalError::require_positive("time_step", time_step)?;
        if n_delay == 0 {
            return Err(VariationalError::invalid("n_delay", "must be at least 1"));
        }
        if !initial_energy.is_finite() {
            return Err(VariationalError::invalid(
                "initial_energy",
                "must be finite",
            ));
        }
        VariationalError::require_positive("initial_spread", initial_spread)?;
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
            Self::validate_walker(&wave_function, &cfg)?;
            walkers.push(Self::build_walker(&wave_function, &cfg, 0, 0));
        }

        Ok(Self {
            wave_function,
            hamiltonian,
            walkers,
            time_step,
            target_population: n_walkers,
            n_delay,
            reference_energy: initial_energy,
            trial_energy: initial_energy,
            pending: VecDeque::with_capacity(n_delay + 1),
            next_id: n_walkers as u64,
            stats: DmcStats::default(),
            grad_scratch: GradBuffer::new(n_particles),
            proposal_scratch: vec![0.0; n_particles * DIM],
        })
    }

    fn validate_walker(wave_function: &W, cfg: &Positions) -> Result<(), VariationalError> {
        let log_psi = wave_function.log_psi(cfg);
        if !log_psi.is_finite() {
            return Err(VariationalError::invalid(
                "walkers",
                "initial configuration has non-finite ln|psi| (singular determinant \
                 or ansatz/particle-count mismatch)",
            ));
        }
        Ok(())
    }

    fn build_walker(wave_function: &W, cfg: &Positions, id: u64, depth: usize) -> DmcWalker {
        let log_psi = wave_function.log_psi(cfg);
        let mut grad_scratch = GradBuffer::new(cfg.n_particles());
        grad_scratch.clear();
        wave_function.log_grad(cfg, &mut grad_scratch);
        DmcWalker {
            cfg: cfg.clone(),
            log_psi,
            drift: grad_scratch.as_slice().to_vec(),
            last_energy: 0.0,
            id_history: VecDeque::with_capacity(depth + 1),
            current_id: id,
        }
    }

    /// Number of live walkers.
    #[inline]
    pub fn population(&self) -> usize {
        self.walkers.len()
    }

    /// The live walker configurations (driver diagnostics and correlated
    /// sampling; mirrors `VmcKernel::walkers`).
    pub fn walker_configurations(&self) -> impl Iterator<Item = &Positions> {
        self.walkers.iter().map(|walker| &walker.cfg)
    }

    /// Time step in use.
    #[inline]
    pub const fn time_step(&self) -> f64 {
        self.time_step
    }

    /// Current reference and trial energies (feedback diagnostics).
    pub const fn energies(&self) -> (f64, f64) {
        (self.reference_energy, self.trial_energy)
    }

    /// Accumulated run statistics (mixed and pure estimators).
    pub const fn stats(&self) -> &DmcStats {
        &self.stats
    }

    /// One imaginary-time step: drift-diffuse every walker, branch,
    /// apply population control, advance the forward-walking ledger
    /// (module docs §1–4). Deterministic given `rng` — per-walker
    /// streams are keyed by (step salt, walker index, DMC domain tag).
    pub fn step(&mut self, rng: &mut impl Rng) -> Result<(), VariationalError> {
        let salt = rng.random::<u64>();
        let tau = self.time_step;
        let sqrt_tau = tau.sqrt();
        let clamp = CLAMP_COEFF / sqrt_tau;

        // 1. Drift-diffusion move + local energies (mixed measurement
        //    happens pre-branching on the moved walkers).
        let mut local_energies = Vec::with_capacity(self.walkers.len());
        for (index, walker) in self.walkers.iter_mut().enumerate() {
            let mut stream: Xoshiro256PlusPlus = RngStreamKey::new(salt)
                .with_replica(index as u64)
                .with_phase(carlo_rs::RngPhase::Measurement)
                .seeded();
            let n_flat = walker.cfg.as_slice().len();
            for coordinate in 0..n_flat {
                // Box-Muller on the walker's dedicated stream (crate
                // practice; see `determinant::random_unit_vector`).
                let u1 = stream.random::<f64>().max(f64::MIN_POSITIVE);
                let u2 = stream.random::<f64>();
                let radius = (-2.0 * u1.ln()).sqrt();
                let angle = std::f64::consts::TAU * u2;
                self.proposal_scratch[coordinate] = walker.cfg.as_slice()[coordinate]
                    + tau * walker.drift[coordinate].clamp(-clamp, clamp)
                    + sqrt_tau * radius * angle.cos();
            }
            let proposal =
                Positions::from_flat(self.proposal_scratch.clone()).expect("finite proposal");
            // Green-function-ratio Metropolis (module docs formula).
            let proposal_log_psi = self.wave_function.log_psi(&proposal);
            let mut proposal_drift = vec![0.0; n_flat];
            self.grad_scratch.clear();
            self.wave_function
                .log_grad(&proposal, &mut self.grad_scratch);
            proposal_drift.copy_from_slice(self.grad_scratch.as_slice());

            // Metropolis for the drift-diffusion Green function:
            // A = [psi^2(R') G(R'->R)] / [psi^2(R) G(R->R')], so
            // ln A = 2 (ln|psi'| - ln|psi|) + (|D_fwd|^2 - |D_bwd|^2) / (2 tau)
            // with D_fwd = R'-R - tau b(R) and D_bwd = R-R' - tau b(R')
            //      = -(R'-R) - tau b(R'), whose square is
            // (|R'-R| + tau b(R'))^2 — the PLUS sign on the backward
            // drift is load-bearing (the backward displacement is
            // R-R', not R'-R); flipping it was caught by the
            // separable-anchor test drifting the population onto a
            // wrong stationary manifold.
            let mut ln_accept = 2.0 * (proposal_log_psi - walker.log_psi);
            for ((&x_new, &x_old), (&b_old, &b_new)) in proposal
                .as_slice()
                .iter()
                .zip(walker.cfg.as_slice())
                .zip(walker.drift.iter().zip(&proposal_drift))
            {
                let displacement = x_new - x_old;
                let forward = displacement - tau * b_old;
                let backward = displacement + tau * b_new;
                ln_accept += (forward * forward - backward * backward) / (2.0 * tau);
            }
            self.stats.attempted_moves += 1;
            if stream.random::<f64>().ln() < ln_accept && proposal_log_psi.is_finite() {
                walker.cfg = proposal;
                walker.log_psi = proposal_log_psi;
                walker.drift = proposal_drift;
                self.stats.accepted_moves += 1;
            }

            let mut grad = std::mem::replace(
                &mut self.grad_scratch,
                GradBuffer::new(walker.cfg.n_particles()),
            );
            let sample = local_energy(
                &self.wave_function,
                &self.hamiltonian,
                &walker.cfg,
                &mut grad,
            );
            local_energies.push(sample.value);
            // `local_energy` leaves the fresh drift in the buffer (it is
            // the same field `log_grad` filled) — keep the walker's copy
            // in sync for the next step's proposal.
            walker.drift.copy_from_slice(grad.as_slice());
            walker.last_energy = sample.value;
            self.grad_scratch = grad;
        }

        // Mixed estimator on the moved, pre-branching population.
        let population_mean =
            local_energies.iter().sum::<f64>() / local_energies.len().max(1) as f64;
        self.stats.energy_sum += population_mean;
        self.stats.energy_squared_sum += population_mean * population_mean;
        self.stats.n_energy_samples += 1;

        // 2. Branching: floor(g + u) copies, g = exp(-tau (E_L - E_T)),
        // with a population-safety cap (`MAX_BRANCH`; the E_T feedback
        // keeps growth O(1) in the validated regime, where the cap is
        // never reached — asserted by the birth/death statistics).
        let mut branched: Vec<DmcWalker> = Vec::with_capacity(self.walkers.len() * 2);
        let mut branch_stream: Xoshiro256PlusPlus =
            RngStreamKey::new(salt.rotate_left(32) ^ 0xD0C0_FFEE).seeded();
        for (walker, &energy) in self.walkers.iter().zip(&local_energies) {
            if !energy.is_finite() {
                return Err(VariationalError::invalid(
                    "local_energy",
                    "a walker has non-finite E_L (ansatz left its physical domain)",
                ));
            }
            let growth = (-tau * (energy - self.trial_energy)).exp();
            let copies = growth + branch_stream.random::<f64>();
            let copies = if copies < 1.0 {
                0
            } else {
                (copies as u64).min(MAX_BRANCH)
            };
            match copies {
                0 => self.stats.total_deaths += 1,
                1 => branched.push(walker.clone()),
                n => {
                    self.stats.total_births += n - 1;
                    for _ in 0..n {
                        branched.push(walker.clone());
                    }
                }
            }
        }
        self.walkers = branched;
        if self.walkers.is_empty() {
            return Err(VariationalError::invalid(
                "population",
                "the walker population died out (E_T feedback diverged; check tau and \
                 the initial energy estimate)",
            ));
        }

        // 3. Population control feedback on the trial energy.
        self.reference_energy += EMA_RATE * (population_mean - self.reference_energy);
        let ratio = self.walkers.len() as f64 / self.target_population as f64;
        self.trial_energy = self.reference_energy - (ratio.ln() / tau).clamp(-50.0, 50.0);

        // 4. Forward-walking ledger. Fresh ids for every surviving
        // lineage; the value measured this step (the moved parent's E_L)
        // rides with each child under its fresh id.
        let mut this_step: Vec<PendingMeasurement> = Vec::with_capacity(self.walkers.len());
        for walker in &mut self.walkers {
            walker.current_id = self.next_id;
            self.next_id += 1;
            walker.id_history.push_back(walker.current_id);
            this_step.push((walker.current_id, walker.last_energy));
        }
        // Credit matured measurements: the record pushed exactly
        // `n_delay` steps ago leaves the ring, and every lineage whose
        // ring now exceeds `n_delay` ids pops the matching ancestor id —
        // once per surviving descendant line, which is descendant
        // weighting.
        let matured = if self.pending.len() >= self.n_delay {
            self.pending.pop_front()
        } else {
            None
        };
        if let Some(record) = matured {
            for walker in &mut self.walkers {
                while walker.id_history.len() > self.n_delay {
                    let popped = walker.id_history.pop_front().expect("lineage nonempty");
                    if let Some(&(_, value)) = record.iter().find(|(id, _)| *id == popped) {
                        self.stats.pure_energy_sum += value;
                        self.stats.pure_energy_squared_sum += value * value;
                        self.stats.n_pure_samples += 1;
                    }
                }
            }
        }
        self.pending.push_back(this_step);

        self.stats.steps += 1;
        Ok(())
    }
}

impl<W: WaveFunctionParams<Config = Positions>> DmcKernel<W> {
    /// Serialize the full DMC state (walkers, lineage rings, pending
    /// measurement ring, energy feedback, statistics) for exact replay.
    pub fn save_snapshot(&self) -> serde_json::Value {
        serde_json::json!({
            "format": DMC_CHECKPOINT_FORMAT,
            "time_step": self.time_step,
            "target_population": self.target_population,
            "n_delay": self.n_delay,
            "reference_energy": self.reference_energy,
            "trial_energy": self.trial_energy,
            "next_id": self.next_id,
            "walkers": self.walkers.iter().map(|walker| serde_json::json!({
                "configuration": walker.cfg.as_slice(),
                "id_history": walker.id_history.iter().collect::<Vec<_>>(),
                "current_id": walker.current_id,
            })).collect::<Vec<_>>(),
            "pending": self.pending.iter().map(|record| record.iter()
                .map(|(id, value)| serde_json::json!({"id": id, "value": value}))
                .collect::<Vec<_>>()).collect::<Vec<_>>(),
            "stats": serde_json::json!({
                "steps": self.stats.steps,
                "energy_sum": self.stats.energy_sum,
                "energy_squared_sum": self.stats.energy_squared_sum,
                "n_energy_samples": self.stats.n_energy_samples,
                "pure_energy_sum": self.stats.pure_energy_sum,
                "pure_energy_squared_sum": self.stats.pure_energy_squared_sum,
                "n_pure_samples": self.stats.n_pure_samples,
                "total_births": self.stats.total_births,
                "total_deaths": self.stats.total_deaths,
                "attempted_moves": self.stats.attempted_moves,
                "accepted_moves": self.stats.accepted_moves,
            }),
        })
    }

    /// Restore from a snapshot written by [`save_snapshot`]. Rejects
    /// unknown formats, structural corruption, and particle-count
    /// mismatches loudly (`CheckpointCorrupted`), never panics on
    /// snapshot data.
    pub fn load_snapshot(&mut self, snapshot: &serde_json::Value) -> Result<(), VariationalError> {
        let corrupt = |detail: &str| VariationalError::CheckpointCorrupted {
            detail: detail.to_string(),
        };
        if snapshot.get("format").and_then(|f| f.as_str()) != Some(DMC_CHECKPOINT_FORMAT) {
            return Err(corrupt("unknown or missing format tag"));
        }
        let n_flat_expected = self.grad_scratch.as_slice().len();
        let walkers = snapshot
            .get("walkers")
            .and_then(|w| w.as_array())
            .ok_or_else(|| corrupt("missing walkers array"))?;
        let mut restored = Vec::with_capacity(walkers.len());
        for walker in walkers {
            let configuration = walker
                .get("configuration")
                .and_then(|c| c.as_array())
                .ok_or_else(|| corrupt("missing walker configuration"))?;
            let coords: Result<Vec<f64>, _> = configuration
                .iter()
                .map(|x| x.as_f64().ok_or_else(|| corrupt("non-numeric coordinate")))
                .collect();
            let cfg =
                Positions::from_flat(coords?).map_err(|_| corrupt("invalid configuration"))?;
            if cfg.as_slice().len() != n_flat_expected {
                return Err(corrupt("particle-count mismatch with this kernel"));
            }
            Self::validate_walker(&self.wave_function, &cfg)?;
            let mut built = Self::build_walker(&self.wave_function, &cfg, 0, 0);
            let history = walker
                .get("id_history")
                .and_then(|h| h.as_array())
                .ok_or_else(|| corrupt("missing id history"))?;
            for id in history {
                built.id_history.push_back(
                    id.as_u64()
                        .ok_or_else(|| corrupt("non-numeric lineage id"))?,
                );
            }
            built.current_id = walker
                .get("current_id")
                .and_then(|i| i.as_u64())
                .ok_or_else(|| corrupt("missing current id"))?;
            restored.push(built);
        }
        if restored.is_empty() {
            return Err(corrupt("empty walker population"));
        }
        let mut pending = VecDeque::with_capacity(self.n_delay + 1);
        let records = snapshot
            .get("pending")
            .and_then(|p| p.as_array())
            .ok_or_else(|| corrupt("missing pending ring"))?;
        for record in records {
            let entries = record
                .as_array()
                .ok_or_else(|| corrupt("malformed pending record"))?;
            let mut decoded = Vec::with_capacity(entries.len());
            for entry in entries {
                decoded.push((
                    entry
                        .get("id")
                        .and_then(|i| i.as_u64())
                        .ok_or_else(|| corrupt("malformed pending id"))?,
                    entry
                        .get("value")
                        .and_then(|v| v.as_f64())
                        .ok_or_else(|| corrupt("malformed pending value"))?,
                ));
            }
            pending.push_back(decoded);
        }

        let f64_of = |key: &str| -> Result<f64, VariationalError> {
            snapshot
                .get(key)
                .and_then(|v| v.as_f64())
                .ok_or_else(|| corrupt("missing numeric field"))
        };
        let stats = snapshot
            .get("stats")
            .ok_or_else(|| corrupt("missing stats"))?;
        let stat_of = |key: &str| -> Result<f64, VariationalError> {
            stats
                .get(key)
                .and_then(|v| v.as_f64())
                .ok_or_else(|| corrupt("missing stats field"))
        };

        self.time_step = f64_of("time_step")?;
        self.target_population = snapshot
            .get("target_population")
            .and_then(|t| t.as_u64())
            .ok_or_else(|| corrupt("missing target population"))?
            as usize;
        self.n_delay = snapshot
            .get("n_delay")
            .and_then(|t| t.as_u64())
            .ok_or_else(|| corrupt("missing n_delay"))? as usize;
        self.reference_energy = f64_of("reference_energy")?;
        self.trial_energy = f64_of("trial_energy")?;
        self.next_id = snapshot
            .get("next_id")
            .and_then(|t| t.as_u64())
            .ok_or_else(|| corrupt("missing next_id"))?;
        self.walkers = restored;
        self.pending = pending;
        self.stats = DmcStats {
            steps: stat_of("steps")? as u64,
            energy_sum: stat_of("energy_sum")?,
            energy_squared_sum: stat_of("energy_squared_sum")?,
            n_energy_samples: stat_of("n_energy_samples")? as u64,
            pure_energy_sum: stat_of("pure_energy_sum")?,
            pure_energy_squared_sum: stat_of("pure_energy_squared_sum")?,
            n_pure_samples: stat_of("n_pure_samples")? as u64,
            total_births: stat_of("total_births")? as u64,
            total_deaths: stat_of("total_deaths")? as u64,
            attempted_moves: stat_of("attempted_moves")? as u64,
            accepted_moves: stat_of("accepted_moves")? as u64,
        };
        Ok(())
    }
}
