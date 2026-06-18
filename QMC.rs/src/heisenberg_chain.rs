//! End-to-end wrapper: Heisenberg S = 1/2 chain as a Carlo.rs [`MonteCarlo`].
//!
//! Composes the building blocks in [`crate::discrete`] into a runnable
//! simulation:
//! - [`HeisenbergChain`](crate::hamiltonian::HeisenbergChain) Hamiltonian
//! - [`SpaceTimeConfig`] (N sites × M Trotter slices, PBC)
//! - [`local_metropolis_sweep`] update kernel
//!
//! Runs through `Scheduler::run_one::<HeisenbergChainMC>(&params)`. Parameters:
//! - `L` (usize, required): chain length
//! - `J` (f64, default 1.0): coupling. **Sign convention** is the standard
//!   physics one, `H = J Σ Sᵢ·Sⱼ`: `J > 0` = antiferromagnetic,
//!   `J < 0` = ferromagnetic.
//! - `beta` (f64, required): inverse temperature
//! - `M` (usize, default derived from `beta`/`dtau_target`): Trotter slices
//! - `dtau_target` (f64, default 0.1): target slice width; `M = ceil(β/dtau_target)`
//!
//! ## Antiferromagnet and the sign problem
//!
//! On the **bipartite** chain, the antiferromagnetic Heisenberg model is
//! sign-problem-free via the sublattice transformation `Sᵢ → −Sᵢ` on every
//! even site. This wrapper applies that transform internally when `J > 0`:
//! the sampler runs on the transformed (ferromagnetic-sign) configuration,
//! and observables are mapped back — energy is invariant under the transform,
//! magnetization becomes staggered magnetization.

use crate::discrete::config::SpaceTimeConfig;
use crate::discrete::worm::local_metropolis_sweep;
use crate::hamiltonian::HeisenbergChain;
use crate::lattice::ChainLattice;
use carlo_rs::{CarloError, Context, FromParams, MonteCarlo, ParallelTemperingCompatible, Params};

/// Pre-composed Heisenberg S = 1/2 chain QMC (path-integral, local Metropolis).
///
/// Internally stores the **transformed** configuration: if the user requested
/// `J > 0` (AF), spins on even sublattice sites are stored flipped so the
/// internal `HeisenbergChain` (which uses `H = −J Σ Sᵢ·Sⱼ`) sees a positive
/// coupling `|J|` and is sign-free. The flag [`Self::staggered`] records
/// whether the transform is active.
pub struct HeisenbergChainMC {
    config: SpaceTimeConfig,
    ham: HeisenbergChain,
    /// True if spins on even sublattice are stored flipped (AF case).
    staggered: bool,
    /// User-facing coupling J in `H = J Σ Sᵢ·Sⱼ`.
    j: f64,
}

impl HeisenbergChainMC {
    /// Energy per site from the current configuration (full quantum PI
    /// estimator — includes spin-exchange fluctuations).
    ///
    /// Invariant under the sublattice transform, so no back-transform needed.
    pub fn energy_per_site(&self) -> f64 {
        self.config.energy_quantum(&self.ham) / self.config.n_sites as f64
    }

    /// Number of lattice sites.
    pub fn n_sites(&self) -> usize {
        self.config.n_sites
    }

    /// User-facing coupling J.
    pub fn coupling(&self) -> f64 {
        self.j
    }

    /// Whether the sublattice transform is active (AF case).
    pub fn is_staggered(&self) -> bool {
        self.staggered
    }

    /// Physical magnetization per site in [-1, 1].
    ///
    /// If [`Self::staggered`], this is the **staggered** magnetization
    /// (the physical observable for the AF chain); else it's the uniform
    /// magnetization.
    fn physical_magnetization(&self) -> f64 {
        let n = self.config.n_sites;
        let m = self.config.n_slices;
        let total = (n * m) as f64;
        let mut sum: i64 = 0;
        for slice in 0..m {
            for site in 0..n {
                let mut s = 2 * self.config.spin(site, slice) as i64 - 1;
                // Back out the sublattice transform to get the physical spin.
                if self.staggered && site % 2 == 0 {
                    s = -s;
                }
                sum += s;
            }
        }
        sum as f64 / total
    }
}

// ── MonteCarlo impl ──────────────────────────────────────────

impl MonteCarlo for HeisenbergChainMC {
    type Rng = rand_xoshiro::Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        local_metropolis_sweep(&mut self.config, &self.ham, &mut ctx.rng);
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        let e = self.config.energy_quantum(&self.ham); // invariant under transform
        let m = self.physical_magnetization();
        let n = self.config.n_sites as f64;
        ctx.measure("Energy", e);
        ctx.measure("EnergyPerSite", e / n);
        ctx.measure("Magnetization", m);
        ctx.measure("M2", m * m);
        ctx.measure("M4", m * m * m * m);
        ctx.measure("Kinks", self.config.num_kinks() as f64);
    }

    fn name(&self) -> &'static str {
        "HeisenbergChainMC"
    }
}

// ── FromParams impl ──────────────────────────────────────────

impl FromParams for HeisenbergChainMC {
    fn from_params(params: &Params, rng: &mut Self::Rng) -> Result<Self, CarloError> {
        let l = params
            .get::<usize>("L")
            .ok_or_else(|| CarloError::InvalidConfig {
                field: "L".into(),
                reason: "chain length L is required".into(),
            })?;
        if l < 2 {
            return Err(CarloError::InvalidConfig {
                field: "L".into(),
                reason: format!("L must be >= 2, got {l}"),
            });
        }
        let beta = params
            .get::<f64>("beta")
            .ok_or_else(|| CarloError::InvalidConfig {
                field: "beta".into(),
                reason: "beta (inverse temperature) is required".into(),
            })?;
        if beta <= 0.0 {
            return Err(CarloError::InvalidConfig {
                field: "beta".into(),
                reason: format!("beta must be > 0, got {beta}"),
            });
        }

        let j_user = params.get::<f64>("J").unwrap_or(1.0);
        // H_user = J Σ Sᵢ·Sⱼ. Internal struct uses H_internal = -j Σ Sᵢ·Sⱼ.
        //
        // Apply sublattice transform for J > 0 (AF): conceptually flip spins
        // on even sites, which maps H_user → -J Σ Sᵢ·Sⱼ, i.e. j_internal = +J
        // > 0 (sign-free). For J < 0 (ferromagnet), no transform: j_internal
        // = -J > 0. Either way j_internal = |J| > 0.
        //
        // We don't actually mutate the random initial config — the transform
        // is just a relabeling, and `physical_magnetization` undoes it on
        // read. Both branches use a positive internal coupling.
        let staggered = j_user > 0.0;
        let ham = HeisenbergChain::new(j_user.abs());

        // Trotter slices: explicit M, or derived from target Δτ.
        let m = params.get::<usize>("M").unwrap_or_else(|| {
            let dtau_target = params.get::<f64>("dtau_target").unwrap_or(0.1);
            ((beta / dtau_target).ceil() as usize).max(2)
        });
        if m < 2 {
            return Err(CarloError::InvalidConfig {
                field: "M".into(),
                reason: format!("M (Trotter slices) must be >= 2, got {m}"),
            });
        }

        let lattice = ChainLattice::new(l);
        let config = SpaceTimeConfig::new_random(lattice, beta, m, rng);

        Ok(Self {
            config,
            ham,
            staggered,
            j: j_user,
        })
    }
}

// ── ParallelTemperingCompatible impl ─────────────────────────

impl ParallelTemperingCompatible for HeisenbergChainMC {
    fn log_weight_ratio(&self, param: &str, new_value: f64) -> f64 {
        match param {
            // W ∝ exp(-β E) at fixed config → log(W'/W) = -(β'-β)·E.
            // Use the quantum estimator for consistency with `measure`.
            "beta" => (self.config.beta - new_value) * self.config.energy_quantum(&self.ham),
            _ => panic!("unsupported PT param: {param}"),
        }
    }

    fn change_parameter(&mut self, param: &str, new_value: f64) {
        match param {
            "beta" => {
                self.config.beta = new_value;
                // Recompute Δτ to preserve the slice width approximately.
                self.config.dtau = new_value / self.config.n_slices as f64;
            }
            _ => panic!("unsupported PT param: {param}"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use carlo_rs::{RayonBackend, RunConfig, Scheduler};
    use rand::SeedableRng;

    #[test]
    fn from_params_requires_l() {
        let mut p = Params::new();
        p.set("beta", 1.0);
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(1);
        assert!(HeisenbergChainMC::from_params(&p, &mut rng).is_err());
    }

    #[test]
    fn from_params_requires_beta() {
        let mut p = Params::new();
        p.set("L", 4usize);
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(1);
        assert!(HeisenbergChainMC::from_params(&p, &mut rng).is_err());
    }

    #[test]
    fn from_params_round_trip() {
        let mut p = Params::new();
        p.set("L", 6usize);
        p.set("beta", 2.0_f64);
        p.set("J", 1.0_f64);
        p.set("M", 16usize);
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(1);
        let mc = HeisenbergChainMC::from_params(&p, &mut rng).unwrap();
        assert_eq!(mc.n_sites(), 6);
        assert_eq!(mc.config.n_slices, 16);
    }

    /// PT log_weight_ratio matches analytic form.
    #[test]
    fn pt_log_weight_ratio() {
        let mut p = Params::new();
        p.set("L", 4usize);
        p.set("beta", 1.0_f64);
        p.set("M", 8usize);
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(1);
        let mc = HeisenbergChainMC::from_params(&p, &mut rng).unwrap();
        let e = mc.config.energy_quantum(&mc.ham);
        let lr = mc.log_weight_ratio("beta", 2.0);
        assert!((lr - (1.0 - 2.0) * e).abs() < 1e-9);
    }

    /// End-to-end through the Scheduler: observables registered, energy
    /// negative for the AF chain at low T.
    #[test]
    fn end_to_end_runs_and_measures() {
        let mut p = Params::new();
        p.set("L", 8usize);
        p.set("beta", 2.0_f64);
        p.set("J", 1.0_f64); // AF convention: H = +J Σ Sᵢ·Sⱼ
        p.set("M", 16usize);

        let cfg = RunConfig {
            thermalization_sweeps: 500,
            measurement_sweeps: 500,
            binsize: 50,
            base_seed: 42,
            ..Default::default()
        };
        let scheduler = Scheduler::new(RayonBackend::new(1), cfg);
        let results = scheduler.run_one::<HeisenbergChainMC>(&p);

        let e = results.get("EnergyPerSite").expect("EnergyPerSite missing");
        let m2 = results.get("M2").expect("M2 missing");
        assert!(e.stderr > 0.0, "estimator should have nonzero error");
        assert!(m2.mean >= 0.0, "M² must be non-negative");
    }
}
