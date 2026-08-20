//! Carlo.rs `MonteCarlo` integration for [`VmcKernel`].
//!
//! Integration route: a direct `impl MonteCarlo for VmcKernel<W>` driven via
//! [`Run::from_parts`](carlo_rs::Run::from_parts), not `FromParams`. A
//! generic `W: WaveFunctionParams` cannot be constructed from a `Params`
//! map without an ansatz registry, and Carlo.rs provides `from_parts`
//! precisely for "models with closures, shared datasets or other
//! non-`Params` construction paths". The worm/lattice solvers keep
//! `FromParams` because their models are string-addressable; the
//! variational family is a typed API.
//!
//! Observable contract (weighted quantities follow the DESIGN.md
//! convention — accumulate `(x·w, w)` pairs as separate observables, form
//! ratios in postprocessing; L0 VMC weights are uniform, so plain means
//! suffice):
//!
//! - `LocalEnergy`, `LocalEnergySquared`, `EnergyPerParticle`,
//!   `LogGradSquared` — one sample per walker per measurement sweep;
//! - acceptance statistics flow through the context attempt clocks.

use carlo_rs::{Context, MonteCarlo, RngPhase, RunPhase};

use super::super::kernel::VmcKernel;
use super::super::wavefunction::{Positions, WaveFunctionParams};

/// Carlo.rs-visible name of the variational Metropolis kernel.
pub const VARIATIONAL_METROPOLIS_NAME: &str = "VariationalMetropolis";

impl<W: WaveFunctionParams<Config = Positions>> MonteCarlo for VmcKernel<W> {
    type Rng = rand_xoshiro::Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        let phase = rng_phase(ctx.phase());
        self.sweep_with_phase(&mut ctx.rng, phase);
        let (attempts, accepts) = self.drain_move_counters();
        ctx.record_attempts(attempts);
        ctx.record_accepted_moves(accepts);
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        let n_particles = self.n_particles() as f64;
        self.measure_population(|sample| {
            ctx.measure("LocalEnergy", sample.value);
            ctx.measure("LocalEnergySquared", sample.value * sample.value);
            ctx.measure("EnergyPerParticle", sample.value / n_particles);
            ctx.measure("LogGradSquared", sample.log_grad_squared);
        });
    }

    fn name(&self) -> &'static str {
        VARIATIONAL_METROPOLIS_NAME
    }
}

/// Map the scheduler's run phase onto the RNG-stream domain phase.
fn rng_phase(phase: RunPhase) -> RngPhase {
    match phase {
        RunPhase::Initialization => RngPhase::Initialization,
        RunPhase::Thermalization => RngPhase::Thermalization,
        RunPhase::Measurement => RngPhase::Measurement,
        RunPhase::Finished => RngPhase::Finished,
    }
}
