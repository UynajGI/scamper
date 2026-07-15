//! Carlo.rs adapter for the explicit occupation worldline solver.

use carlo_rs::{Context, MonteCarlo};

use crate::impurity::spin_boson::occupation::worldline::OccupationWorldlineSampler;

pub struct OccupationWorldlineQmc {
    sampler: OccupationWorldlineSampler,
}

impl OccupationWorldlineQmc {
    pub const fn new(sampler: OccupationWorldlineSampler) -> Self {
        Self { sampler }
    }
    pub const fn sampler(&self) -> &OccupationWorldlineSampler {
        &self.sampler
    }
    pub fn sampler_mut(&mut self) -> &mut OccupationWorldlineSampler {
        &mut self.sampler
    }
}

impl MonteCarlo for OccupationWorldlineQmc {
    type Rng = rand_xoshiro::Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        self.sampler
            .sweep(&mut ctx.rng)
            .unwrap_or_else(|error| panic!("occupation worldline sweep failed: {error}"));
        ctx.record_attempts(self.sampler.slices() as u64);
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        let obs = self
            .sampler
            .measure()
            .unwrap_or_else(|error| panic!("occupation measurement failed: {error}"));
        ctx.measure("OccupationEnergy", obs.energy);
        ctx.measure("OccupationSigmaZ", obs.sigma_z);
        ctx.measure("OccupationSigmaX", obs.sigma_x);
        ctx.measure("OccupationBosonNumber", obs.total_boson_number);
        ctx.measure("OccupationParity", obs.parity);
        ctx.measure("OccupationReducedSpinPurity", obs.reduced_spin_purity);
        ctx.measure_array("OccupationModeNumber", &obs.mode_occupations);
        ctx.measure_array("OccupationModeNumberSquared", &obs.mode_number_squared);
        ctx.measure_array("OccupationModeFactorialMoment", &obs.mode_factorial_moments);
        ctx.measure_array("OccupationModeG2Zero", &obs.mode_g2_zero);
        ctx.measure_array(
            "OccupationModeNumberCorrelations",
            &obs.mode_cross_correlations,
        );
        ctx.measure_array(
            "OccupationSpinBosonCovarianceZn",
            &obs.spin_boson_covariance_z_n,
        );
        ctx.measure(
            "OccupationWorldlineChangeFraction",
            self.sampler.acceptance_fraction(),
        );
    }

    fn name(&self) -> &'static str {
        "OccupationWorldlineQmc"
    }
}
