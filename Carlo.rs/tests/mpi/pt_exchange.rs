//! MPI parallel-tempering exchange protocol test.
//!
//! Run under mpirun:
//! ```bash
//! mpirun -np 2 cargo test --features mpi --test mpi_pt_exchange -- --ignored --nocapture
//! ```

#![cfg(feature = "mpi")]

use carlo_rs::{
    CarloError, Context, FromParams, MonteCarlo, ParallelTemperingCompatible,
    ParallelTemperingConfig, Params, Results,
};
use rand_xoshiro::Xoshiro256PlusPlus;

/// Simple model: measures a parameter-dependent observable.
/// The "energy" is just the parameter value itself, so exchange
/// acceptance is deterministic and easy to reason about.
struct PtTestMc {
    param_value: f64,
    sweep_count: u64,
}

impl MonteCarlo for PtTestMc {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, _ctx: &mut Context<Self::Rng>) {
        self.sweep_count += 1;
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        ctx.measure("ParamValue", self.param_value);
        ctx.measure("SweepCount", self.sweep_count as f64);
    }

    fn name(&self) -> &'static str {
        "PtTestMc"
    }
}

impl FromParams for PtTestMc {
    fn from_params(params: &Params, _rng: &mut Self::Rng) -> Result<Self, CarloError> {
        let param_value = params.get::<f64>("beta").unwrap_or(1.0);
        Ok(Self {
            param_value,
            sweep_count: 0,
        })
    }
}

impl ParallelTemperingCompatible for PtTestMc {
    fn log_weight_ratio(&self, _param: &str, new_value: f64) -> f64 {
        -(new_value - self.param_value) * self.param_value
    }

    fn change_parameter(&mut self, _param: &str, new_value: f64) {
        self.param_value = new_value;
    }
}

#[test]
#[ignore = "requires mpirun -np 2"]
fn pt_exchange_completes_and_returns_results() {
    let config = ParallelTemperingConfig {
        parameter: "beta".to_string(),
        values: vec![1.0, 2.0],
        interval: 5,
    };

    let mut params = Params::new();
    params.set("beta", "1.0");

    let result = carlo_rs::run_parallel_tempering::<PtTestMc, Xoshiro256PlusPlus>(
        &config, &params, 42, 2, 20, 5,
    );

    match result {
        Ok(Some(results)) => {
            // Controller rank gets aggregated results
            assert!(
                results.get("ParamValue").is_some(),
                "ParamValue should be in results"
            );
        }
        Ok(None) => {
            // Worker rank — no aggregate
        }
        Err(e) => {
            // Non-mpirun single process should give a clear error
            panic!("PT exchange failed: {e}");
        }
    }
}
