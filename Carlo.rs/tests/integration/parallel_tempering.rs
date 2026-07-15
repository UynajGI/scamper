use carlo_rs::parallel_tempering::{ParallelTemperingConfig, ParallelTemperingMC};
use carlo_rs::{Context, MonteCarlo, ParallelTemperingCompatible};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

#[test]
fn test_ptmc_config() {
    let config = ParallelTemperingConfig {
        parameter: "T".to_string(),
        values: vec![0.1, 0.5, 1.0, 2.0],
        interval: 100,
    };
    assert_eq!(config.values.len(), 4);
}

// ── PT-compatible test model ──────────────────────────────────────────────

struct PtTestMC {
    beta: f64,
    energy: f64,
}

impl MonteCarlo for PtTestMC {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        ctx.measure("Energy", self.energy);
    }

    fn name(&self) -> &'static str {
        "PtTestMC"
    }
}

impl ParallelTemperingCompatible for PtTestMC {
    fn log_weight_ratio(&self, param: &str, new_value: f64) -> f64 {
        if param == "beta" {
            -(new_value - self.beta) * self.energy
        } else {
            0.0
        }
    }

    fn change_parameter(&mut self, param: &str, new_value: f64) {
        if param == "beta" {
            self.beta = new_value;
        }
    }
}

#[test]
fn test_ptmc_wrapper_creation() {
    let config = ParallelTemperingConfig {
        parameter: "beta".to_string(),
        values: vec![0.5, 1.0, 2.0],
        interval: 10,
    };
    let mc = PtTestMC {
        beta: 0.5,
        energy: -1.0,
    };
    let wrapper = ParallelTemperingMC::new(&config, 0, mc);

    assert_eq!(wrapper.chain_idx(), 0);
    assert!((wrapper.current_value() - 0.5).abs() < 1e-10);
    assert_eq!(wrapper.parameter_name, "beta");
    assert_eq!(wrapper.parameter_values.len(), 3);
    assert_eq!(wrapper.tempering_interval, 10);
}

#[test]
fn test_ptmc_wrapper_current_value() {
    let config = ParallelTemperingConfig {
        parameter: "beta".to_string(),
        values: vec![0.5, 1.0, 2.0],
        interval: 10,
    };
    let mc = PtTestMC {
        beta: 1.0,
        energy: -1.0,
    };
    let wrapper = ParallelTemperingMC::new(&config, 1, mc);
    assert!((wrapper.current_value() - 1.0).abs() < 1e-10);
}

#[test]
fn test_ptmc_set_chain_idx() {
    let config = ParallelTemperingConfig {
        parameter: "beta".to_string(),
        values: vec![0.5, 1.0, 2.0],
        interval: 10,
    };
    let mc = PtTestMC {
        beta: 0.5,
        energy: -1.0,
    };
    let mut wrapper = ParallelTemperingMC::new(&config, 0, mc);

    wrapper.set_chain_idx(2);
    assert_eq!(wrapper.chain_idx(), 2);
    assert!((wrapper.current_value() - 2.0).abs() < 1e-10);
    assert!((wrapper.child_mc.beta - 2.0).abs() < 1e-10);
}

#[test]
fn test_ptmc_pt_measure_and_finalize() {
    let config = ParallelTemperingConfig {
        parameter: "beta".to_string(),
        values: vec![0.5, 1.0],
        interval: 10,
    };
    let mc = PtTestMC {
        beta: 0.5,
        energy: -1.0,
    };
    let mut wrapper = ParallelTemperingMC::new(&config, 0, mc);

    wrapper.pt_measure("CustomObs", 42.0);
    wrapper.pt_measure("CustomObs", 43.0);

    let rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(1);
    let mut ctx = Context::new_with_binsize(rng, 0, 10);
    wrapper.finalize_pt_measurements(&mut ctx);

    let estimates = ctx.finalize_measurements();
    assert!(estimates.contains_key("CustomObs"));
    let est = &estimates["CustomObs"];
    assert!((est.mean - 42.5).abs() < 1e-10);
}

#[test]
fn test_ptmc_log_weight_ratio() {
    let mc = PtTestMC {
        beta: 1.0,
        energy: -2.0,
    };
    // log_weight_ratio = -(new_beta - current_beta) * energy
    // For new_beta=2.0, energy=-2.0: -(2.0 - 1.0) * (-2.0) = 2.0
    let r = mc.log_weight_ratio("beta", 2.0);
    assert!((r - 2.0).abs() < 1e-10);
}

#[test]
fn test_ptmc_log_weight_ratio_unknown_param() {
    let mc = PtTestMC {
        beta: 1.0,
        energy: -2.0,
    };
    let r = mc.log_weight_ratio("temperature", 999.0);
    assert!((r - 0.0).abs() < 1e-10);
}

#[test]
fn test_ptmc_config_clone_debug() {
    let config = ParallelTemperingConfig {
        parameter: "T".to_string(),
        values: vec![0.1, 0.5],
        interval: 50,
    };
    let cloned = config.clone();
    assert_eq!(cloned.values, config.values);
    assert_eq!(cloned.interval, config.interval);

    let debug = format!("{:?}", config);
    assert!(debug.contains("ParallelTemperingConfig"));
}
