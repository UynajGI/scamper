use mcmc_rs::{
    run_parallel_tempering, LogDensity, RandomWalkMetropolis, TemperingConfig, TraceStore,
};

struct TemperedGaussian {
    beta: f64,
}

impl LogDensity<[f64]> for TemperedGaussian {
    fn log_density(&mut self, state: &[f64]) -> f64 {
        -0.5 * self.beta * state[0] * state[0]
    }
}

fn run_once() -> mcmc_rs::TemperingOutput {
    let config = TemperingConfig {
        ladder: vec![1.0, 0.4, 0.15],
        warmup: 24,
        samples: 40,
        thinning: 2,
        exchange_interval: 4,
        base_seed: 1234,
        parameter_names: vec!["x".to_string()],
    };
    run_parallel_tempering(
        |_slot, beta| TemperedGaussian { beta },
        |_slot, _beta| RandomWalkMetropolis::isotropic(1, 0.7).unwrap(),
        vec![vec![-4.0], vec![0.0], vec![4.0]],
        config,
    )
    .unwrap()
}

#[test]
fn replica_exchange_records_fixed_slot_traces_and_edge_statistics() {
    let output = run_once();
    assert_eq!(output.parameter_names, vec!["x".to_string()]);
    assert_eq!(output.replicas.len(), 3);
    assert_eq!(output.exchanges.len(), 2);
    for replica in &output.replicas {
        assert_eq!(replica.trace.len(), 20);
        assert_eq!(replica.trace.dimension(), 1);
        assert!(replica.final_log_density.is_finite());
    }
    assert!(output.exchanges.iter().all(|edge| edge.attempts > 0));
    assert!(output
        .exchanges
        .iter()
        .all(|edge| edge.acceptances <= edge.attempts));
}

#[test]
fn replica_exchange_is_reproducible_for_a_fixed_seed() {
    assert_eq!(run_once(), run_once());
}
