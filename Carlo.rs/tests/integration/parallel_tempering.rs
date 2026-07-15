use carlo_rs::parallel_tempering::ParallelTemperingConfig;

#[test]
fn test_ptmc_config() {
    let config = ParallelTemperingConfig {
        parameter: "T".to_string(),
        values: vec![0.1, 0.5, 1.0, 2.0],
        interval: 100,
    };
    assert_eq!(config.values.len(), 4);
}
