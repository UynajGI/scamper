pub(crate) fn effective_sample_size(chains: &[Vec<f64>]) -> f64 {
    let chain_count = chains.len();
    let draws = chains.first().map_or(0, Vec::len);
    if chain_count < 2 || draws < 2 {
        return f64::NAN;
    }
    let means = chains
        .iter()
        .map(|chain| chain.iter().sum::<f64>() / draws as f64)
        .collect::<Vec<_>>();
    let variances = chains
        .iter()
        .zip(&means)
        .map(|(chain, mean)| autocovariance(chain, *mean, 0) * draws as f64 / (draws - 1) as f64)
        .collect::<Vec<_>>();
    let within = variances.iter().sum::<f64>() / chain_count as f64;
    let mean_of_means = means.iter().sum::<f64>() / chain_count as f64;
    let between = draws as f64
        * means
            .iter()
            .map(|mean| {
                let delta = mean - mean_of_means;
                delta * delta
            })
            .sum::<f64>()
        / (chain_count - 1) as f64;
    let variance_plus = ((draws - 1) as f64 / draws as f64) * within + between / draws as f64;
    if variance_plus <= 0.0 || !variance_plus.is_finite() {
        return (chain_count * draws) as f64;
    }

    let mut correlations = Vec::with_capacity(draws);
    correlations.push(1.0);
    for lag in 1..draws {
        let mean_autocovariance = chains
            .iter()
            .zip(&means)
            .map(|(chain, mean)| autocovariance(chain, *mean, lag))
            .sum::<f64>()
            / chain_count as f64;
        correlations.push(1.0 - (within - mean_autocovariance) / variance_plus);
    }

    let mut paired = Vec::new();
    let mut lag = 1;
    while lag + 1 < correlations.len() {
        let pair_sum = correlations[lag] + correlations[lag + 1];
        if pair_sum < 0.0 {
            break;
        }
        paired.push(pair_sum);
        lag += 2;
    }
    for index in 1..paired.len() {
        if paired[index] > paired[index - 1] {
            paired[index] = paired[index - 1];
        }
    }
    let tau = (-1.0 + 2.0 * (1.0 + paired.iter().sum::<f64>())).max(1.0);
    let total = (chain_count * draws) as f64;
    (total / tau).clamp(1.0, total)
}

fn autocovariance(values: &[f64], mean: f64, lag: usize) -> f64 {
    let count = values.len() - lag;
    values[..count]
        .iter()
        .zip(&values[lag..])
        .map(|(left, right)| (left - mean) * (right - mean))
        .sum::<f64>()
        / values.len() as f64
}
