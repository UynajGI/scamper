pub(crate) fn split_chains(chains: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let draws = chains.iter().map(Vec::len).min().unwrap_or(0);
    let half = draws / 2;
    let mut split = Vec::with_capacity(chains.len() * 2);
    for chain in chains {
        split.push(chain[..half].to_vec());
        split.push(chain[draws - half..draws].to_vec());
    }
    split
}

pub(crate) fn rank_normalize(chains: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let chain_count = chains.len();
    let draws = chains.first().map_or(0, Vec::len);
    let mut indexed = Vec::with_capacity(chain_count * draws);
    for (chain_index, chain) in chains.iter().enumerate() {
        for (draw_index, value) in chain.iter().copied().enumerate() {
            indexed.push((value, chain_index, draw_index));
        }
    }
    indexed.sort_by(|left, right| left.0.total_cmp(&right.0));
    let total = indexed.len();
    let mut normalized = vec![vec![0.0; draws]; chain_count];
    let mut start = 0;
    while start < total {
        let mut end = start + 1;
        while end < total && indexed[end].0 == indexed[start].0 {
            end += 1;
        }
        let average_rank = (start + 1 + end) as f64 / 2.0;
        let probability = (average_rank - 0.375) / (total as f64 + 0.25);
        let normal_score = inverse_standard_normal(probability);
        for &(_, chain_index, draw_index) in &indexed[start..end] {
            normalized[chain_index][draw_index] = normal_score;
        }
        start = end;
    }
    normalized
}

pub(crate) fn rhat(chains: &[Vec<f64>]) -> f64 {
    let chain_count = chains.len();
    let draws = chains.first().map_or(0, Vec::len);
    if chain_count < 2 || draws < 2 {
        return f64::NAN;
    }
    let means = chains
        .iter()
        .map(|chain| chain.iter().sum::<f64>() / draws as f64)
        .collect::<Vec<_>>();
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
    let within = chains
        .iter()
        .zip(&means)
        .map(|(chain, mean)| {
            chain
                .iter()
                .map(|value| {
                    let delta = value - mean;
                    delta * delta
                })
                .sum::<f64>()
                / (draws - 1) as f64
        })
        .sum::<f64>()
        / chain_count as f64;
    if within == 0.0 {
        return if between == 0.0 { 1.0 } else { f64::INFINITY };
    }
    let variance_plus = ((draws - 1) as f64 / draws as f64) * within + between / draws as f64;
    (variance_plus / within).sqrt()
}

pub(crate) fn quantile_sorted(sorted: &[f64], probability: f64) -> f64 {
    if sorted.len() == 1 {
        return sorted[0];
    }
    let location = probability * (sorted.len() - 1) as f64;
    let lower = location.floor() as usize;
    let upper = location.ceil() as usize;
    if lower == upper {
        sorted[lower]
    } else {
        let weight = location - lower as f64;
        sorted[lower] * (1.0 - weight) + sorted[upper] * weight
    }
}

// Peter J. Acklam's rational approximation.
fn inverse_standard_normal(probability: f64) -> f64 {
    const A: [f64; 6] = [
        -3.969_683_028_665_376e1,
        2.209_460_984_245_205e2,
        -2.759_285_104_469_687e2,
        1.383_577_518_672_69e2,
        -3.066_479_806_614_716e1,
        2.506_628_277_459_239,
    ];
    const B: [f64; 5] = [
        -5.447_609_879_822_406e1,
        1.615_858_368_580_409e2,
        -1.556_989_798_598_866e2,
        6.680_131_188_771_972e1,
        -1.328_068_155_288_572e1,
    ];
    const C: [f64; 6] = [
        -7.784_894_002_430_293e-3,
        -3.223_964_580_411_365e-1,
        -2.400_758_277_161_838,
        -2.549_732_539_343_734,
        4.374_664_141_464_968,
        2.938_163_982_698_783,
    ];
    const D: [f64; 4] = [
        7.784_695_709_041_462e-3,
        3.224_671_290_700_398e-1,
        2.445_134_137_142_996,
        3.754_408_661_907_416,
    ];
    const LOWER: f64 = 0.024_25;
    const UPPER: f64 = 1.0 - LOWER;

    if probability < LOWER {
        let q = (-2.0 * probability.ln()).sqrt();
        return (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0);
    }
    if probability <= UPPER {
        let q = probability - 0.5;
        let r = q * q;
        return (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q
            / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0);
    }
    let q = (-2.0 * (1.0 - probability).ln()).sqrt();
    -(((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
        / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
}
