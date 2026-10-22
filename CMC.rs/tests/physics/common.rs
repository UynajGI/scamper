use cmc_rs::{Bond, BondType, CsrLattice, Hamiltonian, IsingModel};

/// Environment variable overriding the seed count of the multi-seed
/// z-score tests (nightly z-score monitoring).
pub const ZSCORE_SEEDS_ENV: &str = "SCUTTLE_ZSCORE_SEEDS";

/// Upper bound accepted for [`ZSCORE_SEEDS_ENV`] — guards against typos that
/// would silently turn a nightly run into a multi-day marathon.
const ZSCORE_SEEDS_MAX: usize = 4096;

/// Seed count for the multi-seed z-score tests.
///
/// Reads `SCUTTLE_ZSCORE_SEEDS` (see [`ZSCORE_SEEDS_ENV`]):
/// - unset or empty → `default` (the documented per-test seed count; this is
///   the regular-CI path and must stay byte-for-byte identical);
/// - an integer in `1..=4096` → that count (nightly high-power monitoring);
/// - anything else → panic with a clear message, because a silently degraded
///   monitoring run is worse than a loudly failed one.
pub fn zscore_seed_count(default: usize) -> usize {
    let raw = match std::env::var(ZSCORE_SEEDS_ENV) {
        Ok(raw) => raw,
        Err(std::env::VarError::NotPresent) => return default,
        Err(std::env::VarError::NotUnicode(raw)) => {
            panic!("{ZSCORE_SEEDS_ENV} is not valid Unicode: {raw:?}")
        }
    };
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return default;
    }
    match trimmed.parse::<usize>() {
        Ok(n) if (1..=ZSCORE_SEEDS_MAX).contains(&n) => {
            eprintln!("[zscore] {ZSCORE_SEEDS_ENV}={n} (default {default})");
            n
        }
        _ => panic!(
            "{ZSCORE_SEEDS_ENV}={raw:?}: expected an integer in 1..={ZSCORE_SEEDS_MAX} \
             (unset the variable to use the default of {default} seeds)"
        ),
    }
}

pub fn assert_close(left: f64, right: f64, tolerance: f64) {
    assert!(
        (left - right).abs() <= tolerance,
        "{left:.17e} != {right:.17e}; |Δ|={:.3e}, tolerance={tolerance:.3e}",
        (left - right).abs()
    );
}

pub fn enumerate_ising(n_sites: usize) -> Vec<Vec<f64>> {
    (0..1usize << n_sites)
        .map(|mask| {
            (0..n_sites)
                .map(|site| if mask & (1 << site) == 0 { -1.0 } else { 1.0 })
                .collect()
        })
        .collect()
}

pub fn direct_ising_energy(spins: &[f64], lattice: &CsrLattice, coupling: f64) -> f64 {
    lattice
        .edges
        .iter()
        .map(|edge| -coupling * edge.weight * spins[edge.source] * spins[edge.target])
        .sum()
}

pub fn exact_ising_moments(lattice: &CsrLattice, coupling: f64, beta: f64) -> (f64, f64, f64, f64) {
    let model = IsingModel::new(coupling);
    let mut z = 0.0;
    let mut e = 0.0;
    let mut e2 = 0.0;
    let mut m2 = 0.0;
    for spins in enumerate_ising(lattice.n_sites) {
        let energy = model.compute_total_energy(&spins, lattice, 1.0);
        let magnetization = spins.iter().sum::<f64>() / lattice.n_sites as f64;
        let weight = (-beta * energy).exp();
        z += weight;
        e += weight * energy;
        e2 += weight * energy * energy;
        m2 += weight * magnetization * magnetization;
    }
    (z, e / z, e2 / z, m2 / z)
}

/// Exact (⟨E⟩, ⟨m²⟩, ⟨cos(θ1−θ3)⟩) of the ferromagnetic XY 4-ring (PBC chain
/// of 4 sites, J = `coupling`) at inverse temperature `beta`.
///
/// The global-rotation zero mode is factored out exactly (θ₄ = 0), leaving a
/// smooth, strictly periodic 3D integrand; periodic trapezoidal quadrature
/// converges on it spectrally, so `grid = 96` per dimension is already far
/// beyond what β ≤ 1.5 features require (verified by grid-doubling in the
/// callers).
pub fn exact_xy_ring4_moments(beta: f64, coupling: f64, grid: usize) -> (f64, f64, f64) {
    let step = std::f64::consts::TAU / grid as f64;
    let mut z = 0.0;
    let mut energy_sum = 0.0;
    let mut m2_sum = 0.0;
    let mut c13_sum = 0.0;
    for i in 0..grid {
        let t1 = i as f64 * step;
        for j in 0..grid {
            let t2 = j as f64 * step;
            for k in 0..grid {
                let t3 = k as f64 * step;
                let energy = -coupling * ((t1 - t2).cos() + (t2 - t3).cos() + t3.cos() + t1.cos());
                let weight = (-beta * energy).exp();
                let mx = t1.cos() + t2.cos() + t3.cos() + 1.0;
                let my = t1.sin() + t2.sin() + t3.sin();
                let m2 = (mx * mx + my * my) / 16.0;
                z += weight;
                energy_sum += weight * energy;
                m2_sum += weight * m2;
                c13_sum += weight * (t1 - t3).cos();
            }
        }
    }
    (energy_sum / z, m2_sum / z, c13_sum / z)
}

// ── N2/N3 network graph builders ─────────────────────────────────────────

/// Star graph: center 0 with `leaves` leaves.
pub fn star_lattice(leaves: usize) -> CsrLattice {
    CsrLattice::from_edges(
        leaves + 1,
        (1..=leaves)
            .map(|leaf| Bond::new(0, leaf, BondType::Generic, 1.0))
            .collect(),
    )
}

/// Sparse G(n, m) Erdos-Renyi graph with distinct non-loop edges.
///
/// Membership is tracked in a HashSet: a Vec-based scan is O(E^2) and makes
/// 60k+-edge graphs take minutes to build in debug test builds.
pub fn erdos_renyi_graph(n: usize, edge_count: usize, seed: u64) -> CsrLattice {
    use rand::{RngExt, SeedableRng};
    use std::collections::HashSet;
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut seen: HashSet<(usize, usize)> = HashSet::with_capacity(edge_count);
    let mut edges = Vec::with_capacity(edge_count);
    while edges.len() < edge_count {
        let left = rng.random_range(0..n);
        let right = rng.random_range(0..n);
        if left != right {
            let pair = if left < right {
                (left, right)
            } else {
                (right, left)
            };
            if seen.insert(pair) {
                edges.push(Bond::new(pair.0, pair.1, BondType::Generic, 1.0));
            }
        }
    }
    CsrLattice::from_edges(n, edges)
}

/// Configuration-model graph with a truncated power-law degree sequence
/// `p(k) ~ k^(-alpha)`, `m <= k <= cutoff` (exact empirical distribution;
/// stub pairing, self-loops repaired by re-pairing, parallel edges allowed).
pub fn config_model_power_law(
    n: usize,
    alpha: f64,
    m: usize,
    cutoff: usize,
    seed: u64,
) -> CsrLattice {
    use rand::{RngExt, SeedableRng};
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(seed);
    // Inverse-CDF degrees on a discrete support m..=cutoff.
    let support: Vec<usize> = (m..=cutoff).collect();
    let mut weights = Vec::with_capacity(support.len());
    let mut total = 0.0_f64;
    for &k in &support {
        let w = (k as f64).powf(-alpha);
        weights.push(w);
        total += w;
    }
    let cumulative: Vec<f64> = weights
        .iter()
        .scan(0.0, |acc, &w| {
            *acc += w / total;
            Some(*acc)
        })
        .collect();
    let mut degrees = Vec::with_capacity(n);
    for _ in 0..n {
        let u: f64 = rng.random();
        let pick = cumulative
            .iter()
            .position(|&c| u <= c)
            .unwrap_or(support.len() - 1);
        degrees.push(support[pick]);
    }
    // Force an even stub count without changing the distribution shape
    // materially: bump one degree by one if needed.
    let stubs: usize = degrees.iter().sum();
    if !stubs.is_multiple_of(2) {
        let bump = degrees
            .iter_mut()
            .find(|d| **d < cutoff)
            .expect("bumpable degree");
        *bump += 1;
    }
    // Stub pairing. Broad degree sequences (kappa >> 2) produce ~kappa/2
    // self-loops per random pairing, so wholesale reshuffling would almost
    // never succeed; instead repair each self-loop locally by swapping
    // partners with a random disjoint pair (each repair strictly removes
    // one self-loop and cannot create a new one). Parallel edges stay
    // allowed, matching the configuration-model ensemble.
    let mut endpoints: Vec<usize> = Vec::with_capacity(stubs);
    for (vertex, &degree) in degrees.iter().enumerate() {
        for _ in 0..degree {
            endpoints.push(vertex);
        }
    }
    shuffle(&mut endpoints, &mut rng);
    let mut pairs: Vec<(usize, usize)> = (0..stubs / 2)
        .map(|pair_index| {
            let left = endpoints[2 * pair_index];
            let right = endpoints[2 * pair_index + 1];
            if left < right {
                (left, right)
            } else {
                (right, left)
            }
        })
        .collect();
    let pair_count = pairs.len();
    for index in 0..pair_count {
        while pairs[index].0 == pairs[index].1 {
            let vertex = pairs[index].0;
            // Find a disjoint partner pair (neither endpoint equals the
            // self-loop vertex, and it is not a self-loop itself).
            let mut other = rng.random_range(0..pair_count);
            let mut attempts = 0;
            while (other == index
                || pairs[other].0 == vertex
                || pairs[other].1 == vertex
                || pairs[other].0 == pairs[other].1)
                && attempts < 4 * pair_count
            {
                other = rng.random_range(0..pair_count);
                attempts += 1;
            }
            let (a, b) = pairs[other];
            pairs[index] = (vertex.min(a), vertex.max(a));
            pairs[other] = (vertex.min(b), vertex.max(b));
        }
    }
    let edges = pairs
        .into_iter()
        .map(|(left, right)| Bond::new(left, right, BondType::Generic, 1.0))
        .collect();
    CsrLattice::from_edges(n, edges)
}

fn shuffle(values: &mut [usize], rng: &mut rand_xoshiro::Xoshiro256PlusPlus) {
    use rand::RngExt;
    for position in (1..values.len()).rev() {
        let pick = rng.random_range(0..=position);
        values.swap(position, pick);
    }
}
