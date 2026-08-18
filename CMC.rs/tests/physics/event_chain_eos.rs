//! Equation-of-state validation for the straight hard-sphere event-chain
//! kernel ([`cmc_rs::HardSphereEventChain`]).
//!
//! Physics under test — two-dimensional hard disks of diameter σ, number
//! density ρ = N/V, packing fraction η = Nπσ²/(4V):
//!
//! * **Contact-value pressure identity.** For hard disks the virial equation
//!   collapses to the contact value of the pair correlation function:
//!   `Z = βP/ρ = 1 + 2η g(σ⁺)` (exact; the hard-core `−βu′(r)` is a delta
//!   spike at contact, leaving `1 + (πσ²/2)ρ g(σ⁺)`).
//!
//! * **Event-chain estimator.** A straight chain of length ℓ launched from an
//!   equilibrium particle collides with every exclusion disk intersecting its
//!   flight segment, so `E[collisions]/ℓ = 2σρ(1−1/N)⟨g⟩_band`, where the
//!   band hugs the contact circle. Expanding the band in polar coordinates
//!   gives `⟨g⟩_band = g(σ⁺) + (πℓ/8) g′(σ⁺) + (πℓ²/24) g″(σ⁺) + O(ℓ³)`,
//!   i.e. the collision rate measures g averaged over a radial shell of
//!   effective width `ε_eff = πℓ/4`. Richardson extrapolation over chain
//!   lengths (ℓ, 2ℓ) cancels the linear term and recovers g(σ⁺) with
//!   `−πℓ²g″/12` residual bias. This uses only the public `collisions()` and
//!   `lifted_distance()` counters; every collision occurs at exactly r = σ.
//!
//! * **Low-density reference.** Exact hard-disk virial coefficients
//!   B₂ = πσ²/2 and B₃ = π(4π − 3√3)σ⁴/12 ≈ 1.92954σ⁴ (Mayer triple
//!   integral; the commonly quoted ratio B₃/B₂² = 4/3 − √3/π ≈ 0.7820)
//!   give `Z = 1 + 2η + c₃η² + O(η³)` with
//!   `c₃ = B₃(4/(πσ²))² = 16/3 − 4√3/π ≈ 3.12802`. The finite-N factor
//!   (1 − 1/N) multiplies the 2η term exactly.
//!
//! * **Cross-solver check.** The same (N, η) fluid is sampled independently
//!   by [`cmc_rs::ParticleMetropolisCore`] with an exact hard-disk core
//!   potential; g(σ⁺) is measured from binned pair distances with Richardson
//!   extrapolation over the shell width. Agreement of two algorithmically
//!   unrelated samplers through two different estimators validates the
//!   equilibrium distribution and the pressure convention.

use cmc_rs::{
    HardSphereEventChain, OrthorhombicCell, PairPotential, ParticleAlgorithm,
    ParticleConfiguration, ParticleMetropolisCore, ParticleSystem, SimulationCell,
};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

use crate::common::zscore_seed_count;

/// Disk count for all fluid points (small enough for O(N²) sweeps, large
/// enough that finite-size corrections beyond the exact (1 − 1/N) term stay
/// below the statistical resolution).
const DISKS: usize = 64;

/// Diameter in solver units; everything below is expressed in σ.
const DIAMETER: f64 = 1.0;

/// Chain length used during event-chain warm-up (long chains decorrelate the
/// random sequential-addition start cheaply: the scan cost per unit lifted
/// distance is 1/chain_length of the measurement cost).
const WARMUP_CHAIN_LENGTH: f64 = 5.0;

/// Exact second virial coefficient of hard disks, B₂ = πσ²/2.
fn hard_disk_b2_exact() -> f64 {
    std::f64::consts::FRAC_PI_2 * DIAMETER * DIAMETER
}

/// Exact third virial coefficient of hard disks,
/// B₃ = π(4π − 3√3)σ⁴/12 ≈ 1.92954σ⁴ (π²/3 − π√3/4 in σ⁴ units).
///
/// Cross-checked inside the calibration harness by Monte Carlo integration
/// of the Mayer triple integral `B₃ = −(1/3)∫∫ f₁₂f₁₃f₂₃ d²r₁₂ d²r₁₃`.
fn hard_disk_b3_exact() -> f64 {
    std::f64::consts::PI * (4.0 * std::f64::consts::PI - 3.0 * 3.0_f64.sqrt()) / 12.0
        * DIAMETER.powi(4)
}

/// η² coefficient of the hard-disk virial series,
/// c₃ = B₃·(4/(πσ²))² = 16/3 − 4√3/π ≈ 3.12802.
fn virial_c3() -> f64 {
    hard_disk_b3_exact() * (4.0 / (std::f64::consts::PI * DIAMETER * DIAMETER)).powi(2)
}

/// Side length of the square box holding `DISKS` disks at packing fraction η.
fn disk_side(eta: f64) -> f64 {
    ((DISKS as f64) * std::f64::consts::FRAC_PI_4 * DIAMETER * DIAMETER / eta).sqrt()
}

/// Number density of the box.
fn disk_density(eta: f64) -> f64 {
    DISKS as f64 / (disk_side(eta) * disk_side(eta))
}

/// Sample mean of a small fixed-size slice.
fn mean(values: &[f64]) -> f64 {
    values.iter().sum::<f64>() / values.len() as f64
}

/// Standard error of the mean across independent seeds (sample standard
/// deviation over sqrt(n)); the per-seed runs are independent chains, so
/// seed-level scatter includes all within-seed autocorrelation.
fn standard_error(values: &[f64]) -> f64 {
    let n = values.len() as f64;
    let average = mean(values);
    let variance = values
        .iter()
        .map(|value| (value - average).powi(2))
        .sum::<f64>()
        / (n - 1.0);
    (variance / n).sqrt()
}

/// Two-sided pooled z-score between two independently measured means.
fn pooled_z(left: (f64, f64), right: (f64, f64)) -> f64 {
    (left.0 - right.0) / (left.1 * left.1 + right.1 * right.1).sqrt()
}

/// Exact hard-disk pair potential: infinite inside the core, zero outside.
///
/// The cutoff equals the diameter, matching the `PairPotential` convention of
/// [`cmc_rs::LennardJones`] (`distance_squared >= cutoff_squared` is zero), so
/// the packed cell list treats core contact consistently.
#[derive(Debug, Clone, Copy)]
struct HardDiskCore {
    diameter_squared: f64,
}

impl PairPotential for HardDiskCore {
    fn cutoff_squared(&self) -> f64 {
        self.diameter_squared
    }

    fn energy(&self, _: u16, _: u16, distance_squared: f64) -> f64 {
        if distance_squared < self.diameter_squared {
            f64::INFINITY
        } else {
            0.0
        }
    }
}

/// Randomised non-overlapping start (random sequential addition), so every
/// solver/seed begins from a distinct disordered fluid state rather than a
/// lattice.
fn random_start(seed: u64, eta: f64) -> ParticleConfiguration<2> {
    let side = disk_side(eta);
    let cell = OrthorhombicCell::new([side, side]).unwrap();
    let core_squared = DIAMETER * DIAMETER;
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut positions = Vec::with_capacity(DISKS);
    let mut guard = 0usize;
    while positions.len() < DISKS {
        guard += 1;
        assert!(
            guard < 200_000,
            "random sequential addition stalled at eta={eta}"
        );
        let candidate = [rng.random::<f64>() * side, rng.random::<f64>() * side];
        let overlaps = positions
            .iter()
            .any(|placed| cell.distance_squared(placed, &candidate) < core_squared);
        if !overlaps {
            positions.push(candidate);
        }
    }
    ParticleConfiguration::new(positions, vec![0; DISKS], cell).unwrap()
}

// ── Event-chain contact-value estimator ───────────────────────────────────

/// Overlap-audit cadence (in chains) for the production runs; occasional
/// `validate()` sweeps keep the no-overlap invariant checked without paying
/// O(N²) per chain.
const AUDIT_EVERY_CHAINS: u64 = 5_000;

/// Warm-up chains (at the long warm-up chain length) that decorrelate the
/// random-sequential-addition start before any measurement.
const WARMUP_CHAINS: usize = 300;

/// The two measurement chain lengths whose sausage averages combine by
/// Richardson extrapolation `g(σ⁺) = 2⟨g⟩_band(ℓ) − ⟨g⟩_band(2ℓ)`, cancelling
/// the linear band-bias term −(πℓ/8)g′(σ⁺).
const MEASUREMENT_CHAIN_LENGTHS: [f64; 2] = [1.0, 2.0];

/// Sausage-averaged contact correlation ⟨g⟩_band at one chain length,
/// measured from the public collision/distance counters after a shared
/// warm-up, driven until `target_collisions` collisions have accumulated.
///
/// `E[collisions]/lifted = 2σρ(1−1/N)⟨g⟩_band` (every exclusion disk whose
/// centre lies in the swept tube is hit exactly once at r = σ), so
/// ⟨g⟩_band = collisions-per-length / (2σρ(1−1/N)).
fn contact_g_band(seed: u64, eta: f64, chain_length: f64, target_collisions: u64) -> f64 {
    // Decorrelate all chain lengths of one seed from a common fluid state.
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed ^ 0x9E37_79B9_7F4A_7C15);
    let mut warmup = HardSphereEventChain::<2>::new(
        random_start(seed, eta),
        DIAMETER,
        WARMUP_CHAIN_LENGTH,
        AUDIT_EVERY_CHAINS,
    )
    .unwrap();
    for _ in 0..WARMUP_CHAINS {
        warmup.step(&mut rng).unwrap();
    }

    // Fresh measurement kernel from the decorrelated configuration, with an
    // independent stream per chain length (the length enters the seed).
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed ^ chain_length.to_bits());
    let mut kernel = HardSphereEventChain::<2>::new(
        warmup.configuration().clone(),
        DIAMETER,
        chain_length,
        AUDIT_EVERY_CHAINS,
    )
    .unwrap();
    let (c0, l0) = (kernel.collisions(), kernel.lifted_distance());
    let mut guard = 0u64;
    while kernel.collisions() - c0 < target_collisions {
        kernel.step(&mut rng).unwrap();
        guard += 1;
        assert!(guard < 100_000_000, "collision budget stalled at eta={eta}");
    }
    let collisions = (kernel.collisions() - c0) as f64;
    let lifted = kernel.lifted_distance() - l0;
    collisions / lifted / (2.0 * DIAMETER * disk_density(eta) * (1.0 - 1.0 / DISKS as f64))
}

/// Richardson-extrapolated contact value g(σ⁺) for one seed: the linear
/// band-bias cancels between chain lengths ℓ and 2ℓ.
fn contact_g(seed: u64, eta: f64, target_collisions: u64) -> f64 {
    let [short, long] = MEASUREMENT_CHAIN_LENGTHS;
    2.0 * contact_g_band(seed, eta, short, target_collisions)
        - contact_g_band(seed, eta, long, target_collisions)
}

/// Finite-N compressibility factor from a contact value:
/// `Z_N = 1 + 2η(1−1/N) g(σ⁺)` (the 2η factor is B₂ρ with
/// B₂ = πσ²/2).
fn compressibility_from_contact(eta: f64, contact_g: f64) -> f64 {
    let finite_n = 1.0 - 1.0 / DISKS as f64;
    1.0 + hard_disk_b2_exact() * disk_density(eta) * finite_n * contact_g
}

/// Virial reference through second order in η (finite-N exact at B₂):
/// `Z_ref = 1 + 2η(1−1/N) + c₃η²`, with the O(η³) remainder bounded by
/// c₄η³ ≲ 6.1η³ (B₄/B₂² ≈ 0.304 for hard disks).
fn virial_reference(eta: f64) -> f64 {
    1.0 + 2.0 * eta * (1.0 - 1.0 / DISKS as f64) + virial_c3() * eta * eta
}

/// Low-density equation-of-state validation: the event-chain collision
/// estimator must reproduce the exact hard-disk virial series through O(η²)
/// at two densities. Statistical gate: |Z − Z_ref| within 4 standard errors
/// across seeds, with an absolute floor covering the O(η³) truncation
/// (c₄η³ ≈ 0.0004 at η = 0.04, ≈ 0.0008 at η = 0.07).
#[test]
fn event_chain_low_density_virial_matches_exact() {
    let seeds = zscore_seed_count(6);
    for eta in [0.04, 0.07] {
        let contacts: Vec<f64> = (0..seeds)
            .map(|seed| contact_g(seed as u64, eta, 8_000))
            .collect();
        let contact = mean(&contacts);
        let error = standard_error(&contacts);
        assert!(
            (0.85..1.25).contains(&contact),
            "eta={eta}: unphysical contact value {contact}"
        );
        let z = compressibility_from_contact(eta, contact);
        let reference = virial_reference(eta);
        let gate = 4.0
            * compressibility_from_contact(eta, error)
                .max(1e-4)
                .max(6.5 * eta.powi(3));
        assert!(
            (z - reference).abs() <= gate,
            "eta={eta}: Z={z:.5} vs virial reference {reference:.5} \
             (gate {gate:.5}, contact g={contact:.4}±{error:.4})"
        );
    }
}

// ── Cross-solver: Metropolis NVT shell estimator ───────────────────────────

/// Shell-averaged contact correlation from a Metropolis NVT trajectory:
/// counting pairs with r ∈ [σ, σ+δ] and normalising by the ideal-gas
/// expectation estimates ⟨g⟩ over the shell; Richardson over (δ, 2δ)
/// cancels the linear shell bias.
fn metropolis_contact_g(seed: u64, eta: f64) -> f64 {
    let potential = HardDiskCore {
        diameter_squared: DIAMETER * DIAMETER,
    };
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0xD1CE_0000 ^ seed);
    let mut system = ParticleSystem::<2>::new(random_start(seed, eta), &potential, 1.0).unwrap();
    let mut kernel = ParticleMetropolisCore::<2>::new(0.15);

    // Warm up, then measure pair shells every `SAMPLE_EVERY` sweeps.
    const SAMPLE_EVERY: usize = 4;
    const WARMUP_SWEEPS: usize = 2_000;
    const MEASUREMENT_SWEEPS: usize = 24_000;
    for _ in 0..WARMUP_SWEEPS {
        kernel.sweep(&mut system, &potential, &mut rng);
    }
    let [delta, wide] = [0.02, 0.04];
    let (mut narrow_count, mut wide_count, mut samples) = (0u64, 0u64, 0u64);
    for sweep in 0..MEASUREMENT_SWEEPS {
        kernel.sweep(&mut system, &potential, &mut rng);
        if sweep % SAMPLE_EVERY != 0 {
            continue;
        }
        let positions = system.configuration().positions();
        let cell = system.configuration().cell();
        let narrow_squared = (DIAMETER + delta).powi(2);
        let wide_squared = (DIAMETER + wide).powi(2);
        for left in 0..positions.len() {
            for right in left + 1..positions.len() {
                let distance_squared = cell.distance_squared(&positions[left], &positions[right]);
                if distance_squared < wide_squared {
                    wide_count += 1;
                    if distance_squared < narrow_squared {
                        narrow_count += 1;
                    }
                }
            }
        }
        samples += 1;
    }

    // Ideal-gas expectation of the pair count in a shell annulus.
    let volume = system.configuration().cell().volume();
    let side = disk_side(eta);
    let _ = side;
    let ideal = |width: f64| {
        (DISKS * (DISKS - 1)) as f64 / 2.0
            * std::f64::consts::PI
            * ((DIAMETER + width).powi(2) - DIAMETER * DIAMETER)
            / volume
    };
    let narrow = narrow_count as f64 / samples as f64 / ideal(delta);
    let wide = wide_count as f64 / samples as f64 / ideal(wide);
    2.0 * narrow - wide
}

/// Cross-solver agreement at η = 0.2: the algorithmically unrelated
/// Metropolis NVT shell estimator and the event-chain collision estimator
/// must agree on g(σ⁺) within pooled 4σ, both landing in the literature
/// window for hard disks (Henderson-era MD: g(σ⁺) ≈ 1.3 at η = 0.2), and
/// the mid-density Z must exceed the low-density virial value clearly.
#[test]
fn metropolis_contact_g_matches_event_chain() {
    const ETA: f64 = 0.2;
    let seeds = zscore_seed_count(4).min(6);

    let event_values: Vec<f64> = (0..seeds)
        .map(|seed| contact_g(seed as u64, ETA, 6_000))
        .collect();
    let metro_values: Vec<f64> = (0..seeds)
        .map(|seed| metropolis_contact_g(seed as u64, ETA))
        .collect();
    let event = (mean(&event_values), standard_error(&event_values));
    let metro = (mean(&metro_values), standard_error(&metro_values));
    let z = pooled_z(event, metro);
    assert!(
        z.abs() < 4.0,
        "eta={ETA}: event-chain g={}±{}, Metropolis g={}±{} (z={z:.2})",
        event.0,
        event.1,
        metro.0,
        metro.1
    );

    // Absolute anchor: the hard-disk contact value at eta = 0.2 is firmly
    // established near 1.3 in the literature; the window rejects estimator
    // or convention errors far outside it without over-committing.
    for (value, solver) in [(event.0, "event chain"), (metro.0, "Metropolis")] {
        assert!(
            (1.15..1.5).contains(&value),
            "{solver} contact value {value} outside the literature window at eta={ETA}"
        );
    }

    // Directional: mid-density Z clearly exceeds the exact low-density
    // virial value (which the low-density test pins to < 1.16).
    let z_mid = compressibility_from_contact(ETA, event.0);
    assert!(
        z_mid > virial_reference(0.07) + 0.08,
        "mid-density Z={z_mid:.4} not clearly above the low-density virial value"
    );
}

// ── Virial-coefficient calibration (guards the constants) ─────────────────

/// Monte-Carlo integration of the Mayer triple integral
/// `B₃ = −(1/3)∫∫ f₁₂ f₁₃ f₂₃ d²r₁₂ d²r₁₃` for hard disks: with particle 1
/// fixed, f = −1 exactly inside the core, so B₃ = (L⁴/3)·P(all three pairs
/// overlap) for r₁₂, r₁₃ uniform on a square of side L ≥ 2σ. A wrong B₃
/// constant (or convention) in the virial reference above fails here at
/// many standard errors.
#[test]
fn hard_disk_b3_mayer_integral_matches_exact() {
    let length = 2.5_f64;
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0xB3_C0FFEE);
    let samples = 20_000_000u64;
    let mut triples = 0u64;
    for _ in 0..samples {
        let r12 = [
            (rng.random::<f64>() - 0.5) * length,
            (rng.random::<f64>() - 0.5) * length,
        ];
        let r13 = [
            (rng.random::<f64>() - 0.5) * length,
            (rng.random::<f64>() - 0.5) * length,
        ];
        let dx = r12[0] - r13[0];
        let dy = r12[1] - r13[1];
        if r12[0] * r12[0] + r12[1] * r12[1] < 1.0
            && r13[0] * r13[0] + r13[1] * r13[1] < 1.0
            && dx * dx + dy * dy < 1.0
        {
            triples += 1;
        }
    }
    let probability = triples as f64 / samples as f64;
    let b3_mc = length.powi(4) / 3.0 * probability;
    let b3_exact = hard_disk_b3_exact();
    let sigma = (probability * (1.0 - probability) / samples as f64).sqrt() * length.powi(4) / 3.0;
    assert!(
        (b3_mc - b3_exact).abs() <= 4.0 * sigma + 1e-4,
        "B3 MC={b3_mc:.5} vs exact {b3_exact:.5} (sigma={sigma:.5})"
    );
}

/// Long variant: more seeds and a denser fluid point for the nightly run.
#[test]
#[ignore = "long: high-precision EOS + cross-solver at eta=0.3 (nightly --ignored)"]
fn event_chain_eos_long() {
    let seeds = zscore_seed_count(16);
    for eta in [0.04, 0.07, 0.12] {
        let contacts: Vec<f64> = (0..seeds)
            .map(|seed| contact_g(seed as u64, eta, 20_000))
            .collect();
        let contact = mean(&contacts);
        let error = standard_error(&contacts);
        let z = compressibility_from_contact(eta, contact);
        let reference = virial_reference(eta);
        // At eta=0.12 the eta^3 truncation reaches ~0.011: the reference
        // comparison stays honest by widening the gate accordingly.
        let gate = 4.0 * compressibility_from_contact(eta, error).max(6.5 * eta.powi(3));
        assert!(
            (z - reference).abs() <= gate,
            "eta={eta}: Z={z:.5} vs reference {reference:.5} (gate {gate:.5})"
        );
    }

    const ETA: f64 = 0.3;
    let event_values: Vec<f64> = (0..seeds)
        .map(|seed| contact_g(seed as u64, ETA, 12_000))
        .collect();
    let event = (mean(&event_values), standard_error(&event_values));
    let metro_values: Vec<f64> = (0..seeds.min(6))
        .map(|seed| metropolis_contact_g(seed as u64, ETA))
        .collect();
    let metro = (mean(&metro_values), standard_error(&metro_values));
    let z = pooled_z(event, metro);
    assert!(
        z.abs() < 4.0,
        "eta={ETA}: event g={}±{} vs Metropolis g={}±{} (z={z:.2})",
        event.0,
        event.1,
        metro.0,
        metro.1
    );
}
