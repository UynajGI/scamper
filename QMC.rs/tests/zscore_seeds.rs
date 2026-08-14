//! Shared seed management for the multi-seed z-score tests.
//!
//! Setting `SCUTTLE_ZSCORE_SEEDS=<n>` raises the seed count of every z-score
//! test for nightly high-power monitoring. When the variable is unset (the
//! regular-CI path) the per-test default seed list is returned unchanged, so
//! default runs are byte-for-byte identical to before.

/// Environment variable overriding the seed count of the multi-seed z-score
/// tests (nightly z-score monitoring).
pub const ZSCORE_SEEDS_ENV: &str = "SCUTTLE_ZSCORE_SEEDS";

/// Upper bound accepted for [`ZSCORE_SEEDS_ENV`] — guards against typos that
/// would silently turn a nightly run into a multi-day marathon.
const ZSCORE_SEEDS_MAX: usize = 4096;

/// Seed count for the multi-seed z-score tests.
///
/// Reads `SCUTTLE_ZSCORE_SEEDS` (see [`ZSCORE_SEEDS_ENV`]):
/// - unset or empty → `default` (the documented per-test seed count);
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

/// Seed list for a multi-seed z-score test.
///
/// Returns `default_seeds` unchanged when the seed count is not overridden.
/// When more seeds are requested, the defaults come first (preserving the
/// default-run statistics) and the list is extended with deterministic
/// splitmix64-mixed values so every additional seed is independent.
pub fn zscore_seeds(default_seeds: &[u64]) -> Vec<u64> {
    let n = zscore_seed_count(default_seeds.len());
    let mut seeds = default_seeds.to_vec();
    while seeds.len() < n {
        let k = seeds.len() as u64;
        let mut z = k.wrapping_add(0x9E37_79B9_7F4A_7C15);
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        seeds.push((z ^ (z >> 31)) ^ default_seeds[0]);
    }
    seeds.truncate(n);
    seeds
}
