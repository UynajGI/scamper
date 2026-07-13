# Cache audit mode

CMC.rs treats accepted-state caches as part of the sampling contract. Auditing is enabled automatically in debug and test builds, and can be enabled in optimized builds with:

```bash
cargo test -p cmc-rs --features cache-audit
cargo run --release -p your-driver --features cmc-rs/cache-audit
```

The default automatic cadence is one audit every 1024 completed sweeps. A kernel-specific non-zero `energy_check_interval` overrides that cadence. Release builds without the feature and without an explicit interval perform no automatic audit.

Audits are detection-only. They never recompute and overwrite a bad cache before comparison.

- Lattice kernels validate graph/configuration shape and compare cached physical energy with a full Hamiltonian evaluation.
- Particle kernels validate finite coordinates, supported species, cached energy, packed cell membership, particle-to-cell and particle-to-slot reverse indices.
- Generalized-ensemble kernels additionally validate histogram/DOS dimensions and the cached macrostate bin against the accepted physical energy.
- Classical worm kernels recompute graph occupations, vertex parity, reduced log weight and endpoint-sector constraints; optional endpoint histogram dimensions are also checked.

The public helpers in `cmc_rs::audit` can also be called by custom kernels:

```rust,ignore
use cmc_rs::{audit_lattice_cache, should_audit_cache};

if should_audit_cache(completed_sweeps, configured_interval) {
    audit_lattice_cache(&system, &model)?;
}
```

A failed built-in periodic audit panics immediately because continuing a Markov chain with inconsistent accepted state would silently corrupt statistics. Direct helper calls return `Result` so applications can choose their own failure policy.
