# Scuttle

Monte Carlo simulation framework in Rust, ported from [Carlo.jl](https://github.com/lukas-weber/Carlo.jl).

## Crates

| Crate | Purpose |
|-------|---------|
| [Carlo.rs](Carlo.rs/README.md) | Core framework — scheduler, measurements, MPI backend, checkpointing |
| QMC.rs | Quantum Monte Carlo — general continuous-time lattice QMC (arbitrary graph/spin) + quantum impurity (retarded-interaction) wormhole QMC |
| CMC.rs | Classical Monte Carlo — lattice (CSR graph, Ising/Potts/O(N), Metropolis/Wolff/SW/Heat-bath/Microcanonical, multi-spin coding, JSON snapshot v2) + particle (periodic cells, Lennard-Jones NVT/NPT/μVT, rigid molecules with dipolar external fields, cell lists, Metropolis-Hastings) + generalized ensembles (Wang-Landau, multicanonical, umbrella sampling) + classical worm (Ising HT graph, persistent physical/worm sectors, multi-component lattices) + classical dynamics (Kawasaki, Gillespie/BKL, hard-sphere event-chain, explicit event time) |
| [MCMC.rs](MCMC.rs/README.md) | Statistical MCMC — RW/slice/Gibbs/composed kernels, static HMC, NUTS, adaptive unit/diagonal/dense metrics, differentiable constrained transforms, replica exchange, traces and multi-chain diagnostics (E-BFMI) |

## Quick Start

```bash
cd Carlo.rs
cargo build --release --features "hdf5 mpi"
```

## Development

```bash
just hooks     # install git hooks (requires lefthook: brew install lefthook / go install ...)
just check     # fmt + clippy + test (all crates)
just deny      # cargo deny (advisories + licenses) — requires cargo-deny
just typos     # spelling check — requires typos
```

Git hooks (via [lefthook](https://github.com/evilmartians/lefthook)) enforce
`cargo fmt` + clippy + typos on commit, Conventional Commits on `commit-msg`,
and `cargo deny` on push. Tests run in CI, not the pre-push hook — use
`just test` locally before pushing. Lint level (including `unsafe_code = "deny"`)
is set in `[workspace.lints]` (`Cargo.toml`). Skip with `LEFTHOOK=0`. CI runs
fmt + clippy + test + deny as parallel jobs with `--all-features`, and a
[nightly workflow](.github/workflows/nightly.yml) runs the long (`--ignored`)
physics tests plus multi-seed z-score monitoring (`just nightly-zscore`
reproduces it locally), and a `carlo-framework` job exercises the Carlo.rs
HDF5 suite and every MPI test under `mpirun` at 1/2/4 ranks
(`just mpi-test [np]` reproduces it locally).

## Acknowledgments

This project is a Rust port of [Carlo.jl](https://github.com/lukas-weber/Carlo.jl)
by Lukas Weber (Max Planck Institute). The core architecture (scheduler,
measurements, MPI backend, checkpointing, parallel tempering) follows the
original design. See the accompanying paper: [arXiv:2408.03386](https://arxiv.org/abs/2408.03386).

## License

Apache-2.0