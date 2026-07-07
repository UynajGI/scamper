# Scuttle

Monte Carlo simulation framework in Rust, ported from [Carlo.jl](https://github.com/lukas-weber/Carlo.jl).

## Crates

| Crate | Purpose |
|-------|---------|
| [Carlo.rs](Carlo.rs/README.md) | Core framework — scheduler, measurements, MPI backend, checkpointing |
| QMC.rs | Quantum Monte Carlo toolbox — worldline objects |
| CMC.rs | Classical Monte Carlo toolbox — CSR lattice, Ising/Potts/XY/Heisenberg, Metropolis/Wolff/SW/Heat-bath algorithms, multi-spin coding, Parallel Tempering |

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
fmt + clippy + test + deny as parallel jobs with `--all-features`.

## Acknowledgments

This project is a Rust port of [Carlo.jl](https://github.com/lukas-weber/Carlo.jl)
by Lukas Weber (Max Planck Institute). The core architecture (scheduler,
measurements, MPI backend, checkpointing, parallel tempering) follows the
original design. See the accompanying paper: [arXiv:2408.03386](https://arxiv.org/abs/2408.03386).

## License

Apache-2.0