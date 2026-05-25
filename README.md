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
just check    # fmt + clippy + test (all crates)
just test     # test all crates
```

## License

Apache-2.0