# Carlo.rs

Monte Carlo simulation framework in Rust with **100% core feature parity** with [Carlo.jl](https://github.com/lukas-weber/Carlo.jl).

See [Carlo.rs/README.md](Carlo.rs/README.md) for details.

## Features

- Autocorrelation analysis (variance ratio, decorrelated mode, covariance estimation)
- MPI distributed execution
- HDF5 checkpointing
- Complex observables
- Parallel tempering
- Performance monitoring

## Build

```bash
cd Carlo.rs
cargo build --release --features "hdf5 mpi"
```

## License

Apache-2.0