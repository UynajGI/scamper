# scamper

Scamper is a Rust Monte Carlo simulation stack. This facade crate re-exports
the registry-published parts of the workspace so applications can depend on a
single package while retaining access to each layer.

## Status

`0.1.0-dev2` is a prerelease. It has no API or behavioral stability guarantee.

## Install

```bash
cargo add scamper@0.1.0-dev2
```

The facade re-exports these crates:

- `scamper::carlo_rs`: scheduling, deterministic random-number setup,
  measurement, error analysis, checkpointing, Rayon, and MPI backends.
- `scamper::cmc_rs`: classical lattice, particle, generalized-ensemble, worm,
  and dynamics kernels.
- `scamper::qmc_rs`: continuous-time lattice QMC, impurity QMC, and
  variational Monte Carlo kernels.

Each crate can also be added directly when an application needs only one
layer.

## Features

`hdf5` forwards HDF5 support to every re-exported crate. `mpi` similarly
forwards MPI support and requires a usable native MPI installation.

```bash
cargo add scamper@0.1.0-dev2 --features hdf5,mpi
```

`mcmc-rs` remains repository-only because that crate name is already occupied
on crates.io by an unrelated package.

## License

Apache-2.0
