# Complete Remaining 10% of MPI Backend

## 1. MonteCarlo Comm-Aware Sweep

Add `sweep_with_comm` default method to MonteCarlo trait, following Julia's `abstract_mc.jl`:

```julia
sweep!(mc::AbstractMC, ctx::MCContext, comm::MPI.Comm) = sweep!(mc, ctx)
```

In Rust: add an optional method with default fallback to `sweep(ctx)`. When ranks_per_run > 1 and the model overrides it, use communicator for inter-rank coordination during sweep.

## 2. Checkpoint MPI Coordination

When ranks_per_run > 1, checkpoint writes/reads need multi-rank coordination:
- Write: MPI.gather contexts from all ranks → rank 0 writes HDF5
- Read: rank 0 reads HDF5 → MPI.scatter/Bcast to all ranks
- Follow Julia's `run.jl` write_checkpoint/read_checkpoint patterns
