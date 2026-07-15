# MPI implementation

## Process topology

`run_distributed` reserves world rank 0 for the controller. The remaining
ranks are divided into equal run groups of `DistributedConfig::ranks_per_run`.
Each group has:

- a group communicator used by `sweep_with_comm` and `measure_with_comm`;
- one leader (group rank 0) in a leader-only communicator with the controller;
- zero or more followers that receive commands by group broadcast.

The required world size is:

```text
1 + number_of_groups * ranks_per_run
```

## Scheduling protocol

The controller receives `Ready`, `Completed`, `Interrupted`, or `Failed`
reports from group leaders and responds with `Assign` or `Stop`. Sweep budgets
are reserved when assigned, so multiple groups cannot overschedule the same
task. Tasks are selected round-robin and split into chunks based on available
worker groups and bin size.

Only world rank 0 returns aggregated results. Worker ranks return an empty
vector after the final world barrier.

## Checkpoint and restart

The controller persists every assignment before it is sent:

```text
job_dir/task_XXXX/runXXXX/mpi-assignment.json
```

A completed run additionally has:

```text
job_dir/task_XXXX/runXXXX/result.json
```

With `hdf5`, each rank writes its own checkpoint file for a multi-rank run.
A two-phase generation protocol first writes `*.next.h5` staging files, then
publishes every rank file, and finally writes `mpi-checkpoint.json` as the
commit marker. A restart only accepts a generation with a matching commit
marker, so a batch-system kill cannot mix checkpoints from different sweeps.
Checkpoint writes use group-wide success consensus before any rank leaves the
collective path. Restart validation rejects changes to task parameters, seed,
bin size, target sweeps, thermalization, or ranks per run while persisted work
exists.

## Model contract

Single-rank models need no MPI-specific implementation. Multi-rank models may
override:

```rust,ignore
fn sweep_with_comm(
    &mut self,
    ctx: &mut Context<Self::Rng>,
    comm: &mpi::topology::SimpleCommunicator,
);

fn measure_with_comm(
    &mut self,
    ctx: &mut Context<Self::Rng>,
    comm: &mpi::topology::SimpleCommunicator,
);
```

All ranks in a run group must execute collectives in the same order. Global
observables should be reduced inside `measure_with_comm` and recorded on group
rank 0. Parallel tempering stores production observables under stable
`pt_chain_XXXX/<observable>` namespaces and gathers all rank-local `Results` on
world rank 0.

## Tests

Run MPI integration binaries separately because each process owns one MPI
initialization lifetime:

```bash
mpirun -np 4 cargo test --features mpi --test mpi_test -- --nocapture
mpirun -np 4 cargo test --features mpi --test mpi_distributed_test -- --nocapture
```

Checkpoint tests require both features:

```bash
cargo test --features "mpi hdf5"
```
