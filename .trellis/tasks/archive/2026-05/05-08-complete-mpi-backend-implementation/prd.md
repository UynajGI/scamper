# Complete MPI Backend Implementation

## Goal
Bring MPI backend from skeleton (~70%) to fully functional, referencing Carlo.jl's `scheduler_mpi.jl`.

## Current Gaps (priority order)

### Critical
1. **Results transfer** — Controller's `is_complete()` handler is empty. Workers accumulate results locally but never send to controller. Need MPI-based results transfer.
2. **run_id tracking** — Controller hardcodes `run_id = 1`. Need per-task run_id counter.
3. **run_comm usage** — `run_comm` is created/split but unused in worker (`_run_comm`). Multi-rank-per-run coordination missing (Bcast timeup, leader-only controller comm).
4. **run_distributed_compat** — Creates dummy TaskSpec (empty Params). Must use actual task params from MpiRunConfig.

### Important
5. **Worker results send** — Worker pushes results to local Vec but never transmits to controller. Must send after task completion.
6. **Controller aggregates received results** — ResultsAggregator exists but never fed real data.
7. **Timeup coordination** — Timeup only checked locally. Julia broadcasts timeup within run group via Bcast.

### Nice-to-have
8. Controller slow warning (from Julia)
9. Channel abstraction (MpiSender/MpiReceiver) per design spec
10. Typed message enums instead of integer magic numbers
