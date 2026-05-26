# Journal - Yuunagi (Part 1)

> AI development session journal
> Started: 2026-05-08

---



## Session 1: Carlo.rs 对齐 Carlo.jl — ResultTools + .githooks + 文档整理

**Date**: 2026-05-08
**Task**: Carlo.rs 对齐 Carlo.jl — ResultTools + .githooks + 文档整理
**Package**: carlo-rs
**Branch**: `master`

### Summary

P1: 新建 output/resulttools.rs，实现 dataframe、measurement_from_obs、recursive_stack、make_scalar。P2: 新增 merge_task_results 接入 register_evaluables 回调。P3/P4: concatenate_results 和 read_progress 已存在。对齐 cli_merge 中 results.json 格式为 {task, parameters, results}。新增 .githooks/pre-commit (fmt+clippy+test)。整理 CLAUDE.md、README.md 及记忆系统。

### Main Changes

(Add details)

### Git Commits

| Hash | Message |
|------|---------|
| `98e643c` | (see git log) |
| `979d96e` | (see git log) |
| `8a67ae0` | (see git log) |

### Testing

- [OK] (Add test results)

### Status

[OK] **Completed**

### Next Steps

- None - task complete


## Session 2: Carlo.rs JSON 格式对齐 Carlo.jl + CLAUDE.md 精简

**Date**: 2026-05-08
**Task**: Carlo.rs JSON 格式对齐 Carlo.jl + CLAUDE.md 精简
**Package**: carlo-rs
**Branch**: `master`

### Summary

ResultObservable 自定义 Serialize 对齐 Carlo.jl JSON.lower()（rebin_len/autocorr_time scalar/rebin_count/internal_bin_len）。修复 pre-commit hook 去重 bug，去掉 hook 中的 cargo fmt。CLAUDE.md 从 152 行精简到 50 行，修正项目描述。记忆系统初始化。

### Main Changes

(Add details)

### Git Commits

| Hash | Message |
|------|---------|
| `56c6f1d` | (see git log) |
| `6c490bd` | (see git log) |
| `6f7ae36` | (see git log) |

### Testing

- [OK] (Add test results)

### Status

[OK] **Completed**

### Next Steps

- None - task complete


## Session 3: Complete MPI backend: results transfer, multi-rank coordination, comm-aware sweep, checkpoint sync

**Date**: 2026-05-08
**Task**: Complete MPI backend: results transfer, multi-rank coordination, comm-aware sweep, checkpoint sync
**Package**: carlo-rs
**Branch**: `master`

### Summary

Brought MPI backend from ~70% skeleton to ~95% parity with Carlo.jl scheduler_mpi.jl. Added results transfer protocol over MPI, controller result aggregation, per-task run_id tracking, multi-rank run_comm coordination (run_follower with broadcast), MonteCarlo sweep_with_comm/measure_with_comm trait methods, Run::step_with_comm, MPI-coordinated checkpoint (rank-specific paths, broadcast existence). Created quality-gate-commit-sync skill combining trellis-check → commit → neat-freak.

### Main Changes

(Add details)

### Git Commits

| Hash | Message |
|------|---------|
| `9b59517` | (see git log) |
| `5964049` | (see git log) |

### Testing

- [OK] (Add test results)

### Status

[OK] **Completed**

### Next Steps

- None - task complete


## Session 4: CMC.rs Phase 1: full rewrite with layered architecture

**Date**: 2026-05-18
**Task**: CMC.rs Phase 1: full rewrite with layered architecture
**Package**: carlo-rs
**Branch**: `main`

### Summary

Designed and implemented CMC.rs rewrite. 6-layer architecture: Lattice (adjacency list), System (pub fields), Model (stateless physics), Algorithm (sweep directly mutates energy), ProposalStrategy (independent, OPSS adaptive), ClassicalMC wrapper (impl MonteCarlo + FromParams). 4 models (Ising, Potts, XY, Heisenberg), 3 algorithms (Metropolis, Wolff, SW), 2 proposals (Standard, OPSS). Wrote Carlo.rs lib.rs API docs, filled code-specs, synced CLAUDE.md/README.md/memory. 166 tests, 0 warnings, 0 unsafe.

### Main Changes

(Add details)

### Git Commits

| Hash | Message |
|------|---------|
| `0932ddc` | (see git log) |
| `de925ed` | (see git log) |

### Testing

- [OK] (Add test results)

### Status

[OK] **Completed**

### Next Steps

- None - task complete


## Session 5: CMC.rs high-priority: derived observables, JSON checkpoint, MultiSpinIsing MonteCarlo integration

**Date**: 2026-05-25
**Task**: CMC.rs high-priority: derived observables, JSON checkpoint, MultiSpinIsing MonteCarlo integration
**Package**: carlo-rs
**Branch**: `main`

### Summary

TDD from easiest to hardest: (1) postprocess.rs — susceptibility, specific_heat, binder_cumulant from E²/M²/M⁴ moments in ClassicalMC::measure(). (2) JSON snapshot checkpoint for ClassicalMC. (3) MultiSpinIsing refactor — fixed bit-plane counting bug, added system/model owned fields, impl MonteCarlo+FromParams+ParallelTemperingCompatible for Scheduler integration. 68 CMC + 128 Carlo tests pass, clippy clean.

### Main Changes

(Add details)

### Git Commits

| Hash | Message |
|------|---------|
| `67b9ab3` | (see git log) |
| `f1edd78` | (see git log) |

### Testing

- [OK] (Add test results)

### Status

[OK] **Completed**

### Next Steps

- None - task complete


## Session 6: CMC.rs medium-priority: lattice_type param, non-square lattice validation

**Date**: 2026-05-25
**Task**: CMC.rs medium-priority: lattice_type param, non-square lattice validation
**Package**: carlo-rs
**Branch**: `main`

### Summary

TDD from easiest to hardest: (1) lattice_type param in build_lattice_from_params() routes to triangular/honeycomb/kagome builders. (2) Non-square lattice smoke tests: triangular ferro Ising + kagome AF Ising end-to-end validation. 73 CMC + 128 Carlo tests, clippy clean. Skipped OPSS benchmark (criterion setup non-critical).

### Main Changes

(Add details)

### Git Commits

| Hash | Message |
|------|---------|
| `d7a7830` | (see git log) |
| `c512e49` | (see git log) |

### Testing

- [OK] (Add test results)

### Status

[OK] **Completed**

### Next Steps

- None - task complete


## Session 7: CMC.rs low-priority: continuous heat-bath, over-relaxation, spatial observables

**Date**: 2026-05-26
**Task**: CMC.rs low-priority: continuous heat-bath, over-relaxation, spatial observables
**Package**: carlo-rs
**Branch**: `main`

### Summary

TDD from easiest to hardest: (1) MicrocanonicalCore — over-relaxation via reflect_spin(), energy exactly preserved. (2) ContinuousHeatBathable trait + ContinuousHeatBathCore — Heisenberg vMF inverse-CDF (no Bessel), XY Best-Fisher von Mises rejection. (3) compute_correlation_1d() for spatially-resolved G(r). 83 CMC + 128 Carlo tests, zero new dependencies, clippy clean.

### Main Changes

(Add details)

### Git Commits

| Hash | Message |
|------|---------|
| `fc97611` | (see git log) |
| `af307e7` | (see git log) |
| `d06a4a2` | (see git log) |
| `819f4ae` | (see git log) |

### Testing

- [OK] (Add test results)

### Status

[OK] **Completed**

### Next Steps

- None - task complete
