# brainstorm: QMC.rs 世界线 QMC 架构设计

## Goal

设计 QMC.rs 的 crate 架构，从世界线 QMC（path integral + worm algorithm）开始实现，预留其他 QMC 方法扩展点。

## What I already know

- QMC.rs 已清空，只剩 `Cargo.toml`（依赖 carlo-rs）
- Carlo.rs 提供核心框架：`MonteCarlo` trait、`Context`、`Scheduler`、`FromParams`
- 世界线 QMC ≠ SSE（两种不同方法）
- Continuous 是主力实现路径，Discrete 作为完整性补充

## Decision (ADR-lite) — 单格点世界线对象

**Context**: 需要在虚时 [0, β) 上表示单个格点的世界线，支持离散和连续虚时。

**Decision**:
- 泛型 trait `Worldline` — 零成本抽象，编译期单态化
- State 用 `u8` 存状态索引，`dim: u8` 记录维度，物理映射由上层模型负责
- Continuous 直接存储 `kinks: Vec<(tau, from, to)>`，零拷贝 `kinks()`
- Discrete 存 `states: Box<[u8]>`，kinks 每次分配 Vec
- trait 用迭代器返回 kinks，不用 `&[..]`，兼容两种实现

**Consequences**:
- 优点：热路径（worm update）零虚表开销，Continuous 零拷贝遍历 kinks
- 缺点：上层算法需要对 `W: Worldline` 泛型化

## Requirements

- QMC.rs 作为 carlo-rs 的下游 crate
- 世界线 QMC = worldline 配置 + worm algorithm
- 架构预留 QMC 类型扩展点（worldline、SSE、DQMC 等）

## Acceptance Criteria

* [ ] `Worldline` trait + `ContinuousWorldline` + `DiscreteWorldline` 实现 + 单元测试
* [ ] Continuous 的 insert_kink / remove_kink / state_at / diagonal 正确

## Out of Scope

- 多格点、格子、Hamiltonian
- Worm algorithm
- SSE、DQMC、VMC
