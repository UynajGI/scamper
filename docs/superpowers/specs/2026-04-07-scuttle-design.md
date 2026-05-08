# Scuttle 设计规格

> 生产级 Rust 蒙特卡洛框架设计文档
> 日期：2026-04-07

---

## I. 项目定位

**Scuttle** 是一个 Rust 实现的生产级蒙特卡洛计算框架，包含三个 Cargo workspace 包：

| 包 | 角色 | 算法范围 |
|----|------|----------|
| **Carlo.rs** | 核心框架 | 抽象层、调度、测量、误差分析、输出 |
| **CMC.rs** | 经典蒙卡 | Metropolis、Wolff、Swendsen-Wang（Ising/Potts/XY） |
| **QMC.rs** | 量子蒙卡 | SSE、Path Integral、Worldline QMC（Heisenberg/XXZ/Hubbard） |

**依赖方向**：`CMC.rs → Carlo.rs`，`QMC.rs → Carlo.rs`，禁止循环/横向依赖。

---

## II. 核心设计决策

| 决策点 | 选择 | 理由 |
|--------|------|------|
| Carlo.rs 核心 trait | 方法钩子风格 + 默认空实现 | 用户实现 `sweep()` 一个方法即可；编译器内联优化充分 |
| 并行策略 | Phase 1 `rayon`，预留 `Backend` trait | 覆盖 90% 单机多核场景；为 MPI/GPU 扩展留接口 |
| 输出格式 | HDF5 默认 + JSON Lines 可选 | 与 Carlo.jl 兼容；物理社区标准 |
| Python 绑定 | Phase 1 不做，API 预留兼容性 | 单人开发资源有限；Phase 2 按需求实现 |
| RNG 管理 | 默认分块种子 + 可选 `--strict-repro` 跳跃序列 | 平衡简单高效与严格确定性 |
| 误差分析 | 均值 ± stderr + 自相关时间 τ | 物理计算主流需求；类似 Carlo.jl |

---

## III. Carlo.rs 核心 API

### MonteCarlo Trait

```rust
pub trait MonteCarlo: Sized {
    fn sweep(&mut self, ctx: &mut Context);
    fn measure(&mut self, ctx: &mut Context) {}  // 默认空
    fn save(&self, _out: &mut Hdf5Group) {}
    fn load(&mut self, _in: &Hdf5Group) {}
    fn name(&self) -> &'static str { "UnnamedMC" }
}
```

### Context

```rust
pub struct Context<R: RngCore + SeedableRng> {
    pub rng: R,
    measurements: Measurements,
    sweep_count: u64,
    thermalization_sweeps: u64,
    thermalized: bool,
}
```

### Backend Trait（并行抽象）

```rust
pub trait Backend: Clone + Send + Sync {
    type Rng: RngCore + SeedableRng + Send;
    fn spawn_tasks<F>(&self, n: usize, seed: u64, f: F)
        where F: Fn(usize, &mut Self::Rng) + Sync;
    fn barrier(&self);
}
```

Phase 1 实现：`RayonBackend`；Phase 2+ 预留：`MpiBackend`、`GpuBackend`。

---

## IV. Workspace 结构

```
scuttle/
├── Cargo.toml          # workspace 定义
├── justfile            # 唯一命令入口
├── rust-toolchain.toml
│
├── Carlo.rs/src/
│   ├── lib.rs
│   ├── monte_carlo.rs  # 核心 trait
│   ├── context.rs      # RNG + 状态
│   ├── measurements.rs # Observable 收集 + binning
│   ├── estimate.rs     # 误差分析
│   ├── backend/        # 并行抽象
│   ├── results.rs      # HDF5/JSON 输出
│   ├── scheduler.rs    # 任务调度
│   └── error.rs        # thiserror 错误定义
│
├── CMC.rs/src/
│   ├── models/         # Ising, Potts, XY (impl MonteCarlo)
│   ├── updates/        # Metropolis, Wolff, Swendsen-Wang
│   └── lattice.rs      # 2D/3D 晶格结构
│
├── QMC.rs/src/
│   ├── models/         # Heisenberg, XXZ, Hubbard
│   ├── sse/            # SSE state, directed loop
│   └── lattice.rs
│
└── examples/
```

---

## V. 执行流程

```
Params → MonteCarlo::from_params()
         ↓
     Scheduler.run()
         ↓
  [热化阶段] sweep × thermalization_sweeps
         ↓
  [测量阶段] sweep + measure × measurement_sweeps
         ↓
     Results.save_hdf5() / to_json_lines()
```

---

## VI. Phase 1 范围

**Carlo.rs**：
- MonteCarlo trait + Context
- RayonBackend
- Measurements + Estimate（binning 分析）
- Results（HDF5 + JSON 输出）
- Scheduler（单任务 + 并行任务）
- 单元测试 + 确定性验证

**CMC.rs**：
- IsingModel（2D/3D）
- Metropolis 单点更新
- Wolff 集群更新
- 与解析解对比测试

**QMC.rs**：
- SSEState + HeisenbergModel
- Directed Loop 更新
- 与 Carlo.jl/StochasticSeriesExpansion.jl 对比测试

---

## VII. Phase 2+ 扩展点

- MPI Backend（跨节点并行）
- GPU Backend（wgpu/cudarc）
- Python 绑定（PyO3）
- 更多模型（Hubbard、XY、Potts）
- 更多更新策略（Swendsen-Wang、worldline）

---

## VIII. 约束与边界

**不包含**：
- ❌ Phase 1 不做 MPI/GPU
- ❌ Phase 1 不做 Python 绑定
- ❌ 不做完整 MCMC 诊断（ESS、$\hat{R}$）
- ❌ 不做贝叶斯推断算法（HMC、NUTS）

**必须满足**：
- ✅ 确定性复现：同种子 → 同输出（`--strict-repro`）
- ✅ 与解析解误差 < 3σ
- ✅ `just check` 全绿才能提交
- ✅ HDF5 输出与 Carlo.jl 格式兼容

---

*文档版本：0.1 · 日期：2026-04-07*