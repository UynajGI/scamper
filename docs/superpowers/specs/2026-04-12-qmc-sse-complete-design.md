# QMC.rs SSE 完整实现设计

**日期**: 2026-04-12
**状态**: 设计中

## 概述

在 Carlo.rs 框架之上完成 QMC.rs 的 SSE（随机级数展开）算法实现。依托 Carlo.rs 的调度器、测量合并、MPI 后端等基础设施，构建一个支持任意格点拓扑、多种量子模型的 QMC 模拟库。

设计哲学：格点统一用邻接表表示，可表达任意维度、任意形状的模型。

## 架构总览

```
Carlo.rs (框架层 - 已完成)
├── MonteCarlo trait, FromParams, Scheduler, Backend
├── Context (RNG, measurements, sweep counter)
└── Results + Merge (rebinning, autocorrelation, jackknife)

QMC.rs (QMC 算法层)
├── lattice/          邻接表拓扑 (已完成)
├── hilbert/          希尔伯特空间 (部分完成)
├── sse/              SSE 算法引擎 (进行中)
│   ├── engine.rs         SSEEngine, Vertex, OperatorSequence
│   ├── diagonal.rs       对角更新
│   ├── loop_.rs          Worm 遍历
│   ├── vertex_data.rs    运行时散射表 (LP 求解)
│   ├── vertex_list.rs    世界线拓扑
│   ├── measurements.rs   标准测量
│   └── improved.rs       Improved estimators (新增)
├── models/           物理模型
│   ├── heisenberg.rs       S=1/2 Heisenberg (已完成)
│   └── xxz.rs              XXZ 模型 (新增)
├── ed/               Lanczos 精确对角化 (新增)
│   ├── lanczos.rs
│   └── hamiltonian.rs
└── tests/
    ├── sse_algorithm_test.rs   单元测试
    ├── physics_test.rs         文献基准测试
    └── ed_test.rs              ED 自动验证
```

## 参考依据

### 算法参考
- **StochasticSeriesExpansion.jl** — Julia 参考实现，使用 LP 求解散射表、worm 遍历、vertex list 世界线拓扑
- **Sandvik 1991** (PhysRevB.43.5950) — 原始 SSE 方法，S=1/S=3/2 链的基态性质
- **Sandvik 1997** (9703200v1) — Directed loop 方程，连续时间世界线
- **Evertz 1997** (9707221v3) — Loop 算法综述，XXZ/Heisenberg 散射规则，性能基准
- **Beard & Wiese 1996** (9602164v1) — 连续时间 loop 算法，2D Heisenberg 手征微扰理论验证

### 可验证数值基准
| 来源 | 系统 | 可观测值 | 期望值 |
|------|------|----------|--------|
| Bethe ansatz | 1D Heisenberg S=1/2 | E/N (T=0) | 1/4 - ln(2) ≈ -0.443147 |
| Sandvik 1991 | 1D S=1, N=64 | χ(π) (T=0) | 20.0 ± 1.5 |
| Beard & Wiese | 2D Heisenberg | ρ_s | 0.185(2) |
| Beard & Wiese | 2D Heisenberg | ℏc | 1.68(1) |
| Beard & Wiese | 2D Heisenberg | M_s | 0.3083(2) |
| Evertz 1997 | 六顶点模型 | z_MC | ~0 (无临界减速) |

---

## 阶段 1: 修复 SSE Worm 算法

**目标**: Worm 算法正确采样反铁磁构型空间，Bethe ansatz 测试通过。

### 1.1 VertexData — 运行时散射表

用 LP 求解器（类似 Julia 的 HiGHS）在模型初始化时计算散射概率表，而非硬编码。

**核心数据结构**（对标 Julia `VertexData`）:

```rust
pub struct Transition {
    pub offset: usize,   // transition_cumprobs 中的起始位置
    pub length: usize,   // 可能结果的数量
}

pub struct VertexData {
    /// 每个顶点的腿状态 [leg, vertex]
    pub leg_states: Vec<u8>,
    /// 顶点权重
    pub weights: Vec<f64>,
    /// 散射表 [leg_in, worm_in, vertex] → Transition
    pub transitions: Vec<Transition>,
    /// 累积概率（每个 transition 对应一段）
    pub cumprobs: Vec<f64>,
    /// 散射目标顶点索引
    pub targets: Vec<u8>,
    /// 出口 (leg_out, worm_out)
    pub step_outs: Vec<(usize, usize)>,
}

impl VertexData {
    /// 从顶点列表和权重构建散射表。
    /// 对于 S=1/2 Heisenberg，LP 会找到确定性解（所有概率 = 1）。
    pub fn build(vertices: &[VertexInfo], dims: (usize, usize)) -> Self;

    /// 执行散射：给定入射腿、worm 类型、当前顶点，返回出射腿和新顶点。
    pub fn scatter(&self, vertex_idx: usize, leg_in: usize, worm_in: usize, rng: &mut impl Rng)
        -> (usize, usize, usize); // (leg_out, worm_out, new_vertex_idx)
}
```

**LP 求解策略**:
- 优先使用纯 Rust LP 求解器（如 `good_lp` + `lpsolve` 或 `minilp`），避免 C 依赖
- 如果 LP 依赖太重，对 S=1/2 Heisenberg/XXZ 直接硬编码 LP 的解析解（文献已知）
- 保留 LP 接口，未来可支持任意模型

### 1.2 Vertex 结构

Vertex 不再用 `vertex_idx: u8` 硬编码语义，而是作为散射表中的索引：

```rust
pub struct Vertex {
    pub bond_idx: usize,
    pub vertex_idx: usize,  // 索引到 VertexData 的 leg_states/transitions
}
```

### 1.3 Worm 遍历修复

关键修正（基于对 Julia 代码的详细分析）:

1. **出口腿配对**: `leg_out = leg_in ^ 1`（0↔1, 2↔3），不是 `leg_in ^ 2`
2. **自旋翻转**: 散射时同时翻转 `leg_in` 和 `leg_out` 的自旋
3. **状态重建**: 从 `v_first` 读取输入腿的自旋状态
4. **Worm 计数**: 自适应（类似 Julia），目标 worm 总长度 ∝ 算符数

### 1.4 对角更新修正

已完成的部分:
- `diagonal_element` 改为 anti-aligned = 1.0, aligned = 0.0（匹配 Julia shift）
- 初始自旋改为 AFM 排列
- `diagonal_shift` 改为 `J * N_bonds / 4`

### 验证标准

```bash
cargo test test_heisenberg_chain_ground_state
# 期望: E/N = -0.443147 ± 3σ (16-site chain, β=10)
```

---

## 阶段 2: Improved Estimators

**目标**: 实现 loop-cluster improved estimators，显著提升测量精度。

### 2.1 原理

在 worm 遍历过程中，worldline 被分割成 cluster。每个 cluster 可以独立翻转。Improved estimator 利用这一结构：

- **关联函数** G(i,j) = 同一 cluster 中 i 和 j 相连的概率
- **均匀磁化率** χ = Σ_i G(0,i) / N
- **交错磁化率** χ_s = Σ_i (-1)^i G(0,i) / N

### 2.2 接口

```rust
pub struct ImprovedEstimators {
    /// 关联函数累加器 G(r) for r = 0..N/2
    correlation: Vec<f64>,
    /// 计数
    n_samples: usize,
}

impl ImprovedEstimators {
    /// 在一次 worm 遍历后更新。
    /// worm_vertices 是该 worm 访问的所有顶点。
    pub fn update_from_worm(&mut self, worm_vertices: &[(usize, usize)]);

    /// 在一次完整 sweep 后归一化。
    pub fn finalize(&self) -> HashMap<String, f64>;
}
```

在 `SSECore::measure()` 中调用，通过 `ctx.measure()` 记录：
- `Chi_uniform` — 均匀磁化率
- `Chi_staggered` — 交错磁化率
- `Correlation_r{r}` — 距离 r 的自旋关联

### 2.3 验证标准

- 磁化率的统计误差应比普通估计器小 5-10 倍
- 1D S=1/2 链的 χ(π) 在低温极限下应与 Sandvik 1991 的结果一致

---

## 阶段 3: XXZ 模型

**目标**: 支持各向异性 XXZ Hamiltonian

### 3.1 模型定义

```
H = J Σ_{<i,j>} [ S^x_i S^x_j + S^y_i S^y_j + Δ S^z_i S^z_j ]
```

Δ=1 时还原为 Heisenberg，Δ=0 时为 XY model，Δ→∞ 时为 Ising。

### 3.2 实现

```rust
pub struct XxzModel {
    lattice: Lattice,
    beta: f64,
    j: f64,
    delta: f64,  // 各向异性参数
}

impl SSEMonteCarlo for XxzModel {
    type HilbertSpace = SpinHalfHS;

    fn bond_operators(&self, _bond_type: BondType) -> Vec<(OpType, f64)> {
        vec![
            (OpType::Diagonal, self.j * self.delta * 0.5),
            (OpType::OffDiagonal, self.j * 0.5),
        ]
    }

    fn hilbert_space(&self) -> &SpinHalfHS { &HS }
    fn beta(&self) -> f64 { self.beta }
}
```

### 3.3 散射表重计算

XXZ 的散射表不同于 Heisenberg（Δ ≠ 1 时 bounce 概率非零）。VertexData 需要从 XXZ 的 bond Hamiltonian 重新计算。

### 3.4 验证标准

| Δ | 系统 | 期望行为 |
|---|------|----------|
| 0.0 | 1D 链 | XY model, E/N = -1/π ≈ -0.3183 |
| 0.5 | 1D 链 | 临界 XY 相 |
| 1.0 | 1D 链 | Heisenberg, E/N = -0.443147 |
| 2.0 | 1D 链 | Ising 反铁磁相 |
| 5.0 | 1D 链 | 强 Ising极限 |

---

## 阶段 4: Lanczos 精确对角化

**目标**: 内置 ED 验证器，自动验证 SSE 结果。

### 4.1 架构

```rust
pub struct SparseHamiltonian {
    /// CSR 格式的稀疏矩阵
    row_ptr: Vec<usize>,
    col_idx: Vec<usize>,
    values: Vec<f64>,
    dim: usize,
}

impl SparseHamiltonian {
    /// 从模型和格点构建 Hamiltonian 矩阵。
    pub fn from_model<M: SSEMonteCarlo>(model: &M) -> Self;

    /// Lanczos 迭代求基态能量。
    pub fn ground_state_energy(&self, tol: f64, max_iter: usize) -> f64;
}
```

### 4.2 测试策略

```rust
#[test]
fn test_ed_vs_sse_heisenberg() {
    // 小系统 ED
    let ed_energy = compute_ed_energy(&lattice, &model);

    // SSE 低温模拟
    let sse_energy = run_sse_simulation(&params);

    // z-score 比较 (应 ≤ 4σ)
    let z_score = (sse_energy.mean - ed_energy) / sse_energy.stderr;
    assert!(z_score.abs() < 4.0);
}
```

### 4.3 验证范围

| 系统 | N | 方法 |
|------|---|------|
| 1D Heisenberg PBC | 4, 6, 8, 10, 12 | ED vs SSE |
| 1D XXZ PBC | 4, 6, 8 | ED vs SSE (Δ=0.5, 1.0, 2.0) |
| 2D Heisenberg | 4, 9 | ED vs SSE |

---

## 阶段 5: 文献基准测试

**目标**: 实现参考文献中的具体数值验证。

### 5.1 Bethe Ansatz 测试（阶段 1 完成后即可运行）

```rust
#[test]
fn test_bethe_ansatz_energy() {
    // 1D Heisenberg S=1/2, N=16, β=10
    // 期望: E/N = 1/4 - ln(2) = -0.443147...
}
```

### 5.2 Beard & Wiese 2D 测试（需要 improved estimators）

```rust
#[test]
fn test_2d_spin_stiffness() {
    // 2D Heisenberg, L=6..20, β=1..100
    // 手征微扰理论拟合:
    // ρ_s = 0.185(2), c = 1.68(1), M_s = 0.3083(2)
}
```

### 5.3 Sandvik S=1 测试（需要 S=1 HilbertSpace）

```rust
#[test]
fn test_s1_staggered_susceptibility() {
    // 1D S=1 chain, N=64, T→0
    // χ(π) = 20.0 ± 1.5
}
```

---

## 依赖关系

```
阶段 1 (Worm 修复) ────────── 基础，所有后续依赖此项
         │
         ▼
阶段 2 (Improved Estimators) ─ 依赖正确的 worm 遍历
         │                            │
         ▼                            ▼
阶段 3 (XXZ)                    阶段 4 (ED 验证)
         │                            │
         ▼                            ▼
         └─────────→ 阶段 5 (文献基准) ← 需要以上全部
```

## 风险和缓解

| 风险 | 缓解 |
|------|------|
| Rust LP 求解器依赖复杂 | 先用已知解析解硬编码 S=1/2 Heisenberg/XXZ 的散射表，LP 作为未来扩展 |
| Improved estimators 实现复杂 | 先做最简单的关联函数估计器，逐步扩展 |
| Lanczos 收敛慢 | 限制在 N≤16 小系统，max_iter=100 足够 |
