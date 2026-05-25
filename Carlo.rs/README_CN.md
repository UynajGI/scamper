# Carlo.rs 用户手册

Monte Carlo 模拟框架 - Rust 实现，与 [Carlo.jl](https://github.com/lukas-weber/Carlo.jl) **100% 核心功能对齐**

---

## 目录

1. [概述](#概述)
2. [快速开始](#快速开始)
3. [核心概念](#核心概念)
4. [实现你的第一个模型](#实现你的第一个模型)
5. [运行模拟](#运行模拟)
6. [结果分析](#结果分析)
7. [高级功能](#高级功能)
8. [API 参考](#api-参考)
9. [常见问题](#常见问题)

---

## 概述

### Carlo.rs 是什么？

Carlo.rs 是一个用于开发高性能蒙特卡洛模拟的 Rust 框架。它处理所有与模型无关的任务：

| 功能 | 描述 |
|------|------|
| **误差分析** | 自动 binning、jackknife 重采样、自相关时间估计、decorrelated 模式、协方差矩阵 |
| **并行执行** | Rayon 多线程 / MPI 分布式 |
| **检查点** | 保存/恢复模拟状态，支持长时间运行 |
| **结果合并** | 多运行结果聚合与统计 |
| **复数支持** | 原生复数观测量支持，re/im 分离存储 |
| **性能监控** | 扫描速率跟踪、耗时显示 |

你需要做的只是实现蒙特卡洛更新和观测量 —— 框架处理剩下的一切。

### 设计理念

```
┌─────────────────────────────────────────────────────────┐
│                    你的代码                              │
│  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐     │
│  │   IsingMC   │  │   SSEMC     │  │   YourMC    │     │
│  │  sweep()    │  │  sweep()    │  │  sweep()    │     │
│  │  measure()  │  │  measure()  │  │  measure()  │     │
│  └──────┬──────┘  └──────┬──────┘  └──────┬──────┘     │
└─────────┼────────────────┼────────────────┼─────────────┘
          │                │                │
          └────────────────┼────────────────┘
                           ▼
┌─────────────────────────────────────────────────────────┐
│                    Carlo.rs 框架                         │
│  ┌─────────── ┐  ┌─────────── ┐  ┌─────────── ┐         │
│  │ Scheduler │  │  Backend  │  │  Context  │           │
│  │ (调度器)   │  │ (并行后端) │  │ (运行时)   │           │
│  └─────────── ┘  └─────────── ┘  └─────────── ┘         │
│  ┌─────────── ┐  ┌─────────── ┐  ┌─────────── ┐         │
│  │Measurements│  │   Merge   │  │Checkpoint │          │
│  │ (测量累积) │  │ (结果合并) │  │ (检查点)   │          │
│  └─────────── ┘  └─────────── ┘  └─────────── ┘         │
└─────────────────────────────────────────────────────────┘
```

---

## 快速开始

### 安装依赖

```bash
# Ubuntu/Debian
sudo apt-get install libhdf5-dev openmpi-bin libopenmpi-dev

# macOS
brew install hdf5 open-mpi
```

### 添加到项目

```toml
# Cargo.toml
[dependencies]
carlo-rs = { path = "path/to/Carlo.rs" }
rand_xoshiro = "0.8"
rand = "0.8"
```

### 最小示例

```rust
use carlo_rs::{MonteCarlo, Context, CarloError, FromParams, Params, Scheduler, RunConfig, RayonBackend};
use rand_xoshiro::Xoshiro256PlusPlus;
use rand::Rng;

// 1. 定义你的模型
struct MyModel {
    // 你的状态变量
}

// 2. 实现 MonteCarlo trait
impl MonteCarlo for MyModel {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        // 执行一次蒙特卡洛扫描（更新构型）
        // 使用 ctx.rng 获取随机数
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        // 测量观测量
        ctx.measure("Energy", 1.0);
    }
}

// 3. 实现 FromParams trait（从参数构建模型）
impl FromParams for MyModel {
    fn from_params(params: &Params, rng: &mut Self::Rng) -> Result<Self, CarloError> {
        Ok(MyModel { /* 初始化 */ })
    }
}

// 4. 运行模拟
fn main() {
    let backend = RayonBackend::new(4); // 4 线程
    let config = RunConfig {
        thermalization_sweeps: 1000,
        measurement_sweeps: 10000,
        binsize: 100,
        base_seed: 42,
        progress_interval: 1000,
        checkpoint_interval: 0,
    };
    let scheduler = Scheduler::new(backend, config);

    let params = Params::new();
    let results = scheduler.run_one::<MyModel>(&params);

    println!("Results: {:?}", results);
}
```

---

## 核心概念

### MonteCarlo Trait

这是你需要实现的唯一核心 trait：

```rust
pub trait MonteCarlo: Sized {
    /// 指定 RNG 类型（必须实现 Rng + SeedableRng + Send）
    type Rng: Rng + SeedableRng + Send;

    /// 执行一次蒙特卡洛扫描
    /// 这是更新构型的地方
    fn sweep(&mut self, ctx: &mut Context<Self::Rng>);

    /// 测量观测量（可选，默认为空）
    /// 在热化完成后每次扫描调用
    fn measure(&mut self, _ctx: &mut Context<Self::Rng>) {}

    /// 算法名称（可选）
    fn name(&self) -> &'static str { "UnnamedMC" }
}
```

### Context

运行时上下文，持有：
- **RNG**: 随机数生成器
- **Measurements**: 观测量累积器
- **Sweep counter**: 扫描计数器

```rust
// 使用 RNG
let r: f64 = ctx.rng.random();
let idx: usize = ctx.rng.random_range(0..n);

// 记录测量
ctx.measure("Energy", energy);
ctx.measure("Magnetization", mag);

// 数组观测量
ctx.measure_array("SpinCorrelation", &correlation_data);

// 复数观测量
ctx.measure_complex("OrderParameter", re, im);

// 检查状态
if ctx.is_thermalized() { /* ... */ }
println!("Sweep: {}", ctx.sweep_count());
```

### Params

参数字典，从 JSON 或代码创建：

```rust
// 从代码创建
let mut params = Params::new();
params.set("L", 100);
params.set("beta", 0.5);

// 获取参数
let l: usize = params.get("L").unwrap_or(100);
let beta: f64 = params.get("beta").unwrap_or(1.0);

// 从 JSON 文件加载
// params.json:
// {"L": 100, "beta": 0.5, "J": 1.0}
```

### RunConfig

模拟配置：

```rust
let config = RunConfig {
    thermalization_sweeps: 1000,  // 热化扫描数
    measurement_sweeps: 10000,    // 测量扫描数
    binsize: 100,                  // bin 大小
    base_seed: 42,                 // 基础种子
    progress_interval: 1000,       // 进度报告间隔
    checkpoint_interval: 0,        // 检查点间隔（0 = 禁用）
};
```

---

## 实现你的第一个模型

### 完整示例：一维 Ising 模型

```rust
use carlo_rs::{
    MonteCarlo, Context, CarloError, FromParams, Params,
    Scheduler, RunConfig, RayonBackend, Results
};
use rand_xoshiro::Xoshiro256PlusPlus;
use rand::Rng;
use std::collections::HashMap;

/// 一维 Ising 模型
struct Ising1D {
    spins: Vec<i8>,     // 自旋数组 (+1 或 -1)
    beta: f64,          // 逆温度
    j: f64,             // 耦合常数
}

impl MonteCarlo for Ising1D {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        let n = self.spins.len();

        // 执行 N 次 Metropolis 更新
        for _ in 0..n {
            // 随机选择一个自旋
            let i = ctx.rng.random_range(0..n);

            // 计算能量变化
            let left = self.spins[(i + n - 1) % n];
            let right = self.spins[(i + 1) % n];
            let dE = 2.0 * self.j * self.spins[i] as f64 * (left + right) as f64;

            // Metropolis 准则
            if dE <= 0.0 || ctx.rng.random::<f64>() < (-self.beta * dE).exp() {
                self.spins[i] *= -1;
            }
        }
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        let n = self.spins.len() as f64;

        // 磁化强度
        let m: f64 = self.spins.iter().sum::<i8>() as f64 / n;
        ctx.measure("Magnetization", m.abs());

        // 能量
        let mut e = 0.0;
        for i in 0..self.spins.len() {
            let j = (i + 1) % self.spins.len();
            e -= self.j * self.spins[i] as f64 * self.spins[j] as f64;
        }
        ctx.measure("Energy", e / n);
    }

    fn name(&self) -> &'static str { "Ising1D" }
}

impl FromParams for Ising1D {
    fn from_params(params: &Params, _rng: &mut Self::Rng) -> Result<Self, CarloError> {
        let l = params.get::<usize>("L").ok_or_else(|| CarloError::InvalidConfig {
            field: "L".into(),
            reason: "Missing system size L".into(),
        })?;
        let beta = params.get::<f64>("beta").unwrap_or(1.0);
        let j = params.get::<f64>("J").unwrap_or(1.0);

        // 随机初始化（或有序初始化）
        Ok(Self {
            spins: vec![1; l], // 全向上
            beta,
            j,
        })
    }
}

// 运行模拟
fn main() -> Result<(), CarloError> {
    // 创建后端
    let backend = RayonBackend::new(4);

    // 配置
    let config = RunConfig {
        thermalization_sweeps: 10000,
        measurement_sweeps: 100000,
        binsize: 1000,
        base_seed: 42,
        progress_interval: 10000,
        checkpoint_interval: 0,
    };

    // 参数
    let mut params = Params::new();
    params.set("L", 100);
    params.set("beta", 0.5);
    params.set("J", 1.0);

    // 运行
    let scheduler = Scheduler::new(backend, config);
    let results = scheduler.run_one::<Ising1D>(&params)?;

    // 输出结果
    for (name, est) in results.estimates() {
        println!("{}: {:.6} ± {:.6}", name, est.mean, est.stderr);
    }

    Ok(())
}
```

### 编译运行

```bash
# 编译
cargo build --release --features hdf5

# 运行
./target/release/your_simulation
```

---

## 运行模拟

### 单任务运行

```rust
let scheduler = Scheduler::new(backend, config);
let results = scheduler.run_one::<YourMC>(&params)?;
```

### 多任务并行

```rust
// 创建多个参数集
let params_list: Vec<Params> = vec![
    params_with("beta", 0.1),
    params_with("beta", 0.2),
    params_with("beta", 0.3),
];

// 并行运行
let results = scheduler.run_parallel::<YourMC>(&params_list);

for (i, r) in results.iter().enumerate() {
    println!("Task {}: Energy = {:.4} ± {:.4}",
        i,
        r.get("Energy").unwrap().mean,
        r.get("Energy").unwrap().stderr,
    );
}
```

### MPI 分布式运行

需要 `mpi` feature：

```bash
# 编译
cargo build --release --features "hdf5 mpi"

# 运行（16 核）
mpirun -np 16 ./target/release/your_simulation
```

代码中：

```rust
use carlo_rs::{MpiBackend, run_distributed, DistributedConfig, TaskSpec};

fn main() -> Result<(), CarloError> {
    let backend = MpiBackend::new()?;

    if backend.is_controller() {
        // 控制器逻辑
        let config = DistributedConfig {
            run_config: RunConfig::default(),
            ranks_per_run: 1,
            run_time: None,
            checkpoint_time: None,
            job_dir: PathBuf::from("./results"),
            tasks: vec![
                TaskSpec { id: 0, target_sweeps: 10000, thermalization: 1000, params: params1 },
                TaskSpec { id: 1, target_sweeps: 10000, thermalization: 1000, params: params2 },
            ],
        };

        let results = run_distributed::<YourMC, _>(config)?;
    }

    Ok(())
}
```

---

## 结果分析

### 获取结果

```rust
let results = scheduler.run_one::<YourMC>(&params)?;

// 获取单个观测量
if let Some(energy) = results.get("Energy") {
    println!("Energy: {:.6} ± {:.6}", energy.mean, energy.stderr);
    println!("Autocorrelation time: {:.2}", energy.autocorr_time);
    println!("Number of bins: {}", energy.n_bins);
}

// 遍历所有观测量
for (name, estimate) in results.estimates() {
    println!("{}: {}", name, estimate.format()); // "1.234567 ± 0.012345"
}
```

### 合并多个结果

```rust
// 多次运行的结果
let results1 = scheduler.run_one::<YourMC>(&params)?;
let results2 = scheduler.run_one::<YourMC>(&params)?;

// 合并（加权平均）
let merged = Results::merge(&[results1, results2]);
```

### Jackknife 分析

用于计算衍生观测量：

```rust
use carlo_rs::{jackknife, Evaluator};

// 从 HDF5 文件读取
let observables = merge_results(&task_dir, &MergeOptions::default())?;

// 创建评估器
let mut evaluator = Evaluator::new(observables, true);

// 定义衍生观测量
// 例如：比热 C = beta^2 * (⟨E²⟩ - ⟨E⟩²)
let result = jackknife(
    |samples| {
        // samples[0] = Energy, samples[1] = Energy^2
        let e = &samples[0];
        let e2 = &samples[1];
        e2 - e * e // variance
    },
    &[energy_samples, energy_sq_samples],
    false,
)?;
```

---

## 高级功能

### 数组观测量

支持记录数组类型的观测量：

```rust
// 记录一维数组
let correlation: Vec<f64> = vec![1.0, 0.8, 0.6, 0.4];
ctx.measure_array("Correlation", &correlation);

// 框架自动处理 shape 推断和统计计算
```

### 复数观测量

原生支持复数观测量（如序参数、格林函数）：

```rust
// 记录复数观测量
ctx.measure_complex("OrderParameter", re, im);

// 结果存储为 {re, im} 格式，匹配 Carlo.jl
```

### 内部性能计时

框架自动记录内部操作耗时（以 `_ll_` 前缀标识）：

- `_ll_sweep_time` - 单次扫描耗时
- `_ll_measure_time` - 测量耗时
- `_ll_checkpoint_read_time` - 检查点读取耗时
- `_ll_checkpoint_write_time` - 检查点写入耗时

```rust
// 结果中包含计时信息
if let Some(sweep_time) = results.get("_ll_sweep_time") {
    println!("Average sweep time: {:.3} s", sweep_time.mean);
}
```

### MultiplexEvaluator（并行回火链求值）

用于并行回火模拟中对所有链同时求值：

```rust
use carlo_rs::MultiplexEvaluator;

// 为 4 个 PT 链创建求值器
let mut multi_eval = MultiplexEvaluator::new(4);

// 为每个链注册求值函数
for chain_idx in 0..4 {
    multi_eval.evaluate("OrderParameter", &["Magnetization"], move |args| {
        // 链特定的计算
        args[0].clone() * temperature_factors[chain_idx]
    });
}

// 运行求值并堆叠结果
multi_eval.run_evaluations(&mut evaluator)?;
```

### 检查点（Checkpoint）

保存模拟状态以恢复：

```rust
// 启用检查点
let config = RunConfig {
    checkpoint_interval: 10000, // 每 10000 扫描保存一次
    ..
};

// Run 结构体支持检查点
#[cfg(feature = "hdf5")]
let run = Run::read_checkpoint(&path, &params, &config, seed)?;

// 写入检查点
run.write_checkpoint(&path)?;
```

检查点包含：
- RNG 状态
- 测量累积器（包括部分 bin）
- 扫描计数
- 模型特定状态（需要实现 `MonteCarloCheckpoint`）

### 并行回火（Parallel Tempering）

用于改进采样：

```rust
use carlo_rs::{
    ParallelTemperingConfig, ParallelTemperingCompatible, ParallelTemperingMC
};

// 实现兼容 trait
impl ParallelTemperingCompatible for YourMC {
    fn log_weight_ratio(&self, param: &str, new_value: f64) -> f64 {
        // 返回 log(w(new)/w(old))
        match param {
            "beta" => /* ... */,
            _ => 0.0,
        }
    }

    fn change_parameter(&mut self, param: &str, new_value: f64) {
        // 更新参数
    }
}

// 配置
let config = ParallelTemperingConfig {
    parameter: "beta".to_string(),
    values: vec![0.1, 0.5, 1.0, 2.0, 5.0],
    interval: 100, // 每 100 扫描尝试交换
};

// 创建回火 MC
let pt_mc = ParallelTemperingMC::new(&config, chain_idx, base_mc);
```

### 自定义 RNG

```rust
use rand_xoshiro::Xoshiro256PlusPlus;  // 默认
use rand_pcg::Pcg64;                    // 替代选项

impl MonteCarlo for YourMC {
    type Rng = Pcg64;  // 使用不同的 RNG
    // ...
}
```

---

## API 参考

### 核心 Traits

| Trait | 必须实现 | 描述 |
|-------|---------|------|
| `MonteCarlo` | `sweep()`, `Rng` | 核心蒙特卡洛算法 |
| `FromParams` | `from_params()` | 从参数构建模型 |
| `MonteCarloCheckpoint` | `write_checkpoint()`, `read_checkpoint()` | 检查点支持（可选） |

### 主要结构体

| 结构体 | 用途 |
|--------|------|
| `Context<R>` | 运行时上下文 |
| `Params` | 参数字典 |
| `RunConfig` | 模拟配置 |
| `Results` | 模拟结果 |
| `Estimate` | 统计估计（均值、误差） |
| `Run<MC, R>` | 单次运行生命周期 |
| `Scheduler` | 调度器 |

### 后端

| 后端 | 特点 |
|------|------|
| `RayonBackend` | 多线程，适合单节点 |
| `MpiBackend` | 分布式，适合多节点 |

---

## 常见问题

### Q: 如何设置不同的随机种子？

```rust
let config = RunConfig {
    base_seed: 12345,  // 每个任务使用 base_seed + task_id
    ..
};
```

### Q: 如何处理复杂的观测量？

```rust
// 在 sweep() 中累积
struct MyMC {
    energy_accum: f64,
    count: usize,
}

impl MonteCarlo for MyMC {
    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        // ... 更新
        self.energy_accum += current_energy;
        self.count += 1;
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        // 计算并记录平均值
        ctx.measure("Energy", self.energy_accum / self.count as f64);
        // 重置或继续累积
    }
}
```

### Q: 如何实现数组观测量？

```rust
// 直接使用 measure_array
let correlation: Vec<f64> = vec![1.0, 0.8, 0.6, 0.4];
ctx.measure_array("Correlation", &correlation);

// 结果会自动计算每个分量的均值、误差和协方差矩阵
```

### Q: 如何实现复数观测量？

```rust
// 使用 measure_complex
ctx.measure_complex("OrderParameter", re, im);

// 结果存储为 {re, im} 格式
if let Some(complex_est) = results.get_complex("OrderParameter") {
    println!("Re: {:.6} ± {:.6}", complex_est.re.mean, complex_est.re.stderr);
    println!("Im: {:.6} ± {:.6}", complex_est.im.mean, complex_est.im.stderr);
}
```

### Q: MPI 运行时如何调试？

```rust
// 在 worker 代码中
eprintln!("[Rank {}] Task {}, sweep {}", rank, task_id, sweep);

// 使用小规模测试
mpirun -np 2 cargo test --features mpi
```

### Q: 如何选择 binsize？

经验法则：
- binsize ≈ 自相关时间 × 10
- 如果不确定，从小开始（100-1000）
- 检查 autocorr_time 估计

---

## 从 Carlo.jl 迁移

如果你熟悉 Carlo.jl，这里是映射关系：

| Carlo.jl | Carlo.rs |
|----------|----------|
| `AbstractMC` | `MonteCarlo` trait |
| `monte_carlo_sweep!` | `MonteCarlo::sweep()` |
| `measure!` | `ctx.measure()` |
| `MCContext` | `Context` |
| `Params` | `Params` |
| `run_distributed` | `run_distributed()` |
| `merge_results` | `merge_results()` |
| `Evaluator` | `Evaluator` |
| `MultiplexEvaluator` | `MultiplexEvaluator` |
| `_ll_sweep_time` 等 | 相同（自动记录） |
| 复数（ResultTools） | 原生 `measure_complex()` |

主要区别：
1. Rust 需要显式类型声明（`type Rng = ...`）
2. 错误处理使用 `Result<T, CarloError>` 而非异常
3. HDF5 功能需要 `hdf5` feature
4. 数组观测量通过 `measure_array()` 原生支持
5. 复数观测量通过 `measure_complex()` 原生支持
6. 进度条和性能监控是内置的

---

## 许可证

Apache-2.0