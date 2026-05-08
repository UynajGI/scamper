# Quantum Monte Carlo Methods: A Discipline-Based Classification

## 引言：为什么QMC需要按学科起源分类

量子蒙特卡洛（Quantum Monte Carlo, QMC）这个名称下覆盖了来自**四个不同学科起源**的方法家族：

| 学科起源 | 空间结构 | 位置变量 | 基本对象 | 代表方法 |
|----------|----------|----------|----------|----------|
| **量子化学** | 连续空间 $\mathbb{R}^{3N}$ | 有坐标 $\mathbf{R}$ | 波函数 $\Psi(\mathbf{R})$ | VMC, DMC, GFMC |
| **凝聚态格点模型** | 离散格点 $i=1,\ldots,L$ | 无连续坐标 | 自旋/占据数 $n_i$, $\sigma_i$ | SSE, Worldline, Worm |
| **量子杂质物理** | 杂质+浴（场） | 杂质无空间扩展 | 自旋+玻色浴耦合 | CT-HYB, 连续虚时聚类 |
| **格点场论（高能物理）** | 格点场 $\phi(x,\tau)$ | 场变量，离散化 | 规范场构型 | HMC, Lattice QCD |

**核心问题**：现有文献和教材往往将这些方法混为一谈，导致学习者困惑。例如：

- 格点模型QMC中的"世界线"是自旋在虚时方向的轨迹，**不是粒子在实空间的位置**
- VMC/DMC中的坐标 $\mathbf{R}$ 是电子在连续空间的位置，**与格点指标 $i$ 完全不同**
- 格点QCD中的"格点"是场论离散化的产物，**与凝聚态格点模型的物理意义不同**

本综述的目标是**明确区分这些学科起源**，帮助读者理解：
1. 每类方法的适用场景和物理对象
2. 跨学科方法迁移时的注意事项
3. 为什么某些方法（如VMC）看似通用，实则隐含学科背景假设

---

## 第一部分：量子化学起源的QMC（连续空间，有坐标 $\mathbf{R}$）

### 1.1 物理背景

量子化学QMC起源于分子和固体的**电子结构计算**。基本对象是：

- **波函数** $\Psi(\mathbf{R})$，其中 $\mathbf{R} = (\mathbf{r}_1, \ldots, \mathbf{r}_N)$ 是 $N$ 个电子在连续空间的坐标
- **哈密顿量** $H = \sum_i \frac{-\hbar^2}{2m}\nabla_i^2 + V(\mathbf{R})$，包含动能和势能

**关键特征**：位置 $\mathbf{r}_i$ 是连续变量 $\in \mathbb{R}^3$，采样在实空间进行。

### 1.2 变分蒙特卡洛（Variational Monte Carlo, VMC）

#### 基本思想

用参数化试探波函数 $\Psi_T(\mathbf{R}, \boldsymbol{\alpha})$ 近似基态，能量期望值：

$$E(\boldsymbol{\alpha}) = \int d\mathbf{R} |\Psi_T(\mathbf{R}, \boldsymbol{\alpha})|^2 E_L(\mathbf{R}, \boldsymbol{\alpha})$$

其中**局域能量** $E_L = H\Psi_T/\Psi_T$。

**蒙特卡洛积分**：按 $|\Psi_T|^2$ 分布采样位形 $\{\mathbf{R}_i\}$，计算平均。

#### 试探波函数形式

典型的 **Slater-Jastrow 形式**：

$$\Psi_{SJ}(\mathbf{R}) = D^\uparrow(\mathbf{R}) D^\downarrow(\mathbf{R}) \times \exp\left[\sum_{i<j} u(r_{ij})\right]$$

其中：
- $D^{\uparrow,\downarrow}$ 是 Slater 行列式（来自 Hartree-Fock 或 DFT）
- Jastrow 因子引入电子-电子关联

#### 现代发展

**神经网络波函数**（2020年后）：

$$\Psi_{NN}(\mathbf{R}) = e^{\mathcal{F}_{\theta}(\mathbf{R})} \cdot D^\uparrow(\mathbf{R}) D^\downarrow(\mathbf{R})$$

其中 $\mathcal{F}_{\theta}$ 由神经网络参数化，可实现更复杂的关联。

#### VMC优化方法

**能量最小化**：调整参数 $\boldsymbol{\alpha}$ 使能量最小化。

梯度公式：

$$\frac{\partial E}{\partial \alpha_k} = 2\left(\langle E_L \frac{\partial \ln \Psi_T}{\partial \alpha_k} \rangle - \langle E_L \rangle \langle \frac{\partial \ln \Psi_T}{\partial \alpha_k} \rangle\right)$$

**优化算法**：
- **随机重构（Stochastic Reconfiguration）**：利用参数间关联矩阵，自然梯度下降
- **线性方法（Linear Method）**：在 $\Psi_T$ 附近构造基态的线性表示，对角化小矩阵
- **最小化方差**：最小化 $\langle (E_L - E)^2 \rangle$，避免能量导数噪声

#### 代表软件

| 软件 | 语言 | 特点 | 适用场景 |
|------|------|------|----------|
| **QMCPACK** | C++/CUDA | 大规模并行，支持GPU | 固体电子结构 |
| **PyQMC** | Python | 集成PySCF，易扩展 | 分子体系、教学 |
| **CASINO** | Fortran | 周期系统优化 | 固体物理 |
| **QWalk** | C++ | 专注于DMC | 分子体系 |
| **DeepQMC** | Python | 神经网络波函数 | 分子、小分子团簇 |

### 1.3 扩散蒙特卡洛（Diffusion Monte Carlo, DMC）

#### 投影原理

利用虚时演化 $e^{-\tau H}$ 的投影性质：

$$|\Psi_0\rangle = \lim_{\tau\to\infty} e^{-\tau H} |\Psi_T\rangle$$

将试探波函数投影到基态。

#### Schrödinger 方程的随机解释

虚时方程：

$$\frac{\partial \Psi}{\partial \tau} = \frac{1}{2}\sum_i \nabla_i^2 \Psi - [V(\mathbf{R}) - E_T]\Psi$$

分解为：
- **扩散项** $\frac{1}{2}\nabla^2\Psi$：对应行走者的随机游走
- **分支项** $-[V-E_T]\Psi$：对应行走者的产生/湮灭

#### 固定节点近似

对于费米子，波函数节点导致符号问题。固定节点近似限制 $\Psi$ 在试探波函数节点处为零，给出给定节点条件下的能量上界。

**节点定理**：对于费米子，试探波函数的能量满足：

$$E_{FN} \geq E_{exact}$$

其中 $E_{FN}$ 是固定节点近似下的能量。节点质量越好，能量越接近精确值。

**释放节点方法**：允许行走者穿过节点，但引入符号，最终恢复精确结果（如果符号问题可控）。

#### DMC算法细节

```
Algorithm: Diffusion Monte Carlo

1. 初始化：从 |Ψ_T|² 采样 N_w 个行走者 {R_i}
2. 时间步演化（每步 Δt）：
   a) 扩散：R_i → R_i + χ，其中 χ 是高斯随机位移
   b) 分支：计算权重 w_i = exp[-Δt(E_L(R_i) - E_T)]
   c) 复制/删除：按权重调整行走者数目
3. 更新参考能量：E_T → E_T + (1/N_w - 1/N_target)/τ
4. 重复直到平衡，测量能量
```

**重要采样**：引入试探波函数 $\Psi_T$，演化方程变为：

$$\frac{\partial f}{\partial \tau} = \frac{1}{2}\nabla^2 f - \nabla \cdot (f \mathbf{F}) - (E_L - E_T)f$$

其中 $f = \Psi\Psi_T$，漂移力 $\mathbf{F} = \nabla \ln |\Psi_T|$。

### 1.4 辅助场量子蒙特卡洛（Auxiliary Field QMC, AFQMC）

#### Hubbard-Stratonovich 变换

将电子-电子相互作用转化为单粒子问题：

$$e^{-\Delta\tau U n_\uparrow n_\downarrow} = \frac{1}{2}\sum_{x=\pm 1} e^{\gamma x (n_\uparrow - n_\downarrow)}$$

其中辅助场 $x$ 是离散随机变量。

#### Slater 行列式演化

配分函数写为 Slater 行列式的路径积分：

$$Z = \sum_{\{x_i\}} \text{Tr}\left[\prod_\tau B_\tau(x_\tau)\right]$$

每个时间步演化 Slater 行列式 $|\Phi(\tau)\rangle$。

#### 约束路径近似

限制行列式满足 $\langle \Phi_T | \Phi \rangle > 0$，解决相位问题（近似）。

**相位问题来源**：Slater行列式的复数权重导致：

$$\langle O \rangle = \frac{\sum_\phi w(\phi) O(\phi)}{\sum_\phi w(\phi)} \to \text{信号淹没在噪声中}$$

当 $w(\phi)$ 可正可负时，分母的方差指数增长。

**约束路径AFQMC（CP-AFQMC）**：
- 强制 $\langle \Phi_T | \Phi \rangle > 0$（Gauge条件）
- 给出能量上界：$E_{CP} \geq E_{exact}$
- 与固定节点DMC类似，受试探波函数质量影响

**无约束AFQMC**：
- 使用 walkers 携带复数权重
- 相位松弛技术控制符号问题
- 可达精确结果，但计算代价更高

#### AFQMC算法流程

```
Algorithm: Auxiliary Field QMC

1. 初始化：N_w 个 Slater 行列式 {|Φ_i⟩}
2. 时间步演化：
   a) 对每个时间片 τ：
      - 采样辅助场 {x_τ}
      - 演化行列式：|Φ_i⟩ → B_τ(x_τ)|Φ_i⟩
   b) 正交化：|Φ_i⟩ → |Φ_i⟩/⟨Φ_i|Φ_i⟩^{1/2}
   c) 重要性采样：按重叠 ⟨Φ_T|Φ_i⟩ 加权
3. 投影收敛后测量能量
```

**计算复杂度**：$O(N^3 M)$ 其中 $N$ 是电子数，$M$ 是时间片数。GPU并行可显著加速。

### 1.5 路径积分蒙特卡洛（Path Integral Monte Carlo, PIMC）

#### 基本原理：量子→经典映射

PIMC的核心思想是将量子统计力学问题转化为经典统计力学问题。

**量子配分函数**：

$$Z_Q = \text{Tr}(e^{-\beta H}) = \sum_n e^{-\beta E_n}$$

无法直接采样（量子态叠加）。

**路径积分变换**：利用 Trotter 分解：

$$Z_Q = \lim_{M\to\infty} \int d\mathbf{R}_1 \cdots d\mathbf{R}_M \, \exp\left[-\sum_{m=1}^{M} \frac{M}{\beta} \cdot \frac{(\mathbf{R}_{m+1}-\mathbf{R}_m)^2}{2\hbar^2} - \frac{\beta}{M} \sum_{m=1}^{M} V(\mathbf{R}_m)\right]$$

这等价于**经典弹性环链系统**：
- $M$ 个"珠子"（beads）$\mathbf{R}_1, \ldots, \mathbf{R}_M$ 串成环
- 相邻珠子间有"弹簧"连接（来自动能项）
- 每个珠子受势场 $V(\mathbf{R}_m)$ 作用

**量子弦的物理图像**：

```
经典粒子（T=∞）：    ●          单个点，无量子涨落

量子粒子（T有限）：  ●─●─●─●─●─●   M个珠子形成环链
                     ↑  "弹簧"连接
                     
玻色子：同一粒子的M个珠子是同一环链
费米子：不同粒子的环链可"交换"，形成交织拓扑
```

**环链周长**：量子涨落的度量

$$\langle R^2 \rangle = \langle |\mathbf{R}_m - \mathbf{R}_{m'}|^2 \rangle \propto \frac{\hbar^2}{mk_B T}$$

温度越低，环链越长（量子涨落越大）。

#### PIMC算法流程

```
Algorithm: Path Integral Monte Carlo

输入：粒子数 N，温度 T，时间片数 M
输出：热力学观测量（能量、密度分布等）

1. 初始化：生成 N 个环链，每个环链有 M 个珠子
   - 玻色子：各环链独立初始化
   - 费米子：需处理交换拓扑

2. 珠子移动更新：
   a) 单珠子移动：随机选择珠子 m，提议位移 δ
   b) 环链整体移动：所有珠子同步位移
   c) 部分段移动：选择连续段进行形状更新
   
3. 交换更新（费米子）：
   - 选择两个环链，尝试交换末端珠子
   - 接受概率包含交换权重符号
   
4. 测量：能量、密度分布、对关联函数等

5. 重复 2-4 直到收敛
```

#### Trotter离散化误差

**误差来源**：Trotter分解是近似：

$$e^{-\beta H} = \lim_{M\to\infty} \left(e^{-\beta H_0/M} e^{-\beta V/M}\right)^M$$

有限 $M$ 时，误差：

$$\Delta Z \propto \frac{\beta^2}{M^2} \langle [H_0, V] \rangle$$

**消除方法**：
- 增大 $M$（代价增加）
- 高阶Trotter公式（更复杂）

**典型 $M$ 取值**：$M \sim 10-100$ 对于分子；$M \sim 100-1000$ 对于强关联系统。

#### 玻色子与费米子的区别

| 特征 | 玻色子PIMC | 费米子PIMC |
|------|------------|------------|
| 环链拓扑 | 各环链独立 | 环链可交换交织 |
| 交换权重 | +1（无符号问题） | $(-1)^{N_{ex}}$（符号问题） |
| 采样难度 | 低（直接采样） | 高（符号问题） |
| 典型应用 | 超流氦、玻色凝聚 | 有限温电子气（受符号问题限制） |

**费米子的符号问题**：交换环链时权重符号改变：

$$w(\text{交换拓扑}) = (-1)^{\text{交换次数}} \times |w|$$

平均符号 $\langle s \rangle$ 随温度降低和粒子数增加指数衰减。

#### PIMC与其他量子化学QMC的关系

| 方法 | 目标态 | 温度 | 路径积分 |
|------|--------|------|----------|
| VMC | 基态（近似） | $T=0$ | 无 |
| DMC | 基态（投影） | $T=0$ | 虚时投影 |
| PIMC | 有限温平衡态 | $T>0$ | 完整路径积分 |
| AFQMC | 基态或有限温 | $T=0$ 或有限温 | 辅助场路径积分 |

**关键区别**：
- DMC：虚时演化到 $\tau\to\infty$，投影到基态
- PIMC：虚时周期 $\tau\in[0, \beta]$，采样有限温配分函数

#### PIMC的代表应用

| 应用领域 | 系统 | 观测量 |
|----------|------|--------|
| 超流氦-4 | $^4$He液体 | 超流密度、结构因子 |
| 固态氢 | H$_2$分子固体 | 相变、结构性质 |
| 电子气 | 有限温Jellium模型 | 相关能、结构因子（符号问题限制） |
| 团簇 | 原子团簇 | 热力学稳定性 |

### 1.6 符号问题的深入讨论

#### 符号问题的本质

费米子波函数的反对称性导致：

$$\Psi(\ldots, \mathbf{r}_i, \mathbf{r}_j, \ldots) = -\Psi(\ldots, \mathbf{r}_j, \mathbf{r}_i, \ldots)$$

在路径积分中，这表现为**负权重**：

$$Z = \int \mathcal{D}[\text{paths}] \, (-1)^{N_{ex}} e^{-S[\text{paths}]}$$

其中 $N_{ex}$ 是交换次数。

**"符号问题"**：正负权重相互抵消，导致信号噪声比指数下降：

$$\frac{\text{Signal}}{\text{Noise}} \sim e^{-\beta N \Delta E}$$

其中 $\Delta E$ 是最低激发能与基态能之差。

#### 符号问题的处理策略

| 方法 | 原理 | 适用场景 | 局限 |
|------|------|----------|------|
| **固定节点** | 限制波函数在试探节点内 | 量子化学 | 能量上界，受试探质量影响 |
| **约束路径** | 限制行列式重叠为正 | AFQMC | 同上 |
| **永久近似** | 用永久代替行列式 | 玻色子 | 不适用于费米子 |
| **解析延拓** | 从无符号问题参数外推 | 某些模型 | 外推误差 |
| **重新加权** | 以正权重为参考分布 | 弱符号问题 | 指数代价 |

### 1.6 量子化学QMC的特点总结

| 特点 | 说明 |
|------|------|
| **空间结构** | 连续 $\mathbb{R}^{3N}$ |
| **位置变量** | $\mathbf{R} = (\mathbf{r}_1, \ldots, \mathbf{r}_N)$，电子坐标 |
| **采样对象** | 位形（电子位置配置） |
| **试探波函数** | Slater-Jastrow 或神经网络形式 |
| **符号问题** | 费米子节点问题，需固定节点近似 |
| **典型系统** | 分子、固体电子结构 |

---

## 第二部分：凝聚态格点模型QMC（离散格点，无连续坐标）

### 2.1 物理背景

凝聚态格点模型QMC起源于**量子磁体和玻色-费米Hubbard模型**。基本对象是：

- **自旋** $\sigma_i \in \{\uparrow, \downarrow\}$ 或 **占据数** $n_i \in \{0, 1, 2, \ldots\}$
- **格点指标** $i \in \{1, \ldots, L\}$（离散标签，**不是空间坐标**）
- **哈密顿量** 如 Heisenberg 模型或 Bose-Hubbard 模型

**关键区别**：
- 格点 $i$ 是离散标签，没有空间距离的概念
- "世界线"是自旋/粒子在**虚时方向**的轨迹，不是实空间轨迹
- **无位置坐标 $\mathbf{R}$**

### 2.2 路径积分与世界线表示

#### Suzuki-Trotter 分解

配分函数通过 Trotter 分解：

$$Z = \text{Tr}(e^{-\beta H}) = \lim_{M\to\infty} \text{Tr}\left[\left(e^{-\beta H_0/M} e^{-\beta V/M}\right)^M\right]$$

虚时 $\tau \in [0, \beta]$ 被离散化为 $M$ 个时间片。

#### 世界线（Worldline）

对于自旋模型，自旋状态 $\sigma_i(\tau)$ 在虚时方向形成轨迹：

```
τ = β  ───────────────────────
       |  ↑  |  ↓  |  ↑  |  ↑  |
       |     |     |     |     |
       |  ↑  |  ↓  |  ↑  |  ↑  |
τ = 0  ───────────────────────
        i=1   i=2   i=3   i=4
```

**注意**：横轴是格点指标 $i$（离散），纵轴是虚时 $\tau$。"世界线"是沿虚时的轨迹，**不是粒子在空间移动的轨迹**。

### 2.3 随机序列展开（Stochastic Series Expansion, SSE）

#### 基本思想 [Sandvik 1999, 2019]

SSE 直接展开配分函数：

$$Z = \sum_{n=0}^\infty \sum_{\alpha} \sum_{S_n} \frac{\beta^n}{n!} \langle \alpha | H_{a_n} \cdots H_{a_1} | \alpha \rangle$$

其中算符序列 $S_n = (a_1, \ldots, a_n)$ 来自算符基 $\{H_a\}$。

**算符基构造**：对于Heisenberg模型 $H = J\sum_{\langle ij \rangle} \mathbf{S}_i \cdot \mathbf{S}_j$：

| 算符类型 | 形式 | 作用 |
|----------|------|------|
| $H_{ij}^{(0)}$ | $J S_i^z S_j^z$ | 对角，保持自旋态 |
| $H_{ij}^{(1)}$ | $\frac{J}{2}(S_i^+ S_j^- + S_i^- S_j^+)$ | 非对角，翻转自旋对 |

**截断**：配分函数截断到 $n \leq L_{max}$，其中 $L_{max}$ 自适应调整以覆盖重要构型。

**关键创新**：
- 无 Trotter 误差（直接展开，不离散化虚时）
- 算符循环（operator-loop）更新，克服临界慢化

#### SSE详细算法

```
Algorithm: SSE Operator-Loop Update

输入：自旋构型 |α⟩，算符序列 S_L = {H_a}
输出：更新后的 |α'⟩, S_L'

1. 扩展算符序列到长度 L
2. 构造算符-自旋图：
   - 对角算符：两条腿连接相同的自旋
   - 非对角算符：腿连接不同自旋，形成翻转边
3. 选择入口腿：
   - 从所有腿中随机选择
4. 构造循环：
   - 按有向规则遍历：进入顶点后选择出口腿
   - 规则：保持图连通性
5. 翻转循环上的算符：
   - 对角 H^{(0)} ↔ 非对角 H^{(1)}
   - 同时翻转相关自旋
6. 接受概率：p_accept = min(1, W_new/W_old)
```

**计算复杂度**：$O(L \cdot N)$ 其中 $L$ 是算符序列长度，$N$ 是格点数。临界点 $L \propto N$，总复杂度 $O(N^2)$。

#### 临界慢化与加速

**临界慢化**：接近临界点时，构型关联时间发散：

$$\tau \sim \xi^z$$

其中 $\xi$ 是关联长度，$z$ 是动力学临界指数。

**SSE的循环更新**：构建大范围循环，实现非局域更新：

| 更新类型 | 动力学指数 $z$ | 效率 |
|----------|----------------|------|
| 单自旋翻转 | $z \approx 2$ | 临界慢化严重 |
| 小循环 | $z \approx 1$ | 部分改善 |
| 大循环（SSE） | $z \approx 0$ | 克服临界慢化 |

### 2.4 有向循环算法（Directed Loop）[Syljuåsen & Sandvik 2002, Alet et al. 2005]

传统循环算法在某些参数区域出现**回溯（backtracking）**，效率下降。有向循环引入方向概念，设计转移概率减少回溯。

**回溯问题**：在强场或阻挫系统中，循环可能多次通过同一条边，降低更新效率。

**有向循环解决方案**：
- 引入边的方向标记
- 设计转移矩阵 $P(exit|entry)$ 使得回溯概率最小
- 最优转移概率由矩阵方程解出

### 2.5 虫子算法（Worm Algorithm）[Prokof'ev & Svistunov 1998, 2001]

#### 基本思想

扩展构型空间，添加"虫子"（一对开放端点）：

1. 从闭合构型插入虫子
2. 虫子头局域移动
3. 头尾相遇时虫子消失，产生新闭合构型

#### 优势

- 动力学临界指数接近零
- 可直接计算格林函数：虫子端点关联给出 $G(i, \tau; j, \tau')$

**虫子构型的物理意义**：
- 虫子头位置 $(i, \tau)$ 代表粒子产生
- 虫子尾位置 $(j, \tau')$ 代表粒子湮灭
- 虫子传播给出粒子格林函数

#### 虫子算法在玻色子系统的应用

对于Bose-Hubbard模型：

$$H = -t\sum_{\langle ij \rangle}(b_i^\dagger b_j + h.c.) + \frac{U}{2}\sum_i n_i(n_i-1) - \mu\sum_i n_i$$

虫子更新：
```
1. 随机选择格点 i 和虚时 τ
2. 插入虫子：|n⟩ → |n+1⟩，创建头在 (i, τ)
3. 虫子移动：沿世界线或跳跃到相邻格点
4. 头尾相遇：虫子湮灭，构型更新完成
```

### 2.6 格点模型符号问题的条件

#### Marshall符号规则 [Marshall 1955]

对于**双线性自旋模型**：

$$H = \sum_{\langle ij \rangle} J_{ij} \mathbf{S}_i \cdot \mathbf{S}_j$$

**定理**：若所有 $J_{ij} < 0$（铁磁耦合），则基态波函数在标准基下系数同号，无符号问题。

**证明**：作变换 $\sigma_i \to -\sigma_i$ 在一个子格点上，可将所有 $J_{ij}$ 变为正（反铁磁）。但基态仍是正定波函数。

#### 反铁磁Heisenberg模型的符号问题

对于**无阻挫反铁磁**：

$$H = J\sum_{\langle ij \rangle} \mathbf{S}_i \cdot \mathbf{S}_j, \quad J > 0$$

作子格变换 $S_i^x \to -S_i^x$, $S_i^y \to -S_i^y$ 在子格 $B$ 上，哈密顿量变为：

$$H' = -J\sum_{\langle ij \rangle} (S_i^z S_j^z - S_i^x S_j^x - S_i^y S_j^y)$$

**关键**：所有相互作用项变号，但可通过基矢变换使权重为正。

**无符号问题的条件**：
1. 双线性相互作用
2. 可二分格点（无奇数环）
3. 无阻挫（所有环满足乘积为正）

#### 阻挫导致符号问题

**阻挫定义**：格点存在奇数环，无法同时满足所有反铁磁键。

**例子**：三角格点反铁磁Heisenberg模型

$$H = J\sum_{\langle ij \rangle \in \triangle} \mathbf{S}_i \cdot \mathbf{S}_j$$

对于三自旋环，无法同时让所有自旋对反平行。每个三角形有6种基态构型，权重有正有负。

**符号问题强度**：平均符号 $\langle s \rangle$ 随系统尺寸指数衰减：

$$\langle s \rangle \sim e^{-c N}$$

其中 $c \propto \text{阻挫度}$。

#### 费米子Hubbard模型的符号问题

**Hubbard模型**：

$$H = -t\sum_{\langle ij \rangle, \sigma}(c_{i\sigma}^\dagger c_{j\sigma} + h.c.) + U\sum_i n_{i\uparrow} n_{i\downarrow}$$

**半满条件下的符号问题**：
- 无阻挫格点：通过变换可消除符号问题
- 阻挫格点：符号问题存在

**粒子-空穴对称性**：半满时，通过变换：

$$c_{i\uparrow} \to c_{i\uparrow}, \quad c_{i\downarrow} \to (-1)^i c_{i\downarrow}^\dagger$$

可将跳跃项变号，部分消除符号问题。

| 模型 | 格点类型 | 掺杂 | 符号问题 |
|------|----------|------|----------|
| Hubbard | 可二分 | 半满 | 无 |
| Hubbard | 可二分 | 掺杂 | 有 |
| Hubbard | 三角/阻挫 | 任意 | 有 |
| t-J模型 | 任意 | 任意 | 有 |

### 2.7 格点模型QMC的特点总结

| 特点 | 说明 |
|------|------|
| **空间结构** | 离散格点 $i \in \{1, \ldots, L\}$ |
| **位置变量** | **无连续坐标**，只有格点指标 |
| **采样对象** | 自旋构型、算符序列、世界线构型 |
| **"世界线"含义** | 虚时方向的轨迹，不是空间轨迹 |
| **Trotter误差** | SSE无误差；Worldline方法需离散化 |
| **符号问题** | 无阻挫系统无符号问题；阻挫系统有 |
| **典型系统** | 量子磁体、Bose-Hubbard、费米Hubbard |

---

## 第三部分：量子杂质物理QMC（杂质+浴耦合）

### 3.0 连续时间QMC的一般框架

#### 为什么需要连续时间？

**离散虚时的局限**：
- Trotter误差：$\Delta Z \propto (\Delta\tau)^2$
- 存储代价：需存储所有时间片状态
- 效率问题：稀疏事件（如自旋翻转）浪费大量存储

**连续时间的优势**：
- **无Trotter误差**：精确结果
- **稀疏表示**：只记录"事件"发生的时间点
- **效率提升**：事件密度低时计算量大幅减少

#### CT-QMC的核心思想

**配分函数展开**：以微扰展开替代时间离散化

$$Z = \text{Tr}(e^{-\beta H}) = \sum_{n=0}^\infty \frac{1}{n!} \int_0^\beta d\tau_1 \cdots \int_0^\beta d\tau_n \langle T_\tau H(\tau_1) \cdots H(\tau_n) \rangle$$

关键：积分变量 $\tau_i$ 是**连续的**，而非离散格点。

**事件表示**：配分函数由"事件"序列描述

$$Z = \sum_{\{\mathcal{C}\}} w(\mathcal{C})$$

其中 $\mathcal{C} = \{(\tau_1, op_1), (\tau_2, op_2), \ldots\}$ 是事件列表（时间+操作类型）。

#### CT-QMC与离散时间QMC的对比

| 特征 | 离散时间QMC | 连续时间QMC |
|------|-------------|-------------|
| 时间表示 | $\tau_m = m\Delta\tau$（格点） | $\tau \in [0, \beta]$（连续） |
| 存储对象 | 所有时间片状态 | 事件时间序列 $\{\tau_i\}$ |
| Trotter误差 | **存在** | **不存在** |
| 计算复杂度 | $O(M \cdot N)$ 每步 | $O(k)$ 每步，$k$ 是事件数 |
| 适用场景 | 连续空间粒子 | 格点模型、量子杂质 |
| 代表方法 | PIMC、Worldline | SSE、CT-HYB、CT-AUX |

**事件密度**：

$$\langle k \rangle \sim \beta \cdot \text{相互作用强度}$$

强耦合时事件多，计算量增加；弱耦合时事件少，效率高。

#### CT-QMC的采样策略

**Metropolis更新**：提议新构型 $\mathcal{C}'$，接受概率

$$p_{accept} = \min\left(1, \frac{w(\mathcal{C}')}{w(\mathcal{C})} \times \frac{q(\mathcal{C}|\mathcal{C}')}{q(\mathcal{C}'|\mathcal{C})}\right)$$

其中 $q$ 是提议概率。

**典型更新操作**：
1. **插入事件**：在随机时间 $\tau$ 插入新事件
2. **删除事件**：移除已有事件
3. **改变事件类型**：修改事件的操作类型
4. **移动事件时间**：调整事件发生时间

**测量**：事件关联给出物理量

$$\langle O \rangle = \sum_{\mathcal{C}} w(\mathcal{C}) O(\mathcal{C}) / Z$$

#### CT-QMC的三种主要变体

| 方法 | 展开方式 | 采样对象 | 典型应用 |
|------|----------|----------|----------|
| **CT-HYB** | 杂化项展开 | 杂化片段 $[\tau_{start}, \tau_{end}]$ | DMFT杂质求解 |
| **CT-AUX** | 辅助场展开 | 辅助场构型 $\{x(\tau)\}$ | 强关联杂质 |
| **CT-INT** | 相互作用项展开 | 相互作用事件 | Hubbard模型 |

**CT-HYB**（最常用）：

配分函数展开为杂化事件的路径积分：

$$Z = \sum_k \int_0^\beta d\tau_1 \cdots d\tau_k \, \text{Tr}\left[e^{-\beta H_{imp}} T_\tau H_{hyb}(\tau_1) \cdots H_{hyb}(\tau_k)\right]$$

每个杂化事件产生一个"片段"，占据虚时区间。

**CT-INT**：

展开相互作用项 $H_{int} = Un_\uparrow n_\downarrow$：

$$Z = \sum_k \frac{(-U)^k}{k!} \int_0^\beta d\tau_1 \cdots d\tau_k \, \langle n_\uparrow(\tau_1)n_\downarrow(\tau_1) \cdots \rangle$$

适用于Hubbard模型，但费米子行列式产生符号问题。

### 3.1 CT-HYB：杂化展开连续时间QMC

#### Anderson杂质模型

$$H = H_{imp} + H_{bath} + H_{hyb}$$

其中：
- $H_{imp} = \sum_{m} \epsilon_m n_m + \sum_{mm'} U_{mm'} n_m n_{m'}$：杂质局域能级+相互作用
- $H_{bath} = \sum_{p,m} \epsilon_{pm} c_{pm}^\dagger c_{pm}$：费米浴
- $H_{hyb} = \sum_{p,m} (V_{pm} c_{pm}^\dagger d_m + V_{pm}^* d_m^\dagger c_{pm})$：杂化耦合

**核心思想**：对杂化项 $H_{hyb}$ 进行微扰展开，浴自由度积分消去。

#### 杂化展开的数学形式

配分函数：

$$Z = \text{Tr}(e^{-\beta H}) = \text{Tr}\left(e^{-\beta(H_{imp}+H_{bath})} T_\tau e^{-\int_0^\beta H_{hyb}(\tau) d\tau}\right)$$

展开指数：

$$Z = \sum_{k=0}^\infty \int_0^\beta d\tau_1 \cdots \int_{\tau_{k-1}}^\beta d\tau_k \, \text{Tr}\left[e^{-\beta(H_{imp}+H_{bath})} T_\tau H_{hyb}(\tau_1) \cdots H_{hyb}(\tau_k)\right]$$

**关键步骤**：对浴求迹，得到仅含杂质自由度的表达式。

每个杂化事件 $H_{hyb}(\tau)$ 包含一对操作：
- $c_{pm}^\dagger d_m$：电子从浴跳到杂质（产生杂质电子）
- $d_m^\dagger c_{pm}$：电子从杂质跳回浴（湮灭杂质电子）

#### Segment表示

**物理图像**：每个轨道的杂质占据状态在虚时方向形成"片段"（segment）。

```
轨道 m 的占据历史：

τ = β  ───────────────────────────────────
       │     ████████          ██████████  │ ← 片段（杂质被占据）
       │                                      │ ← 空隙（杂质空）
τ = 0  ───────────────────────────────────
       ↑τ_start    ↑τ_end    ↑τ_start'   ↑τ_end'
```

**片段参数**：每个片段由时间区间 $[\tau_{start}, \tau_{end}]$ 定义。

**构型表示**：完整构型 $\mathcal{C} = \{(m, \tau_{start}, \tau_{end})_i\}$ 包含所有轨道的所有片段。

#### 权重计算

**杂质部分**：片段间的相互作用贡献权重。

$$w_{imp}(\mathcal{C}) = e^{-\sum_m \epsilon_m L_m} \times \prod_{m<m'} e^{-U_{mm'} d_{mm'}}$$

其中 $L_m$ 是轨道 $m$ 的总占据时长，$d_{mm'}$ 是轨道 $m$ 和 $m'$ 片段的重叠时长。

**杂化部分**：由浴的Green函数决定。

$$w_{hyb}(\mathcal{C}) = \det\left[\mathcal{G}_{ij}\right]$$

其中 $\mathcal{G}_{ij} = \Delta(\tau_i - \tau_j)$ 是杂化函数矩阵，$\Delta(\tau)$ 是Weiss函数。

**总权重**：

$$w(\mathcal{C}) = w_{imp}(\mathcal{C}) \cdot w_{hyb}(\mathcal{C})$$

#### CT-HYB算法流程

```
Algorithm: CT-HYB for Anderson Impurity Model

输入：Weiss函数 Δ(τ)，相互作用 U，温度 T
输出：杂质Green函数 G(τ)，自能 Σ(iω)

1. 初始化：空构型或随机构型

2. 更新操作（单个轨道）：
   a) 插入片段：
      - 随机选择轨道 m
      - 随机选择时间段 [τ_start, τ_end]
      - 计算新权重 w_new（含行列式更新）
      - 接受概率：p = min(1, w_new/w_old × β²Δτ²/(k+1))
      
   b) 删除片段：
      - 随机选择已有片段
      - 计算删除后的权重
      - 接受概率：p = min(1, w_new/w_old × k/β²Δτ²)
      
   c) 移动片段边界：
      - 调整 τ_start 或 τ_end
      - 局部权重更新

3. 全局更新（多轨道系统）：
   - 同时插入/删除多个轨道的片段
   - 处理轨道间相互作用约束

4. 测量：
   a) 杂质Green函数：G(τ) = ⟨d(τ)d†(0)⟩ 由片段端点计算
   b) 双占据：⟨n↑n↓⟩ 由片段重叠时长计算
   c) 其他观测量

5. 重复 2-4 直到收敛
```

#### 行列式更新技巧

**挑战**：每次插入/删除需要重新计算行列式，代价 $O(k^3)$。

**快速更新**：利用行列式的递推关系，代价降低到 $O(k^2)$。

插入新时间点 $\tau_{new}$：

$$\det(\mathcal{G}') = \det(\mathcal{G}) \times \left(\mathcal{G}_{new,new} - \sum_{i,j} \mathcal{G}_{new,i} (\mathcal{G}^{-1})_{i,j} \mathcal{G}_{j,new}\right)$$

**Lazy Skip Lists** [Sémon 2014]：进一步优化，达到 $O(k \log k)$。

#### CT-HYB的优势与局限

| 优势 | 局限 |
|------|------|
| 无Trotter误差（连续虚时） | 需计算行列式，强耦合时片段多 |
| 精确处理多轨道相互作用 | 解析延拓获得实频信息（ill-posed） |
| 无符号问题（玻色浴积分后） | 复杂相互作用时效率下降 |
| DMFT自洽循环高效 | 需要Weiss函数输入 |

**适用场景**：
- DMFT杂质求解（最主要应用）
- 多轨道强关联系统
- 密度矩阵嵌入理论（DMET）

### 3.2 CT-AUX：辅助场连续时间QMC

#### 核心思想

不同于CT-HYB展开杂化项，CT-AUX展开相互作用项，引入**离散辅助场**。

对于Hubbard相互作用 $Un_\uparrow n_\downarrow$，使用Hubbard-Stratonovich变换：

$$e^{-\Delta\tau U n_\uparrow n_\downarrow} = \frac{1}{2}\sum_{s=\pm 1} e^{\gamma s (n_\uparrow - n_\downarrow - 1/2)}$$

其中 $\gamma = \text{arccosh}(e^{\Delta\tau U/2})$，辅助场 $s \in \{+1, -1\}$ 是离散随机变量。

#### 连续时间辅助场

将离散时间推广到连续时间：辅助场在随机时间点 $\{\tau_i\}$ 取值 $s_i \in \{+1, -1\}$。

配分函数：

$$Z = \sum_{k=0}^\infty \frac{U^k}{k!} \int_0^\beta d\tau_1 \cdots d\tau_k \sum_{\{s_i\}} \text{Tr}\left[e^{-\beta H_0} T_\tau \prod_{i=1}^k V_{s_i}(\tau_i)\right]$$

其中 $V_s = e^{\gamma s (n_\uparrow - n_\downarrow - 1/2)}$。

#### CT-AUX vs CT-HYB

| 特征 | CT-AUX | CT-HYB |
|------|--------|--------|
| 展开对象 | 相互作用项 $H_{int}$ | 杂化项 $H_{hyb}$ |
| 辅助变量 | 离散场 $s_i \in \pm 1$ | 连续时间片段 |
| 权重形式 | Slater行列式 | 杂化行列式 |
| 符号问题 | 可能存在 | 通常无 |
| 适用模型 | Hubbard、强关联杂质 | Anderson杂质、DMFT |

**CT-AUX的优势**：
- 可处理任意相互作用形式（通过选择合适的辅助场分解）
- 与格点QMC的AFQMC方法类似，便于联合计算

**CT-AUX的局限**：
- 离散辅助场可能引入符号问题
- 多轨道系统辅助场组合复杂

### 3.3 CT-INT：相互作用展开连续时间QMC

#### 核心思想

直接展开相互作用项，不引入辅助场。

对于Hubbard模型 $H = H_0 + Un_\uparrow n_\downarrow$：

$$Z = \text{Tr}(e^{-\beta H_0} e^{-\beta U n_\uparrow n_\downarrow})$$

展开相互作用：

$$Z = \sum_{k=0}^\infty \frac{(-U)^k}{k!} \int_0^\beta d\tau_1 \cdots d\tau_k \, \text{Tr}\left[e^{-\beta H_0} T_\tau \prod_{i=1}^k n_\uparrow(\tau_i) n_\downarrow(\tau_i)\right]$$

#### 权重计算

对于自由部分 $H_0$，迹可以解析计算：

$$\text{Tr}\left[e^{-\beta H_0} T_\tau n_\uparrow(\tau_1) n_\downarrow(\tau_1) \cdots\right] = \det(\mathcal{G}_\uparrow) \times \det(\mathcal{G}_\downarrow)$$

其中 $\mathcal{G}_\sigma$ 是自旋 $\sigma$ 的自由Green函数矩阵。

**符号问题**：$(-U)^k$ 项在 $U > 0$ 时正负交替，导致符号问题！

$$w(\mathcal{C}) = (-U)^k \det(\mathcal{G}_\uparrow) \det(\mathcal{G}_\downarrow)$$

当 $k$ 增加时，权重符号震荡。

#### CT-INT的适用场景

| 模型 | 符号问题 | 可用性 |
|------|----------|--------|
| 吸引Hubbard（$U < 0$） | 无（$(-U)^k > 0$） | 可用 |
| 排斥Hubbard（$U > 0$） | 有（$(-U)^k$ 震荡） | 受限 |
| 半满排斥Hubbard | 可通过变换消除 | 可用 |
| 阻挫系统 | 有 | 受限 |

**应用**：
- 吸引Hubbard模型（超导研究）
- 半满无阻挫系统
- 作为基准测试方法

### 3.4 三种CT-QMC方法的综合比较

| 方法 | 展开项 | 采样对象 | 权重形式 | 符号问题 | 主要应用 |
|------|--------|----------|----------|----------|----------|
| **CT-HYB** | $H_{hyb}$ | 片段 $[\tau_s, \tau_e]$ | 杂化行列式 | 通常无 | DMFT杂质求解 |
| **CT-AUX** | $H_{int}$（HS变换） | 辅助场 $(s, \tau)$ | Slater行列式 | 可能存在 | 强关联杂质 |
| **CT-INT** | $H_{int}$（直接展开） | 相互作用事件 | 双行列式 | 存在（$U>0$） | 吸引Hubbard |

**选择建议**：

```
Anderson杂质模型 → CT-HYB（首选）
                  ↘ CT-AUX（复杂相互作用）

Hubbard格点模型 → CT-INT（吸引U）
                  ↘ CT-AUX（排斥U）
                  
DMFT自洽循环 → CT-HYB（标准）
```

### 3.5 模型家族与Spin-Boson模型

量子杂质物理研究**局域自由度与大环境（浴）的耦合**，是一类重要的多体问题。典型模型包括：

| 模型 | 杂质 | 浴 | 耦合形式 | 典型应用 |
|------|------|-----|----------|----------|
| **Spin-Boson** | 二能级系统 | 玻色浴 | $\sigma^z \sum_i \lambda_i (a_i + a_i^\dagger)$ | 量子退相干、耗散相变 |
| **Anderson杂质** | 局域电子 | 费米浴 | $c^\dagger d + h.c.$ | DMFT、重费米子 |
| **Rabi模型** | 二能级 | 单模腔场 | $g(a + a^\dagger)\sigma^x$ | Cavity QED、电路QED |
| **Dicke模型** | 多个二能级 | 单模腔场 | $g(a + a^\dagger)\sum_i \sigma_i^x$ | 超辐射相变 |
| **Jaynes-Cummings-Hubbard** | 格点自旋+腔 | 各格点腔场 | $g(\sigma_i^+ a_i + h.c.)$ | 极化子Mott绝缘体 |

**关键特征**：
- 杂质本身无空间扩展（单自旋/轨道或有限格点）
- 浴是**场变量**（玻色场或费米场），不是实空间坐标
- QMC采样在**虚时空间**进行，采样对象是杂质构型

### 3.6 Spin-Boson模型的连续虚时QMC

#### 模型定义

$$H = \frac{\Delta}{2}\sigma^x + \frac{\sigma^z}{2}\sum_i \lambda_i (a_i + a_i^\dagger) + \sum_i \omega_i a_i^\dagger a_i$$

浴的性质由谱函数描述：

$$J(\omega) = \pi \sum_i \lambda_i^2 \delta(\omega - \omega_i) = 2\pi \alpha \omega_c^{1-s} \omega^s$$

其中 $\alpha$ 是耦合强度，$\omega_c$ 是截止频率，$s$ 是谱指数：
- $s = 1$：Ohmic浴
- $s < 1$：Sub-Ohmic浴（可能导致局域化相变）
- $s > 1$：Super-Ohmic浴

#### 路径积分与浴积分消去

**关键技巧**：对玻色浴求迹，得到仅含杂质的有效作用量。

配分函数：

$$Z = \text{Tr}_\sigma \text{Tr}_B e^{-\beta H} = \sum_{\sigma(\tau)} e^{-S_{eff}[\sigma(\tau)]}$$

有效作用量：

$$S_{eff}[\sigma(\tau)] = -\int_0^\beta d\tau \int_0^\tau d\tau' \sigma(\tau) K_\beta(\tau-\tau') \sigma(\tau') + \frac{\Delta}{2}\int_0^\beta d\tau \sigma^x(\tau)$$

核函数：

$$K_\beta(\tau) = \int_0^\infty d\omega \frac{J(\omega)}{\pi} \frac{\cosh\left(\frac{\beta\omega}{2}[1-\frac{2\tau}{\beta}]\right)}{\sinh(\frac{\beta\omega}{2})}$$

渐近行为（对于 $\omega_c^{-1} \ll \tau \ll \beta$）：

$$K(\tau) \propto \frac{\alpha}{\tau^{1+s}}$$

**物理意义**：浴在虚时方向诱导**长程相互作用**，强度由 $\alpha$ 和 $s$ 决定。

#### 连续虚时聚类算法 [Winter et al. 2009, PRL 102, 030601]

**算法核心**：直接在连续虚时中采样自旋构型，避免Trotter离散化误差。

**采样对象**：自旋翻转时间序列 $\{\tau_1, \tau_2, \ldots, \tau_n\}$

**更新步骤**：

```
Algorithm: Continuous-Time Cluster Update for Spin-Boson

1. 表示：自旋构型 σ(τ) 由翻转时间 {τ_i} 完全确定
2. 插入：按 Poisson 分布（率 Γ）插入候选翻转点
3. 聚类构造：
   - 对于自旋取向相同的片段对 (I, II)
   - 以概率 p(s_I, s_II) 连接
   - p = 1 - exp(-2∫∫ dτ dτ' K_β(τ-τ'))
4. 聚类翻转：每个聚类独立以 1/2 概率翻转
5. 清理：移除未被翻转的候选点
```

**关键优势**：
- 无Trotter误差（连续虚时）
- 聚类更新克服临界慢化
- 精确处理长程虚时相互作用

#### 量子相变与临界行为

Sub-Ohmic Spin-Boson模型（$s < 1$）展现量子相变：

| 相 | 耦合强度 | 物理特征 |
|-----|----------|----------|
| **退局域相** | $\alpha < \alpha_c$ | 自旋在两态间隧穿 |
| **局域相** | $\alpha > \alpha_c$ | 自旋被冻结在某一态 |

**临界指数争议** [Winter 2009 vs NRG]：
- NRG预言非平均场指数
- QMC揭示**危险无关变量**效应，实际为平均场指数

$$y_t^* = \frac{1}{2}, \quad y_h^* = \frac{3}{4} \quad (\text{for } s < 1/2)$$

#### 路径积分与Berry相位 [Kirchner 2010]

自旋路径积分中存在**Berry相位**项：

$$S_{Berry} = i\pi \int_0^\beta d\tau \frac{d\phi}{d\tau} \cos\theta$$

Berry相位影响基态选择和临界行为，需在QMC中正确处理。

### 3.7 腔-QED系统的QMC方法

#### 从Spin-Boson到Cavity-QED的方法差异

| 特征 | Spin-Boson（连续浴） | Cavity-QED（离散模） |
|------|----------------------|---------------------|
| 浴谱 | 连续谱 $J(\omega)$ | 离散模 $\omega_1, \omega_2, \ldots$ |
| 光子数 | 无限（热浴） | 有限（需截断 $n_{max}$） |
| 路径积分 | 浴积分消去 → 有效作用量 | 浴积分消去 → **有限维矩阵** |
| QMC采样 | 自旋翻转时间 $\{\tau_i\}$ | 自旋翻转时间 + 光子数变化时间 |

**关键区别**：腔系统光子数有限，浴积分消去后不产生长程虚时相互作用，而是有限维问题。

#### Jaynes-Cummings模型的QMC处理

**模型**：

$$H = \omega_c a^\dagger a + \omega_a \sigma^z/2 + g(\sigma^+ a + \sigma^- a^\dagger)$$

**光子数截断**：设 $n_{max}$ 为最大光子数，Hilbert空间维度 $= 2(n_{max}+1)$。

**路径积分**：
- 对光子模求迹后，有效作用量包含有限个关联函数
- 可用精确对角化或小规模QMC处理

**QMC策略**：
1. **截断SSE**：将光子模视为有限玻色子格点，应用SSE算法
2. **混合方法**：光子部分精确处理，自旋部分路径积分

#### 多模腔与腔阵列的QMC

**Jaynes-Cummings-Hubbard模型**：

$$H = \sum_i \left[\omega_c a_i^\dagger a_i + g(\sigma_i^+ a_i + \sigma_i^- a_i^\dagger)\right] - J\sum_{\langle ij \rangle} a_i^\dagger a_j$$

**QMC方法选择**：

| 模型规模 | 推荐方法 |
|----------|----------|
| 单腔（1个JC单元） | 精确对角化（无需QMC） |
| 小阵列（≤8腔） | SSE（光子截断）+ DMRG |
| 大阵列 | SSE/Worm（将极化子视为格点玻色子） |

**光子截断准则**：$n_{max}$ 需覆盖基态和低激发态的典型光子数。

$$n_{max} \sim \frac{g^2}{\omega_c^2} + \text{thermal contribution}$$

#### Cavity-Bose-Hubbard模型的QMC

**模型**：

$$H = -t\sum_{\langle ij \rangle}(b_i^\dagger b_j + h.c.) + \frac{U}{2}\sum_i n_i(n_i-1) + g(a^\dagger + a)\sum_i(b_i + b_i^\dagger)$$

**采样对象**：
- 玻色子构型 $\{n_i(\tau)\}$：格点QMC（SSE/Worm）
- 腔场构型：可积分消去或视为额外格点

**方法**：
1. 将腔场视为第0号格点的玻色子，腔光子数截断
2. 应用标准Worm算法，测量腔-格点关联

#### 腔-QED QMC的符号问题

**单腔JC模型**：无符号问题（有限维，精确可解）

**腔阵列**：
- 无阻挫时：通常无符号问题（类比Bose-Hubbard）
- 强耦合极限：可能需检查腔-格点耦合是否引入负权重

**开放腔系统**：量子轨迹方法无符号问题（纯态演化）

### 3.8 Jaynes-Cummings-Hubbard模型的QMC方法

#### 模型定义

Jaynes-Cummings-Hubbard (JCH) 模型描述腔阵列中光-物质耦合：

$$H = \sum_i H_i^{JC} - t\sum_{\langle ij \rangle}(a_i^\dagger a_j + h.c.) + \sum_{\langle ij \rangle} V n_i n_j$$

其中局域JC哈密顿量：

$$H_i^{JC} = \omega_c a_i^\dagger a_i + \omega_a \sigma_i^+ \sigma_i^- + g(\sigma_i^+ a_i + \sigma_i^- a_i^\dagger)$$

**与Rabi-Hubbard模型的区别**：

| 模型 | 相互作用形式 | 守恒量 | 对称性 |
|------|--------------|--------|--------|
| **JCH** | $\sigma^+ a + \sigma^- a^\dagger$（旋转波近似） | 激发数 $N = N_\gamma + N_s$ | U(1) |
| **Rabi-Hubbard** | $(\sigma^+ + \sigma^-)(a + a^\dagger)$（含反旋转项） | 无 | $Z_2$ |

#### 方法一：双层映射 + Worm QMC [Wei et al. 2021, PRB 103, 045115]

**核心思想**：将JCH模型映射为双层Bose-Hubbard模型。

**映射方案** [Wei 2021]：

```
顶层：光子层（z=1）
      - 光子跳跃 t
      - 无排斥相互作用
      
底层：原子层（z=2）
      - 无跳跃
      - 原子激发间排斥 V
      
层间耦合：
      - 原子-光子耦合 g（跳跃）
```

**双层有效哈密顿量**：

$$H = \sum_i \left[\omega_c n_i^a + \omega_a n_i^\sigma + g(a_i^\dagger \sigma_i + a_i \sigma_i^\dagger)\right] - t\sum_{\langle ij \rangle}(a_i^\dagger a_j + h.c.) + V\sum_{\langle ij \rangle} n_i^\sigma n_j^\sigma$$

**Worm QMC算法** [Prokof'ev 1998]：

```
Algorithm: Worm QMC for JCH Model

1. 双层格点初始化：
   - 光子层：玻色子占据数 {n_i^a}
   - 原子层：自旋状态 {σ_i}

2. Worm更新：
   a) 插入虫子：在光子层或原子层创建头-尾对
   b) 虫子移动：
      - 光子层：沿跳跃边移动
      - 跨层：通过g耦合翻转自旋
   c) 头尾相遇：虫子湮灭

3. 测量 [Wei 2021]：
   - 光子密度 ρ_a = ⟨n_i^a⟩
   - 原子激发密度 ρ_σ = ⟨n_i^σ⟩
   - 超流刚度：ρ_s = L^{2-d} W^2 / (2dβt)
     （W是光子层的winding number）
   - 结构因子：S(Q) = ⟨ρ_Q ρ_Q†⟩ / N

4. 光子截断：检查结果对n_max收敛
```

**关键发现** [Wei 2021]：
- 双分格点（一维链、正方格点）：无稳定的超辐射固态相
- 三角格点：存在稳定的超辐射固态相

#### 方法二：Stochastic Green Function QMC [Flottat et al. 2016, EPJD]

**适用场景**：Rabi-Hubbard模型（含反旋转项，激发数不守恒）。

**SGF算法优势**：
- 可处理激发数不守恒的情况
- 可处理大光子数构型

**Flottat 2016的关键发现**：
- RH模型不存在Mott绝缘体相
- 相图只有：相干相 ↔ 非相干压缩相
- 深入相干相时，光子数发散

#### 方法三：SSE方法 [Hohenadler 2011, PRA]

**SSE算符基构造**：

| 算符类型 | 形式 | 作用 |
|----------|------|------|
| $H_{hop}^{(1)}$ | $-t(a_i^\dagger a_j + h.c.)$ | 光子跳跃 |
| $H_{JC}^{(1)}$ | $g(\sigma_i^+ a_i + h.c.)$ | 原子-光子耦合 |
| $H_{diag}^{(0)}$ | $\omega_c n_i^a + \omega_a \sigma_i^+ \sigma_i^-$ | 对角能量 |

**SSE更新**：循环更新需要处理光子数变化，确保在截断范围内。

**临界行为研究** [Hohenadler 2011]：
- 测量动力学临界指数 $z$
- JCH模型在超流-Mott转变处 $z = 1$

#### QMC可测量的物理量

| 观测量 | 公式 | 物理意义 |
|--------|------|----------|
| 光子数 | $\langle a^\dagger a \rangle$ | 腔场激发 |
| 原子激发 | $\langle \sigma^+ \sigma^- \rangle$ | 集体激发 |
| 超流刚度 | $\rho_s = L^{2-d} W^2 / (2d\beta t)$ | 相干相序参量 |
| 结构因子 | $S(Q) = \langle \rho_Q \rho_Q^\dagger \rangle$ | 固态序 |
| 光子Green函数 | $G_{a^\dagger a}(R)$ | 光子关联 |
| 光子凝聚分数 | $C_l = \sum_R G_{a^\dagger a}(R) / N_l$ | 相干序参量 |

#### 光子截断收敛性

**截断准则** [Wei 2021]：

$$n_{max} \geq \langle n^a \rangle + 3\sqrt{\langle (n^a)^2 \rangle - \langle n^a \rangle^2}$$

**收敛检验**：增大 $n_{max}$ 直到观测量变化小于统计误差。

#### 文献综述

| 文献 | 模型 | 方法 | 主要发现 |
|------|------|------|----------|
| **Wei 2021** [PRB 103, 045115] | 扩展JCH | Worm QMC | 三角格点存在超辐射固态 |
| **Flottat 2016** [EPJD] | Rabi-Hubbard | SGF QMC | 无Mott相，只有相干/非相干转变 |
| **Hohenadler 2011** [PRA 84, 041608] | JCH | SSE | 动力学临界指数 $z=1$ |

### 3.9 Tavis-Cummings模型与多原子腔系统

#### 从JCH到TC

| 模型 | 每腔原子数 | 原子间耦合 | 典型应用 |
|------|------------|------------|----------|
| **Jaynes-Cummings** | 1 | 无 | 单腔量子电动力学 |
| **Tavis-Cummings** | $N > 1$ | 通过腔场间接耦合 | 集体效应、超辐射 |
| **Jaynes-Cummings-Hubbard** | 1 | 腔间光子跳跃 | 腔阵列量子相变 |

**Tavis-Cummings模型**：单腔内$N$个原子与同一腔场耦合。

$$H_{TC} = \omega_c a^\dagger a + \frac{\omega_a}{2}\sum_{i=1}^{N} \sigma_i^z + g\sum_{i=1}^{N}(\sigma_i^+ a + \sigma_i^- a^\dagger)$$

**对称性**：总激发数守恒 $N_{ex} = a^\dagger a + \sum_i \sigma_i^+ \sigma_i^-$

#### TC模型的QMC方法

**基于JCH方法的推广**：将TC模型视为单腔的特例，应用Worm QMC或SSE方法。

**关键区别**：
- 单腔TC：有限维系统，可用精确对角化
- TC晶格（多腔）：应用JCH的QMC方法，每腔内含多个原子

**计算复杂度**：单腔Hilbert空间维度

$$\dim(\mathcal{H}) = 2^N \times (n_{max} + 1)$$

当 $N$ 较大时，精确对角化受限，QMC方法更高效。

#### 超辐射相变的QMC研究

**Dicke模型**：TC模型在$N \to \infty$极限的经典化，展现超辐射相变。

**QMC挑战**：
- 大 $N$ 极限：原子自由度多
- 热力学极限：需研究晶格系统

**方法选择**：

| 系统规模 | 推荐方法 |
|----------|----------|
| 单腔，$N \leq 5$ | 精确对角化 |
| 单腔，$N > 5$ | 截断SSE或量子轨迹 |
| TC晶格 | Worm QMC（推广JCH方法） |

### 3.10 非平衡与耗散QMC方法

#### 开放量子系统的Monte Carlo方法

**量子轨迹方法**（Quantum Trajectory / Monte Carlo Wave Function）：

$$|\psi(t+dt)\rangle = \frac{M_k |\psi(t)\rangle}{\|M_k |\psi(t)\rangle\|}$$

其中 $M_k$ 是量子跳跃算符。

**Hierarchy of Pure States (HOPS)** [Hartmann 2017, Suess 2014]：

利用层级结构展开密度矩阵，用随机态演化模拟开放系统。

#### 耗散相变的QMC研究 [Liu 2025, PRA]

驱动-耗散Kerr腔：

$$H = -\Delta a^\dagger a + \frac{U}{2}a^{\dagger 2}a^2 + F(a + a^\dagger)$$

耗散：$\dot{\rho} = -i[H, \rho] + \kappa\mathcal{D}[a]\rho$

QMC可研究：
- 非平衡稳态
- 多临界点
- 临界慢化与动力学指数

### 3.10 NRG与QMC方法对比

#### 数值重整化群（NRG）原理

NRG由Wilson为Kondo问题开发[Wilson 1975]，核心思想是将浴对数离散化，构造Wilson链：

$$H_{NRG} = H_{imp} + \sum_{n=0}^{N} \left[ \epsilon_n f_n^\dagger f_n + t_n (f_n^\dagger f_{n+1} + h.c.) \right]$$

沿链方向能量尺度按 $\Lambda^{-n}$ 衰减，实现多尺度重整化。

#### NRG的质量流误差 [Vojta et al. 2010, PRB 81, 075122]

**关键发现**：NRG迭代过程中，浴传播子的实部产生杂质参数重整化：

$$\text{Re}\,\Gamma_n(\omega=0) - \text{Re}\,\Gamma_\infty(\omega=0) \propto T_n^s$$

这导致**序参量质量的虚假温度依赖**，在量子临界点附近产生定性错误结果。

**对Spin-Boson模型的影响**：
- 标准NRG预言 $s < 1/2$ 时非平均场临界指数（错误）
- 修正NRG（考虑质量流）给出平均场指数 $\nu = 1/2$, $\eta = 0$（正确）

| 方法 | 优势 | 局限 | 适用场景 |
|------|------|------|----------|
| **NRG** | 能解析极低能标；RG流图直观 | 质量流误差；玻色子截断误差 | Kondo问题、谱函数 |
| **CT-QMC** | 无Trotter误差；精确处理长程相互作用 | 虚时到实频需解析延拓 | Spin-Boson临界行为、DMFT |
| **精确对角化** | 无随机误差 | Hilbert空间受限 | 小系统基准验证 |

### 3.11 开放量子系统的随机方法

#### 量子轨迹方法（Quantum Trajectory / Monte Carlo Wave Function）

开放系统Lindblad主方程：

$$\dot{\rho} = -i[H, \rho] + \sum_k \gamma_k \left( L_k \rho L_k^\dagger - \frac{1}{2}\{L_k^\dagger L_k, \rho\} \right)$$

**量子跳跃算法**：
1. 按概率 $dp_k = \gamma_k \langle \psi|L_k^\dagger L_k|\psi\rangle dt$ 执行跳跃 $|\psi'\rangle = L_k|\psi\rangle/\|L_k|\psi\rangle\|$
2. 否则执行非幺正演化 $|\psi'\rangle = e^{-iH_{eff}dt}|\psi\rangle$，其中 $H_{eff} = H - \frac{i}{2}\sum_k \gamma_k L_k^\dagger L_k$

#### 纯态层级方法（HOPS）[Suess et al. 2014, PRL; Hartmann et al. 2017, JTCC]

**核心思想**：将非马尔可夫量子态扩散（NMQSD）中的泛函导数层级化：

$$\psi_t^{(n)} = \int_0^t ds_1 \cdots \int_0^{s_{n-1}} ds_n \, \alpha(t-s_1)\cdots\alpha(t-s_n) \frac{\delta^n \psi_t}{\delta z_{s_1}^* \cdots \delta z_{s_n}^*}$$

**层级方程**（指数关联浴）：

$$\partial_t \psi_t^{(n)} = (-iH - n w + L z_t^*) \psi_t^{(n)} + \alpha(0) L \psi_t^{(n-1)} - L^\dagger \psi_t^{(n+1)}$$

**优势**：
- 处理纯态而非密度矩阵，维度更低
- 易于并行化（独立轨迹）
- 可系统检验收敛性（增加层级数）

#### 层级运动方程（HEOM）[Yan 2016, Frontiers Phys.]

HEOM直接演化密度矩阵层级：

$$\dot{\rho}_{\vec{n}} = -\left(i\mathcal{L} + \sum_k n_k \gamma_k\right) \rho_{\vec{n}} + \sum_k \mathcal{V} \rho_{\vec{n}+\vec{e}_k} + \sum_k n_k \mathcal{V}^\dagger \rho_{\vec{n}-\vec{e}_k}$$

其中 $\vec{n} = (n_1, n_2, \ldots)$ 是层级指标，$\mathcal{L}\rho = [H, \rho]$，$\mathcal{V}\rho = [L, \rho]$。

| 方法 | 态空间 | 温度处理 | 谱密度要求 |
|------|--------|----------|------------|
| **量子轨迹** | 纯态 | Lindblad（马尔可夫） | 无特定要求 |
| **HOPS** | 纯态层级 | 可处理非马尔可夫 | 指数展开 |
| **HEOM** | 密度矩阵层级 | 可处理非马尔可夫 | 指数展开 |
| **CT-QMC** | 虚时路径 | 平衡态 | 任意谱密度 |

### 3.12 DMFT与杂质求解器生态

#### DMFT自洽循环

动力学平均场理论（DMFT）将格点问题映射为自洽杂质问题：

$$G_{loc}(i\omega_n) = \sum_k \frac{1}{i\omega_n + \mu - \epsilon_k - \Sigma(i\omega_n)}$$

$$\mathcal{G}_0^{-1}(i\omega_n) = G_{loc}^{-1}(i\omega_n) + \Sigma(i\omega_n)$$

杂质求解器计算 $\Sigma(i\omega_n) = \mathcal{G}_0^{-1} - G_{imp}^{-1}$。

#### CT-HYB求解器软件生态

| 软件 | 语言 | 特点 | 引用 |
|------|------|------|------|
| **TRIQS/CTHYB** | C++/Python | 模块化，Python接口，支持复杂相互作用 | [Seth 2016, CPC] |
| **iQIST** | C++/Fortran | 开源工具包，多轨道支持，包含多种求解器 | [Huang 2015, CPC] |
| **ALPS/CT-HYB** | C++ | 与ALPS框架集成 | ALPS项目 |
| **Wien2k+DMFT** | 混合 | 第一性原理+DMFT联合计算 | 各种实现 |

#### 实时动力学方法

**问题**：虚时QMC需解析延拓获得实频信息，ill-posed问题。

**解决方案**：
1. **单自旋杂化展开** [Kubiczek 2019]：直接计算实时Green函数
2. **张量网络时间演化** [Thoenniss 2023]：MPS在时间域演化
3. **CT-AFQMC** [Gull 2008]：辅助场连续时间方法

### 3.13 量子杂质/腔系统QMC的特点总结

| 特点 | 说明 |
|------|------|
| **空间结构** | 杂质单点/有限格点 + 浴场（腔场） |
| **位置变量** | 杂质无空间扩展；浴是场变量 $a, a^\dagger$ |
| **采样对象** | 自旋翻转时间、杂化片段、光子数构型 |
| **虚时表示** | 连续虚时，无 Trotter 误差 |
| **浴积分消去** | 核心技巧，将浴影响编码为有效作用量 |
| **符号问题** | 量子杂质通常无符号问题；腔系统取决于模型 |
| **典型系统** | Spin-Boson, Anderson杂质, DMFT, Cavity-QED, JCH模型 |
| **关联方法** | NRG、HOPS、HEOM、量子轨迹 |

---

## 第四部分：格点场论QMC（高能物理背景）

### 4.1 物理背景

格点场论起源于**高能物理的QCD（量子色动力学）**。基本对象是：

- **规范场** $U_\mu(x)$（连接，Link variable）
- **场变量** $\phi(x)$ 或费米场
- **格点** 是连续场论的离散化，**物理意义与凝聚态格点不同**

**关键区别**：
- 格点QCD的格点是空间-时间离散化，格点间距 $a$ 有物理意义
- 凝聚态格点模型的格点是物理实体（原子位置），不需要连续极限

**QCD拉格朗日量**：

$$\mathcal{L}_{QCD} = -\frac{1}{4}F_{\mu\nu}^a F^{a\mu\nu} + \bar{\psi}(i\gamma^\mu D_\mu - m)\psi$$

其中 $F_{\mu\nu}^a$ 是规范场强，$D_\mu$ 是协变导数，$\psi$ 是夸克场。

**格点离散化**：

$$S_G[U] = \beta \sum_{\text{plaquettes}} \left(1 - \frac{1}{N_c}\text{Re Tr }U_{\square}\right)$$

其中 $U_{\square}$ 是四连接组成的plaquette，$\beta = 6/g^2$ 是耦合参数。

### 4.2 Hybrid Monte Carlo（HMC）

#### 基本思想

将场构型视为经典力学系统：

- 场变量 $\phi$ 对应"坐标"
- 引入共轭"动量" $\pi$
- 构造哈密顿量 $H = \frac{1}{2}\pi^2 + S[\phi]$

按哈密顿动力学演化，再 Metropolis 接受/拒绝。

#### HMC详细算法

```
Algorithm: Hybrid Monte Carlo

1. 初始化场构型 φ，采样动量 π ~ N(0, 1)
2. 分子动力学演化（Leapfrog积分）：
   for n = 1 to N_md do
     π ← π - (Δτ/2) ∂S/∂φ
     φ ← φ + Δτ π
     π ← π - (Δτ/2) ∂S/∂φ
   end for
3. Metropolis 检验：
   接受概率 p = min(1, exp(-ΔH))
4. 重复步骤 1-3
```

**接受率优化**：
- 选择步数 $N_{md}$ 和步长 $\Delta\tau$ 使得接受率 ~ 70-80%
- 自适应调整：根据接受率调整步长

**优势**：
- 全局更新：整个场构型同时改变
- 无临界慢化：大步长保持高接受率
- 精确性：Metropolis步骤保证正确平衡分布

### 4.3 费米子行列式问题

#### 费米子积分

对于费米场 $\psi$，积分给出行列式：

$$Z = \int \mathcal{D}U \det M[U] e^{-S_G[U]}$$

其中 $M[U]$ 是费米子矩阵（如Wilson费米子或staggered费米子）。

**问题**：
- $\det M[U]$ 可能是复数（符号问题）
- 计算行列式代价 $O(V^3)$，其中 $V$ 是格点体积

#### 解决方案

**1. 虚拟分子动力学（Pseudofermion）**：

引入辅助玻色场 $\phi$：

$$\det M^\dagger M = \int \mathcal{D}\phi \mathcal{D}\phi^* e^{-\phi^\dagger (M^\dagger M)^{-1} \phi}$$

作用量变为：

$$S_{eff} = S_G[U] + \phi^\dagger (M^\dagger M)^{-1} \phi$$

求导需要解线性方程组 $(M^\dagger M)x = \phi$。

**2. 多质量求解器**：

利用 $M^\dagger M$ 的稀疏性，用共轭梯度法求解。对于多个质量参数，可用多重质量共轭梯度（multi-mass CG）一次求解多个方程。

**3. 行列式估计**：

对于大系统，直接计算行列式不可行。使用：
- 随机估计：$\det M \approx \frac{1}{N}\sum_i \eta_i^\dagger M^{-1} \eta_i$
- 特征值截断：保留低模贡献

### 4.4 费米子符号问题

#### 符号问题的来源

对于有限化学势 $\mu$，费米子矩阵 $M[U; \mu]$ 非厄米：

$$M[U; \mu] \neq M^\dagger[U; \mu]$$

导致 $\det M$ 为复数：

$$\det M = |\det M| e^{i\theta}$$

**平均符号**：

$$\langle s \rangle = \frac{Z}{Z_{|s|}} = \langle e^{i\theta} \rangle_{|s|}$$

在QCD相变点附近，$\langle s \rangle \sim e^{-V f(\mu)}$ 指数衰减。

#### 处理方法

| 方法 | 原理 | 适用场景 | 局限 |
|------|------|----------|------|
| **重加权** | 以 $\mu=0$ 配分函数采样 | 小 $\mu$ | 指数代价 |
| **Taylor展开** | 展开 $\ln Z(\mu)$ 到有限阶 | 小 $\mu$ | 截断误差 |
| **解析延拓** | 从虚化学势延拓 | 中等 $\mu$ | ill-posed |
| **复Langevin** | 复数场演化 | 大 $\mu$ | 收敛性问题 |
| **Lefschetz thimble** | 在复流形上积分 | 理论研究 | 计算复杂 |

### 4.5 Heatbath 和 Overrelaxation [Kennedy 1985]

**Heatbath**：按局部条件概率更新单个连接。

对于 $SU(2)$ 规范场，连接 $U$ 的条件概率为：

$$P(U|\{U'\}) \propto \exp\left(\beta \text{Re Tr}(U J^\dagger)\right)$$

其中 $J$ 是邻接连接的乘积。

**Overrelaxation**：确定性地变换场构型，保持能量。

对于 $SU(2)$，作反射变换：

$$U \to U' = J^\dagger U J^\dagger / \det J$$

这保持作用量不变，但改变构型，加速遍历。

**组合使用**：
- Heatbath：遍历构型空间
- Overrelaxation：加速混合，减少自相关时间

### 4.6 连续极限外推

**格点物理量**：

$$O_{lattice}(a) = O_{continuum} + c_1 a + c_2 a^2 + \ldots$$

**外推步骤**：
1. 计算多个格点间距 $a_i$ 下的观测量
2. 拟合到多项式形式
3. 取 $a \to 0$ 极限

**改进作用量**：设计作用量使高阶修正项消失，加速收敛。

### 4.7 格点场论QMC的特点总结

| 特点 | 说明 |
|------|------|
| **空间结构** | 格点场论，离散化连续场 |
| **位置变量** | 场变量 $\phi(x)$，格点间距 $a$ 有物理意义 |
| **采样对象** | 规范场构型、费米场构型 |
| **算法** | HMC、Heatbath、Overrelaxation |
| **连续极限** | 需取 $a \to 0$，外推物理量 |
| **符号问题** | 费米子行列式可为负（问题） |
| **典型系统** | QCD、格点Higgs模型 |

---

## 第五部分：跨学科方法迁移的注意事项

### 5.1 "世界线"的含义差异

| 领域 | "世界线"含义 | 图示 |
|------|-------------|------|
| 量子化学（PIMC） | 粒子在实空间的路径 $\mathbf{R}(\tau)$ | 粒子轨迹（空间-虚时） |
| 格点模型 | 自旋在虚时的轨迹 $\sigma_i(\tau)$ | 格点-虚时图上的线 |
| 格点场论 | 场构型的时空路径 | 场值 $\phi(x,t)$ |

**常见误区**：将格点模型的"世界线"理解为粒子轨迹。

**正确理解**：
- 格点模型中，"世界线"是**单个格点**上自旋/占据数随虚时变化的记录
- 横轴是格点指标（离散），纵轴是虚时
- 世界线"穿过"某个格点表示该格点的自旋发生翻转

### 5.2 "格点"的物理意义差异

| 领域 | 格点意义 | 是否有连续极限 | 格点间距物理意义 |
|------|----------|----------------|------------------|
| 凝聚态格点模型 | 物理实体（原子位置） | 无 | 不存在（格点就是物理） |
| 格点场论 | 计算工具（离散化） | 有，需外推 | 有（截断效应） |

**常见误区**：在格点模型中引入"格点间距"概念。

**案例**：Bose-Hubbard模型的临界行为
- 凝聚态：临界指数是物理量，与格点间距无关
- 格点场论：需外推到连续极限得到物理临界指数

### 5.3 符号问题的处理差异

| 领域 | 符号问题来源 | 处理方法 | 能否精确解决 |
|------|-------------|----------|--------------|
| 量子化学 | 费米子节点 | 固定节点近似 | 近似（依赖试探波函数） |
| 格点模型 | 阻挫/费米子 | 特定对称性可避免 | 部分模型可精确解决 |
| 格点场论 | 费米子行列式 | 重加权、部分淬火 | 近似（有限化学势） |

**案例：符号问题在不同领域的表现**

1. **量子化学**：水分子（H₂O）的DMC计算
   - 固定节点近似：能量精度取决于Slater-Jastrow试探波函数的节点质量
   - 典型误差：0.1-1 eV（取决于试探波函数质量）

2. **凝聚态**：三角格点Heisenberg模型
   - 阻挫导致符号问题
   - 平均符号 $\langle s \rangle \sim e^{-cN}$，大系统无法处理

3. **格点QCD**：有限化学势QCD
   - 相图上半部分（高密度）无法用传统QMC访问
   - 需要其他方法（如复Langevin、Lefschetz thimble）

### 5.4 具体案例研究

#### 案例1：将DMC方法应用于格点模型？

**问题**：想用DMC计算Bose-Hubbard模型的基态能量。

**分析**：
- DMC是连续空间方法，采样电子坐标 $\mathbf{R}$
- Bose-Hubbard模型是格点模型，无连续坐标
- **不兼容**

**正确方法**：使用SSE或Worm算法。

#### 案例2：将SSE方法应用于分子体系？

**问题**：想用SSE计算水分子的基态能量。

**分析**：
- SSE要求离散格点指标 $i$
- 分子中电子位置是连续变量
- **需要先离散化空间**（但这引入误差）

**替代方案**：
- 使用VMC/DMC（连续空间方法）
- 或使用格点DFT（将分子放到格点上计算）

#### 案例3：CT-HYB与DMRG的比较

**问题**：Anderson杂质模型，应该用CT-HYB还是DMRG？

| 方法 | 优势 | 劣势 |
|------|------|------|
| CT-HYB | 无系统误差（连续虚时）；适应强关联 | 需解析延拓；有限温 |
| DMRG | 可达零温；实频信息直接 | 截断误差；一维效率最高 |

**决策**：
- 需要实频谱函数 → DMRG（或NRG）
- DMFT自洽循环 → CT-HYB
- 强关联极限 → 两者都可用，需比较

### 5.5 方法迁移建议

1. **明确学科背景**：阅读QMC文献时，首先判断作者的学科背景（量子化学/凝聚态/高能物理）

2. **检查空间结构**：方法是否隐含连续空间/离散格点假设？

3. **理解采样对象**：采样的是电子坐标、自旋构型，还是场构型？

4. **注意术语差异**：同一术语（世界线、格点）在不同领域含义不同

5. **验证符号问题**：检查目标系统是否有符号问题，及现有方法的处理能力

---

## 第六部分：方法选择指南

### 6.1 按物理系统选择QMC方法

```
物理系统 → 推荐QMC方法

分子/固体电子结构
  ├─ 基态能量 → VMC/DMC（固定节点）
  ├─ 激发态   → AFQMC或时间依赖DMC
  └─ 强关联   → AFQMC（约束路径）

量子磁体
  ├─ 无阻挫   → SSE/Directed Loop（无符号问题）
  ├─ 有阻挫   → DMRG或变分方法
  └─ 阻挫+费米子 → 需处理符号问题

玻色子系统
  ├─ Bose-Hubbard → SSE/Worm
  └─ 连续玻色子   → PIMC

量子杂质/DMFT
  ├─ Anderson杂质 → CT-HYB
  ├─ Spin-Boson   → 连续虚时聚类
  └─ 实时动力学   → HOPS/HEOM或CT-HYB实时间

格点QCD
  ├─ 零化学势   → HMC（无符号问题）
  └─ 有限化学势 → 重加权/Taylor展开/复Langevin
```

### 6.2 方法比较矩阵

| 方法 | 空间类型 | 符号问题 | 温度 | 实频信息 | 计算代价 |
|------|----------|----------|------|----------|----------|
| VMC | 连续 | 有（固定节点） | 基态 | 无 | 中等 |
| DMC | 连续 | 有（固定节点） | 基态 | 无 | 中等 |
| AFQMC | 连续/格点 | 有（约束路径） | 有限温 | 无 | 高 |
| SSE | 格点 | 条件性 | 有限温 | 无 | 中等 |
| Worm | 格点 | 无（玻色子） | 有限温 | 可通过解析延拓 | 中等 |
| CT-HYB | 杂质 | 通常无 | 有限温 | 需解析延拓 | 中等-高 |
| HMC | 场论 | 有（费米子） | 有限温 | 无 | 高 |
| DMRG | 格点/杂质 | 无 | 零温/有限温 | 直接 | 低-中等 |
| NRG | 杂质 | 无 | 有限温 | 直接 | 中等 |

### 6.3 计算复杂度对比

| 方法 | 时间复杂度 | 空间复杂度 | 并行化潜力 |
|------|------------|------------|------------|
| VMC | $O(N^3)$ | $O(N^2)$ | 高（独立样本） |
| DMC | $O(N^3 N_w)$ | $O(N^2 N_w)$ | 高 |
| AFQMC | $O(N^3 M N_w)$ | $O(N^2 N_w)$ | 中等 |
| SSE | $O(N^2 \beta)$ | $O(N \beta)$ | 低（串行更新） |
| CT-HYB | $O(N_{orb}^3 \beta)$ | $O(N_{orb}^2 \beta)$ | 中等 |
| HMC | $O(V M_{md})$ | $O(V)$ | 中等 |

其中 $N$：粒子数/电子数，$N_w$：行走者数，$M$：时间片数，$V$：格点体积，$N_{orb}$：轨道数，$M_{md}$：分子动力学步数。

### 6.4 软件生态概览

| 领域 | 软件 | 语言 | 特点 |
|------|------|------|------|
| **量子化学QMC** | QMCPACK | C++/CUDA | 大规模并行，GPU加速 |
| | PyQMC | Python | 集成PySCF，易扩展 |
| | CASINO | Fortran | 周期系统优化 |
| **格点模型QMC** | ALPS | C++ | 多种算法集成 |
| | SSE代码（各研究组） | Fortran/C++ | 定制化实现 |
| | Quantum TEA | Python | 张量网络+QMC |
| **量子杂质/DMFT** | TRIQS | C++/Python | 模块化DMFT框架 |
| | iQIST | C++/Fortran | 开源杂质求解器 |
| | ALPS/CT-HYB | C++ | 集成框架 |
| **格点场论** | MILC | C | QCD专用 |
| | Chroma | C++ | 高性能QCD |
| | Grid | C++ | 现代架构 |

---

## 参考文献

### 量子化学QMC

1. Kim, J. et al. (2018). QMCPACK: an open source ab initio quantum Monte Carlo package for the electronic structure of atoms, molecules and solids. *J. Phys.: Condens. Matter*, 30, 195901. [被引: 270]

2. Kent, P. R. C. et al. (2020). QMCPACK: Advances in the development, efficiency, and application of auxiliary field and real-space variational and diffusion quantum Monte Carlo. *J. Chem. Phys.*, 152, 174105. [被引: 136]

3. Wheeler, W. A. et al. (2023). PyQMC: an all-Python real-space quantum Monte Carlo module in PySCF. *J. Chem. Phys.*, 159, 234108. [被引: 19]

4. Barker, J. A. (1979). A quantum-statistical Monte Carlo method; path integrals with boundary conditions. *J. Chem. Phys.*, 70, 2914. [被引: 405]

### 格点模型QMC

5. Sandvik, A. W. (1999). Stochastic series expansion method with operator-loop update. *Phys. Rev. B*, 59, R14157. [被引: 753]

6. Sandvik, A. W. (2003). Stochastic series expansion method for quantum Ising models with arbitrary interactions. *Phys. Rev. E*, 68, 056701. [被引: 171]

7. Sandvik, A. W. (2019). Stochastic Series Expansion Methods. *arXiv:1909.10585*. [被引: 19]

8. Sandvik, A. W. (2010). Computational Studies of Quantum Spin Systems. *AIP Conf. Proc.*, 1297, 135. [被引: 570]

9. Melko, R. G. (2013). Stochastic Series Expansion Quantum Monte Carlo. *Springer Series in Solid-State Sciences*.

10. Syljuåsen, O. F. & Sandvik, A. W. (2002). Quantum Monte Carlo with directed loops. *Phys. Rev. E*, 66, 046701. [被引: 830]

11. Alet, F. et al. (2005). Generalized directed loop method for quantum Monte Carlo simulations. *Phys. Rev. E*, 71, 036706. [被引: 248]

12. Prokof'ev, N. V. & Svistunov, B. V. (1998). "Worm" algorithm in quantum Monte Carlo simulations. *Phys. Lett. A*, 241, 259. [被引: 272]

13. Prokof'ev, N. & Svistunov, B. (2001). Worm Algorithms for Classical Statistical Models. *Phys. Rev. Lett.*, 87, 160601. [被引: 369]

14. Blöte, H. W. J. et al. (2002). Cluster Monte Carlo simulation of the transverse Ising model. *Phys. Rev. E*, 66, 066701. [被引: 282]

#### DMRG方法（相关张量网络方法）

15. White, S. R. (1992). Density matrix formulation for quantum renormalization groups. *Phys. Rev. Lett.*, 69, 2863. [被引: 7586]

16. White, S. R. (1993). Density-matrix algorithms for quantum renormalization groups. *Phys. Rev. B*, 48, 10345. [被引: 3198]

17. Schollwöck, U. (2004). The density-matrix renormalization group. *Rev. Mod. Phys.*, 77, 259. [被引: 3280]

### 量子杂质QMC

18. Winter, A. et al. (2009). Quantum Phase Transition in the Sub-Ohmic Spin-Boson Model: Quantum Monte Carlo Study with a Continuous Imaginary Time Cluster Algorithm. *Phys. Rev. Lett.*, 102, 030601. [被引: 160]

19. Kirchner, S. (2010). Spin Path Integrals, Berry Phase, and the Quantum Phase Transition in the Sub-Ohmic Spin-Boson Model. *J. Low Temp. Phys.*, 161, 244. [被引: 23]

20. Sperstad, I. B. et al. (2012). Quantum criticality in spin chains with non-Ohmic dissipation. *Phys. Rev. B*, 85, 014504. [被引: 29]

21. Thoenniss, J. et al. (2023). Efficient method for quantum impurity problems out of equilibrium. *Phys. Rev. B*, 107, 155120. [被引: 58]

22. Casiano-Diaz, E. et al. (2023). A path integral ground state Monte Carlo algorithm for entanglement of lattice bosons. *SciPost Phys.*, 14, 054. [被引: 8]

#### CT-HYB方法与算法进展

23. Gull, E. et al. (2008). Continuous-time auxiliary field Monte Carlo for quantum impurity models. *EPL*, 82, 57003. [被引: 100+]

24. Hafermann, H. et al. (2013). Efficient implementation of the continuous-time hybridization expansion quantum impurity solver. *Comput. Phys. Commun.*, 184, 1280. [被引: 150+]

25. Shinaoka, H. et al. (2017). Continuous-time hybridization expansion quantum impurity solver for multi-orbital systems with complex hybridizations. *Comput. Phys. Commun.*, 215, 128. [被引: 80+]

26. Sémon, P. et al. (2014). Lazy skip lists: a new algorithm for fast hybridization-expansion quantum Monte Carlo. *Phys. Rev. B*, 90, 075149. [被引: 50+]

27. Kubiczek, P. et al. (2019). Exact real-time dynamics of single-impurity Anderson model from a single-spin hybridization-expansion. *SciPost Phys.*, 7, 016. [被引: 20+]

#### NRG方法

28. Vojta, M. et al. (2010). The mass-flow error in the Numerical Renormalization Group method and the critical behavior of the sub-Ohmic spin-boson model. *Phys. Rev. B*, 81, 075122. [被引: 160+]

29. Pižorn, I. et al. (2012). Variational Numerical Renormalization Group: Bridging the gap between NRG and Density Matrix Renormalization Group. *Phys. Rev. Lett.*, 108, 067202. [被引: 40+]

#### 软件实现

30. Seth, P. et al. (2016). TRIQS/CTHYB: A continuous-time quantum Monte Carlo hybridisation expansion solver for quantum impurity problems. *Comput. Phys. Commun.*, 200, 274. [被引: 200+]

31. Huang, L. et al. (2015). iQIST: An open source continuous-time quantum Monte Carlo impurity solver toolkit. *Comput. Phys. Commun.*, 195, 140. [被引: 80+]

#### 开放量子系统方法

32. Süß, D. et al. (2014). Hierarchy of Stochastic Pure States for Open Quantum System Dynamics. *Phys. Rev. Lett.*, 113, 150403. [被引: 211]

33. Hartmann, R. et al. (2017). Exact Open Quantum System Dynamics Using the Hierarchy of Pure States (HOPS). *J. Chem. Theory Comput.*, 13, 5618. [被引: 68]

34. Werner, A. H. et al. (2016). Positive Tensor Network Approach for Simulating Open Quantum Many-Body Systems. *Phys. Rev. Lett.*, 116, 237201. [被引: 181]

35. Yan, Y. et al. (2016). Dissipation equation of motion approach to open quantum systems. *Frontiers Phys.*, 11, 110306. [被引: 85]

36. Thoenniss, J. et al. (2023). Nonequilibrium quantum impurity problems via matrix-product states in the temporal domain. *Phys. Rev. B*, 107, 155121. [被引: 48]

### 腔-QED系统与耗散相变

37. Hwang, M.-J. et al. (2015). Quantum Phase Transition and Universal Dynamics in the Rabi Model. *Phys. Rev. Lett.*, 115, 180404. [被引: 440]

38. Liu, M. et al. (2017). Universal Scaling and Critical Exponents of the Anisotropic Quantum Rabi Model. *Phys. Rev. Lett.*, 119, 180601. [被引: 163]

39. Hwang, M.-J. et al. (2016). Quantum Phase Transition in the Finite Jaynes-Cummings Lattice Systems. *Phys. Rev. Lett.*, 116, 153601. [被引: 139]

40. Li, B.-W. et al. (2022). Observation of Non-Markovian Spin Dynamics in a Jaynes-Cummings-Hubbard Model Using a Trapped-Ion Quantum Simulator. *Phys. Rev. Lett.*, 128, 023601. [被引: 27]

41. De Filippis, G. et al. (2023). Signatures of Dissipation Driven Quantum Phase Transition in Rabi Model. *Phys. Rev. Lett.*, 130, 210404. [被引: 24]

42. Puebla, R. et al. (2017). Probing the Dynamics of a Superradiant Quantum Phase Transition with a Single Trapped Ion. *Phys. Rev. Lett.*, 118, 073601. [被引: 118]

43. Zheng, R.-H. et al. (2023). Observation of a Superradiant Phase Transition with Emergent Cat States. *Phys. Rev. Lett.*, 130, 263601. [被引: 49]

44. Nataf, P. et al. (2010). No-go theorem for superradiant quantum phase transitions in cavity QED and counter-example in circuit QED. *Nat. Commun.*, 1, 117. [被引: 337]

45. Liu, J.-Y. et al. (2025). Universal criticality of nonequilibrium quantum phase transition in a driven-dissipative Kerr cavity. *Phys. Rev. A*, 111, 013714. [被引: 2]

46. Wang, K. et al. (2026). Superradiant strongly correlated quantum states in cavity Hubbard model. *arXiv:2603.13657*.

47. Liang, Y. et al. (2026). Frustrated Rydberg Atom Arrays Meet Cavity QED: Emergence of the Superradiant Clock Phase. *Phys. Rev. Lett.*, 136, 113601. [被引: 1]

### 格点场论QMC

48. Joseph, A. (2020). Markov Chain Monte Carlo Methods in Quantum Field Theories. *SpringerBriefs in Physics*. [被引: 16]

49. Kennedy, A. D. et al. (1985). Improved heatbath method for Monte Carlo calculations in lattice gauge theories. *Phys. Lett. B*, 156, 393. [被引: 305]

50. Knechtli, F. (2017). Lattice Quantum Chromodynamics. *arXiv:1705.06990*.

51. Halliday, I. G. (1984). Lattice field theories. *Rep. Prog. Phys.*, 47, 891. [被引: 9]

### 综述与教材

52. Becca, F. & Sorella, S. (2017). Quantum Monte Carlo Approaches for Correlated Systems. *Cambridge University Press*. [被引: 506]

53. Benedict, K. A. (2019). Quantum Monte Carlo methods: algorithms for lattice models. *Contemp. Phys.*, 60, 83. [被引: 50]

54. Gubernatis, J. et al. (2016). Quantum Monte Carlo Methods. *Cambridge University Press*. [被引: 208]

---

## 附录：术语对照表

### A. 核心术语跨领域对照

| 术语 | 量子化学含义 | 格点模型含义 | 格点场论含义 |
|------|-------------|-------------|-------------|
| **格点** | 无格点概念 | 物理实体（原子） | 计算工具（离散化） |
| **世界线** | 粒子路径 $\mathbf{R}(\tau)$ | 自旋虚时轨迹 | 场构型时空路径 |
| **坐标** | 电子位置 $\mathbf{r}_i$ | 无连续坐标 | 场变量 $\phi(x)$ |
| **采样** | 位形 $\{\mathbf{R}\}$ | 自旋构型/算符序列 | 规范场构型 |
| **Trotter分解** | PIMC需要 | Worldline需要；SSE不需要 | HMC不需要 |
| **符号问题** | 费米子节点 | 阻挫/费米子 | 费米子行列式 |

### B. 方法缩写对照

| 缩写 | 全称 | 领域 |
|------|------|------|
| VMC | Variational Monte Carlo | 量子化学 |
| DMC | Diffusion Monte Carlo | 量子化学 |
| GFMC | Green's Function Monte Carlo | 量子化学 |
| AFQMC | Auxiliary Field QMC | 量子化学/凝聚态 |
| PIMC | Path Integral Monte Carlo | 量子化学/凝聚态 |
| SSE | Stochastic Series Expansion | 凝聚态格点模型 |
| CT-HYB | Continuous-Time Hybridization | 量子杂质/DMFT |
| CT-AUX | Continuous-Time Auxiliary Field | 量子杂质 |
| HMC | Hybrid Monte Carlo | 格点场论 |
| NRG | Numerical Renormalization Group | 量子杂质 |
| DMRG | Density Matrix Renormalization Group | 凝聚态格点模型 |
| DMFT | Dynamical Mean Field Theory | 强关联电子系统 |
| HOPS | Hierarchy of Pure States | 开放量子系统 |
| HEOM | Hierarchical Equations of Motion | 开放量子系统 |

### C. 哈密顿量符号约定

| 符号 | 含义 | 典型模型 |
|------|------|----------|
| $\sigma_i^{\alpha}$ | 自旋算符（$\alpha = x, y, z$） | Heisenberg, Ising |
| $S_i^{\alpha}$ | 自旋算符（自旋>1/2时用） | 高自旋模型 |
| $c_{i\sigma}^\dagger, c_{i\sigma}$ | 费米子产生/湮灭算符 | Hubbard, Anderson |
| $b_i^\dagger, b_i$ | 玻色子产生/湮灭算符 | Bose-Hubbard |
| $a^\dagger, a$ | 光子/声子产生/湮灭算符 | Rabi, Spin-Boson |
| $n_i$ | 占据数算符 | Hubbard, Bose-Hubbard |
| $U_\mu(x)$ | 规范场连接（Link） | 格点QCD |

### D. 路径积分表示

| 领域 | 路径积分变量 | 配分函数形式 |
|------|--------------|--------------|
| 量子化学 | 电子坐标 $\mathbf{R}(\tau)$ | $Z = \int \mathcal{D}\mathbf{R} \, e^{-S[\mathbf{R}]}$ |
| 格点模型 | 自旋构型 $\{\sigma_i(\tau)\}$ 或算符序列 | $Z = \sum_{\{S_n\}} \frac{\beta^n}{n!} \langle S_n \rangle$ |
| 量子杂质 | 自旋翻转时间 $\{\tau_i\}$ | $Z = \sum_{\{\tau_i\}} e^{-S_{eff}[\{\tau_i\}]}$ |
| 格点场论 | 场构型 $\phi(x)$ | $Z = \int \mathcal{D}\phi \, e^{-S[\phi]}$ |

### E. 符号问题判据

| 系统 | 无符号问题条件 | 原因 |
|------|----------------|------|
| 自旋模型（可二分格点） | 所有 $J_{ij} < 0$ 或所有 $J_{ij} > 0$（变换后） | Marshall规则 |
| Hubbard模型（可二分格点） | 半满 + 无阻挫 | 粒子-空穴对称性 |
| 玻色子系统 | 总是 | 权重为正 |
| 量子杂质（Spin-Boson） | 总是 | 浴积分后无符号问题 |
| QCD（零化学势） | $\mu = 0$ | $\det M = \det M^\dagger$ |

### F. 常见错误与纠正

| 错误 | 正确理解 |
|------|----------|
| "格点模型的格点间距$a$趋近于零" | 格点模型的格点是物理实体，没有"间距趋近于零"的概念 |
| "世界线是粒子在空间运动的轨迹" | 格点模型的世界线是虚时方向的轨迹，与空间运动无关 |
| "SSE有Trotter误差" | SSE直接展开配分函数，无Trotter误差；Worldline方法才有 |
| "量子杂质QMC有符号问题" | 通常无符号问题（浴积分消去后）；腔系统可能有 |
| "DMRG可以处理任何维度的系统" | DMRG对一维系统最高效；二维需要大面积，三维几乎不可行 |

---

## 致谢

本文综述基于以下开源软件和文献资源：
- ScholarAIO知识管理系统
- arXiv预印本服务器
- 各开源QMC软件项目（QMCPACK, TRIQS, iQIST, ALPS等）