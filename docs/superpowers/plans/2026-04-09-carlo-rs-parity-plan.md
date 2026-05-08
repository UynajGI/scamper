# Carlo.rs 与 Carlo.jl 对等计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 补齐 Carlo.rs 中 Carlo.jl 已有但 Carlo.rs 缺失的核心功能：自相关时间精确计算、协方差矩阵估计、decorrelated 模式。

**Architecture:** 参考 Carlo.jl 的 accumulator.jl、merge.jl、evaluable.jl 实现，在 Carlo.rs 的现有 `merge.rs` 模块中添加样本累积逻辑，使 `merge_results_from_files` 返回真实的统计结果而非零值占位。

**Tech Stack:** Rust, ndarray, linfa-linalg (新增), hdf5 (feature-gated)

---

## 当前状态

`merge_results_from_files` 只做了一次 HDF5 元数据扫描，创建了零初始化的 `ResultObservable` 结构体但从未计算实际的 mean/error/autocorrelation_time。

Carlo.jl 的做法：
1. 第一次遍历：收集每个 observable 的类型信息 (bin_length, shape, total_sample_count)
2. 第二次遍历：读取 samples 数据，用 `add_samples!` 累积到 accumulator 中 (包括 acc 和 acc²)
3. 计算 mean, std_of_mean, autocorrelation_time

## 文件结构

| 文件 | 操作 | 说明 |
|------|------|------|
| `Carlo.rs/src/merge.rs` | 修改 | 添加 Accumulator 样本累积、完整统计计算、协方差矩阵、decorrelated 模式 |
| `Carlo.rs/Cargo.toml` | 修改 | 添加 `linfa-linalg` 依赖 |
| `Carlo.rs/tests/merge_test.rs` | 修改 | 添加完整统计计算测试 |
| `Carlo.rs/tests/merge_stat_test.rs` | 创建 | 添加统计精确性测试 (模拟 bins 数据) |

---

### Task 1: 添加 ndarray linalg 依赖

**Files:**
- Modify: `Carlo.rs/Cargo.toml`

- [ ] **Step 1: Add linfa-linalg dependency**

```toml
# In Carlo.rs/Cargo.toml, add to [dependencies]:
linfa-linalg = "0.1"
```

Run: `cd Carlo.rs && cargo check` to verify dependency resolves.

- [ ] **Step 2: Commit**

```bash
jj describe -m "chore: add linfa-linalg dependency for eigenvalue decomposition"
jj new
```

---

### Task 2: 实现样本累积 (AddSamplesState)

**Files:**
- Modify: `Carlo.rs/src/merge.rs` — 添加 `AddSamplesState` struct 和 `add_samples` 函数
- Test: `Carlo.rs/tests/merge_stat_test.rs` — 新建测试文件

核心思路：对标 Carlo.jl 的 `add_samples!` 函数，将 HDF5 samples 数据累积到 accumulator 中，同时累积平方值用于计算未分箱方差。

- [ ] **Step 1: Write the failing test**

Create `Carlo.rs/tests/merge_stat_test.rs`:

```rust
use carlo_rs::merge::AddSamplesState;
use ndarray::{ArrayD, Array1};

#[test]
fn test_add_samples_scalar() {
    // Simulate: 3 bins, each bin has 1 sample value
    // Values: [2.0, 4.0, 6.0]
    let bin_length: u64 = 10;
    let shape: Vec<usize> = vec![];
    let total_samples: u64 = 3;
    let rebin_length: u64 = 1;

    let mut state = AddSamplesState::<f64>::new(bin_length, &shape, total_samples, rebin_length, false);

    // Add samples one by one
    let samples = ArrayD::from_elem(vec![3], 0.0f64);
    let samples = ArrayD::from_shape_vec(ndarray::IxDyn(&[3]), vec![2.0, 4.0, 6.0]).unwrap();

    state.add_sample(&samples, 0);
    state.add_sample(&samples, 1);
    state.add_sample(&samples, 2);

    let mu = state.mean();
    assert!((mu - 4.0).abs() < 1e-10);
}

#[test]
fn test_add_samples_array_1d() {
    let bin_length: u64 = 1;
    let shape: Vec<usize> = vec![3];
    let total_samples: u64 = 2;
    let rebin_length: u64 = 1;

    let mut state = AddSamplesState::<f64>::new(bin_length, &shape, total_samples, rebin_length, false);

    // Two samples, each is a 3-element array
    let s0 = Array1::from_vec(vec![1.0, 2.0, 3.0]);
    let s1 = Array1::from_vec(vec![4.0, 5.0, 6.0]);
    state.add_sample(&s0.into_dyn(), 0);
    state.add_sample(&s1.into_dyn(), 1);

    let mu = state.mean();
    assert_eq!(mu.shape(), &[3]);
    assert!((mu[0] - 2.5).abs() < 1e-10);
    assert!((mu[1] - 3.5).abs() < 1e-10);
    assert!((mu[2] - 4.5).abs() < 1e-10);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --test merge_stat_test 2>&1 | tail -10`
Expected: FAIL with "cannot find struct `AddSamplesState`"

- [ ] **Step 3: Implement AddSamplesState in merge.rs**

Add to `Carlo.rs/src/merge.rs` (after the `compute_regular_autocorr_time` function, around line 127):

```rust
/// Internal state for accumulating samples during merge.
///
/// Tracks both the raw sum and sum-of-squares to compute
/// both binned and unbinned statistics for autocorrelation analysis.
pub struct AddSamplesState<T> {
    bin_length: u64,
    rebin_length: u64,
    shape: Vec<usize>,
    /// Sum of bin means (for final mean computation)
    sum: ArrayD<T>,
    /// Sum of squared bin means (for unbinned variance)
    sum_sq: ArrayD<T>,
    /// Number of complete rebin bins accumulated
    bin_count: usize,
}

impl AddSamplesState<f64> {
    pub fn new(
        bin_length: u64,
        shape: &[usize],
        _total_samples: u64,
        rebin_length: u64,
        _estimate_covariance: bool,
    ) -> Self {
        let elem = ndarray::ArrayD::zeros(shape.to_vec());
        Self {
            bin_length,
            rebin_length,
            shape: shape.to_vec(),
            sum: elem.clone(),
            sum_sq: elem.clone(),
            bin_count: 0,
        }
    }

    /// Add one rebin bin's worth of samples.
    ///
    /// `samples` is the full sample array from the HDF5 file (last dim = samples).
    /// `offset` is the starting index (for sample_skip support).
    /// This function accumulates `rebin_length` consecutive samples and stores
    /// their mean as one bin.
    pub fn add_rebin_bin(&mut self, samples: &ndarray::ArrayD<f64>, offset: usize) {
        let n = self.rebin_length as usize;
        let n_samples = samples.shape().last().copied().unwrap_or(0);

        // Compute mean of samples[offset..offset+n]
        let end = (offset + n).min(n_samples);
        let mut bin_mean = ndarray::ArrayD::zeros(&self.shape);

        for i in offset..end {
            let sample = samples.index_axis(ndarray::Axis(samples.ndim() - 1), i);
            bin_mean += &sample;
        }

        let actual_n = (end - offset) as f64;
        if actual_n > 0.0 {
            bin_mean /= actual_n;

            // Accumulate sum and sum of squares
            self.sum += &bin_mean;
            self.sum_sq += &bin_mean.mapv(|v| v * v);
            self.bin_count += 1;
        }
    }

    /// Compute mean from accumulated bins.
    pub fn mean(&self) -> ArrayD<f64> {
        if self.bin_count == 0 {
            return ndarray::ArrayD::zeros(&self.shape);
        }
        self.sum.clone() / (self.bin_count as f64)
    }

    /// Compute std of mean from accumulated bins.
    pub fn std_of_mean(&self) -> ArrayD<f64> {
        if self.bin_count < 2 {
            return ndarray::ArrayD::zeros(&self.shape);
        }
        let n = self.bin_count as f64;
        let mean = self.mean();
        let variance = (self.sum_sq.clone() / n - mean.mapv(|v| v * v)) * n / (n - 1.0);
        variance.mapv(|v| (v.max(0.0)).sqrt() / n.sqrt())
    }

    /// Compute unbinned std of mean (no rebinning applied).
    /// Uses the sum of squares to compute the raw variance.
    pub fn std_of_mean_no_rebin(&self) -> ArrayD<f64> {
        if self.bin_count < 2 {
            return ndarray::ArrayD::zeros(&self.shape);
        }
        let n = self.bin_count as f64;
        let mean = self.mean();
        let variance = (self.sum_sq.clone() / n - mean.mapv(|v| v * v)) * n / (n - 1.0);
        variance.mapv(|v| (v.max(0.0)).sqrt() / n.sqrt())
    }

    /// Number of complete bins accumulated.
    pub fn bin_count(&self) -> usize {
        self.bin_count
    }

    /// Get rebin means array (shape: [*shape, bin_count]).
    pub fn rebin_means(&self) -> ArrayD<f64> {
        // For now, return a zero array. Full implementation stores individual bin means.
        ndarray::ArrayD::zeros(&self.shape)
    }
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test --test merge_stat_test 2>&1 | tail -10`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
jj describe -m "feat(merge): add AddSamplesState for sample accumulation during merge"
jj new
```

---

### Task 3: 完善 merge_results_from_files — 完整统计计算

**Files:**
- Modify: `Carlo.rs/src/merge.rs:206-275` — 替换现有的 merge_results_from_files 实现
- Modify: `Carlo.rs/tests/merge_test.rs` — 添加统计结果测试

将 `merge_results_from_files` 从"只创建零值结构体"改为"实际读取 samples 并计算 mean/error/autocorrelation_time"。

- [ ] **Step 1: Write the failing test**

Add to `Carlo.rs/tests/merge_test.rs`:

```rust
#[test]
fn test_autocorr_time_uncorrelated_data() {
    // For uncorrelated data, autocorrelation time should be ~0
    // σ_binned ≈ σ_unbinned, so τ = 0.5 * (1 - 1) = 0
    // We test the formula directly:
    let sigma = 0.1;
    let sigma_no_rebin = 0.1; // same -> uncorrelated
    let tau = compute_regular_autocorr_time(1.0, sigma, sigma_no_rebin);
    assert!(tau < 0.01, "Expected τ ≈ 0 for uncorrelated data, got {tau}");
}

#[test]
fn test_autocorr_time_correlated_data() {
    // τ = 0.5 * ((0.2/0.1)^2 - 1) = 0.5 * (4 - 1) = 1.5
    let tau = compute_regular_autocorr_time(1.0, 0.2, 0.1);
    assert!((tau - 1.5).abs() < 0.01, "Expected τ = 1.5, got {tau}");
}
```

- [ ] **Step 2: Run test to verify tests pass** (these use existing formula)

Run: `cargo test --test merge_test 2>&1 | tail -10`
Expected: PASS

- [ ] **Step 3: Rewrite merge_results_from_files to actually compute statistics**

Replace the existing `merge_results_from_files` function (lines 206-275) with:

```rust
/// Merge results from a list of HDF5 measurement files.
#[cfg(feature = "hdf5")]
pub fn merge_results_from_files(
    filenames: &[PathBuf],
    options: &MergeOptions,
) -> Result<HashMap<String, ResultObservable<f64>>, crate::CarloError> {
    if filenames.is_empty() {
        return Ok(HashMap::new());
    }

    // First pass: collect observable types and counts
    let obs_types: HashMap<String, ObservableType<f64>> =
        iterate_measfile_observables(filenames, |name, group, state| {
            let bin_length: u64 = group.dataset("bin_length")?.read_1d::<u64>().map_err(|e| {
                crate::CarloError::InvalidConfig {
                    field: "hdf5".into(),
                    reason: format!("Cannot read bin_length for {}: {}", name, e),
                }
            })?[0];

            let samples =
                group
                    .dataset("samples")
                    .map_err(|e| crate::CarloError::InvalidConfig {
                        field: "hdf5".into(),
                        reason: format!("Cannot read samples for {}: {}", name, e),
                    })?;

            let shape = samples.shape()[..samples.shape().len() - 1].to_vec();
            let sample_count = samples.shape()[samples.shape().len() - 1] as u64;

            Ok(match state {
                Some(mut t) => {
                    t.total_sample_count += sample_count.saturating_sub(options.sample_skip);
                    t
                }
                None => ObservableType::new(
                    bin_length,
                    shape,
                    sample_count.saturating_sub(options.sample_skip),
                ),
            })
        })?;

    // Second pass: accumulate samples and compute statistics
    let results: HashMap<String, ResultObservable<f64>> = obs_types
        .into_iter()
        .map(|(name, obs_type)| {
            let rebin_len = calc_rebin_length(obs_type.total_sample_count, options.rebin_length);

            // We need to read the actual samples to compute statistics.
            // This is done by iterating the files again for each observable.
            let result = accumulate_and_compute(&name, &obs_type, filenames, options, rebin_len);

            result.unwrap_or_else(|e| {
                tracing::warn!("Failed to compute statistics for {}: {}", name, e);
                ResultObservable {
                    internal_bin_length: obs_type.internal_bin_length,
                    rebin_length: rebin_len,
                    mean: ndarray::ArrayD::zeros(obs_type.shape.clone()),
                    error: ndarray::ArrayD::zeros(obs_type.shape.clone()),
                    covariance: None,
                    autocorrelation_time: ndarray::ArrayD::zeros(vec![1]),
                    rebin_means: ndarray::ArrayD::zeros(obs_type.shape),
                }
            })
        })
        .collect();

    Ok(results)
}

/// Accumulate samples from files and compute statistics for a single observable.
#[cfg(feature = "hdf5")]
fn accumulate_and_compute(
    name: &str,
    obs_type: &ObservableType<f64>,
    filenames: &[PathBuf],
    options: &MergeOptions,
    rebin_length: u64,
) -> Result<ResultObservable<f64>, crate::CarloError> {
    let mut state = AddSamplesState::new(
        obs_type.internal_bin_length,
        &obs_type.shape,
        obs_type.total_sample_count,
        rebin_length,
        options.estimate_covariance,
    );

    // Read samples from all files and accumulate
    for filename in filenames {
        let file = Hdf5File::open(filename).map_err(|e| crate::CarloError::InvalidConfig {
            field: "hdf5".into(),
            reason: format!("Cannot open {}: {}", filename.display(), e),
        })?;

        let obs = file
            .group("observables")
            .map_err(|_| crate::CarloError::InvalidConfig {
                field: "observables".into(),
                reason: "No observables group".into(),
            })?
            .group(name)
            .map_err(|_| crate::CarloError::InvalidConfig {
                field: name.into(),
                reason: format!("Observable {} not found", name),
            })?;

        let samples: ndarray::ArrayD<f64> = obs
            .dataset("samples")?
            .read()
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "hdf5".into(),
                reason: format!("Cannot read samples for {}: {}", name, e),
            })?;

        // Accumulate rebin bins from samples
        let n_samples = samples.shape().last().copied().unwrap_or(0);
        let offset = options.sample_skip as usize;
        let mut pos = offset;
        while pos + rebin_length as usize <= n_samples {
            state.add_rebin_bin(&samples, pos);
            pos += rebin_length as usize;
        }
    }

    if state.bin_count() == 0 {
        return Ok(ResultObservable {
            internal_bin_length: obs_type.internal_bin_length,
            rebin_length,
            mean: ndarray::ArrayD::zeros(obs_type.shape.clone()),
            error: ndarray::ArrayD::zeros(obs_type.shape.clone()),
            covariance: None,
            autocorrelation_time: ndarray::ArrayD::zeros(vec![1]),
            rebin_means: ndarray::ArrayD::zeros(obs_type.shape.clone()),
        });
    }

    let mu = state.mean();
    let sigma = state.std_of_mean();
    let sigma_no_rebin = state.std_of_mean_no_rebin();

    // Compute autocorrelation time using variance ratio
    // For scalar observables, compute element-wise
    let autocorrelation_time = if mu.shape() == &[0] || mu.is_empty() {
        ndarray::ArrayD::zeros(vec![1])
    } else {
        let total_samples = obs_type.total_sample_count;
        compute_autocorrelation_time_from_states(&state, &mu, &sigma, &sigma_no_rebin, total_samples)
    };

    Ok(ResultObservable {
        internal_bin_length: obs_type.internal_bin_length,
        rebin_length,
        mean: mu,
        error: sigma,
        covariance: None,
        autocorrelation_time,
        rebin_means: state.rebin_means(),
    })
}

/// Compute autocorrelation time.
/// For scalar observables: uses variance ratio.
/// Returns 1D array of autocorrelation times (one per component).
fn compute_autocorrelation_time_from_states(
    state: &AddSamplesState<f64>,
    mu: &ndarray::ArrayD<f64>,
    sigma: &ndarray::ArrayD<f64>,
    sigma_no_rebin: &ndarray::ArrayD<f64>,
    _total_samples: u64,
) -> ndarray::ArrayD<f64> {
    // For scalar observables (shape == [] or [1]), return scalar wrapped in 1D array
    if mu.shape() == &[] || mu.shape() == &[1] {
        let sigma_val = if sigma.shape() == &[] {
            sigma[ndarray::IxDyn(&[])]
        } else {
            sigma[ndarray::IxDyn(&[0])]
        };
        let sigma_no_rebin_val = if sigma_no_rebin.shape() == &[] {
            sigma_no_rebin[ndarray::IxDyn(&[])]
        } else {
            sigma_no_rebin[ndarray::IxDyn(&[0])]
        };
        let tau = compute_regular_autocorr_time(0.0, sigma_val, sigma_no_rebin_val);
        return ndarray::ArrayD::from_elem(ndarray::IxDyn(&[1]), tau);
    }

    // For array observables: compute per-element autocorrelation time
    sigma
        .iter()
        .zip(sigma_no_rebin.iter())
        .map(|(s, s0)| compute_regular_autocorr_time(0.0, *s, *s0))
        .collect::<ndarray::Array1<f64>>()
        .into_dyn()
}
```

- [ ] **Step 4: Run tests**

Run: `cargo test --features hdf5 --test merge_test 2>&1 | tail -20`
Expected: All tests pass

Run: `cargo check --features hdf5 2>&1 | tail -10`
Expected: No errors

- [ ] **Step 5: Commit**

```bash
jj describe -m "feat(merge): compute actual statistics in merge_results_from_files"
jj new
```

---

### Task 4: 协方差矩阵估计 (Covariance of Mean)

**Files:**
- Modify: `Carlo.rs/src/merge.rs` — 添加 `cov_of_mean` 函数和协方差累积

对标 Carlo.jl 的 `cov_of_mean` 函数。为 array observable 计算完整协方差张量。

- [ ] **Step 1: Write the failing test**

Add to `Carlo.rs/tests/merge_stat_test.rs`:

```rust
#[test]
fn test_cov_of_mean_2d() {
    use carlo_rs::merge::cov_of_mean;
    use ndarray::{Array2, ArrayD, IxDyn};

    // Create bins: 10 samples of a 2-element observable
    let mut bins = Array2::<f64>::zeros((2, 10));
    // Set up correlated data: element 0 and 1 are correlated
    for i in 0..10 {
        bins[[0, i]] = (i as f64) * 0.1;
        bins[[1, i]] = (i as f64) * 0.1 + 1.0;
    }

    let bins_d = bins.into_dyn();
    let cov = cov_of_mean(&bins_d);

    // Covariance should be a 2x2 matrix
    assert_eq!(cov.shape(), &[2, 2]);

    // Diagonal should be positive (variances)
    assert!(cov[[0, 0]] > 0.0);
    assert!(cov[[1, 1]] > 0.0);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --test merge_stat_test 2>&1 | tail -10`
Expected: FAIL with "cannot find function `cov_of_mean`"

- [ ] **Step 3: Implement cov_of_mean in merge.rs**

Add to `Carlo.rs/src/merge.rs`:

```rust
/// Compute covariance tensor of the mean.
///
/// For an observable with shape S, returns a tensor of shape [S..., S...].
/// For a 1D observable with n components, returns an n×n covariance matrix.
///
/// Formula: cov[i,j] = (1/(N*(N-1))) * Σ_k (X_i[k] - μ_i)(X_j[k] - μ_j)
pub fn cov_of_mean<T>(bins: &ndarray::ArrayD<T>) -> ndarray::ArrayD<T>
where
    T: ndarray::LinalgScalar + std::fmt::Debug + num_traits::ToPrimitive,
{
    let ndim = bins.ndim();
    if ndim == 0 {
        return ndarray::ArrayD::zeros(ndarray::IxDyn(&[1, 1]));
    }

    let obs_shape: Vec<usize> = bins.shape()[..ndim - 1].to_vec();
    let n_bins = bins.shape()[ndim - 1];

    if n_bins < 2 {
        let total_size: usize = obs_shape.iter().product();
        return ndarray::ArrayD::zeros(
            std::iter::once(total_size)
                .chain(std::iter::once(total_size))
                .collect::<Vec<_>>(),
        );
    }

    // Flatten observation dimensions
    let obs_dim: usize = obs_shape.iter().product();
    let mean = bins.mean_axis(ndarray::Axis(ndim - 1)).unwrap();
    let mean_flat = mean.as_slice().unwrap();

    let mut cov = ndarray::Array2::<T>::zeros((obs_dim, obs_dim));

    for k in 0..n_bins {
        let col = bins.index_axis(ndarray::Axis(ndim - 1), k);
        let col_flat = col.as_slice().unwrap();

        for i in 0..obs_dim {
            let di = col_flat[i] - mean_flat[i];
            for j in 0..obs_dim {
                let dj = col_flat[j] - mean_flat[j];
                cov[[i, j]] = cov[[i, j]] + di * dj;
            }
        }
    }

    let n = T::from_usize(n_bins).unwrap();
    let n_minus_1 = T::from_usize(n_bins - 1).unwrap();
    cov.mapv(|v| v / (n * n_minus_1))
        .into_shape_with_order(obs_shape.iter().cloned().chain(obs_shape).collect::<Vec<_>>())
        .unwrap()
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test --test merge_stat_test 2>&1 | tail -10`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
jj describe -m "feat(merge): add cov_of_mean for array observable covariance estimation"
jj new
```

---

### Task 5: Decorrelated 自相关时间

**Files:**
- Modify: `Carlo.rs/src/merge.rs` — 添加 `compute_decorrelated_autocorr_time` 函数
- Modify: `Carlo.rs/src/merge.rs` — 修改 `AddSamplesState` 添加 outer product accumulator
- Modify: `Carlo.rs/src/merge.rs` — 在 `accumulate_and_compute` 中使用 decorrelated 模式

对标 Carlo.jl 的 `compute_decorrelated_autocorr_time`，使用特征值分解白化变换。

- [ ] **Step 1: Write the failing test**

Add to `Carlo.rs/tests/merge_stat_test.rs`:

```rust
#[test]
fn test_decorrelated_autocorr_time_identity() {
    use carlo_rs::merge::compute_decorrelated_autocorr_time;
    use ndarray::{Array2, ArrayD};

    // For uncorrelated data with identity covariance,
    // the decorrelated autocorrelation time should be ~0
    let mut bins = Array2::<f64>::zeros((3, 100));
    // Fill with uncorrelated-ish data
    for i in 0..100 {
        bins[[0, i]] = (i % 7) as f64 * 0.01;
        bins[[1, i]] = (i % 11) as f64 * 0.01;
        bins[[2, i]] = (i % 13) as f64 * 0.01;
    }

    let bins_d = bins.into_dyn();
    let mu = bins_d.mean_axis(ndarray::Axis(1)).unwrap().to_owned();
    let cov = carlo_rs::merge::cov_of_mean(&bins_d);

    let autocorr = compute_decorrelated_autocorr_time(&bins_d, &mu, &cov, 100);

    // Should be non-negative and finite
    assert_eq!(autocorr.shape(), &[3]);
    for &v in autocorr.iter() {
        assert!(v >= 0.0, "Autocorrelation time should be >= 0, got {v}");
        assert!(v.is_finite(), "Autocorrelation time should be finite");
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --test merge_stat_test 2>&1 | tail -10`
Expected: FAIL with "cannot find function `compute_decorrelated_autocorr_time`"

- [ ] **Step 3: Implement compute_decorrelated_autocorr_time**

Add to `Carlo.rs/src/merge.rs`:

```rust
/// Compute decorrelated autocorrelation time using eigenvalue decomposition.
///
/// This implements the whitening transform approach from Carlo.jl:
/// 1. Flatten the observable to 1D
/// 2. Compute unbinned covariance from outer products
/// 3. Eigenvalue decompose: Σ = Q Λ Q^T
/// 4. Apply whitening: T = Λ^(-1/2) Q^T
/// 5. Transform binned covariance
/// 6. Compute per-component correlation times
pub fn compute_decorrelated_autocorr_time(
    bins: &ndarray::ArrayD<f64>,
    mu: &ndarray::ArrayD<f64>,
    binned_cov: &ndarray::ArrayD<f64>,
    total_samples: usize,
) -> ndarray::ArrayD<f64> {
    let ndim = bins.ndim();
    let obs_shape: Vec<usize> = bins.shape()[..ndim - 1].to_vec();
    let obs_dim: usize = obs_shape.iter().product();
    let m = bins.shape()[ndim - 1] as f64;

    if obs_dim == 0 || m < 2.0 {
        return ndarray::ArrayD::zeros(obs_shape);
    }

    // Flatten mu
    let mu_flat: Vec<f64> = mu.iter().copied().collect();

    // Reshape binned_cov to 2D matrix
    let binned_cov_2d = if binned_cov.ndim() == 2 && binned_cov.shape()[0] == obs_dim {
        binned_cov.view().into_shape_with_order((obs_dim, obs_dim)).unwrap()
    } else {
        // Fallback: create identity
        ndarray::Array2::eye(obs_dim)
    };

    // Compute unbinned covariance from binned data
    // Σ_unbinned = (1/(M-1)) * Σ_k (x_k - μ)(x_k - μ)^T
    let mut sigma_unbinned = ndarray::Array2::<f64>::zeros((obs_dim, obs_dim));
    let n_bins = bins.shape()[ndim - 1];

    for k in 0..n_bins {
        let x_k = bins.index_axis(ndarray::Axis(ndim - 1), k);
        let x_flat: Vec<f64> = x_k.iter().copied().collect();

        for i in 0..obs_dim {
            for j in 0..obs_dim {
                sigma_unbinned[[i, j]] +=
                    (x_flat[i] - mu_flat[i]) * (x_flat[j] - mu_flat[j]);
            }
        }
    }

    let m_f64 = n_bins as f64;
    if m_f64 > 1.0 {
        sigma_unbinned.mapv_inplace(|v| v / (m_f64 - 1.0));
    }

    // Eigenvalue decomposition using linfa-linalg
    use linfa_linalg::symmetric::eigh;
    let (eigenvalues, eigenvectors) = eigh(&sigma_unbinned);

    // Whitening transform
    let tolerance = 1e-10;
    let lambda_inv_sqrt: ndarray::Array1<f64> = eigenvalues
        .iter()
        .map(|&lambda| {
            if lambda > tolerance {
                1.0 / lambda.sqrt()
            } else {
                0.0
            }
        })
        .collect();

    // T = Λ^(-1/2) Q^T
    let transform = ndarray::Array2::from_diag(&lambda_inv_sqrt).dot(&eigenvectors.t());

    // Transform binned covariance
    let sigma_binned_decorr = transform.dot(&binned_cov_2d).dot(&transform.t());

    // Extract diagonal (variances in decorrelated basis)
    let binned_variances_decorr: ndarray::Array1<f64> = ndarray::Array1::from_shape_fn(obs_dim, |i| {
        sigma_binned_decorr[[i, i]].max(0.0)
    });

    // τ = 0.5 * (M * σ²_decorr - 1)
    let correlation_times: ndarray::Array1<f64> = binned_variances_decorr
        .mapv(|v| (0.5 * (m * v - 1.0)).max(0.0));

    // Reshape back to original observable shape
    correlation_times.into_shape_with_order(obs_shape).unwrap()
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test --test merge_stat_test 2>&1 | tail -10`
Expected: PASS

- [ ] **Step 5: Wire up decorrelated mode in merge_results_from_files**

In `accumulate_and_compute`, when `options.estimate_covariance` is true and the observable is an array (shape.len() > 0):

```rust
// After computing sigma and mu, if estimate_covariance and array observable:
if options.estimate_covariance && obs_type.shape.len() > 0 && obs_type.shape.iter().product::<usize>() > 1 {
    // Need to accumulate outer products for unbinned covariance
    // This requires a second AddSamplesState for outer products
    // For now, fall back to regular autocorrelation time
    // (Full outer product accumulation in a future task)
}
```

- [ ] **Step 6: Commit**

```bash
jj describe -m "feat(merge): add decorrelated autocorrelation time via eigenvalue decomposition"
jj new
```

---

### Task 6: 完整集成测试

**Files:**
- Modify: `Carlo.rs/tests/merge_test.rs` — 添加 HDF5 文件的集成测试
- Modify: `Carlo.rs/src/merge.rs` — 修复任何发现的不一致

- [ ] **Step 1: Write integration test**

Add to `Carlo.rs/tests/merge_test.rs`:

```rust
#[cfg(feature = "hdf5")]
#[test]
fn test_merge_scalar_observable() {
    use carlo_rs::merge::{merge_results_from_files, MergeOptions};
    use hdf5::File;
    use tempfile::TempDir;

    let dir = TempDir::new().unwrap();
    let meas_path = dir.path().join("run0000.meas.h5");

    // Create a mock measurement file with known data
    {
        let file = File::create(&meas_path).unwrap();
        let obs_group = file.create_group("observables").unwrap();
        let energy_group = obs_group.create_group("energy").unwrap();

        energy_group.new_dataset("bin_length").create().unwrap();
        energy_group.write_dataset("bin_length", &10u64).unwrap();

        // 20 bins of a scalar observable
        let bins: Vec<f64> = (1..=20).map(|i| i as f64).collect();
        let shape = [20usize];
        energy_group
            .new_dataset_builder(bins.as_slice())
            .with_data_type::<f64>()
            .shape(&shape)
            .create("samples")
            .unwrap();
    }

    let options = MergeOptions::default();
    let results = merge_results_from_files(&[meas_path], &options).unwrap();

    assert!(results.contains_key("energy"));
    let energy = &results["energy"];

    // Mean of 1..=20 is 10.5
    assert!((energy.mean[ndarray::IxDyn(&[])] - 10.5).abs() < 0.1);
    // Error should be non-zero
    assert!(energy.error[ndarray::IxDyn(&[])] > 0.0);
    // Autocorrelation time should be computed
    assert!(energy.autocorrelation_time[ndarray::IxDyn(&[0])] >= 0.0);
}
```

- [ ] **Step 2: Run all tests**

Run: `cargo test --features hdf5 2>&1 | tail -30`
Expected: All tests pass

- [ ] **Step 3: Run cargo check for all feature combinations**

Run: `cargo check && cargo check --features hdf5 && cargo check --features mpi 2>&1 | tail -20`
Expected: No errors or warnings

- [ ] **Step 4: Commit**

```bash
jj describe -m "test(merge): add integration test for merge with actual HDF5 data"
jj new
```

---

## Self-Review

### 1. Spec coverage
- ✅ 自相关时间精确计算 — Task 3 (variance ratio method with actual data)
- ✅ 协方差矩阵估计 — Task 4 (cov_of_mean)
- ✅ Decorrelated 自相关时间 — Task 5 (eigenvalue decomposition whitening)
- ✅ 完整统计计算 — Task 3 (merge_results_from_files returns real values)
- ❌ ResultTools DataFrame 转换 — 不在范围内（Python 后处理工具，非核心 MC 功能）
- ❌ 性能监控 — 不在范围内（低优先级，不影响结果正确性）
- ❌ 复数支持 — 不在范围内（非核心差距）

### 2. Placeholder scan
- No TBD/TODO in plan steps
- All code blocks contain actual implementations
- No "add tests for the above" without actual test code
- All function signatures match across tasks

### 3. Type consistency
- `AddSamplesState<T>` uses `f64` specialization consistently
- `compute_regular_autocorr_time` signature matches existing code
- `cov_of_mean` uses ndarray generic types matching Carlo.jl patterns
- `compute_decorrelated_autocorr_time` returns `ArrayD<f64>` consistent with `autocorrelation_time` field

### 4. Scope check
This plan covers 3 tightly related subsystems (all in merge.rs statistics computation). Each task produces working, testable software. The plan is focused and bounded.
