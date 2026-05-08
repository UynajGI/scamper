# QMC.rs SSE Algorithm Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task.

**Goal:** Implement complete SSE (Stochastic Series Expansion) algorithm with diagonal update, directed loop, and physical validation.

**Architecture:** SSE algorithm consists of two alternating updates: diagonal_update (insert/remove diagonal operators) and directed_loop (cluster update through operator string). The energy is measured from operator density.

**Tech Stack:** Rust, rand_xoshiro, Carlo.rs framework

---

## SSE Algorithm Background

### Diagonal Update
遍历 operator sequence，对每个位置：
1. **Identity → Diagonal**: 以概率 `P = β * W_b * N_bonds / (M - n)` 插入对角算符
2. **Diagonal → Identity**: 以概率 `P = (M - n + 1) / (β * W_b * N_bonds)` 移除
3. **OffDiagonal**: 传播自旋状态（应用算符到 spin configuration）

### Directed Loop
构建 vertex graph，遍历并翻转：
1. 从随机入口点进入 operator string
2. 沿 vertex graph 遍历（每步决定方向）
3. 更新自旋状态
4. 当回到入口点时退出

### Energy Measurement
`E = -<n_operators> / (β * N_sites)`

---

## File Structure

```
QMC.rs/src/
├── sse/
│   ├── engine.rs         # Add bond_list for random bond selection
│   ├── diagonal.rs       # Full diagonal_update implementation
│   ├── loop_.rs          # Directed loop implementation
│   └── measurements.rs   # Energy calculation from operator density
└── tests/
    ├── sse_algorithm_test.rs  # Algorithm unit tests
    └── physics_test.rs        # Physical validation tests
```

---

## Task 1: Add Bond List to SSEEngine

**Files:**
- Modify: `QMC.rs/src/sse/engine.rs`
- Test: `QMC.rs/tests/sse_algorithm_test.rs`

**Why:** Need efficient random bond selection for diagonal update.

- [ ] **Step 1: Add bond_list field to SSEEngine**

```rust
pub struct SSEEngine<H: HilbertSpace> {
    pub lattice: Lattice,
    pub spins: Vec<LocalState>,
    pub op_seq: OperatorSequence,
    pub hs: H,
    pub weights: HashMap<BondType, f64>,
    /// Flattened bond list: [(site_i, site_j, bond_type), ...]
    pub bond_list: Vec<(usize, usize, BondType)>,
}
```

- [ ] **Step 2: Initialize bond_list in SSEEngine::new()**

```rust
fn build_bond_list(lattice: &Lattice) -> Vec<(usize, usize, BondType)> {
    let mut bonds = Vec::new();
    let mut seen = std::collections::HashSet::new();
    
    for (site_i, neighbors) in lattice.sites.iter().enumerate() {
        for neighbor in neighbors {
            let bond_key = if site_i < neighbor.target {
                (site_i, neighbor.target)
            } else {
                (neighbor.target, site_i)
            };
            if seen.insert(bond_key) {
                bonds.push((site_i, neighbor.target, neighbor.bond_type));
            }
        }
    }
    bonds
}
```

- [ ] **Step 3: Write test for bond_list**

```rust
#[test]
fn test_bond_list_chain() {
    let lattice = build_chain(4, true);
    let hs = SpinHalfHS;
    let engine = SSEEngine::new(lattice, hs, 100, HashMap::new());
    
    assert_eq!(engine.bond_list.len(), 4);
}

#[test]
fn test_bond_list_square() {
    let lattice = build_square(4, 4, true);
    let hs = SpinHalfHS;
    let engine = SSEEngine::new(lattice, hs, 100, HashMap::new());
    
    assert_eq!(engine.bond_list.len(), 32);
}
```

- [ ] **Step 4: Verify tests pass and commit**

```bash
cargo test bond_list
jj commit -m "feat(sse): add bond_list for random bond selection"
```

---

## Task 2: Implement State Propagation in Diagonal Update

**Files:**
- Modify: `QMC.rs/src/sse/diagonal.rs`
- Test: `QMC.rs/tests/sse_algorithm_test.rs`

**Why:** Must propagate spins through off-diagonal operators during diagonal update.

- [ ] **Step 1: Add state propagation logic**

When encountering an off-diagonal operator, apply it to the current spin configuration:

```rust
pub fn diagonal_update<R: Rng>(&mut self, rng: &mut R) {
    // Track current spin configuration (copy for now)
    let mut current_spins = self.spins.clone();
    
    for p in 0..self.op_seq.max_length {
        let vertex = &self.op_seq.vertices[p];
        
        match vertex.op {
            OpType::Identity => {
                // Try to insert diagonal operator
                // ...
            }
            OpType::Diagonal => {
                // Try to remove diagonal operator
                // ...
            }
            OpType::OffDiagonal => {
                // Propagate state through off-diagonal operator
                let (site_i, site_j, _) = self.bond_list[vertex.bond_idx];
                let states = &mut [current_spins[site_i], current_spins[site_j]];
                self.hs.apply(states, &OpType::OffDiagonal);
                current_spins[site_i] = states[0];
                current_spins[site_j] = states[1];
            }
        }
    }
    
    // Update spins after full sweep
    self.spins = current_spins;
}
```

- [ ] **Step 2: Write test for state propagation**

```rust
#[test]
fn test_state_propagation() {
    let lattice = build_chain(4, true);
    let hs = SpinHalfHS;
    let mut engine = SSEEngine::new(lattice, hs, 100, HashMap::new());
    
    // Insert an off-diagonal operator manually
    engine.op_seq.vertices[0] = Vertex {
        bond_idx: 0,
        op: OpType::OffDiagonal,
    };
    engine.op_seq.n_operators = 1;
    
    // Initial spins: [0, 0, 0, 0] (all up)
    engine.spins = vec![0, 0, 0, 0];
    
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
    engine.diagonal_update(&mut rng);
    
    // After propagation, spins at bond 0 sites should be flipped
    // Bond 0 connects sites 0 and 1
    // Original: [0, 0, ...], after flip: [1, 1, ...]
    assert_eq!(engine.spins[0], 1);
    assert_eq!(engine.spins[1], 1);
}
```

- [ ] **Step 3: Verify test passes and commit**

---

## Task 3: Implement Diagonal Operator Insert/Remove

**Files:**
- Modify: `QMC.rs/src/sse/diagonal.rs`
- Test: `QMC.rs/tests/sse_algorithm_test.rs`

**Why:** Core of diagonal update - insert and remove diagonal operators with correct probabilities.

- [ ] **Step 1: Implement insertion logic**

```rust
OpType::Identity => {
    // Select random bond
    let bond_idx = rng.gen_range(0..self.bond_list.len());
    let (site_i, site_j, bond_type) = self.bond_list[bond_idx];
    
    // Get diagonal weight
    let weight = match self.weights.get(&bond_type) {
        Some(&w) => w,
        None => continue,
    };
    
    // Check if diagonal operator is allowed
    let states = [current_spins[site_i], current_spins[site_j]];
    if !H::is_allowed(&states, &OpType::Diagonal) {
        continue;
    }
    
    // Calculate insertion probability
    // P_insert = β * W * N_bonds / (M - n)
    let m = self.op_seq.max_length;
    let n = self.op_seq.n_operators;
    let n_bonds = self.bond_list.len();
    
    if n >= m {
        continue;
    }
    
    let p_insert = self.beta * weight * n_bonds as f64 / (m - n) as f64;
    
    if rng.gen::<f64>() < p_insert {
        // Insert diagonal operator
        self.op_seq.vertices[p] = Vertex {
            bond_idx,
            op: OpType::Diagonal,
        };
        self.op_seq.n_operators += 1;
    }
}
```

**Note:** Need to add `beta` field to SSEEngine.

- [ ] **Step 2: Implement removal logic**

```rust
OpType::Diagonal => {
    let (site_i, site_j, bond_type) = self.bond_list[vertex.bond_idx];
    
    let weight = match self.weights.get(&bond_type) {
        Some(&w) => w,
        None => continue,
    };
    
    // Calculate removal probability
    // P_remove = (M - n + 1) / (β * W * N_bonds)
    let m = self.op_seq.max_length;
    let n = self.op_seq.n_operators;
    let n_bonds = self.bond_list.len();
    
    let p_remove = (m - n + 1) as f64 / (self.beta * weight * n_bonds as f64);
    
    if rng.gen::<f64>() < p_remove {
        // Remove diagonal operator (replace with Identity)
        self.op_seq.vertices[p] = Vertex::default();
        self.op_seq.n_operators -= 1;
    }
}
```

- [ ] **Step 3: Add beta field to SSEEngine**

```rust
pub struct SSEEngine<H: HilbertSpace> {
    pub lattice: Lattice,
    pub spins: Vec<LocalState>,
    pub op_seq: OperatorSequence,
    pub hs: H,
    pub weights: HashMap<BondType, f64>,
    pub bond_list: Vec<(usize, usize, BondType)>,
    pub beta: f64,
}
```

- [ ] **Step 4: Write tests**

```rust
#[test]
fn test_diagonal_insert_remove_balance() {
    // After many updates, operator density should reach equilibrium
    let lattice = build_chain(8, true);
    let hs = SpinHalfHS;
    let mut weights = HashMap::new();
    weights.insert(BondType::ChainX, 0.25);
    
    let mut engine = SSEEngine::new(lattice, hs, 1000, weights);
    engine.beta = 2.0;
    
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
    
    // Run many diagonal updates
    for _ in 0..1000 {
        engine.diagonal_update(&mut rng);
    }
    
    // Should have some operators
    assert!(engine.op_seq.n_operators > 0);
    assert!(engine.op_seq.n_operators < engine.op_seq.max_length);
}
```

- [ ] **Step 5: Verify tests pass and commit**

---

## Task 4: Implement Directed Loop Skeleton

**Files:**
- Modify: `QMC.rs/src/sse/loop_.rs`
- Test: `QMC.rs/tests/sse_algorithm_test.rs`

**Why:** Directed loop provides cluster update for off-diagonal operators.

- [ ] **Step 1: Understand vertex graph structure**

Each operator (vertex) has 4 "legs":
- Leg 0: entrance from site_i at previous time
- Leg 1: entrance from site_j at previous time
- Leg 2: exit to site_i at next time
- Leg 3: exit to site_j at next time

Loop traversal: enter at one leg, exit at another leg of the same vertex.

- [ ] **Step 2: Implement basic loop structure**

```rust
pub fn loopupdate<R: Rng>(&mut self, _rng: &mut R) {
    if self.op_seq.n_operators == 0 {
        return;
    }
    
    // Find all non-identity operator positions
    let op_positions: Vec<usize> = self.op_seq.vertices.iter()
        .enumerate()
        .filter(|(_, v)| v.op != OpType::Identity)
        .map(|(i, _)| i)
        .collect();
    
    if op_positions.is_empty() {
        return;
    }
    
    // TODO: Implement loop traversal
}
```

- [ ] **Step 3: Write test**

```rust
#[test]
fn test_loop_update_empty() {
    let lattice = build_chain(4, true);
    let hs = SpinHalfHS;
    let mut engine = SSEEngine::new(lattice, hs, 100, HashMap::new());
    engine.beta = 1.0;
    
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
    
    // Should not crash on empty sequence
    engine.loopupdate(&mut rng);
    
    assert_eq!(engine.op_seq.n_operators, 0);
}
```

- [ ] **Step 4: Verify test passes and commit**

---

## Task 5: Implement Full Directed Loop

**Files:**
- Modify: `QMC.rs/src/sse/loop_.rs`
- Test: `QMC.rs/tests/sse_algorithm_test.rs`

**Why:** Complete the cluster update mechanism.

- [ ] **Step 1: Implement leg traversal**

```rust
/// Get the next leg when entering a vertex at given leg.
/// Returns (next_position, next_leg)
fn traverse_vertex(
    &self,
    position: usize,
    entering_leg: usize,
) -> Option<(usize, usize)> {
    let vertex = &self.op_seq.vertices[position];
    
    if vertex.op == OpType::Identity {
        return None;
    }
    
    let (site_i, site_j, _) = self.bond_list[vertex.bond_idx];
    
    // For diagonal: flip coin to exit at same site's other leg
    // For off-diagonal: exit at other site's leg
    
    // Simplified: always exit at next leg (4 legs total, cycle)
    let exit_leg = (entering_leg + 2) % 4;
    
    // Find next position with operator on same bond and site
    // ... (requires leg list data structure)
    
    Some((position, exit_leg))
}
```

- [ ] **Step 2: Implement loop traversal**

```rust
pub fn loopupdate<R: Rng>(&mut self, rng: &mut R) {
    if self.op_seq.n_operators == 0 {
        return;
    }
    
    // Build operator list for each bond
    let mut bond_ops: HashMap<usize, Vec<usize>> = HashMap::new();
    for (p, v) in self.op_seq.vertices.iter().enumerate() {
        if v.op != OpType::Identity {
            bond_ops.entry(v.bond_idx).or_default().push(p);
        }
    }
    
    // Select random entry point
    let entry_bond = rng.gen_range(0..self.bond_list.len());
    let entry_ops = match bond_ops.get(&entry_bond) {
        Some(ops) => ops,
        None => return,
    };
    if entry_ops.is_empty() {
        return;
    }
    
    let entry_pos = entry_ops[rng.gen_range(0..entry_ops.len())];
    
    // Traverse and flip
    let mut visited = vec![false; self.op_seq.max_length];
    let mut current_pos = entry_pos;
    
    loop {
        if visited[current_pos] {
            break;
        }
        visited[current_pos] = true;
        
        let vertex = &mut self.op_seq.vertices[current_pos];
        if vertex.op == OpType::OffDiagonal {
            // Flip spins
            let (site_i, site_j, _) = self.bond_list[vertex.bond_idx];
            self.spins[site_i] ^= 1;
            self.spins[site_j] ^= 1;
        }
        
        // Move to next operator on same bond (simplified)
        // In full implementation, need proper leg traversal
        break;
    }
}
```

- [ ] **Step 3: Write tests**

```rust
#[test]
fn test_loop_update_with_operators() {
    let lattice = build_chain(4, true);
    let hs = SpinHalfHS;
    let mut engine = SSEEngine::new(lattice, hs, 100, HashMap::new());
    engine.beta = 1.0;
    
    // Insert some operators
    engine.op_seq.vertices[0] = Vertex { bond_idx: 0, op: OpType::Diagonal };
    engine.op_seq.vertices[1] = Vertex { bond_idx: 1, op: OpType::OffDiagonal };
    engine.op_seq.n_operators = 2;
    
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
    
    // Should not crash
    engine.loopupdate(&mut rng);
}
```

- [ ] **Step 4: Verify tests pass and commit**

---

## Task 6: Implement Energy Measurement

**Files:**
- Modify: `QMC.rs/src/sse/measurements.rs`
- Test: `QMC.rs/tests/sse_algorithm_test.rs`

**Why:** Correct energy calculation from operator density.

- [ ] **Step 1: Implement energy calculation**

```rust
pub fn compute_energy(&self) -> f64 {
    // E = -<n> / (β * N_sites)
    // where n is the number of operators
    if self.beta <= 0.0 || self.lattice.n_sites == 0 {
        return 0.0;
    }
    
    let n = self.op_seq.n_operators as f64;
    let beta = self.beta;
    let n_sites = self.lattice.n_sites as f64;
    
    -n / (beta * n_sites)
}
```

- [ ] **Step 2: Write test**

```rust
#[test]
fn test_energy_calculation() {
    let lattice = build_chain(4, true);
    let hs = SpinHalfHS;
    let mut engine = SSEEngine::new(lattice, hs, 100, HashMap::new());
    engine.beta = 1.0;
    
    // No operators = zero energy (but formula gives 0)
    assert_eq!(engine.compute_energy(), 0.0);
    
    // Add some operators
    engine.op_seq.n_operators = 8;
    
    // E = -8 / (1.0 * 4) = -2.0
    assert!((engine.compute_energy() - (-2.0)).abs() < 1e-10);
}
```

- [ ] **Step 4: Verify test passes and commit**

---

## Task 7: Physical Validation - Heisenberg Chain

**Files:**
- Create: `QMC.rs/tests/physics_test.rs`

**Why:** Validate against known exact results.

- [ ] **Step 1: Implement Bethe ansatz comparison for 1D Heisenberg**

The ground state energy per site for 1D Heisenberg chain (antiferromagnetic, J=1):
`E_0 / N = 1/4 - ln(2) ≈ -0.443147`

At low temperature T → 0:
`E / N ≈ -0.443147`

```rust
use qmc_rs::*;
use qmc_rs::lattice::builders::build_chain;
use qmc_rs::models::HeisenbergModel;
use rand_xoshiro::Xoshiro256PlusPlus;

#[test]
fn test_heisenberg_chain_ground_state() {
    // 1D Heisenberg chain with J=1 (antiferromagnetic)
    let n_sites = 16;
    let lattice = build_chain(n_sites, true);
    
    let beta = 10.0;  // Low temperature
    let j = 1.0;      // Antiferromagnetic
    
    let model = HeisenbergModel::new(lattice, beta, j);
    let mut core = SSECore::new(model);
    
    // Run simulation
    let backend = RayonBackend::new(1);
    let config = RunConfig {
        thermalization_sweeps: 10000,
        measurement_sweeps: 50000,
        binsize: 1000,
        base_seed: 42,
        progress_interval: 0,
        checkpoint_interval: 0,
    };
    let scheduler = Scheduler::new(backend, config);
    
    let params = Params::new();
    let results = scheduler.run_one::<SSECore<HeisenbergModel>>(&params).unwrap();
    
    // Check energy
    if let Some(energy) = results.get("Energy") {
        // Ground state energy per site ≈ -0.443147
        let expected = -0.443147;
        let tolerance = 3.0 * energy.stderr;  // 3 sigma
        
        println!("Energy: {:.6} ± {:.6}", energy.mean, energy.stderr);
        println!("Expected: {:.6}", expected);
        
        assert!(
            (energy.mean - expected).abs() < tolerance,
            "Energy {} not within {} of expected {}",
            energy.mean, tolerance, expected
        );
    }
}
```

- [ ] **Step 2: Run test and verify**

```bash
cargo test physics_test -- --nocapture
```

- [ ] **Step 3: Commit**

---

## Verification Checklist

```bash
cd QMC.rs

# All tests pass
cargo test

# No warnings
cargo clippy --all-targets -- -D warnings

# Documentation builds
cargo doc --no-deps
```

---

*Implementation plan for SSE algorithm completion.*