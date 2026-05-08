# QMC.rs Philosophy

> This document defines core principles and design criteria, serving as the judgment basis for all architectural decisions.
>
> Established through Socratic dialogue on 2026-04-09.

---

## I. Ontology: What is QMC.rs?

### Core Definition

QMC.rs is a **Quantum Monte Carlo algorithm toolbox** — a component library between the Carlo.rs framework layer and concrete model implementations. It provides mid-granularity components that users assemble to implement QMC models.

### Position in the Stack

```
┌─────────────────────────────────────────────────────┐
│  User's Model Implementation                        │
│  (implements QMC traits, defines Hamiltonian)       │
├─────────────────────────────────────────────────────┤
│  QMC.rs — Algorithm Toolbox                         │
│  ├─ qmc_lattice::SSEMonteCarlo trait               │
│  ├─ qmc_lattice::SSECore (operator sequence + loop)│
│  ├─ qmc_continuous::VMMonteCarlo trait             │
│  ├─ qmc_impurity::CTHybMonteCarlo trait            │
│  └─ qmc_field::HMCMonteCarlo trait                 │
├─────────────────────────────────────────────────────┤
│  Carlo.rs — Core Framework                          │
│  ├─ MonteCarlo trait (base)                        │
│  ├─ Backend trait (parallel abstraction)           │
│  ├─ Scheduler (thermalization → measurement)       │
│  └─ Measurements + Results                         │
└─────────────────────────────────────────────────────┘
```

### What Users Implement

Users implement **Hamiltonian decomposition**: defining bond operators and their weights via `fn bond_operators(bond: Bond) -> Vec<Operator>`. The framework automatically constructs VertexData and handles loop updates internally.

**User responsibility**: Define what operators exist on each bond.
**Framework responsibility**: Handle operator sequence management, loop traversal, weight computation, measurement accumulation.

### Domain Coverage (Long-term Vision)

Based on the discipline-based classification framework, QMC.rs covers all four disciplinary origins:

| Layer | Trait | Methods |
|-------|-------|---------|
| Lattice QMC | `LatticeQMC` → `SSEMonteCarlo` / `WormMonteCarlo` | SSE, Worldline, Worm, Directed Loop |
| Continuous QMC | `ContinuousQMC` → `VMMonteCarlo` / `DMMonteCarlo` | VMC, DMC, AFQMC |
| Impurity QMC | `ImpurityQMC` → `CTHybMonteCarlo` / `SpinBosonMonteCarlo` | CT-HYB, CT-AUX, Spin-Boson |
| Field QMC | `FieldQMC` → `HMCMonteCarlo` | HMC (Lattice Field Theory) |

### Trait Hierarchy

```rust
// Carlo.rs (base)
trait MonteCarlo { fn sweep(&mut self, ctx: &mut Context); }

// QMC.rs (domain layer)
trait LatticeQMC: MonteCarlo { /* lattice-specific methods */ }
trait ContinuousQMC: MonteCarlo { /* continuous-space methods */ }
trait ImpurityQMC: MonteCarlo { /* impurity-specific methods */ }
trait FieldQMC: MonteCarlo { /* field-theory methods */ }

// QMC.rs (method layer)
trait SSEMonteCarlo: LatticeQMC {
    fn bond_operators(&self, bond: Bond) -> Vec<Operator>;
    fn measure(&mut self, ctx: &mut Context);
}
```

---

## II. Teleology: What is Success?

### Priority Hierarchy

**C > B > A > D**

| Priority | Criterion | Description |
|----------|-----------|-------------|
| **C (Foundation)** | Correctness | Results must match literature and reference implementations. Without correctness, nothing else matters. |
| **B (Necessity)** | Performance | Must achieve competitive performance (not slower than reference Fortran/C++ codes). Without performance, users won't adopt. |
| **A (Goal)** | Efficiency | Users implement new models in days, not weeks. This is the ultimate user value. |
| **D (Optional)** | Education | Code readability as teaching material — nice-to-have, not required. |

### Success Metrics

1. **Correctness**: Results from QMC.rs models match reference values (Carlo.jl, literature) within statistical error (≤ 3σ).
2. **Performance**: Sweep rate ≥ 80% of optimized Fortran implementations for standard models.
3. **Efficiency**: New model implementation time ≤ 3 days for users familiar with the trait system.

---

## III. Methodology: How to Achieve It?

### Core Design Principle: Progressive Complexity

Like `serde` in Rust ecosystem: simple defaults, advanced customization available.

**Principle**: Users can start with minimal implementation, progressively replace components with advanced versions.

```
Level 1: Implement bond_operators() only → default SSECore handles everything
Level 2: Override measure() for custom observables
Level 3: Replace SSECore with custom implementation for special cases
Level 4: Implement full SSEMonteCarlo trait manually
```

### Cross-Domain Architecture: Layered Trait Design

Different disciplinary origins have fundamentally different abstractions (continuous space vs discrete lattice vs impurity+bath). The trait hierarchy mirrors this:

- **Domain traits** (`LatticeQMC`, `ContinuousQMC`, etc.) encode domain-specific constraints
- **Method traits** (`SSEMonteCarlo`, `VMMonteCarlo`, etc.) provide algorithm-specific interfaces
- Each domain can evolve independently without breaking others

### Model Placement: Inside QMC.rs

Concrete models live inside QMC.rs domain modules:
- `qmc_lattice::models::Heisenberg`
- `qmc_lattice::models::BoseHubbard`
- `qmc_continuous::models::HydrogenMolecule` (future)

**Rationale**: Models demonstrate trait usage, provide test coverage, and reduce user friction. A separate Models.rs package would add management overhead without clear benefit.

### Implementation Order

Phase 1: `qmc_lattice` (SSE + Directed Loop)
Phase 2: `qmc_lattice::models` (Heisenberg, XXZ)
Phase 3: `qmc_impurity` (CT-HYB)
Phase 4: `qmc_continuous` (VMC)
Phase 5: `qmc_field` (HMC)

---

## IV. Boundaries: What is QMC.rs NOT?

### Explicit Exclusions

| Exclusion | Reason | Where it belongs |
|-----------|--------|------------------|
| **Visualization/Analysis** | Post-processing is domain-specific | Python scripts, Jupyter notebooks |
| **Parallel Scheduling** | Framework-level responsibility | Carlo.rs Backend trait |
| **MCMC Diagnostics** | Not QMC-specific | Carlo.rs or external tools |
| **Bayesian Inference** | Different problem class | Separate package |

### What it CAN Include

| Included | Reason |
|----------|--------|
| **Concrete models** | Demonstrate usage, provide test coverage |
| **Algorithm components** | Core value proposition |
| **Measurement helpers** | QMC-specific measurements (Green's functions, correlation lengths) |

---

## V. Decision Criteria

When facing design choices, test with these questions:

### Ontology Check
1. Is this a framework-level feature (belongs to Carlo.rs) or QMC-specific?
2. Which disciplinary origin does this feature serve? (Lattice/Continuous/Impurity/Field)
3. What granularity level is appropriate? (User responsibility vs framework responsibility)

### Teleology Check
4. Does this design preserve correctness guarantees?
5. Does this design enable competitive performance?
6. Does this design reduce user implementation effort?

### Methodology Check
7. Can users start simple and progressively customize?
8. Does the trait hierarchy respect domain boundaries?
9. Is the implementation order logical?

### Boundary Check
10. Is this feature within QMC.rs scope (not Carlo.rs, not post-processing)?

---

## VI. Key Architectural Decisions Summary

| Decision Point | Choice | Rationale |
|----------------|--------|-----------|
| **Position** | Algorithm toolbox (component library) | Between framework and models, maximum flexibility |
| **Granularity** | Mid-level (bond operators → framework handles loops) | Balance ease-of-use with customization |
| **Trait design** | Layered hierarchy by domain → method | Respect disciplinary origin differences |
| **Model placement** | Inside QMC.rs domain modules | Reduce friction, provide examples |
| **Success priority** | Correctness → Performance → Efficiency | Foundation before value |
| **Design principle** | Progressive complexity | Low barrier to entry, high ceiling for experts |

---

*This philosophy document was established through Socratic dialogue. It serves as the design criterion for all QMC.rs architectural decisions.*