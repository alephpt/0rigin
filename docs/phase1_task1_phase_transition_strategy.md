# Phase Transition Universality Class Determination Strategy

## Executive Summary

**Current Status**: The SMFT system shows a noise-driven phase transition with critical exponent β = 0.0988 ± 0.0038, which represents a **7σ deviation** from the claimed 2D Ising universality class (β = 0.125). This discrepancy demands rigorous investigation to determine whether SMFT represents:
1. A novel universality class with unique critical behavior
2. Finite-size effects masking true 2D Ising behavior
3. Crossover behavior between universality classes
4. Numerical artifacts from insufficient analysis

## Part 1: Current Status Assessment

### 1.1 What We Have

**Implementation**:
- `PhaseTransitionAnalyzer` class with basic critical point detection
- Power-law fitting for β exponent extraction
- Single grid size (64×64) noise scan with 21 points
- Simple derivative-based critical point location (σ_c = 0.85)

**Results**:
- β = 0.0988 ± 0.0038 (measured)
- High fit quality (R² = 0.9882)
- Consistent reproducibility across runs
- Clear second-order transition behavior

### 1.2 Critical Gaps

1. **No finite-size scaling**: Only tested 64×64 grid
2. **Single exponent**: Only β measured, need ν, γ, η for full classification
3. **No data collapse**: Haven't verified scaling hypothesis
4. **Limited σ resolution**: 21 points may miss critical region details
5. **No Binder cumulant**: Can't precisely locate σ_c(L)
6. **No correlation length**: Missing ξ divergence analysis

## Part 2: Theoretical Background

### 2.1 Universality Classes in 2D

| Class | β | ν | γ | η | α | System Examples |
|-------|---|---|---|---|---|-----------------|
| **2D Ising** | 1/8 = 0.125 | 1 | 7/4 = 1.75 | 1/4 = 0.25 | 0 | Ferromagnets, Z₂ symmetry |
| **2D XY** | ≈0.23 | ≈0.67 | ≈1.32 | ≈0.04 | ≈-0.01 | Superfluids, U(1) symmetry |
| **2D 3-state Potts** | 1/9 ≈ 0.111 | 5/6 ≈ 0.833 | 13/9 ≈ 1.44 | 4/15 ≈ 0.267 | 1/3 | Z₃ symmetry |
| **2D 4-state Potts** | 1/12 ≈ 0.083 | 2/3 ≈ 0.667 | 7/6 ≈ 1.17 | 1/4 = 0.25 | 2/3 | Z₄ symmetry |
| **Mean Field** | 0.5 | 0.5 | 1 | 0 | 0 | High dimensions, long-range |

### 2.2 Physics Determining Universality

**Key Factors**:
1. **Symmetry**: Order parameter symmetry group (discrete Z₂ vs continuous U(1))
2. **Dimensionality**: System dimension (2D in our case)
3. **Range of interactions**: Short-range vs long-range
4. **Conservation laws**: Conserved vs non-conserved order parameter

**SMFT Characteristics**:
- **Order parameter**: R (synchronization amplitude) - real, positive
- **Symmetry**: Phase θ has U(1) symmetry, but R breaks to Z₂-like
- **Interactions**: Kuramoto coupling is local (nearest-neighbor dominant)
- **Dynamics**: Non-conserved order parameter with damping

### 2.3 Finite-Size Scaling Theory

Near criticality, observables follow scaling forms:

```
⟨R⟩(σ, L) = L^(-β/ν) f_R((σ - σ_c)L^(1/ν))
χ(σ, L) = L^(γ/ν) f_χ((σ - σ_c)L^(1/ν))
ξ(σ, L) = L g_ξ((σ - σ_c)L^(1/ν))
```

Where:
- L = system size
- f_R, f_χ, g_ξ = universal scaling functions
- Data collapse occurs when plotted as X = (σ - σ_c)L^(1/ν), Y = ⟨R⟩L^(β/ν)

## Part 3: Implementation Strategy

### 3.1 Phase I: Multi-Size Scaling Analysis

**Objective**: Extract true critical exponents via finite-size scaling

**Grid Sizes**: L = 16, 32, 64, 128, 256 (5 system sizes)

**Config Template**:
```yaml
# config/phase_transition_FSS_L{size}.yaml
test_name: "phase_transition_FSS_L{size}"
grid:
  size_x: {size}
  size_y: {size}
  L_domain: 100.0  # Keep physical size constant
physics:
  noise_scan: [0.70, 0.72, 0.74, ..., 0.98, 1.00]  # 31 points near σ_c
  total_steps: 8000  # Longer for larger systems
validation:
  scenario: "phase_transition_FSS"
```

**New Analysis Methods Needed**:
```cpp
// src/validation/FiniteSizeScaling.h
class FiniteSizeScaling {
public:
    struct ScalingResult {
        float sigma_c;        // True critical point
        float beta;           // Order parameter exponent
        float nu;             // Correlation length exponent
        float gamma;          // Susceptibility exponent
        float eta;            // Anomalous dimension
        float quality;        // Data collapse quality metric
    };

    // Binder cumulant: U_L = 1 - ⟨R⁴⟩/(3⟨R²⟩²)
    static float computeBinderCumulant(const std::vector<float>& R_field);

    // Susceptibility: χ = L²(⟨R²⟩ - ⟨R⟩²)
    static float computeSusceptibility(const std::vector<float>& R_field, int L);

    // Find σ_c from Binder cumulant crossing
    static float findCriticalPointFromBinder(
        const std::map<int, std::vector<DataPoint>>& data_by_L);

    // Extract exponents from FSS collapse
    static ScalingResult performDataCollapse(
        const std::map<int, std::vector<DataPoint>>& data_by_L);
};
```

### 3.2 Phase II: Multiple Observable Analysis

**Objective**: Measure independent exponents for consistency check

**Observables to Track**:
1. **Order Parameter**: ⟨R⟩ → β
2. **Susceptibility**: χ = L²Var(R) → γ
3. **Correlation Length**: ξ from R(x)R(0) decay → ν
4. **Specific Heat**: C = ∂E/∂T (or noise analog) → α
5. **Binder Cumulant**: U_L → σ_c(L) precisely

**Implementation**:
```cpp
// Extend ObservableComputer
void computeCorrelationFunction(const std::vector<float>& R_field,
                                std::vector<float>& correlation) {
    // G(r) = ⟨R(x+r)R(x)⟩ - ⟨R⟩²
    // Fit to G(r) ~ exp(-r/ξ) to extract ξ
}

float extractCorrelationLength(const std::vector<float>& correlation) {
    // Exponential fit to extract ξ
    // ξ ~ |σ - σ_c|^(-ν)
}
```

### 3.3 Phase III: Data Collapse Validation

**Objective**: Verify scaling hypothesis and universality

**Process**:
1. Plot ⟨R⟩L^(β/ν) vs (σ - σ_c)L^(1/ν)
2. Optimize σ_c, β/ν, 1/ν for best collapse
3. Repeat for χL^(-γ/ν) vs (σ - σ_c)L^(1/ν)
4. Verify consistency of ν from different observables

**Python Analysis Script**:
```python
# analyze_FSS_collapse.py
import numpy as np
from scipy.optimize import minimize

def collapse_quality(params, data_dict):
    """Measure quality of data collapse"""
    sigma_c, beta_over_nu, one_over_nu = params

    # Transform data to scaling variables
    X_all = []
    Y_all = []
    for L, data in data_dict.items():
        X = (data['sigma'] - sigma_c) * L**one_over_nu
        Y = data['R_mean'] * L**beta_over_nu
        X_all.extend(X)
        Y_all.extend(Y)

    # Measure spread in Y for similar X values
    # Lower spread = better collapse
    return compute_spread_metric(X_all, Y_all)

# Optimize for best collapse
result = minimize(collapse_quality, x0=[0.85, 0.125, 1.0])
sigma_c, beta_over_nu, one_over_nu = result.x
beta = beta_over_nu / one_over_nu
nu = 1.0 / one_over_nu
```

### 3.4 Phase IV: Universality Class Identification

**Decision Tree**:
```
Measure: β, ν, γ, η from FSS
         ↓
β ≈ 0.125 ± 0.01?
    YES → Check ν ≈ 1.0 ± 0.05?
              YES → Check γ ≈ 1.75 ± 0.1?
                        YES → 2D ISING CLASS ✓
                        NO  → Mixed/Crossover
              NO  → Not 2D Ising
    NO  → β ≈ 0.23 ± 0.02?
              YES → Check ν ≈ 0.67?
                        YES → 2D XY CLASS ✓
                        NO  → Mixed behavior
              NO  → β ≈ 0.111 ± 0.01?
                        YES → 3-STATE POTTS ✓
                        NO  → β ≈ 0.099?
                                  YES → NOVEL CLASS! 🎉
                                  NO  → Check other Potts
```

## Part 4: Expected Results

### 4.1 Scenario A: True 2D Ising

**Expected**:
- FSS reveals β → 0.125 as L → ∞
- Current β = 0.099 due to finite-size corrections
- Data collapse with ν = 1.0 works perfectly
- γ = 1.75, η = 0.25 confirmed

**Implications**: SMFT has emergent Z₂ symmetry despite U(1) phase

### 4.2 Scenario B: 2D XY Class

**Expected**:
- β ≈ 0.23, very different from current measurement
- Would require re-analysis of transition region
- Kosterlitz-Thouless physics possible

**Implications**: Continuous symmetry dominates

### 4.3 Scenario C: Novel Universality Class

**Expected**:
- β remains ≈ 0.099 even with FSS
- Consistent set of exponents not matching known classes
- Unique data collapse with non-standard ν

**Implications**:
- **Major discovery**: New critical behavior in Kuramoto systems
- Requires theoretical explanation
- Publishable result for physics journals

### 4.4 Scenario D: Crossover Behavior

**Expected**:
- Exponents vary with L (non-universal)
- No clean data collapse
- Different observables give inconsistent exponents

**Implications**: System exhibits crossover between universality classes

## Part 5: Success Criteria

### 5.1 Minimum Requirements

1. **Multi-size data**: At least 5 grid sizes with L_max/L_min ≥ 16
2. **Critical point precision**: σ_c determined to ±0.001
3. **Exponent precision**: β to ±0.005, ν to ±0.02
4. **Data collapse quality**: R² > 0.99 for scaling plots
5. **Consistency**: All observables give compatible exponents

### 5.2 Validation Checklist

- [ ] Binder cumulant curves cross at single σ_c
- [ ] χ shows correct L^(γ/ν) scaling at σ_c
- [ ] Correlation length ξ ~ L at criticality
- [ ] Hyperscaling relations satisfied: 2β + γ = 2 - α (2D)
- [ ] FSS collapse works for multiple observables
- [ ] Results stable against:
  - [ ] Timestep variation (dt/2, dt/4)
  - [ ] Equilibration time (2x, 4x longer)
  - [ ] Initial conditions (random vs ordered)

### 5.3 Publication Standards

For claiming "novel universality class":
1. **Statistical significance**: >10σ deviation from all known classes
2. **Systematic checks**: Rule out artifacts, crossover, corrections
3. **Theoretical support**: Identify symmetry/mechanism for new class
4. **Independent verification**: Multiple analysis methods agree

## Part 6: Implementation Roadmap

### Week 1: Infrastructure
- [ ] Implement `FiniteSizeScaling` class
- [ ] Add Binder cumulant calculation
- [ ] Add susceptibility calculation
- [ ] Create multi-size config generator

### Week 2: Data Collection
- [ ] Run L=16 full noise scan (0.7-1.0, 61 points)
- [ ] Run L=32 full noise scan
- [ ] Run L=64 full noise scan
- [ ] Run L=128 focused scan (0.8-0.9, 41 points)
- [ ] Run L=256 focused scan (if computationally feasible)

### Week 3: Analysis
- [ ] Binder cumulant crossing → σ_c
- [ ] FSS fits for β, ν from ⟨R⟩
- [ ] FSS fits for γ from χ
- [ ] Data collapse optimization
- [ ] Error analysis via bootstrap

### Week 4: Validation & Report
- [ ] Consistency checks across observables
- [ ] Systematic error estimation
- [ ] Theoretical interpretation
- [ ] Final report with universality classification

## Appendix A: Computational Requirements

**Estimated Runtime**:
- L=16: 1 min × 61 points = 1 hour
- L=32: 4 min × 61 points = 4 hours
- L=64: 16 min × 61 points = 16 hours
- L=128: 64 min × 41 points = 44 hours
- L=256: 256 min × 41 points = 175 hours

**Total**: ~240 CPU hours (10 days single-core, 1 day with 10 cores)

**Storage**: ~50 GB for all runs with snapshots

## Appendix B: Code Modifications Summary

**New Files**:
1. `src/validation/FiniteSizeScaling.{h,cpp}` - FSS analysis
2. `scripts/analyze_FSS_collapse.py` - Data collapse
3. `scripts/generate_FSS_configs.py` - Config generator
4. `config/phase_transition_FSS_*.yaml` - Size-specific configs

**Modified Files**:
1. `ObservableComputer` - Add Binder, susceptibility
2. `PhaseTransitionAnalyzer` - Integrate FSS methods
3. `SMFTTestRunner` - Handle FSS campaign mode

## Conclusion

The current β = 0.099 result is scientifically significant, representing either:
1. Strong finite-size effects masking 2D Ising behavior
2. A genuinely novel universality class in synchronization transitions

This strategy provides a rigorous path to definitively determine which scenario is correct through comprehensive finite-size scaling analysis. The discovery of a new universality class would be a major scientific contribution, while confirming 2D Ising would validate theoretical predictions about synchronization transitions.

**Recommendation**: Proceed with Phase I (multi-size scaling) immediately to determine whether finite-size effects explain the discrepancy. This requires ~1 week of focused effort but will definitively resolve the universality question.