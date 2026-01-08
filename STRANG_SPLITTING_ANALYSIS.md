# Strang Splitting Implementation & Analysis

**Date**: 2026-01-08
**Objective**: Achieve <0.01% energy conservation for Sine-Gordon particle scattering

## Implementation Summary

###  Strang Splitting (T-V-T)

Implemented proper 2nd-order symplectic Strang splitting for Sine-Gordon equation:

```
∂²θ/∂t² = ∇²θ - sin(θ)
```

**Splitting**:
1. **Half-step kinetic (T/2)**: θ += (dt/2)·π
2. **Full-step potential (V)**: π += dt·(∇²θ - sin(θ))
3. **Half-step kinetic (T/2)**: θ += (dt/2)·π

**Key difference from Velocity Verlet**: Nonlinear force evaluated at midpoint position θ_{n+1/2}, providing better symmetry for strong nonlinearity.

## Comparative Results

**Test Configuration**:
- Grid: 64×64×64
- dx = 1.0
- dt = 0.005
- Initial conditions: Gaussian wave packet (σ=5, A=0.1)
- Duration: 2000 steps (total time = 10)

| Integrator | Energy Drift | Time Reversibility | Status |
|------------|--------------|-------------------|---------|
| **Velocity Verlet** | 0.0679% | 7.45×10⁻⁹ rad | ✓ Symplectic |
| **Strang Splitting** | 0.0686% | 7.45×10⁻⁹ rad | ✓ Symplectic |
| **RK2 Symplectic** | 0.0684% | 3.73×10⁻⁹ rad | ✓ Symplectic |

## Key Finding

**All three integrators perform identically** (±0.001% variation), confirming:

1. ✅ **Temporal integration is correct** - all methods are 2nd-order symplectic
2. ✅ **Time reversibility is excellent** - O(10⁻⁹) phase error
3. ❌ **Energy drift ~0.068% is fundamental limit** of 2nd-order spatial discretization

## Root Cause Analysis

The 0.068% energy drift is **NOT due to time integration** but due to **spatial discretization error**:

### Dispersion Relation for Discrete Wave Equation

**Continuous**:
```
ω² = k²  (exact dispersion)
```

**Discrete (2nd-order centered differences)**:
```
ω²_discrete = (2/dx²)·[2 - cos(k·dx)] ≈ k²·[1 - k²·dx²/12 + O(dx⁴)]
```

**Error**: O(k²·dx²) truncation error causes frequency-dependent phase velocity → energy sloshes between spatial modes at ~0.07% level.

### Energy Oscillation Pattern

Observed energy oscillates with period ~600 steps:
```
E(t): 3.6838 → 3.6863 → 3.6838 → 3.6863
```

This is **resonant beating between discrete modes** - characteristic of grid dispersion.

## Achieving <0.01% Conservation

### ✅ SOLUTION IMPLEMENTED: 4th-Order Spatial Stencils

**Implementation** (2026-01-08):

1. **4th-order Laplacian** in evolution:
   ```cpp
   ∇²θ ≈ [-θ_{i±2} + 16θ_{i±1} - 30θ_i] / (12·dx²)  per direction
   Error: O(dx⁴)
   ```

2. **4th-order gradients** in energy computation:
   ```cpp
   ∂θ/∂x ≈ [θ_{i-2} - 8θ_{i-1} + 8θ_{i+1} - θ_{i+2}] / (12·dx)
   Error: O(dx⁴)
   ```

3. **Critical insight**: BOTH evolution AND energy measurement must use same discretization order

**Results**:

| Spatial Order | Energy Drift | Improvement | Status |
|---------------|--------------|-------------|--------|
| 2nd-order | 0.0686% | Baseline | ❌ 7× above threshold |
| 4th-order | **0.0038%** | **18× better** | ✅ **GO/NO-GO MET** |

**Verification**:
- ✅ Energy drift: 0.0038% < 0.01% threshold
- ✅ Time reversibility: <1e-9 rad phase error
- ✅ Laplacian accuracy: 62× better than 2nd-order (verified on analytic function)

**Cost**: 18-neighbor stencil vs 6-neighbor (3× stencil width, ~1.5× compute cost)

### Alternative Approaches (not pursued)

**Option 2: Finer Grid** - Tested and rejected
- 128³ grid showed worse conservation (0.14%) due to changing physical setup
- 8× memory cost impractical

**Option 3: Spectral Methods (FFT)** - Not needed
- Would achieve machine precision but 4th-order already meets requirement
- Higher implementation complexity unnecessary

## Validation Status

### ✅ FINAL STATUS (4th-order stencils) - 2026-01-08
- ✅ Strang Splitting implemented correctly (T-V-T)
- ✅ 4th-order Laplacian in evolution (18-neighbor stencil)
- ✅ 4th-order gradients in energy computation (consistent order)
- ✅ Time reversibility excellent (10⁻⁹ rad)
- ✅ **Energy drift 0.0038% < 0.01% threshold** ← **GO/NO-GO MET**
- ✅ All integrators symplectic and working

### Implementation Complete
**Target achieved**: <0.01% energy conservation ✓

**Files modified**:
- `include/ConservativeSolver.h` - Added SpatialOrder enum, computeLaplacian4thOrder()
- `src/ConservativeSolver.cpp` - 4th-order Laplacian + gradient implementation
- `CMakeLists.txt` - Added validation tests

**Tests created**:
- `test_spatial_order_comparison.cpp` - Verifies 2nd vs 4th order (PASS)
- `test_laplacian_accuracy.cpp` - Confirms 62× accuracy improvement
- `test_spatial_order_timestep.cpp` - Timestep independence validation
- `test_grid_resolution.cpp` - Grid convergence study

## Code Locations

- **ConservativeSolver.h**: `/home/persist/neotec/0rigin/include/ConservativeSolver.h`
- **ConservativeSolver.cpp**: `/home/persist/neotec/0rigin/src/ConservativeSolver.cpp`
- **Test**: `/home/persist/neotec/0rigin/test_integrator_comparison.cpp`
- **YAML Config**: `/home/persist/neotec/0rigin/config/particle_scattering_sine_gordon.yaml`

## References

- Hairer, E., Lubich, C., & Wanner, G. (2006). *Geometric Numerical Integration*
- Strikwerda, J. C. (2004). *Finite Difference Schemes and Partial Differential Equations*
- Trefethen, L. N. (2000). *Spectral Methods in MATLAB*

---

**Conclusion**: Strang splitting implementation is correct. Energy conservation limited by 2nd-order spatial discretization (~0.07%). Higher-order stencils required for <0.01% target.
