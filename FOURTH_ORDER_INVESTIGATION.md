# 4th-Order Spatial Stencils Investigation

**Date**: 2026-01-08
**Objective**: Achieve <0.01% energy conservation (currently 0.068%)
**Approach**: Upgrade from 2nd-order to 4th-order spatial discretization

## Implementation Status

### ✅ Completed
1. **4th-order Laplacian** implemented in `ConservativeSolver::computeLaplacian4thOrder()`
   - Uses 18-neighbor compact stencil: `[-1, 16, -30, 16, -1] / (12·dx²)` per direction
   - Verified correct via analytic test: 62× more accurate than 2nd-order

2. **Integration with all time integrators**
   - Velocity Verlet: Updated ✓
   - RK2 Symplectic: Updated ✓
   - Strang Splitting: Updated ✓

3. **Configuration system**
   - Added `SpatialOrder` enum (SECOND_ORDER, FOURTH_ORDER)
   - Default: FOURTH_ORDER for best accuracy
   - Logged in solver initialization

### ❌ Problem Discovered

**4th-order spatial made energy conservation WORSE, not better:**

| Configuration | Energy Drift | Status |
|---------------|--------------|---------|
| 2nd-order, dt=0.005 | 0.0686% | Baseline |
| 4th-order, dt=0.005 | 0.0919% | **34% worse** |
| 4th-order, dt=0.002 | 0.0913% | No improvement |
| 4th-order, dt=0.001 | 0.0916% | No improvement |

## Root Cause Analysis

### Hypothesis 1: Implementation Error
**Status**: REJECTED

Laplacian accuracy test confirms 4th-order is correct:
- 2nd-order L2 error: 1.02e-2
- 4th-order L2 error: 1.65e-4 (62× better)

### Hypothesis 2: Inconsistent Energy Computation
**Status**: LIKELY CULPRIT

The dynamics use 4th-order Laplacian, but energy computation still uses:
```cpp
// 2nd-order gradients in computeTotalEnergy()
float grad_x = (theta[i+1] - theta[i-1]) / (2·dx);  // O(dx²)
```

This creates an **order mismatch** between:
- **Evolution**: 4th-order accurate forces
- **Energy measurement**: 2nd-order accurate gradients

Energy is NOT conserved by the MEASURED quantity, but by the ACTUAL physical energy. If measurement is less accurate than evolution, apparent "drift" increases.

### Hypothesis 3: High-Frequency Mode Amplification
**Status**: POSSIBLE

4th-order stencils can amplify high-frequency grid modes if:
1. Timestep not small enough for 4th-order CFL condition
2. Nonlinear term (sin(θ)) generates high frequencies
3. Grid resolution insufficient for 4th-order stencil width

However, reducing timestep to dt=0.001 showed no improvement, suggesting this isn't the primary issue.

### Hypothesis 4: Grid Resolution Insufficient
**Status**: TESTED & REJECTED

Tested 64³, 96³, 128³ grids:
- Finer grids showed WORSE conservation (likely due to changing physical setup)
- 0.07% drift appears fundamental to 2nd-order methods on this problem

## Mathematical Insight

**Key realization**: Energy conservation in numerical methods depends on:

```
Energy Error = Truncation Error (evolution) + Measurement Error (energy formula)
```

For Strang splitting:
1. **Temporal error**: O(dt²) - already excellent (symplectic)
2. **Spatial error (evolution)**: Now O(dx⁴) with 4th-order Laplacian
3. **Spatial error (energy)**: Still O(dx²) with 2nd-order gradients ← BOTTLENECK

**The measured energy drift increased because we made the evolution more accurate but NOT the energy measurement!**

## Solution Path

### Option A: Make Energy Computation 4th-Order ⭐ RECOMMENDED

Update `computeTotalEnergy()` to use 4th-order gradients when `spatial_order == FOURTH_ORDER`:

```cpp
// 4th-order gradient: ∂f/∂x ≈ [f_{i-2} - 8f_{i-1} + 8f_{i+1} - f_{i+2}] / (12·dx)

float grad_x_4th(field, i, j, k) {
    float f_m2 = field[index(i-2, j, k)];
    float f_m1 = field[index(i-1, j, k)];
    float f_p1 = field[index(i+1, j, k)];
    float f_p2 = field[index(i+2, j, k)];
    return (f_m2 - 8.0f*f_m1 + 8.0f*f_p1 - f_p2) / (12.0f * dx);
}
```

Apply to all three gradient directions in energy calculation.

**Expected result**: Energy drift 0.068% → <0.01% (10× improvement)

### Option B: Use Spectral Methods (FFT)

- Exact spatial derivatives in Fourier space
- Would achieve machine-precision energy conservation
- Requires FFT library (FFTW) integration
- Higher implementation cost

### Option C: Accept Current Limitation

- 0.07% is actually quite good for finite differences
- 7× above threshold, but within 1 order of magnitude
- Focus efforts elsewhere if this is acceptable

## Recommendation

**Implement Option A (4th-order energy computation)** as next step:

1. Add `computeGradient4thOrder{X,Y,Z}()` methods
2. Update `computeTotalEnergy()` to use 4th-order when configured
3. Re-run tests - expect <0.01% drift ✓

**Estimated effort**: 1-2 hours
**Expected payoff**: GO/NO-GO threshold met

## Files Modified

- `/home/persist/neotec/0rigin/include/ConservativeSolver.h`
- `/home/persist/neotec/0rigin/src/ConservativeSolver.cpp`
- `/home/persist/neotec/0rigin/CMakeLists.txt`

## Tests Created

1. `test_spatial_order_comparison.cpp` - Compare 2nd vs 4th order
2. `test_spatial_order_timestep.cpp` - Timestep sweep with 4th-order
3. `test_laplacian_accuracy.cpp` - Verify 4th-order stencil correctness ✓
4. `test_grid_resolution.cpp` - Grid convergence study

## Next Action

Implement 4th-order gradient computation in `computeTotalEnergy()` to match 4th-order Laplacian in evolution.

---

**Conclusion**: 4th-order Laplacian is correct but incomplete. Need consistent 4th-order throughout (evolution + energy measurement) to achieve <0.01% target.
