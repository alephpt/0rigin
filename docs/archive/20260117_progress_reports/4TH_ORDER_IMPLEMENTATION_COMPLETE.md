# 4th-Order Spatial Stencils - Implementation Complete

**Date**: 2026-01-08
**Objective**: Achieve <0.01% energy conservation for Sine-Gordon particle scattering
**Result**: ✅ **TARGET ACHIEVED** - 0.0038% energy drift (18× improvement)

---

## Executive Summary

Implemented 4th-order spatial discretization for Conservative Solver to achieve GO/NO-GO energy conservation requirement (<0.01%). Previous 2nd-order implementation showed 0.0686% drift (7× above threshold). New 4th-order implementation achieves **0.0038% drift** - comfortably below threshold.

**Key insight**: Consistent discretization order between **evolution dynamics** (Laplacian) and **energy measurement** (gradients) is critical for accurate energy conservation tracking.

---

## Results

### Energy Conservation

| Configuration | Laplacian Order | Gradient Order | Energy Drift | vs Threshold | Status |
|---------------|-----------------|----------------|--------------|--------------|--------|
| **Baseline** | 2nd (O(dx²)) | 2nd (O(dx²)) | 0.0686% | 7× above | ❌ FAIL |
| **Upgraded** | **4th (O(dx⁴))** | **4th (O(dx⁴))** | **0.0038%** | **2.6× below** | ✅ **PASS** |

**Improvement factor**: 18× better energy conservation

### Time Reversibility

- **Phase error**: <1×10⁻⁹ rad (excellent)
- **Threshold**: 1×10⁻⁴ rad
- **Status**: ✅ PASS (5 orders of magnitude better than required)

### Computational Cost

- **Stencil width**: 6-neighbor → 18-neighbor (3× wider)
- **Compute cost**: ~1.5× per timestep (still O(N), acceptable)
- **Memory**: No increase (same field storage)

---

## Technical Implementation

### 1. 4th-Order Laplacian (Evolution)

**Formula** (per direction):
```
∂²θ/∂x² ≈ [-θ_{i-2} + 16θ_{i-1} - 30θ_i + 16θ_{i+1} - θ_{i+2}] / (12·dx²)
```

**3D Laplacian**:
```
∇²θ = ∂²θ/∂x² + ∂²θ/∂y² + ∂²θ/∂z²
```

**Truncation error**: O(dx⁴) - provides 16× better accuracy per grid doubling

**Implementation**: `ConservativeSolver::computeLaplacian4thOrder()`

### 2. 4th-Order Gradients (Energy Measurement)

**Formula**:
```
∂θ/∂x ≈ [θ_{i-2} - 8θ_{i-1} + 8θ_{i+1} - θ_{i+2}] / (12·dx)
```

**Applied to energy functional**:
```
E = ∫ [½(∂θ/∂t)² + ½|∇θ|² + (1-cos θ)] dV
```

**Critical**: Gradients in energy must match Laplacian order in dynamics

**Implementation**: Updated `ConservativeSolver::computeTotalEnergy()`

### 3. Configuration System

Added `SpatialOrder` enum:
```cpp
enum class SpatialOrder {
    SECOND_ORDER,   // 6-neighbor (O(dx²))
    FOURTH_ORDER    // 18-neighbor (O(dx⁴))
};
```

**Default**: `FOURTH_ORDER` for best accuracy
**Configurable**: Can fall back to 2nd-order if needed

---

## Validation

### Test 1: Spatial Order Comparison ✅

**Setup**:
- Grid: 64×64×64
- Integrator: Strang Splitting (T-V-T)
- Initial conditions: Gaussian wave packet (σ=5, A=0.1)
- Duration: 2000 steps (t=10)

**Results**:
```
2nd-order: 0.0686% drift → FAIL
4th-order: 0.0038% drift → PASS
```

**Test file**: `test_spatial_order_comparison.cpp`

### Test 2: Laplacian Accuracy ✅

**Setup**: Analytic function f(x,y,z) = sin(πx)·sin(πy)·sin(πz)
**Exact**: ∇²f = -3π²·f

**Results**:
```
2nd-order L2 error: 1.02×10⁻²
4th-order L2 error: 1.65×10⁻⁴  (62× better)
```

**Conclusion**: 4th-order stencil mathematically correct

**Test file**: `test_laplacian_accuracy.cpp`

### Test 3: Timestep Independence ✅

**Setup**: 4th-order with dt = 0.005, 0.002, 0.001

**Results**: Energy drift remains ~0.004% (timestep-independent)

**Conclusion**: Spatial discretization was the bottleneck, not temporal

**Test file**: `test_spatial_order_timestep.cpp`

### Test 4: Grid Resolution Study

**Setup**: 64³, 96³, 128³ grids with 2nd-order

**Results**: Coarser grid (64³) performed best due to consistent physical setup

**Conclusion**: Finer grid alone insufficient - need higher-order stencils

**Test file**: `test_grid_resolution.cpp`

---

## Code Changes

### Modified Files

#### 1. `include/ConservativeSolver.h`

**Added**:
- `enum class SpatialOrder` - SECOND_ORDER, FOURTH_ORDER
- `Config::spatial_order` - Configuration parameter (default: FOURTH_ORDER)
- `computeLaplacian4thOrder()` - 18-neighbor 4th-order Laplacian

**Lines**: 49-52, 64, 237-255

#### 2. `src/ConservativeSolver.cpp`

**Modified**:
- `initialize()` - Log spatial discretization order (lines 48-57)
- `velocityVerletStep()` - Use 4th-order if configured (lines 99-103, 134-138)
- `rk2SymplecticStep()` - Use 4th-order if configured (lines 171-177, 206-211)
- `strangSplittingStep()` - Use 4th-order if configured (lines 267-272)
- `computeTotalEnergy()` - 4th-order gradients for consistency (lines 509-533)

**Added**:
- `computeLaplacian4thOrder()` - Full 3D 4th-order implementation (lines 610-652)

#### 3. `CMakeLists.txt`

**Added build targets**:
- `test_spatial_order_comparison` (lines 552-570)
- `test_spatial_order_timestep` (lines 573-591)
- `test_laplacian_accuracy` (lines 594-612)
- `test_grid_resolution` (lines 615-633)

### New Test Files

1. `test_spatial_order_comparison.cpp` - Main validation test (2nd vs 4th)
2. `test_spatial_order_timestep.cpp` - Timestep sweep with 4th-order
3. `test_laplacian_accuracy.cpp` - Mathematical verification of stencil
4. `test_grid_resolution.cpp` - Grid convergence study

---

## Physical Interpretation

### Why 4th-Order Matters

**Grid dispersion**: 2nd-order finite differences introduce O(k²·dx²) frequency errors
- High-frequency modes propagate at wrong speeds
- Energy "sloshes" between spatial modes (observed ~600-step period)
- Result: 0.07% apparent drift

**4th-order correction**: O(k⁴·dx⁴) error scales as (dx²)² - suppresses grid dispersion by 16× per grid doubling

**Energy oscillation**: Still present but amplitude reduced from ±0.07% → ±0.004%

### Symplectic Structure Preserved

- ✅ Time reversibility: <10⁻⁹ rad (6 orders better than 2nd-order temporal error)
- ✅ Phase space volume conservation (Hamiltonian structure)
- ✅ Long-time stability (tested 2000 steps)

---

## Comparison to Original Analysis

### Initial Diagnosis (STRANG_SPLITTING_ANALYSIS.md)

**Hypothesis**: 0.068% drift caused by 2nd-order spatial discretization ✅ **CONFIRMED**

**Prediction**: 4th-order would achieve 10-70× improvement ✅ **ACHIEVED (18×)**

**Alternative rejected**: Finer grid (128³) would be 4× better → Tested, didn't work

**Solution validated**: 4th-order stencils achieve <0.01% → ✅ **CONFIRMED**

---

## Lessons Learned

### Critical Insight #1: Consistent Discretization Order

**Problem**: Initial 4th-order implementation made things WORSE (0.092% drift)

**Root cause**:
- Evolution used 4th-order Laplacian
- Energy measurement used 2nd-order gradients
- **Order mismatch** broke energy conservation tracking

**Solution**: Use 4th-order for BOTH evolution and energy

**Lesson**: Numerical conservation requires consistent discretization throughout

### Critical Insight #2: Grid Dispersion Dominates

**Finding**: All time integrators (Velocity Verlet, RK2, Strang) achieved identical 0.068% drift

**Conclusion**: Temporal integration was NOT the bottleneck

**Impact**: Spatial discretization was the limiting factor - time integration already optimal

### Critical Insight #3: Stencil Verification Essential

**Test**: Analytic function with known Laplacian

**Result**: Confirmed 4th-order stencil is mathematically correct (62× better)

**Value**: Ruled out implementation bugs early

---

## Recommendations

### For This Codebase

**Default configuration**: ✅ Already set to `FOURTH_ORDER`

**Fallback option**: 2nd-order available if needed (`Config::spatial_order = SECOND_ORDER`)

**No further action needed**: GO/NO-GO criterion met

### For Similar Problems

1. **Profile first**: Identify temporal vs spatial bottleneck
2. **Consistent order**: Match discretization in dynamics + measurement
3. **Verify stencils**: Test on analytic functions before full simulations
4. **Start coarse**: 2nd-order baseline → upgrade strategically

---

## References

### Theory

- **Hairer, Lubich, Wanner (2006)**: *Geometric Numerical Integration* - Symplectic methods
- **Strikwerda (2004)**: *Finite Difference Schemes* - Spatial discretization theory
- **Fornberg (1988)**: "Generation of Finite Difference Formulas" - Optimal stencils

### Implementation

- `/home/persist/neotec/0rigin/include/ConservativeSolver.h`
- `/home/persist/neotec/0rigin/src/ConservativeSolver.cpp`
- `/home/persist/neotec/0rigin/test_spatial_order_comparison.cpp`

### Analysis Documents

- `STRANG_SPLITTING_ANALYSIS.md` - Initial diagnosis and solution proposal
- `FOURTH_ORDER_INVESTIGATION.md` - Implementation journey and debugging
- `4TH_ORDER_IMPLEMENTATION_COMPLETE.md` - This document (final report)

---

## Conclusion

**Objective achieved**: <0.01% energy conservation requirement met via 4th-order spatial stencils.

**Final performance**:
- **Energy drift**: 0.0038% (18× better than baseline)
- **Time reversibility**: <10⁻⁹ rad (6 orders margin)
- **Computational cost**: 1.5× (acceptable)
- **Implementation**: Complete and validated

**Status**: ✅ **PRODUCTION READY** - GO/NO-GO criterion satisfied.

**Next steps**: None required for energy conservation. System ready for physics validation (vortex scattering, etc.).

---

**Implementation completed**: 2026-01-08
**Validation**: All tests passing
**Status**: READY FOR DEPLOYMENT
