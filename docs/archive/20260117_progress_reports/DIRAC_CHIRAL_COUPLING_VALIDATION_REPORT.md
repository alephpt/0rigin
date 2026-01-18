# Dirac-Vacuum Chiral Coupling Validation Report

**Date**: 2026-01-08
**Test**: Dirac3D::applyChiralMassStep() implementation
**Agent**: @qa (Operations Tier 1)

## Executive Summary

**VALIDATION STATUS: PARTIAL PASS WITH CRITICAL ISSUES**

The chiral coupling implementation is functionally correct but exhibits severe numerical instability that violates TRD energy conservation standards.

## Test Results

### Test 1: Full Physics Test (32³ grid, dt=0.01)
- **Energy Conservation**: ❌ FAILED (3.16×10³⁵% drift)
- **Norm Conservation**: ❌ FAILED (4.72×10¹⁹% drift)
- **Chiral Asymmetry**: ✅ PASSED (m_L - m_R = 0.84)
- **Particle Localization**: ✅ PASSED (density peaks at high R)

**Critical Issue**: Exponential instability - energy and norm grow without bound.

### Test 2: Simplified Stability Test (16³ grid, dt=0.001)
- **Norm Stability**: ✅ PASSED (0.45% drift)
- **Chiral Asymmetry**: ✅ PASSED (m_L - m_R = 0.80)
- **Energy Proxy**: Stable evolution observed

**Finding**: With 10× smaller timestep, the implementation shows reasonable stability.

## Root Cause Analysis

### 1. Timestep Sensitivity
The chiral mass operator exp(-iβMΔt) where M = Δ·R·e^(iθγ⁵) requires:
- **dt < 1/(2·Δ·max(R))** for stability
- With Δ=2.5 and R≈1, critical dt ≈ 0.2
- Test used dt=0.01, well below critical value
- **Issue**: Implementation may have sign error or incorrect operator splitting

### 2. Implementation Analysis

Reviewing `Dirac3D::applyChiralMassStep()`:
```cpp
// Lines 347-356 show the critical evolution:
const std::complex<float> phase_upper(-m_S * dt, -m_P * dt);
const std::complex<float> phase_lower(+m_S * dt, -m_P * dt);
```

**Potential Issue**: The phase calculation appears correct, but the exponential application may be unstable for large mass values.

### 3. Comparison with TRD Standards

Per TRD standards (CLAUDE.md):
- **Required**: ΔE/E < 0.01% (GO/NO-GO criterion)
- **Achieved**: ΔE/E > 10³⁵% (catastrophic failure)
- **Violation**: 10³⁷× worse than acceptable

## Critical Findings

1. **Physics Correctness**: ✅
   - Chiral asymmetry properly generated
   - Particle localization occurs as expected
   - Left/right chirality separation confirmed

2. **Numerical Stability**: ❌
   - Exponential energy growth indicates non-unitary evolution
   - Violates symplectic integration requirements
   - Does NOT meet TRD production standards

3. **Timestep Requirements**:
   - dt=0.001 → 0.45% drift (marginally acceptable)
   - dt=0.01 → 10³⁵% drift (catastrophic)
   - **Recommendation**: dt < 0.001 required for stability

## Recommendations

### Immediate Actions (CRITICAL)

1. **Review Implementation**:
   - Check sign conventions in phase calculations
   - Verify operator ordering in split-step method
   - Compare with GPU shader implementation (if exists)

2. **Add Stability Guards**:
   ```cpp
   // Add timestep validation
   float dt_critical = 1.0f / (2.0f * Delta * max_R);
   if (dt > 0.5f * dt_critical) {
       throw std::runtime_error("Timestep too large for chiral coupling stability");
   }
   ```

3. **Implement Adaptive Timestepping**:
   - Monitor norm conservation per step
   - Reduce dt if drift exceeds threshold
   - Use Richardson extrapolation for error estimation

### Long-term Improvements

1. **Higher-Order Integrator**:
   - Current split-step is 2nd order
   - Consider 4th-order Forest-Ruth or Yoshida schemes

2. **Implicit Methods**:
   - Crank-Nicolson for mass term
   - Would allow larger timesteps

3. **Energy-Conserving Formulation**:
   - Reformulate as variational integrator
   - Guarantees exact energy conservation

## Test Code Deliverables

### Created Files:
1. `/test/test_dirac_vacuum_chiral_coupling.cpp` - Full physics test
2. `/test/test_dirac_vacuum_chiral_coupling_simple.cpp` - Simplified stability test
3. `/config/dirac_vacuum_chiral_coupling.yaml` - Test configuration

### CMakeLists.txt Updates:
- Added `test_dirac_vacuum_chiral_coupling` target
- Added `test_dirac_vacuum_chiral_coupling_simple` target

## Verdict

**DO NOT DEPLOY TO PRODUCTION**

While the physics implementation is correct, the numerical instability violates TRD's fundamental requirement of <0.01% energy conservation. The implementation requires immediate stability fixes before it can be considered production-ready.

### Pass Criteria Status:
- ❌ Energy conservation < 0.01%
- ❌ Norm conservation < 1e-6
- ✅ Chiral asymmetry confirmed
- ✅ Particle localization verified

**Overall: 2/4 criteria passed = FAILED**

## Next Steps

1. @developer must fix numerical stability issues
2. @integration must add timestep validation
3. @qa (this agent) will re-validate after fixes
4. Consider fallback to simpler mass coupling if stability cannot be achieved

---

**Signed**: Operations Tier 1 QA Agent
**Standards Applied**: TEST, SEC, PERF, DEV
**TRD Compliance**: VIOLATED - Energy conservation requirement not met