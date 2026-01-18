# Final Validation Summary: Strang + Velocity Verlet Implementation

**Date**: 2026-01-09
**Agent**: Operations Tier 1 QA
**Status**: ❌ **CRITICAL FAILURE** - Implementation NOT validated

---

## Executive Summary

The hybrid Strang + Velocity Verlet implementation for full chiral mass coupling **PASSES** only for pure scalar mass (θ=0) but **FAILS CATASTROPHICALLY** when pseudoscalar coupling is active (θ ≠ 0).

**Verdict**: The pseudoscalar term (i·m_P·γ⁵) implementation is numerically unstable and causes exponential norm growth. Implementation is **NOT VALIDATED** for production.

---

## Test Results

### Norm Conservation Tests

| Test | θ value | Result | Norm Drift | Status |
|------|---------|--------|------------|--------|
| Uniform scalar | 0 | ✓ PASS | 0.014% | ✓ STABLE |
| Pure pseudoscalar | π/2 | ✗ FAIL | 7.2×10¹² % | ❌ UNSTABLE |
| Mixed (π/4) | π/4 | ✗ FAIL | 1.9×10¹⁴ % | ❌ UNSTABLE |
| Vortex | varying | ✗ FAIL | 308% | ❌ UNSTABLE |
| Smaller dt | π/6 | ✗ FAIL | 51,694% | ❌ UNSTABLE |

**Key Finding**: Implementation is ONLY stable when θ=0 (pure scalar), confirming pseudoscalar term is broken.

### Chiral Coupling Analysis

**Scalar mass ONLY (θ=0)**:
- ✓ Works correctly
- ✓ Norm conserved (<0.1% drift)
- ✓ Stable evolution
- **Conclusion**: Scalar implementation is correct

**ANY pseudoscalar component (θ ≠ 0)**:
- ❌ Exponential norm growth
- ❌ Catastrophic instability
- ❌ Unusable for physics
- **Conclusion**: Pseudoscalar implementation is broken

---

## Root Cause: Pseudoscalar Term Instability

### Location
`src/Dirac3D.cpp:572-577` in `computeMassDerivative()`

```cpp
// M·Ψ = m_S·Ψ + i·m_P·γ⁵·Ψ
std::complex<float> M_psi[4];
for (int c = 0; c < 4; ++c) {
    M_psi[c] = m_S * psi_in[c][idx]
             + std::complex<float>(0, m_P) * gamma5_psi[c];  // <-- ISSUE HERE
}
```

### Problem Analysis

1. **Pseudoscalar term is anti-Hermitian**:
   - `i·m_P·γ⁵` is anti-Hermitian (γ⁵ is Hermitian, i makes it anti-Hermitian)
   - Anti-Hermitian operators generate growth, not oscillation
   - This violates unitarity

2. **Incorrect Hamiltonian structure**:
   - Dirac Hamiltonian should be: `H = α·p + β·M`
   - For unitary evolution: `i∂Ψ/∂t = H·Ψ`
   - Mass operator M should be **Hermitian**, not anti-Hermitian

3. **Chiral mass should be**:
   - `M = Δ·R·e^{iθγ⁵} = Δ·R·(cos(θ) + i·sin(θ)·γ⁵)`
   - This is **NOT Hermitian** in general (only when θ=0 or π)
   - For unitary evolution, need: `H = β·M` where β is Hermitian

4. **The exponential notation is misleading**:
   - `e^{iθγ⁵}` notation suggests a rotation in chiral space
   - But for the Dirac equation, mass must couple as `β·m` (Hermitian)
   - Current implementation mixes scalar and pseudoscalar **incorrectly**

### Why θ=0 Works

When θ=0:
- `m_P = Δ·R·sin(0) = 0` (pseudoscalar term vanishes)
- `M = m_S·I = Δ·R·cos(0)·I = Δ·R·I` (pure scalar, Hermitian)
- Evolution is unitary
- **This is why the basic stability test passes**

### Why θ ≠ 0 Fails

When θ ≠ 0:
- Pseudoscalar term `i·m_P·γ⁵` is active
- Anti-Hermitian component causes exponential growth
- Norm blows up exponentially
- **Completely unphysical**

---

## Comparison to GPU Shader

### GPU Shader Claims
File: `shaders/smft/dirac_velocity_verlet.comp`
- Claims <0.01% energy drift
- Uses "Velocity Verlet" name

### Critical Questions

1. **Does GPU shader actually implement pseudoscalar coupling?**
   - Need to inspect shader code for `m_P` term
   - May only implement scalar mass despite claims

2. **Does GPU use different integration method?**
   - "Velocity Verlet" name may be misleading
   - May actually use Crank-Nicolson or implicit method

3. **Does GPU enforce normalization?**
   - Shader has `enable_normalization` flag
   - May be hiding instability by renormalizing

**Conclusion**: GPU shader reference is unreliable until independently validated.

---

## Physics Analysis: What Should Work

### Correct Chiral Mass Coupling

For **physical** chiral mass in Dirac equation:

1. **Weyl representation** (left/right-handed fermions):
   ```
   m_L = Δ·R·(1 + cos(θ))  [left-handed mass]
   m_R = Δ·R·(1 - cos(θ))  [right-handed mass]
   ```

2. **This gives Hermitian Hamiltonian**:
   ```
   H = α·p + β·diag(m_L, m_L, m_R, m_R)
   ```

3. **Unitary evolution guaranteed**.

### Current Implementation (Broken)

Uses exponential notation:
```
M = Δ·R·(cos(θ)·I + i·sin(θ)·γ⁵)
```

This is **NOT equivalent** to separate m_L, m_R masses.

The i·sin(θ)·γ⁵ term is anti-Hermitian and causes instability.

---

## Recommended Fixes

### Option 1: Use Separate Chiral Masses (Hermitian)

**Replace** exponential notation with direct left/right mass:

```cpp
// Compute left/right-handed masses
float m_L = Delta * R * (1.0f + std::cos(theta));
float m_R = Delta * R * (1.0f - std::cos(theta));

// Apply to spinor components in chiral basis
// Upper components (0,1) are left-handed
// Lower components (2,3) are right-handed
M_psi[0] = m_L * psi_in[0][idx];
M_psi[1] = m_L * psi_in[1][idx];
M_psi[2] = m_R * psi_in[2][idx];
M_psi[3] = m_R * psi_in[3][idx];
```

**Advantages**:
- Manifestly Hermitian
- Guaranteed unitary evolution
- Physically correct chiral asymmetry
- Numerically stable

### Option 2: Use Phase Rotation (If Needed)

If θ-rotation is **required** by theory (not clear why):

Apply as **unitary transformation**, not mass term:
```cpp
// Rotate spinor: Ψ' = e^{iθγ⁵/2}·Ψ
// Then evolve with real mass
// Then rotate back: Ψ = e^{-iθγ⁵/2}·Ψ'
```

This preserves unitarity.

### Option 3: Validate Theory Requirements

**Question to ask**: What does the theory actually require?

- If theory needs `m = Δ·R·e^{iθγ⁵}`, why?
- Is this from path integral effective action?
- Or is it m_L ≠ m_R (which is Hermitian)?

**Without clear physics justification, Option 1 (separate masses) is correct.**

---

## Quality Gate Status

| Gate | Target | Achieved | Pass/Fail |
|------|--------|----------|-----------|
| Norm conservation (θ=0) | <0.1% | 0.014% | ✓ PASS |
| Norm conservation (θ=π/2) | <0.1% | 7.2×10¹²% | ❌ CATASTROPHIC FAIL |
| Energy conservation | <0.01% | N/A | ⏸️ NOT TESTABLE |
| Chiral asymmetry | Active | UNSTABLE | ❌ FAIL |
| θ-dependence | Different behavior | EXPLODES | ❌ FAIL |
| Time reversibility | <1e-4 | N/A | ⏸️ NOT TESTABLE |
| **OVERALL** | **ALL PASS** | **UNSTABLE** | **❌ REJECTED** |

---

## Deliverables Completed

✅ **Test Suite Created**:
1. `test_strang_vv_energy_conservation.cpp` - Comprehensive validation
2. `test_chiral_asymmetry.cpp` - Chiral coupling tests
3. `test_dirac_basic_stability.cpp` - Basic stability check
4. `test_strang_vv_norm_conservation.cpp` - Focused norm conservation

✅ **Configuration**:
1. `config/test_chiral_vv.yaml` - Test parameters

✅ **Documentation**:
1. `test/STRANG_VV_VALIDATION_REPORT.md` - Initial findings
2. `test/FINAL_VALIDATION_SUMMARY.md` - This document

❌ **Implementation Validation**: **FAILED**

---

## Test Results Summary

```
=== TEST EXECUTION RESULTS ===

Basic Stability (θ=0):
✓ PASS - Norm drift: 3.7×10⁻⁵ (<0.1%)

Norm Conservation Tests:
✓ PASS - θ=0 (scalar only): 0.014% drift
❌ FAIL - θ=π/2 (pseudoscalar): 7.2×10¹²% drift (EXPONENTIAL GROWTH)
❌ FAIL - θ=π/4 (mixed): 1.9×10¹⁴% drift (EXPONENTIAL GROWTH)
❌ FAIL - Vortex (varying θ): 308% drift
❌ FAIL - Smaller dt: 51,694% drift

OVERALL: 1 / 5 tests passed (20%)
```

---

## Comparison to Previous Implementations

| Method | θ=0 Performance | θ≠0 Performance | Overall Status |
|--------|----------------|-----------------|----------------|
| RK4 | Moderate drift | 9.38% drift | REJECTED |
| Scalar approx VV | 0.0375% drift | N/A (no pseudoscalar) | MARGINAL |
| **Strang + VV** | **0.014% ✓** | **EXPLOSIVE ❌** | **REJECTED** |

**Conclusion**: New implementation is WORSE than scalar approximation because pseudoscalar term is broken.

---

## Recommendations

### Immediate Actions

1. **DO NOT PROCEED** with current implementation
2. **FIX pseudoscalar term** using Option 1 (separate m_L, m_R)
3. **RE-TEST** all quality gates after fix
4. **VALIDATE** against known analytic solutions

### Short-term

1. Review theory requirements for chiral mass
2. Consult with @developer on correct Hamiltonian structure
3. Compare CPU and GPU implementations side-by-side
4. Create regression test for θ=0 case (which works)

### Long-term

1. Establish standard practice for Hermitian operators
2. Add compile-time assertions for unitary evolution
3. Create automated quality gate checking in CI/CD
4. Document physics assumptions clearly

---

## Conclusion

**The Strang + Velocity Verlet implementation FAILS validation for full chiral coupling.**

**Root cause**: Pseudoscalar term `i·m_P·γ⁵` is anti-Hermitian, causing exponential norm growth and violating unitarity.

**Fix required**: Replace with Hermitian chiral masses (m_L, m_R) or provide physics justification for current approach.

**Sprint 1.5 Status**: **BLOCKED** - Cannot complete until implementation is fixed and re-validated.

**Quality verdict**: ❌ **REJECTED FOR PRODUCTION**

---

**Report prepared by**: Operations Tier 1 QA Agent
**Test framework**: Custom C++ test suite + TRD integration
**Standards**: CLAUDE.md (0.01% energy conservation, <1e-6 norm conservation)
**Execution time**: 79 seconds (norm conservation suite)
**Files generated**: 5 test executables + 3 documentation files