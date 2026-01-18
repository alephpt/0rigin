# Strang + Velocity Verlet Implementation Validation Report

**Date**: 2026-01-09
**Test Suite**: test_strang_vv_energy_conservation.cpp
**Configuration**: config/test_chiral_vv.yaml
**Status**: ❌ CRITICAL FAILURE - UNSTABLE IMPLEMENTATION

---

## Executive Summary

The hybrid Strang + Velocity Verlet implementation for full chiral mass coupling **FAILS** all quality gates due to catastrophic numerical instability. Energy conservation violates the GO/NO-GO criterion by 5 orders of magnitude.

**Verdict**: Implementation is NOT validated for production. Root cause analysis required.

---

## Test Results

### Test 1: Energy Conservation ❌ FAIL

| Metric | Target | Achieved | Status |
|--------|--------|----------|--------|
| Energy drift (ΔE/E) | <0.01% | >150,000% | ❌ FAIL |
| Norm conservation | <1e-6 | >54 | ❌ FAIL |
| Time to failure | N/A | ~400 steps | ❌ UNSTABLE |

**Energy Evolution**:
```
Step 0:   E = 0.094580  (initial)
Step 100: E = 0.107565  (ΔE/E = 13.7%)
Step 200: E = 0.342900  (ΔE/E = 263%)
Step 300: E = 4.692103  (ΔE/E = 4,861%)
Step 400: E = 142.362   (ΔE/E = 150,420%)
```

**Observation**: Exponential energy growth indicates numerical instability, not physical behavior.

### Test 2: θ-Dependence Verification ⏸️ NOT EXECUTED

Test aborted due to Test 1 failure. Cannot verify pseudoscalar term until stability is achieved.

### Test 3: Chiral Asymmetry Detection ⏸️ NOT EXECUTED

Test aborted. Requires stable energy conservation first.

### Test 4: Time Reversibility ⏸️ NOT EXECUTED

Test aborted. Symplectic property cannot be validated with unstable integrator.

---

## Root Cause Analysis

### Issue 1: Excessive Sub-stepping

**Location**: `src/Dirac3D.cpp:478`

```cpp
const int N_substeps = 100;  // Oppenheimer refinement factor
```

**Problem**: The implementation performs **100 Velocity Verlet iterations** per outer time step `dt`. With `dt = 0.01`, the effective sub-step is `dt_sub = 0.0001`.

**Impact**:
- Each full step performs 100 mass operator applications
- Accumulated numerical errors compound exponentially
- Total iterations for 1000 steps: **100,000 VV sub-steps**

**Expected behavior**: Oppenheimer refinement is for handling rapid phase variations, but 100 substeps is excessive for `dt = 0.01`.

### Issue 2: Velocity Verlet Implementation

**Location**: `src/Dirac3D.cpp:489-534`

The Velocity Verlet implementation uses the pattern:
```cpp
k1 = dΨ/dt at t
Ψ_half = Ψ + (dt/2)·k1
k2 = dΨ/dt at t+dt/2
Ψ_new = Ψ + dt·k2
```

**Concern**: This is NOT standard Velocity Verlet. Standard VV is:
```
Ψ_half = Ψ + (dt/2)·k1
k2 = dΨ/dt(Ψ_half)
Ψ_new = Ψ + dt·k2
```

But for the Schrödinger-type equation (i∂Ψ/∂t = H·Ψ), the correct symplectic integrator is **Strang splitting** or **explicit midpoint**, NOT this pattern.

### Issue 3: Mass Operator Sign

**Location**: `src/Dirac3D.cpp:592`

```cpp
dpsi_dt[c][idx] = std::complex<float>(0, -1) * beta_M_psi[c];
```

The `-i` factor in `dΨ/dt = -i·β·M·Ψ` ensures unitary evolution. However, when combined with Velocity Verlet (designed for real-valued second-order ODEs), this may cause instability.

**Proper approach**: For unitary evolution of `i∂Ψ/∂t = H·Ψ`, use:
- Crank-Nicolson (implicit, unconditionally stable)
- Strang splitting (T-M-T with explicit operators)
- NOT Velocity Verlet (designed for Hamiltonian mechanics, not unitary QM)

---

## Comparison to Reference Implementation

### GPU Shader (dirac_velocity_verlet.comp)

**Claims**: <0.01% energy drift

**Key differences**:
1. **Shader uses SINGLE sub-step per time step** (no 100× sub-stepping)
2. **Different integration pattern** (may be using Crank-Nicolson or implicit method)
3. **Proper normalization enforcement** (`enable_normalization` flag)

**Conclusion**: CPU implementation does NOT match the validated GPU shader.

---

## Physics Validation

### Full Chiral Coupling Verification ❓ UNKNOWN

**Implementation includes**:
```cpp
float m_S = Delta * R * std::cos(theta);  // Scalar mass
float m_P = Delta * R * std::sin(theta);  // Pseudoscalar mass

M_psi[c] = m_S * psi_in[c][idx] + std::complex<float>(0, m_P) * gamma5_psi[c];
```

**Status**: ✓ Code correctly implements both scalar and pseudoscalar terms.

**However**: Cannot validate physically until numerical stability is achieved.

### Chiral Representation Verification ✓ CORRECT

**Implementation**:
```cpp
std::complex<float> gamma5_psi[4] = {
    psi_in[0][idx],      // γ⁵ acts as +1 on upper components
    psi_in[1][idx],
    -psi_in[2][idx],     // γ⁵ acts as -1 on lower components
    -psi_in[3][idx]
};
```

**Status**: ✓ Correct chiral basis where γ⁵ = diag(1, 1, -1, -1)

---

## Comparison to Previous Methods

| Method | Energy Drift | Status | Notes |
|--------|--------------|--------|-------|
| RK4 | 9.38% | REJECTED | High drift, not symplectic |
| Scalar approx VV | 0.0375% | MARGINAL | Only scalar mass, no pseudoscalar |
| **Strang + VV (this)** | **>150,000%** | **REJECTED** | **Catastrophic instability** |
| Target | <0.01% | N/A | GO/NO-GO criterion |

**Regression**: New implementation is **4 million times worse** than RK4 and **400 million times worse** than scalar approximation.

---

## Recommended Fixes

### Priority 1: Remove Excessive Sub-stepping

**Action**: Change `N_substeps` from 100 to **1** and test stability.

**Rationale**:
- Outer Strang splitting already handles operator decomposition
- Inner VV should operate at same time scale as outer loop
- 100× sub-stepping is overkill for `dt = 0.01`

**Test**: Recompile and run test_strang_vv_energy_conservation

### Priority 2: Validate Integration Pattern

**Action**: Review Velocity Verlet implementation against standard symplectic integrators.

**Consider**:
- Is VV the right choice for unitary evolution (i∂Ψ/∂t = H·Ψ)?
- Should we use **Strang splitting alone** for both T and M operators?
- GPU shader may use different method despite name

### Priority 3: Add Normalization

**Action**: Implement periodic normalization as done in GPU shader.

**Code**:
```cpp
if (step % 10 == 0) {
    float norm = getNorm();
    if (norm > 0) {
        for (int c = 0; c < 4; ++c) {
            for (uint32_t i = 0; i < _N_total; ++i) {
                _psi[c][i] /= std::sqrt(norm);
            }
        }
    }
}
```

### Priority 4: Time Step Stability Analysis

**Action**: Implement CFL condition check for Dirac equation.

**Typical constraint**: `dt < dx / c` where `c` is effective wave speed.

For current setup: `dx = 1.0`, so `dt < 1.0`. Current `dt = 0.01` should be stable.

---

## Quality Gate Status

| Gate | Target | Achieved | Pass/Fail |
|------|--------|----------|-----------|
| Energy conservation | <0.01% | >150,000% | ❌ FAIL |
| Norm conservation | <1e-6 | >54 | ❌ FAIL |
| Chiral asymmetry | Detectable | N/A | ⏸️ NOT TESTED |
| θ-dependence | Different for θ=0,π/4,π/2 | N/A | ⏸️ NOT TESTED |
| Time reversibility | <1e-4 | N/A | ⏸️ NOT TESTED |
| **OVERALL** | **ALL PASS** | **ALL FAIL** | **❌ REJECTED** |

---

## Deliverables Status

- ✅ Test suite created: test_strang_vv_energy_conservation.cpp
- ✅ Chiral asymmetry test: test_chiral_asymmetry.cpp
- ✅ YAML configuration: config/test_chiral_vv.yaml
- ❌ Implementation validated: **FAILED - UNSTABLE**

---

## Next Steps

1. **IMMEDIATE**: Fix N_substeps parameter (reduce from 100 to 1)
2. **SHORT-TERM**: Review and correct Velocity Verlet implementation
3. **VALIDATION**: Re-run full test suite after fixes
4. **COMPARISON**: Test against GPU shader to understand discrepancies
5. **DOCUMENTATION**: Update implementation notes in Dirac3D.h

---

## Conclusion

**The hybrid Strang + Velocity Verlet implementation as currently written is NUMERICALLY UNSTABLE and unsuitable for production.**

Root cause appears to be excessive sub-stepping (100× per outer step) combined with potential issues in the Velocity Verlet pattern for unitary quantum evolution.

**Recommendation**:
1. Fix N_substeps
2. Validate against known stable integrators (RK2, Strang)
3. Re-test against quality gates
4. Only proceed to Sprint 1.5 completion after achieving <0.01% energy drift

**Current Status**: Sprint 1.5 **BLOCKED** pending implementation fixes.

---

**Report prepared by**: Operations Tier 1 QA Agent
**Test framework**: TRD test harness (./trd --test)
**Standards**: CLAUDE.md §3 (0.01% energy conservation GO/NO-GO criterion)