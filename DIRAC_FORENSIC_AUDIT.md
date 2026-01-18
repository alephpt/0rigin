# DIRAC IMPLEMENTATION FORENSIC AUDIT

**Date**: 2026-01-09
**Auditor**: Operations Tier 1 Agent

## EXECUTIVE SUMMARY

CRITICAL FINDINGS:
1. **DUPLICATE IMPLEMENTATIONS**: Two completely separate Dirac implementations (Dirac3D vs DiracEvolution)
2. **INCONSISTENT PHYSICS**: Different formulas for chiral mass coupling across implementations
3. **NUMERICAL INSTABILITY**: 35% norm drift in 1000 steps (test_dirac_vacuum_chiral_coupling_simple)
4. **INCOMPLETE IMPLEMENTATIONS**: applyChiralMassStep uses only scalar part, ignores pseudoscalar
5. **UNTESTED CODE PATHS**: Velocity Verlet path in Dirac3D never called in production

## 1. CODE INVENTORY

### CPU Implementations
| File | Status | Purpose |
|------|--------|---------|
| include/Dirac3D.h | ACTIVE | 3D Dirac solver header |
| src/Dirac3D.cpp | ACTIVE | 3D implementation with TWO evolution methods |
| src/DiracEvolution.h | ACTIVE | 2D Dirac solver header (separate implementation!) |
| src/DiracEvolution.cpp | ACTIVE | 2D implementation (different from Dirac3D) |
| src/ConservativeSolver.cpp | ACTIVE | Calls Dirac3D::stepWithChiralMass |
| src/TRDEngine3D.cpp | ACTIVE | Calls ConservativeSolver::evolveDirac |

### GPU Implementations
| File | Status | Purpose |
|------|--------|---------|
| shaders/smft/dirac_velocity_verlet.comp | COMPILED | Full VV implementation with both m_S and m_P |
| shaders/smft/dirac_rk4.comp | COMPILED | RK4 implementation (deprecated per standards) |
| shaders/smft/dirac_stochastic.comp | COMPILED | Stochastic evolution (not used) |

### Test Files
- **Using Dirac3D**: test_dirac_vacuum_chiral_coupling_simple.cpp, test_dirac_vacuum_chiral_coupling.cpp, test_dirac3d_free.cpp, etc.
- **Using DiracEvolution**: test_dirac_em_coupling.cpp, test_trd_em_integration.cpp, test_geodesic_verification.cpp

## 2. METHODS INVENTORY

### Dirac3D Methods (src/Dirac3D.cpp)

| Method | Line | Integrator | m_S | m_P | Status | Called By |
|--------|------|------------|-----|-----|--------|-----------|
| step() | 383 | Strang | YES | NO | ACTIVE | Tests only |
| stepWithChiralMass() | 458 | Strang+VV | YES | YES | ACTIVE | ConservativeSolver |
| applyMassStep() | 306 | Exponential | YES | NO | ACTIVE | step() |
| applyChiralMassStep() | 332 | Exponential | YES (cos θ only) | NO! | **INCOMPLETE** | Not called! |
| applyMassVelocityVerlet() | 489 | Velocity Verlet | YES | YES | ACTIVE | stepWithChiralMass() |
| computeMassDerivative() | 536 | N/A (derivative) | YES | YES | ACTIVE | applyMassVelocityVerlet() |

### Key Findings:
1. **applyChiralMassStep (line 332)**: Only uses cos(θ) term, completely ignores sin(θ) pseudoscalar mass!
2. **stepWithChiralMass (line 458)**: Uses VV sub-stepping (100 substeps) which includes FULL chiral mass
3. **Inconsistency**: Two mass application methods with DIFFERENT physics

## 3. ACTUAL EXECUTION PATH

```
TRDEngine3D::runSimulation()
    ↓
ConservativeSolver::evolveDirac() [line 326]
    ↓
Dirac3D::stepWithChiralMass() [line 458]
    ↓
applyKineticHalfStep() [Strang T/2]
    ↓
for (100 substeps) {
    applyMassVelocityVerlet() [line 489]
        ↓
    computeMassDerivative() [line 536] - FULL m_S + i·m_P·γ⁵
}
    ↓
applyKineticHalfStep() [Strang T/2]
```

**CRITICAL**: Production code uses VV with FULL chiral mass (both scalar and pseudoscalar)

## 4. FORMULA DISCREPANCIES

### applyChiralMassStep (line 369) - INCORRECT
```cpp
// ONLY uses scalar part:
const float m_eff = Delta * R * std::cos(theta);  // Missing sin(θ) term!
```

### computeMassDerivative (line 561) - CORRECT
```cpp
// FULL chiral mass:
float m_S = Delta * R * std::cos(theta);  // Scalar mass
float m_P = Delta * R * std::sin(theta);  // Pseudoscalar mass
M_psi[c] = m_S * psi_in[c][idx] + std::complex<float>(0, m_P) * gamma5_psi[c];
```

### GPU dirac_velocity_verlet.comp (line 218) - CORRECT
```glsl
float m_S = params.Delta * R * cos(params.chiral_angle);  // Scalar mass
float m_P = params.Delta * R * sin(params.chiral_angle);  // Pseudoscalar mass
```

## 5. DUPLICATES IDENTIFIED

1. **Dirac3D vs DiracEvolution**: Two completely separate implementations
   - Dirac3D: 3D solver
   - DiracEvolution: 2D solver with EM coupling
   - NO CODE SHARING!

2. **Mass Evolution Methods**:
   - applyMassStep: Simple exponential for scalar mass only
   - applyChiralMassStep: Incomplete exponential (scalar only, ignores pseudoscalar)
   - applyMassVelocityVerlet: Full VV with both scalar and pseudoscalar

3. **GPU Shaders**:
   - dirac_rk4.comp: RK4 (should be removed per standards)
   - dirac_velocity_verlet.comp: VV implementation
   - dirac_stochastic.comp: Unused stochastic evolution

## 6. ROOT CAUSE OF FAILURES

### Test: test_dirac_vacuum_chiral_coupling_simple
**Result**: 35% norm drift in 1000 steps
**Root Cause**:
1. Velocity Verlet with dt=0.01 is unstable for rapidly oscillating chiral phase
2. 100 substeps may not be enough for convergence
3. No adaptive timestep control

### Energy Conservation Issues
**Problem**: Tests show energy drift > 0.01%
**Root Cause**:
1. Mixing of integrators (Strang for kinetic, VV for mass)
2. Fixed timestep without stability analysis
3. No symplectic structure preservation check

## 7. INCONSISTENCIES FOUND

1. **applyChiralMassStep**: Claims to handle chiral mass but only uses scalar part
2. **Comment vs Code**: Line 339-381 has extensive comment about unitarity issues, then uses scalar-only approximation
3. **GPU vs CPU**: GPU shader implements full physics, CPU applyChiralMassStep doesn't
4. **Test Coverage**: Some tests use DiracEvolution (2D), others use Dirac3D (3D)

## 8. UNUSED/DEAD CODE

1. **applyChiralMassStep()**: Never called in production (line 332-381)
2. **dirac_rk4.comp**: Compiled but not used (violates standards)
3. **dirac_stochastic.comp**: Compiled but not used
4. **DiracEvolution class**: Only used in 3 tests, not in production

## 9. CONSOLIDATION PLAN

### Immediate Actions Required

1. **REMOVE DEAD CODE**:
   - Delete applyChiralMassStep() - it's incorrect and unused
   - Remove dirac_rk4.comp - violates symplectic standards
   - Remove dirac_stochastic.comp - unused

2. **FIX VELOCITY VERLET STABILITY**:
   - Implement adaptive timestep based on max(|m_S|, |m_P|)
   - Add CFL-like condition: dt < π / max_frequency
   - Consider implicit midpoint for mass evolution

3. **UNIFY IMPLEMENTATIONS**:
   - Decide: Keep Dirac3D or DiracEvolution (not both!)
   - Standardize on single integrator (recommend: pure Strang with exponential)
   - Share code between CPU and GPU

4. **CORRECT THE PHYSICS**:
   - Ensure ALL paths include both m_S and m_P
   - Fix matrix exponential for non-Hermitian mass operator
   - Add energy conservation monitoring

5. **TEST CONSOLIDATION**:
   - All tests should use production code path
   - Remove test-specific implementations
   - Add energy conservation criterion < 0.01%

### Recommended Architecture

```
TRDEngine3D
    ↓
ConservativeSolver
    ↓
Dirac3D (single implementation)
    ↓
stepWithSymplectic() - Pure Strang splitting with:
    - Kinetic: exp(-iα·k dt/2) via FFT
    - Mass: Magnus expansion for exp(-iβM dt)
    - Adaptive substeps based on ||M||
```

## 10. CRITICAL BUGS

1. **BUG #1**: applyChiralMassStep ignores pseudoscalar mass (line 369)
   - Impact: Wrong physics if this method were used
   - Fix: Delete the method (it's unused)

2. **BUG #2**: Fixed 100 substeps regardless of parameters (line 478)
   - Impact: Instability for large Delta or rapidly varying theta
   - Fix: Adaptive substeps based on max(Δ·R·|∇θ|)

3. **BUG #3**: No energy conservation check in production
   - Impact: Can't detect numerical drift
   - Fix: Add energy computation and threshold check

## CONCLUSION

The Dirac implementation has serious architectural issues:
- Multiple duplicate implementations with different physics
- Incomplete chiral mass handling in unused code paths
- Numerical instability from fixed timestep
- No energy conservation monitoring

The production code path (stepWithChiralMass → VV) includes correct physics but suffers from stability issues. Immediate consolidation and stability fixes are required to meet the 0.01% energy conservation standard.

**Recommendation**: Complete rewrite using single, well-tested symplectic integrator with adaptive timestepping and energy monitoring.

## VALIDATION RESULTS

### Debug Test Results (test_dirac_debug.cpp)
- **dt = 0.001, Delta = 0.1**: Norm drift = 0.977% ✓ PASS
- **dt = 0.01, Delta = 0.5**: Norm drift = 35% ✗ FAIL

### Root Cause Confirmed
The issue is **timestep stability**, not the integrator itself. The Velocity Verlet integrator is correct but requires:
```
dt < 2π / (20 * ω_max)
where ω_max = Delta * max(R) * max(1 ± cos(θ))
```

For Delta = 0.5, R = 1, this gives dt_critical ≈ 0.003, but tests use dt = 0.01.

### Immediate Fix Required
1. Add adaptive timestep selection based on CFL-like condition
2. Increase substeps from 100 to ceil(dt / dt_critical)
3. Add runtime warning when dt exceeds stability limit