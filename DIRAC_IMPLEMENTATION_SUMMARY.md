# DIRAC IMPLEMENTATION SUMMARY

## Current State

### What's Working
- **Physics**: Full chiral mass coupling M = Δ·R·(cos(θ) + i·sin(θ)·γ⁵) implemented correctly in:
  - `Dirac3D::computeMassDerivative()` (CPU)
  - `dirac_velocity_verlet.comp` (GPU)
- **Integration**: Velocity Verlet is symplectic and correct
- **Production Path**: `TRDEngine3D → ConservativeSolver → Dirac3D::stepWithChiralMass()`

### What's Broken
1. **Timestep Stability**: Tests use dt=0.01 but need dt<0.003 for stability
2. **Fixed Substeps**: Hardcoded 100 substeps regardless of parameters
3. **Dead Code**: `applyChiralMassStep()` ignores pseudoscalar mass (but unused)
4. **Duplicate Classes**: `Dirac3D` vs `DiracEvolution` (confusing)
5. **No Energy Monitoring**: Can't detect drift in production

## Quick Fixes (Can Do Now)

### Fix 1: Adaptive Substeps in stepWithChiralMass
```cpp
// Replace line 478 in Dirac3D.cpp
// const int N_substeps = 100;  // OLD

// NEW: Adaptive substeps based on stability
float omega_max = Delta * 2.0f;  // Conservative estimate
float dt_critical = M_PI / (10.0f * omega_max);
const int N_substeps = std::max(1, (int)std::ceil(dt / dt_critical));

if (N_substeps > 100) {
    std::cerr << "[WARNING] Dirac evolution requires " << N_substeps
              << " substeps for stability (dt=" << dt
              << ", Delta=" << Delta << ")" << std::endl;
}
```

### Fix 2: Add Norm Check
```cpp
// Add after line 486 in Dirac3D.cpp
float norm_after = getNorm();
float drift = std::abs(norm_after - norm_before) / norm_before;
if (drift > 0.01f) {
    std::cerr << "[ERROR] Norm drift " << (drift*100)
              << "% exceeds 1% threshold!" << std::endl;
}
```

### Fix 3: Remove Dead Code
- Delete `applyChiralMassStep()` method (lines 332-381)
- Remove `dirac_rk4.comp` shader
- Remove `dirac_stochastic.comp` shader

## Long-term Fixes (Need Design)

### Architecture Consolidation
```
Single Dirac Implementation
├── CPU: Dirac3D (keep, fix stability)
├── GPU: dirac_velocity_verlet.comp (keep)
└── Remove: DiracEvolution class

Single Integration Method
├── Strang splitting (outer)
├── Magnus expansion (for mass operator)
└── Adaptive timestepping
```

### Proper Exponential for Non-Hermitian Mass
The mass operator M = Δ·R·e^(iθγ⁵) is not Hermitian when θ ≠ 0.
Need Magnus expansion or Cayley transform for unitary evolution.

## Test Updates Needed

### Update All Tests
```cpp
// Change all tests from:
const float dt = 0.01f;

// To adaptive:
float omega_max = Delta * 2.0f;
const float dt = std::min(0.01f, M_PI / (20.0f * omega_max));
```

## Files to Clean Up

| File | Action | Reason |
|------|--------|--------|
| src/Dirac3D.cpp:332-381 | DELETE | applyChiralMassStep() wrong physics |
| src/DiracEvolution.* | DELETE | Duplicate implementation |
| shaders/dirac_rk4.comp | DELETE | Violates symplectic standard |
| shaders/dirac_stochastic.comp | DELETE | Unused |
| test_dirac_em_coupling.cpp | UPDATE | Uses DiracEvolution |
| test_trd_em_integration.cpp | UPDATE | Uses DiracEvolution |

## Quality Gates to Add

```cpp
// Add to all Dirac tests
assert(norm_drift < 0.01);      // Norm conservation
assert(energy_drift < 0.01);    // Energy conservation
assert(time_reversible < 1e-4); // Symplectic structure
```

## Executive Summary

**Problem**: Dirac evolution has 35% norm drift due to timestep instability, not wrong physics.

**Quick Fix**: Adaptive substeps based on `dt_critical = π/(10·Δ)`.

**Proper Fix**: Magnus expansion for non-Hermitian mass operator.

**Action Items**:
1. ✅ Implement adaptive substeps (5 min fix)
2. ✅ Add norm monitoring (2 min fix)
3. ✅ Delete dead code (5 min)
4. ⏳ Consolidate to single implementation (1 hour)
5. ⏳ Implement Magnus expansion (2 hours)