# Operator Splitting Bug Fix - Summary

**Date:** 2025-12-18
**Issue:** Operator splitting parameter N was not being used in evolution loop
**Status:** ✓ FIXED

---

## Critical Bug Identified

### Root Cause
The test runner passed N to the test loop but **never passed it to the evolution function**:

```cpp
// Bug in SMFTTestRunner.cpp:164
for (int N : _config.operator_splitting.substep_ratios) {
    // N is defined here...
    for (int step = 0; step < total_steps; ++step) {
        // ...but NOT used here!
        _engine->stepWithDirac(_config.physics.dt, _config.physics.coupling);
    }
}
```

**Consequence:** All three tests (N=1, 10, 100) ran identically, producing perfect bit-for-bit agreement. This was incorrectly interpreted as "perfect convergence" but was actually a test infrastructure bug.

---

## Fixes Applied

### 1. Added substep_ratio Parameter to stepWithDirac()

**File:** `src/SMFTEngine.h`
```cpp
// Before:
void stepWithDirac(float dt, float lambda_coupling);

// After:
void stepWithDirac(float dt, float lambda_coupling, int substep_ratio = 1,
                   float K = 1.0f, float damping = 0.1f);
```

### 2. Implemented N-fold Kuramoto Substep Loop

**File:** `src/SMFTEngine.cpp`
```cpp
void SMFTEngine::stepWithDirac(float dt, float lambda_coupling, int substep_ratio,
                               float K, float damping) {
    // Operator splitting: Execute N Kuramoto steps per Dirac step
    float substep_dt = dt / static_cast<float>(substep_ratio);

    for (int n = 0; n < substep_ratio; ++n) {
        // Step Kuramoto dynamics (fast timescale)
        step(substep_dt, K, damping);
    }

    // Get current synchronization field for mass coupling
    std::vector<float> R_field = getSyncField();

    // Compute mass field: m(x,y) = Δ·R(x,y)
    std::vector<float> mass_field(_Nx * _Ny);
    for (uint32_t i = 0; i < _Nx * _Ny; i++) {
        mass_field[i] = _Delta * R_field[i];
    }

    // Split-operator evolution step (unitary, preserves norm exactly)
    _dirac_evolution->step(mass_field, dt);

    // ... feedback logic ...
}
```

### 3. Updated Test Runner to Pass N

**File:** `src/simulations/SMFTTestRunner.cpp`
```cpp
// Before:
_engine->stepWithDirac(_config.physics.dt, _config.physics.coupling);

// After:
_engine->stepWithDirac(_config.physics.dt, _config.physics.coupling, N,
                       _config.physics.K, _config.physics.damping);
```

---

## Verification

### Observable Performance Impact
- **N=1:** 100 Dirac steps = 100 Kuramoto steps (baseline)
- **N=10:** 100 Dirac steps = 1,000 Kuramoto steps (10× slower) ✓
- **N=100:** 100 Dirac steps = 10,000 Kuramoto steps (100× slower) ✓

The dramatic slowdown for N=10, 100 confirms the substep loop is executing correctly.

### Why Previous Tests Showed Perfect Convergence

**Initial conditions in Phase 1 validation:**
- Kuramoto phases: All uniform (θ = 0)
- Natural frequencies: All zero (ω = 0)
- Result: R(x,y) = 0 everywhere → m(x) = Δ·0 = 0

**Physical interpretation:**
- With uniform phases and zero frequencies, Kuramoto has NO dynamics
- R_avg remains exactly 0 for all time
- Dirac evolves as a free particle (m=0) regardless of N
- **This is physically correct but a bad test case**

### Why Born-Oppenheimer Still Valid

Despite the bug, the physical interpretation remains true:
- In the limit τ_Kuramoto ≪ τ_Dirac, N→∞ is unnecessary
- For uniform/near-uniform phase distributions, Kuramoto equilibrates instantly
- N=1 is computationally optimal when Born-Oppenheimer holds

However, the **test was invalid** because it never actually tested different N values.

---

## Next Steps

### 1. Rerun Phase 1 Validation with Fixed Code
```bash
bin/smft --test config/phase1_full_validation.yaml
```

Expected outcomes with non-uniform initial conditions:
- N=1, 10, 100 should produce **different** intermediate states
- But final observables should **converge** (< 5% difference)
- Convergence validates Born-Oppenheimer approximation

### 2. Create Coupled Dynamics Test

**New config:** `config/coupled_validation.yaml`
```yaml
initial_conditions:
  kuramoto:
    phase_distribution: "random"  # Non-uniform phases
    omega_distribution: "gaussian"  # Non-zero frequencies
    omega_mean: 0.0
    omega_std: 0.1

physics:
  K: 1.0  # Enable Kuramoto coupling
  coupling: 0.5  # Strong Ψ-θ feedback
```

This will test:
- True operator splitting convergence
- Kuramoto-Dirac coupling dynamics
- GPU-CPU data transfer (R-field readback)

### 3. Generate Complete Validation Report

Include:
- Quantitative convergence metrics vs N
- Timing data (confirm 10×, 100× slowdown)
- R-field evolution plots (not just R=0)
- Coupled wavepacket-sync dynamics

---

## Status

✓ **Bug Fixed:** Operator splitting now functional
⚠️ **Validation Incomplete:** Need to rerun tests with correct implementation
⚠️ **Plots Outdated:** Current visualizations show buggy results (perfect convergence)

**Next Action:** Rerun Phase 1 validation and regenerate plots showing actual convergence behavior.
