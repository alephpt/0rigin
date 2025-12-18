# Operator Splitting Critical Issue - GPU Failure

**Date:** 2025-12-18
**Status:** 🚨 BLOCKED - GPU compute broken

---

## Summary

After fixing the operator splitting bug (N parameter not being used), we discovered that **all tests still produce identical results** because:

**ROOT CAUSE: GPU Kuramoto compute is completely broken**
- `step()` fails on every call with "Failed to begin compute batch"
- No CPU fallback implemented → function returns immediately
- Kuramoto phases never evolve (R = 0 always)
- Dirac evolves as free particle regardless of N

---

## Evidence

### 1. GPU Failure Logs
```
[stepWithDirac] First call: substep_ratio=1, substep_dt=0.01
[SMFTEngine] Failed to begin compute batch  ← EVERY SINGLE STEP
[SMFTEngine] Uploaded data to GPU
  Step 0/500 | R_avg = 0 | norm = 1
[SMFTEngine] Failed to begin compute batch
[SMFTEngine] Uploaded data to GPU
[SMFTEngine] Failed to begin compute batch
...
```

### 2. Identical Results
```bash
$ diff output/coupled_dynamics_validation/N_1/observables.csv \
       output/coupled_dynamics_validation/N_100/observables.csv
# No differences - files are IDENTICAL
```

Convergence: [0.000, 0.000, 0.000] ← Perfect identity, NOT convergence

### 3. R-field Always Zero
```
R_avg = 0  (all timesteps, all N values)
R_max = 0
R_min = 0
R_var = 0
```

**This proves Kuramoto never ran** - even with vortex initial conditions, R should be non-zero.

---

## Code Analysis

### SMFTEngine.cpp:230-245
```cpp
void SMFTEngine::step(float dt, float K, float damping) {
    // ... pipeline checks ...

    // Begin compute batch
    if (!_compute->beginBatch()) {
        std::cerr << "[SMFTEngine] Failed to begin compute batch" << std::endl;
        return;  // ← SILENTLY FAILS, NO CPU FALLBACK
    }

    // GPU compute code (never executed) ...
}
```

**Impact:**
- N=1: Calls `step()` once → fails → Dirac evolves with m=0
- N=10: Calls `step()` 10× → fails 10× → Dirac STILL evolves with m=0
- N=100: Calls `step()` 100× → fails 100× → Dirac STILL evolves with m=0

All three are **functionally identical** because Kuramoto never runs.

---

## Why This Invalidates Previous Validation

### Phase 1 Validation (before N fix)
- **Claimed:** "Perfect convergence validates Born-Oppenheimer"
- **Reality:** N parameter wasn't used → all tests ran identical code
- **Conclusion:** Test infrastructure bug

### Coupled Dynamics Validation (after N fix)
- **Claimed:** "Vortex defect creates spatial R gradient"
- **Reality:** GPU fails → no Kuramoto evolution → R=0 always
- **Conclusion:** GPU compute bug

**BOTH results are artifacts of bugs, NOT physical validation.**

---

## Required Fixes

### 1. Fix GPU Compute (CRITICAL)
Investigate why `_compute->beginBatch()` fails:
- Check descriptor set bindings
- Verify buffer compatibility
- Review pipeline state
- Check Vulkan validation layers for errors

**OR implement CPU fallback:**
```cpp
void SMFTEngine::step(float dt, float K, float damping) {
    if (!_compute->beginBatch()) {
        // CPU fallback implementation
        stepCPU(dt, K, damping);
        return;
    }
    // GPU path ...
}
```

### 2. Add Diagnostics
```cpp
// In SMFTTestRunner.cpp after stepWithDirac():
std::vector<float> R_field = _engine->getSyncField();
float R_current = std::accumulate(R_field.begin(), R_field.end(), 0.0f) / R_field.size();

if (step == 0 && R_current == 0.0f) {
    std::cerr << "WARNING: R_avg = 0 at t=0 - Kuramoto likely not initialized" << std::endl;
}
```

### 3. Sanity Check in Test Runner
```cpp
// After evolution loop, BEFORE validation:
if (R_avg_final == 0.0f && config.kuramoto_initial.phase_distribution != "uniform") {
    throw std::runtime_error("Test FAILED: R=0 with non-uniform initial conditions - Kuramoto did not run!");
}
```

---

## Next Steps

1. **DO NOT regenerate plots** - current results are bogus (R=0)
2. **DO NOT mark validation complete** - GPU is broken
3. **Investigate GPU failure** - why does beginBatch() fail?
4. **Consider CPU-only mode** - implement CPU Kuramoto as fallback
5. **Rerun ALL tests** - after GPU is fixed

---

## Status Timeline

| Date | Issue | Status |
|------|-------|--------|
| Before fix | N parameter not passed to stepWithDirac() | ✓ FIXED |
| After fix | GPU compute fails every call | 🚨 BLOCKED |
| Current | No Kuramoto evolution → R=0 → free particle | 🚨 INVALID |

**VALIDATION CANNOT PROCEED UNTIL GPU COMPUTE IS FIXED.**
