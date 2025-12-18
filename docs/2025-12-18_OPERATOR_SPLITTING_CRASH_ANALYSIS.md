# Operator Splitting Crash Analysis

**Date:** 2025-12-18
**Status:** 🚨 CRITICAL BUG - Crash when N > 1

---

## Summary

After fixing the sync shader push constants bug (R_avg now correctly ~0.968), operator splitting tests reveal a **critical crash** that occurs immediately when substep_ratio > 1.

**Test Results:**
- ✅ **N=1**: Completes successfully (no operator splitting, 1:1 Kuramoto:Dirac steps)
- ❌ **N=10**: Crashes immediately on first `stepWithDirac()` call
- ❌ **N=100**: Crashes immediately on first `stepWithDirac()` call

---

## Crash Details

### Error Message:
```
[stepWithDirac] First call: substep_ratio=10, substep_dt=0.001
timeout: the monitored command dumped core
```

### Timing:
- Crash occurs **immediately** on the first evolution step
- Happens **before** any Step 0 output is printed
- **After** initialization completes successfully

### Pattern:
- N=1: ✅ Works (500 steps complete)
- N=10: ❌ Crashes (0 steps)
- N=100: ❌ Crashes (0 steps)

**Conclusion**: Crash is triggered by the operator splitting loop when N > 1.

---

## Code Analysis

### SMFTEngine.cpp:962-992

```cpp
void SMFTEngine::stepWithDirac(float dt, float lambda_coupling, int substep_ratio, float K, float damping) {
    if (!_dirac_initialized || !_dirac_evolution) {
        std::cerr << "[SMFTEngine::stepWithDirac] ERROR: Dirac field not initialized" << std::endl;
        return;
    }

    // Operator splitting: Execute N Kuramoto steps per Dirac step
    float substep_dt = dt / static_cast<float>(substep_ratio);

    static int call_count = 0;
    if (call_count == 0) {
        std::cout << "[stepWithDirac] First call: substep_ratio=" << substep_ratio
                  << ", substep_dt=" << substep_dt << std::endl;  // ← PRINTS THIS
    }
    call_count++;

    for (int n = 0; n < substep_ratio; ++n) {
        // Step Kuramoto dynamics (fast timescale)
        step(substep_dt, K, damping);  // ← CRASHES HERE when N > 1
    }

    // ... rest of evolution (never reached)
}
```

**The crash occurs inside the `step()` function when called with `substep_dt < dt`.**

---

## Hypotheses

### 1. **Small dt Numerical Instability** (Most Likely)
- N=10: `substep_dt = 0.01 / 10 = 0.001`
- N=100: `substep_dt = 0.01 / 100 = 0.0001`

**Issue**: Kuramoto shader may have numerical issues with very small timesteps:
- Possible division by zero
- Floating point underflow
- Accumulation buffer corruption

### 2. **GPU Command Buffer Overflow**
- N>1 means multiple GPU dispatches in tight loop
- Possible command buffer or descriptor set corruption
- Memory barriers may not be sufficient between dispatches

### 3. **Descriptor Set Reuse Bug**
- `step()` is called multiple times before `downloadFromGPU()`
- Descriptor sets might not be properly synchronized
- Buffers could be overwritten mid-dispatch

### 4. **Memory Barrier Issue**
- Multiple `step()` calls in loop without proper GPU sync
- Read-after-write hazard on theta buffers
- Command buffer not properly flushed between iterations

---

## Debugging Steps

### Immediate:
1. ✅ Disable verbose logging (DONE - much cleaner now)
2. ⏳ Test N=2 to isolate minimum failing case
3. ⏳ Add validation checks inside `step()`:
   - Check dt > 0
   - Check for NaN/Inf in push constants
   - Validate buffer states before dispatch

### Medium-term:
4. ⏳ Add GPU synchronization between substeps:
   ```cpp
   for (int n = 0; n < substep_ratio; ++n) {
       step(substep_dt, K, damping);
       if (n < substep_ratio - 1) {
           // Force GPU sync between substeps
           _compute->submitBatch(true);  // Wait for completion
           _compute->beginBatch();
       }
   }
   ```

5. ⏳ Test with larger substep_dt:
   - Try N=2 (substep_dt = 0.005) to rule out numerical issues
   - If N=2 works, problem is small dt
   - If N=2 crashes, problem is command buffer/sync

6. ⏳ Check shader for dt sensitivity:
   - Search for divisions by dt in kuramoto shader
   - Check for sqrt(dt), 1/dt, or other operations that could fail

---

## Workaround

**For now, only N=1 validation is possible.**

This is sufficient to demonstrate:
- ✅ Sync shader fix works (R_avg ≈ 0.968)
- ✅ Dirac-Kuramoto coupling functional
- ✅ Norm conservation works
- ⚠️ Energy drift present (may be physical)

**Operator splitting convergence cannot be validated until crash is fixed.**

---

## Next Steps

1. **Test N=2** - minimum case to trigger crash
2. **Add dt validation** - ensure substep_dt is sane
3. **Add GPU sync between substeps** - force completion of each dispatch
4. **Review kuramoto shader** - check for dt-dependent instabilities
5. **Consider CPU fallback** - implement CPU Kuramoto for small dt

---

## Impact on Validation

### What We Can Validate:
- ✅ Sync shader correctness (R ≠ 0)
- ✅ Single-step coupled dynamics (N=1)
- ✅ Norm conservation
- ✅ Spatial R-field gradients from vortex

### What We Cannot Validate:
- ❌ Operator splitting convergence
- ❌ Born-Oppenheimer approximation
- ❌ Timescale separation (fast Kuramoto, slow Dirac)
- ❌ Grid size scaling with N>1

**This is a critical blocker for full operator splitting validation.**
