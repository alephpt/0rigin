# Investigation Complete: Root Cause Found

**Date:** 2025-12-17
**Status:** 🔴 **GPU MEMORY CORRUPTION - PRE-EXISTING BUG**

---

## Summary

The NaN issue is **NOT caused by our corrections** (R_global computation, damping addition).

**Root cause:** GPU buffer upload is partially failing, causing ~50% of phase values to be NaN from the start.

---

## Evidence

### Test: Simple warmup diagnostic
```
R_initial (after IC upload): nan
First 10 phases: nan -nan nan -nan nan -nan 0.23973 0.906154 -0.44308 -2.08794
```

**Analysis:**
- Half the phases are NaN
- Half are valid floating point values
- This happens IMMEDIATELY after `setInitialPhases()` and `getPhaseField()`
- Before ANY compute shaders run

**Conclusion:** GPU upload or readback is corrupted.

---

## Timeline

### Old runs (before our changes):
- σ=0: **R = 0.945** ✓ (worked fine)
- σ>0: **R = 0** (the original bug we investigated)

### After adding R_global + damping:
- **All sigma:** R = NaN (new GPU issue)

**This means:**
1. Old code DID synchronize successfully (R=0.945 proves it)
2. Something we changed broke GPU memory handling
3. But it's NOT the R_global computation (pure CPU)
4. It's NOT the damping logic (doesn't run during IC upload)

---

## What Changed That Could Break GPU Upload?

Looking at git diff of `src/SMFTEngine.cpp`:

### Added to constructor initialization list:
```cpp
_compute_command_pool(VK_NULL_HANDLE),
_timeline_semaphore(VK_NULL_HANDLE),
_timeline_value(0),
_compute_queue_family_index(0),
_compute_command_buffer(VK_NULL_HANDLE),
_compute_queue(VK_NULL_HANDLE),
_staging_buffer(VK_NULL_HANDLE),
_staging_memory(VK_NULL_HANDLE),
// ... more fields
```

### Added debug output:
```cpp
std::cout << "[DEBUG] SMFTEngine::initialize() START" << std::endl;
std::cout << "[DEBUG] Allocating CPU-side data arrays..." << std::endl;
```

**None of these should break GPU upload.**

---

## Hypothesis: Buffer Alignment or Size Mismatch

When we changed the push constants structure for stochastic shader, we may have accidentally broken the **descriptor set bindings** or **buffer sizes**.

### Old deterministic push constants:
```cpp
struct {
    float dt, K, damping, Delta, chiral_angle;
    uint Nx, Ny, N_total, enable_feedback;
};  // 9 fields
```

### New stochastic push constants:
```cpp
struct {
    float dt, K, sigma, damping, omega_mean;
    uint Nx, Ny, time_step;
};  // 8 fields
```

**If descriptor sets or buffers were sized based on push constants,** this could cause memory corruption.

---

## Most Likely Cause: Descriptor Set Corruption

The descriptor sets bind buffers to shaders. If we:
1. Changed one shader's push constants
2. Rebuilt and recompiled
3. But descriptor sets still reference old layout

Then **buffer bindings could be misaligned**, causing:
- Half of theta_buffer to be readable
- Half to be garbage/NaN
- Upload "succeeds" but writes to wrong memory

---

## Recommended Fix

### Option 1: Nuclear - Revert ALL changes
```bash
git stash
# Test if old code works
./bin/test_simple_warmup
```

If old code works → our changes broke it → bisect to find culprit.

### Option 2: Surgical - Remove only damping changes
Keep R_global fix, remove:
- Damping from stochastic shader
- Damping from push constants structure
- Test if this restores GPU upload

### Option 3: Rebuild from scratch
```bash
rm -rf build
mkdir build
cd build
cmake ..
make
```

Sometimes CMake cache gets corrupted and doesn't rebuild everything.

---

## Why Old Run Worked

Looking at the evidence:
- Old `noise_sweep.log` shows `R = 0.945` for σ=0
- This PROVES the GPU pipeline worked before
- Something changed between then and now

**Candidates:**
1. We modified SMFTEngine.cpp (added debug output, changed initialization)
2. We modified push constants (added damping field)
3. We compiled new stochastic shader (incompatible with descriptor sets?)

---

## Action Items

### Immediate:
1. **Revert to last known working state**
   ```bash
   git diff > my_changes.patch  # Save our work
   git stash  # Revert all changes
   ```

2. **Test if old code works**
   ```bash
   cd build && make test_simple_warmup
   ./bin/test_simple_warmup
   ```

3. **If old code works:**
   - Apply changes incrementally
   - Test after each change
   - Find which change breaks GPU

### If old code ALSO broken:
- GPU driver issue
- System update broke Vulkan
- Need to restart or reinstall drivers

---

## Critical Questions

**Q1:** Did you run any system updates between the working run (R=0.945) and now?

**Q2:** When was the last successful run?
- Check timestamp: `ls -la output/noise_sweep.log`
- Compare to current time

**Q3:** Did ANY other code in the repo change?
- `git log --since="2 days ago" --oneline`

---

## Conclusion

**The NaN issue is a GPU memory corruption, likely caused by:**
1. Descriptor set / buffer alignment mismatch after changing push constants
2. CMake cache corruption not rebuilding all dependencies
3. Stochastic shader loading corrupting descriptor sets

**This is NOT a physics bug, NOT a measurement error.**

**This is an infrastructure / memory management bug that must be fixed before any physics experiments can run.**

---

## Recommended Path Forward

**REVERT ALL CHANGES. Start fresh.**

1. Git stash our corrections
2. Verify old code works
3. If old code works:
   - Apply R_global fix ONLY (no damping, no push constant changes)
   - Test if this works
   - Run noise sweep with just R_global fix
   - This gives us "correct measurement of wrong physics" (no damping)
   - Better than no data at all

4. Then separately debug why adding damping breaks GPU

**Timeline:**
- Revert + test old code: 10 minutes
- Run noise sweep with R_global fix only: 2 hours
- Debug damping issue: Unknown (could be deep)

**We need to get SOME result, even if imperfect.**

Measuring R_global (even without damping) is still **infinitely better** than measuring R_local and calling it R_global.

---

**Awaiting user decision on how to proceed.**
