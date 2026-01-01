# R_avg Corruption Bug - Root Cause and Fix

## The Bug

**Symptom**: `R_avg` observable computed correctly (0.9998) inside `ObservableComputer::compute()` but appeared as 0 in caller/output.

**All Failed Attempts**:
1. Return by value
2. Output parameter by reference
3. Output parameter by pointer
4. Global pointer hack (`g_result_hack`)
5. Memcpy
6. Inline computation bypass

## Root Cause

**The actual problem was NOT in ObservableComputer or struct passing.**

The root cause was in `SMFTEngine::stepWithDirac()`:

```cpp
// BEFORE (BROKEN):
for (int n = 0; n < substep_ratio; ++n) {
    step(substep_dt, K, damping);  // Updates _R_field_data on CPU
}
std::vector<float> R_field = getSyncField();  // Returns _R_field_data
// ... later ...
downloadFromGPU();  // OVERWRITES _R_field_data with zeros from GPU!
```

### The Issue

1. `step()` uses **CPU fallback** (no GPU shaders loaded)
2. CPU fallback correctly updates `_R_field_data` with values ~0.9998
3. `getSyncField()` reads the good CPU data
4. **THEN** `downloadFromGPU()` is called, which **overwrites** the good CPU data with zeros from uninitialized GPU buffers
5. Subsequent calls to `getSyncField()` in test runner return all zeros

## The Fix

**File**: `/home/persist/neotec/0rigin/src/SMFTEngine.cpp`

```cpp
// AFTER (FIXED):
for (int n = 0; n < substep_ratio; ++n) {
    step(substep_dt, K, damping);
}

// === CRITICAL FIX: Only download if using GPU ===
// CPU fallback updates _R_field_data directly
// GPU path needs download, but we're in CPU fallback mode
if (_kuramoto_pipeline && _sync_pipeline && _gravity_pipeline) {
    downloadFromGPU();
}

std::vector<float> R_field = getSyncField();  // Now returns correct data
```

### Why This Works

- **CPU fallback mode**: Pipelines are NULL, so we skip `downloadFromGPU()`, preserving CPU-computed data
- **GPU mode**: Pipelines exist, so we download GPU data before reading
- **Correct sequencing**: Download (if needed) happens BEFORE getSyncField() is called

## Verification

```bash
cd build && make -j4
./bin/smft --test ../config/stuckelberg_integration_test.yaml 2>&1 | grep "R_avg"
```

**Output**:
```
Step 0/20000 | R_avg = 0.999802 | norm = 1
Step 1000/20000 | R_avg = 0.999803 | norm = 0.999674
```

## Lessons Learned

1. **The bug was NOT where it appeared to be** - symptoms suggested struct corruption, actual cause was data synchronization
2. **Debug the data flow, not just the function** - the value was correct inside the function but became corrupted later
3. **CPU/GPU synchronization is subtle** - blindly calling downloadFromGPU() can overwrite good CPU data
4. **Check your assumptions** - the "nuclear option" (inline computation) revealed the field data itself was corrupted, not the struct passing

## Files Modified

1. `/home/persist/neotec/0rigin/src/SMFTEngine.cpp` - Added conditional GPU download
2. `/home/persist/neotec/0rigin/src/simulations/SMFTTestRunner.cpp` - Cleaned up debug code

## Status

✅ **RESOLVED** - R_avg now displays correct non-zero values throughout simulation.
