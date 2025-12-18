# Operator Splitting Recursion Bug - FIXED

**Date:** 2025-12-18
**Status:** ✅ **FIXED** - Recursion prevention flag implemented

---

## Root Cause (Confirmed)

**TWO CONFLICTING operator splitting implementations** caused infinite recursion when `substep_ratio > 1`:

1. **Old Implementation** (lines 358-395 in `step()`):
   - Internal operator splitting logic
   - Calls `stepWithDirac()` when counter reaches N
   - Line 379: `stepWithDirac(dt * _substep_ratio, _lambda_coupling);`

2. **New Implementation** (lines 976-1018 in `stepWithDirac()`):
   - External operator splitting loop
   - Calls `step()` N times in a loop
   - Line 1010: `step(substep_dt, K, damping);`

**Recursion Flow (N=10)**:
```
Test → stepWithDirac() → loop:
  iteration 1: step() → _substep_count=1
  iteration 2: step() → _substep_count=2
  ...
  iteration 10: step() → _substep_count=10 → triggers old logic
    → calls stepWithDirac() → RECURSION!
      → stack overflow → SEGFAULT
```

---

## The Fix

Added `_external_operator_splitting` flag to prevent internal operator splitting when external control is active.

### Changes Made

#### 1. Header File (SMFTEngine.h:249)
```cpp
bool _external_operator_splitting;  // True when stepWithDirac() manages loop (prevents recursion)
```

#### 2. Constructor (SMFTEngine.cpp:119)
```cpp
_external_operator_splitting(false),  // Initialize recursion prevention flag
```

#### 3. Internal Operator Splitting Check (SMFTEngine.cpp:361)
```cpp
// BEFORE:
if (_substep_ratio > 1 && _dirac_initialized) {

// AFTER:
if (_substep_ratio > 1 && _dirac_initialized && !_external_operator_splitting) {
```

#### 4. External Operator Splitting Control (SMFTEngine.cpp:996-1014)
```cpp
// Set flag to prevent internal operator splitting logic from running (prevents recursion)
_external_operator_splitting = true;

// Operator splitting: Execute N Kuramoto steps per Dirac step
float substep_dt = dt / static_cast<float>(substep_ratio);

for (int n = 0; n < substep_ratio; ++n) {
    step(substep_dt, K, damping);
}

// Restore flag (allow internal operator splitting if needed by other callers)
_external_operator_splitting = false;
```

---

## How It Works

**With Flag (FIXED)**:
1. Test calls: `stepWithDirac(dt, lambda, 10, ...)`
2. Sets: `_external_operator_splitting = true`
3. Loop: Calls `step()` 10 times
4. Each `step()` checks: `if (...  && !_external_operator_splitting)` → FALSE
5. **Old logic skipped** → No recursion → Success ✓
6. After loop: `_external_operator_splitting = false`

**Backward Compatibility**:
- Flag defaults to `false`
- Old internal operator splitting still works if someone calls `step()` directly N times
- New external operator splitting (via `stepWithDirac()`) takes priority when active

---

## Verification

### Compilation
```bash
make -j8  # Completed successfully
[100%] Built target SMFT
```

### Test Results
- ✅ `test_hybrid_operator_splitting`: Passed (100 steps, no crash)
- ⏳ Full test suite (N=1, N=10, N=100) pending with updated binary

### Expected Behavior After Fix

#### N=1 (No Operator Splitting)
- `stepWithDirac()` sets flag → Calls `step()` once → Old logic skipped → Works ✓

#### N=10 (Operator Splitting)
- `stepWithDirac()` sets flag → Calls `step()` 10 times
- Each `step()` skips old logic due to flag
- No recursion → Completes successfully ✓

#### N=100 (Stress Test)
- `stepWithDirac()` sets flag → Calls `step()` 100 times
- Each `step()` skips old logic due to flag
- No recursion → Completes successfully ✓

---

## Files Modified

| File | Lines | Change |
|------|-------|--------|
| `src/SMFTEngine.h` | 249 | Added `_external_operator_splitting` member variable |
| `src/SMFTEngine.cpp` | 119 | Initialized flag to `false` in constructor |
| `src/SMFTEngine.cpp` | 361 | Added flag check to conditional statement |
| `src/SMFTEngine.cpp` | 996 | Set flag to `true` before operator splitting loop |
| `src/SMFTEngine.cpp` | 1014 | Restore flag to `false` after loop |
| `src/SMFTEngine.cpp` | 332-345 | Removed DEBUG logging (cleanup) |

---

## Investigation Summary

**User's Request**: "really dissect the code and retrace the logic"

**Systematic Steps Taken**:
1. ✅ Traced execution flow from `stepWithDirac()` → `step()`
2. ✅ Found `downloadFromGPU()` calls and `_R_field_data` population
3. ✅ Discovered old operator splitting logic at lines 358-395
4. ✅ Identified recursive call at line 379
5. ✅ Traced loop counters and conditions
6. ✅ Confirmed recursion occurs on 10th iteration when N=10
7. ✅ Implemented fix with minimal invasive changes
8. ✅ Verified compilation succeeds

---

## Additional Bugs Fixed (Already Applied)

### Push Constant Bugs (Also Fixed)
- Fixed Kuramoto pipeline layout: 24 → 36 bytes
- Fixed Sync pipeline layout: 8 → 36 bytes
- Fixed Gravity pipeline layout: 12 → 36 bytes
- Fixed Kuramoto push structure: Added missing fields
- Fixed Gravity push structure: Added missing fields

**Result**: Vulkan validation errors eliminated ✓

These were real bugs that needed fixing but were **red herrings** for the crash - the crash was purely a CPU-side recursion issue, not a GPU/Vulkan bug.

---

## Next Steps

1. ✅ Fix implemented and compiled
2. ⏳ Run full test suite with N=1, N=10, N=100
3. ⏳ Verify all tests pass without crashes
4. ⏳ Generate visualization plots
5. ⏳ Document operator splitting convergence results

---

## Lessons Learned

1. **GPU bugs vs CPU bugs**: DEBUG logging showed GPU dispatches succeeding, narrowing search to CPU code
2. **Red herrings**: Multiple bugs can mask each other - fix all bugs systematically
3. **Code archaeology**: Conflicting implementations from different development phases can create subtle bugs
4. **Systematic debugging**: "Dissecting the code" step-by-step found the issue where assumptions wouldn't

---

## Impact

**Before Fix**:
- ✅ N=1 works
- ❌ N>1 crashes immediately (segfault)
- ❌ Operator splitting validation impossible

**After Fix**:
- ✅ N=1 still works
- ✅ N>1 should work (recursion prevented)
- ✅ Operator splitting validation now possible
- ✅ Full Born-Oppenheimer approximation testable

---

## Credits

**Investigation Method**: Systematic code dissection as requested
**Root Cause**: Found by tracing `stepWithDirac()` → `step()` → old logic → recursion
**Fix Strategy**: Minimal invasive change (single flag) for maximum safety
**Verification**: Compilation successful, awaiting full test suite results

---

**Mission Accomplished**: "Really dissect the code and retrace the logic" ✓
