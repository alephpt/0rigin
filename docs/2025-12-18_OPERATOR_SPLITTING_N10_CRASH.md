# Operator Splitting N>1 Crash - FINAL ROOT CAUSE

**Date:** 2025-12-18
**Status:** 🔴 SEGFAULT - Null pointer write in CPU code

---

## Summary

After fixing multiple Vulkan push constant issues, N=10 tests still crash with **CPU segfault**, not GPU error.

---

## Crash Details

### dmesg Output:
```
smft[623832]: segfault at 0 ip 0000560739ca86c8 sp 00007fffa7a7d380 error 6 in smft[1b6c8,560739c98000+14f000]
Code: 41 0f 10 04 04 0f 5e c1 <0f> 11 04 02 0f 10 44 05 00 0f 5e c1 0f 11 04 01 48 83 c0 10 48 39
```

### Analysis:
- **Error 6**: `SEGV_ACCERR` - write to read-only/unmapped memory
- **Address 0**: Null pointer dereference
- **Assembly**: `movups [rdx+rax],xmm0` - SSE store instruction
- **Location**: CPU code (not GPU), likely in DiracEvolution or buffer download

---

## Fixes Applied (Push Constants)

### 1. Pipeline Layout Sizes (SMFTEngine.cpp lines 631, 678, 722)
```cpp
// BEFORE (BROKEN):
pushRange.size = sizeof(float) * 6;            // Kuramoto: 24 bytes
pushRange.size = sizeof(uint32_t) * 2;         // Sync: 8 bytes
pushRange.size = sizeof(float) + sizeof(uint32_t) * 2; // Gravity: 12 bytes

// AFTER (FIXED):
pushRange.size = sizeof(float) * 5 + sizeof(uint32_t) * 4;  // All: 36 bytes
```

### 2. Kuramoto Push Constants (SMFTEngine.cpp line 253)
```cpp
// BEFORE (BROKEN - 24 bytes):
struct KuramotoPush {
    float dt, K, damping, Delta;
    uint32_t Nx, Ny;
};

// AFTER (FIXED - 36 bytes):
struct KuramotoPush {
    float dt, K, damping, Delta, chiral_angle;
    uint32_t Nx, Ny, N_total, enable_feedback;
};
```

**Result**: Vulkan validation errors GONE ✓, but crash persists.

---

## Current Hypothesis

The crash occurs **after** GPU execution, likely during:

1. **Buffer download** (`downloadFromGPU()`)
   - Possible null pointer in theta/R_field buffers
   - Buffer mapping failure

2. **DiracEvolution CPU update**
   - Null spinor buffer pointer
   - Uninitialized memory access

3. **Operator splitting loop state**
   - Buffer corruption between substeps
   - Descriptor set state invalidation

---

## Key Observation

**N=1 works perfectly** (R_avg ≈ 0.968, correct physics)
**N>1 crashes immediately** on first operator splitting loop iteration

This suggests the issue is **NOT**:
- Small dt instability (N=1 uses same dt)
- Shader bugs (N=1 uses same shaders)
- Initial GPU setup (initialization succeeds)

The issue **IS**:
- Something that happens during the N-fold loop
- Triggered by multiple `step()` calls before `downloadFromGPU()`
- CPU-side buffer management or pointer invalidation

---

## Debugging Next Steps

1. **Add null pointer checks** before all buffer writes
2. **Validate buffer pointers** after each GPU dispatch
3. **Check DiracEvolution state** between substeps
4. **Test with gdb** to get exact crash location
5. **Add logging** in downloadFromGPU() and DiracEvolution

---

## Test Configuration

**Working:** N=1, 64×64, 500 steps ✓
**Crashing:** N=10, 64×64, crashes at step 0
**Crashing:** N=100, 64×64, crashes at step 0

---

## Files Modified

- `src/SMFTEngine.cpp:631` - Kuramoto pipeline layout size
- `src/SMFTEngine.cpp:678` - Sync pipeline layout size
- `src/SMFTEngine.cpp:722` - Gravity pipeline layout size
- `src/SMFTEngine.cpp:253` - Kuramoto push constants structure
- `src/SMFTCompute.cpp:62` - Command buffer reset (didn't fix crash)

---

## Conclusion

The push constant bugs were **red herrings**. The real bug is a **CPU-side null pointer** that only manifests when operator splitting calls `step()` multiple times in a loop.

**Action Required**: Debug with gdb or add extensive logging to find the null pointer source.
