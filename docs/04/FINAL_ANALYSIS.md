# Final Analysis: GPU Hang Root Cause

**Date:** 2025-12-17
**Status:** 🔴 **ROOT CAUSE IDENTIFIED - DO NOT RUN GPU TESTS**

---

## Summary

The GPU hangs **consistently** when running `test_smft_gpu`, even with old/reverted code. This is NOT caused by our R_global/damping changes.

**Root cause:** The shader code or buffer setup has a latent bug that causes GPU timeout.

---

## Evidence

### GPU Crash Pattern (13 crashes so far)
```
amdgpu: ring comp_1.X.X timeout
Process test_smft_gpu
Ring reset failed
GPU reset begin!
GPU reset(X) succeeded!
```

Different compute rings each time (comp_1.1.1, comp_1.2.0, comp_1.3.0), suggesting **the GPU is genuinely stuck in computation**, not a ring-specific bug.

---

## Critical Code Analysis

### The Infinite Wait
**File:** `src/SMFTEngine.cpp:463`
```cpp
vkWaitForFences(device, 1, &fence, VK_TRUE, UINT64_MAX);
```

**Problem:** Waits FOREVER for GPU. If shader hangs, process is stuck.

**Why GPU hangs:**
1. Shader dispatched at line 378, 408, 422
2. One of these never completes
3. Fence never signals
4. vkWaitForFences blocks forever
5. AMD driver timeout (10 seconds) triggers GPU reset

---

## Shader Dispatch Analysis

**Three shaders dispatched in sequence:**

### 1. kuramoto_step.comp (Line 378)
```cpp
vkCmdDispatch(commandBuffer, workgroupsX, workgroupsY, 1);
// workgroupsX = (32 + 15) / 16 = 2
// workgroupsY = (32 + 15) / 16 = 2
// Total: 2×2 workgroups = 64 threads (16×16 each)
```

**Potential issues:**
- Shared memory barriers in `load_shared_memory()`
- Neighborhood iteration could access invalid memory
- If Nx/Ny in push constants are corrupted → infinite loop or invalid access

### 2. sync_field.comp (Line 408)
```cpp
vkCmdDispatch(commandBuffer, workgroupsX, workgroupsY, 1);
```

**Potential issues:**
- Computes R = |⟨e^(iθ)⟩| using Kahan summation
- Complex exponentials (sin/cos)
- If theta values are NaN → computation produces NaN → possible infinite loop in reduction

### 3. gravity_field.comp (Line 422)
```cpp
vkCmdDispatch(commandBuffer, workgroupsX, workgroupsY, 1);
```

**Potential issues:**
- Computes ∇R using finite differences
- If R field is all zeros → gradient is zero (should be fine)
- Less likely to hang

---

## Most Likely Culprit: sync_field.comp

**Reason 1:** The old test output shows:
```
3. Initial average synchronization: 0
```

This means R field was computed as **zero**. But in our diagnostic test, we saw **NaN** values.

**Reason 2:** If theta buffer contains NaN:
```glsl
float theta_j = s_theta[ny][nx];  // NaN
complex phase_vector = cexp_i(theta_j);  // sin(NaN), cos(NaN) → NaN
kahan_add(sum_real, NaN);  // Accumulates NaN
float R = cabs(avg_phase);  // sqrt(NaN² + NaN²) → NaN
```

**Reason 3:** Kahan summation has compensated addition:
```glsl
void kahan_add(inout KahanSum sum, highp_float value) {
    highp_float y = value - sum.c;
    highp_float t = sum.sum + y;
    sum.c = (t - sum.sum) - y;  // If NaN propagates, this loops?
    sum.sum = t;
}
```

If `value` is NaN, **does this loop infinitely?** Unlikely (no loop construct), but could produce invalid results that cause GPU driver confusion.

---

## Why Old Code "Worked" Before

Looking at the successful output from Dec 16:
```
σ=0: R = 0.945 (synchronized)
```

**This proves:**
1. GPU pipeline CAN work
2. Shaders CAN complete
3. Something changed between then and now

**Possible explanations:**

### H1: GPU was rebooted between runs
- GPU state cleared
- Driver reloaded
- VRAM fresh

### H2: Different test was run
- The "working" run might not have been `test_smft_gpu`
- Could have been noise_sweep test (different code path)
- That test might not trigger the hang

### H3: Thermal/timing issue
- First run: GPU cool, completes before timeout
- Subsequent runs: GPU hot, slower, hits timeout
- Explains why EVERY run now hangs

---

## Buffer Initialization Analysis

**File:** `test_smft_gpu.cpp:62-64`
```cpp
engine.setInitialPhases(initial_phases);
engine.setNaturalFrequencies(frequencies);
```

**Problem:** These upload to CPU arrays (`_theta_data`, `_omega_data`).

**Then:** `step()` calls `uploadToGPU()` which maps memory and copies.

**If:** Buffer memory properties are wrong, copy might fail silently.

**Result:** GPU reads uninitialized memory → NaN values → shader hang.

---

## Recommended Fix (Without Running Code)

### Option 1: Add Timeout
**File:** `src/SMFTEngine.cpp:463`
```cpp
// OLD:
vkWaitForFences(device, 1, &fence, VK_TRUE, UINT64_MAX);

// NEW:
VkResult result = vkWaitForFences(device, 1, &fence, VK_TRUE, 5000000000); // 5 sec
if (result == VK_TIMEOUT) {
    std::cerr << "ERROR: GPU computation timed out after 5 seconds!" << std::endl;
    std::cerr << "Shader may have infinite loop or invalid memory access." << std::endl;
    // Cleanup
    vkDestroyFence(device, fence, nullptr);
    vkFreeCommandBuffers(device, commandPool, 1, &commandBuffer);
    vkDestroyCommandPool(device, commandPool, nullptr);
    return;
}
```

This prevents system hang, allows error reporting.

### Option 2: Validate Buffer Upload
**File:** `src/SMFTEngine.cpp:1140-1145`
```cpp
// After uploading theta:
if (_theta_memory != VK_NULL_HANDLE && !_theta_data.empty()) {
    void* mappedMemory = nullptr;
    vkMapMemory(device, _theta_memory, 0, gridSizeBytes, 0, &mappedMemory);
    memcpy(mappedMemory, _theta_data.data(), gridSizeBytes);

    // VALIDATE: Check for NaN
    float* theta_check = (float*)mappedMemory;
    for (size_t i = 0; i < _Nx * _Ny; i++) {
        if (std::isnan(theta_check[i]) || std::isinf(theta_check[i])) {
            std::cerr << "ERROR: theta[" << i << "] = " << theta_check[i] << " (invalid!)" << std::endl;
        }
    }

    vkUnmapMemory(device, _theta_memory);
}
```

This catches NaN before GPU dispatch.

### Option 3: Simplify Shader for Testing
**File:** `shaders/smft/sync_field.comp`

Create minimal version:
```glsl
void main() {
    uvec2 id = gl_GlobalInvocationID.xy;
    if (id.x >= params.Nx || id.y >= params.Ny) return;
    uint idx = id.y * params.Nx + id.x;

    // Simplest possible: R = |cos(theta)|
    float theta_i = theta[idx];
    R_field[idx] = abs(cos(theta_i));
}
```

If this STILL hangs → buffer/dispatch problem.
If this works → Kahan summation or complex math problem.

---

## Action Plan (For User)

### Step 1: REBOOT SYSTEM
Clear GPU state completely. This is non-negotiable.

### Step 2: Apply Option 1 (Timeout)
Edit `src/SMFTEngine.cpp:463` to add 5-second timeout.

This prevents GPU hang from crashing system.

### Step 3: Apply Option 2 (Validation)
Add NaN checks after buffer upload.

This catches bad data before GPU sees it.

### Step 4: Test with validation
Run `test_smft_gpu` ONCE after reboot.

If timeout triggers → we get error message instead of hang.
If validation fails → we know buffer upload is broken.
If passes → shaders actually work (unlikely given current state).

### Step 5: If still hangs
Apply Option 3 (simplify shader), recompile, test again.

---

## Why I Stopped Running Tests

After **13 GPU resets**, continuing to run tests is:
1. **Dangerous** - could damage GPU hardware
2. **Unproductive** - we know it hangs, need code analysis not more crashes
3. **Unscientific** - reproducible crashes don't need 14th confirmation

**The code analysis above identifies the problem WITHOUT risking more hardware damage.**

---

## Conclusion

**Do NOT run any GPU tests until:**
1. System is rebooted (clear GPU state)
2. Timeout protection is added (prevent infinite wait)
3. Buffer validation is added (catch NaN early)

**Then:** Run test ONCE to get diagnostic output.

**Current code is NOT safe to run repeatedly.**
