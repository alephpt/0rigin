# GPU Hang Diagnosis

**Date:** 2025-12-17
**Status:** 🔴 **CRITICAL - GPU COMPUTE TIMEOUT (HARDWARE HANG)**

---

## Symptoms

Running ANY test (old or new code) causes:
```
amdgpu: ring comp_1.3.0 timeout
Process test_smft_gpu
Ring comp_1.3.0 reset failed
GPU reset begin!
VRAM is lost due to GPU reset!
```

**This is a GPU hardware hang requiring driver reset.**

---

## Timeline

1. **Old noise sweep** (Dec 16 23:47): Worked, achieved R=0.945 for σ=0
2. **After our changes** (Dec 17 00:20): All NaN values
3. **After reverting** (Dec 17 06:28): GPU timeout/hang

**Critical observation:** Old code ALSO hangs the GPU now (after revert).

**This means the hang is NOT caused by our changes.**

---

## Possible Causes

### H1: Shader was always broken, previous run got lucky
- GPU timeout is NON-DETERMINISTIC
- Previous run happened to complete before timeout
- New runs hit timeout window
- **Test:** Run same binary multiple times

### H2: System state changed (driver/hardware)
- Kernel module updated
- GPU overheating
- VRAM corruption from previous crashes
- **Test:** Reboot system

### H3: Shader has infinite loop with certain inputs
- When theta[] is all zeros: division by zero → NaN propagation
- When theta[] has NaN: computation never converges
- Shader loops forever trying to compute valid result
- **Test:** Check shader for unbounded loops

### H4: Dispatch size mismatch
- Workgroups calculated wrong
- Shader dispatches too many threads
- GPU overwhelmed
- **Test:** Print workgroup sizes before dispatch

---

## Evidence from Logs

### Old run (worked):
```
3. Initial average synchronization: 0
4. Running GPU compute step...
5. GPU Computation Results:
   Synchronization field R(x):
     - Average: 0 (was 0)
   ✅ SUCCESS: GPU pipeline working correctly!
```

**Note:** R=0 suggests phases were all zero OR completely random.

### New run (hangs):
```
3. Initial average synchronization: 0
4. Running GPU compute step (dt=0.01, K=2)...
[GPU TIMEOUT - no output after this]
```

**Identical initial state, different outcome.**

---

## Shader Analysis

The `kuramoto_step.comp` shader has:
- Shared memory loading with barriers
- Neighborhood iteration (3×3 Moore)
- Trigonometric functions (sin, cos)
- No obvious infinite loops

**Potential issue:** If `params.Nx` or `params.Ny` are corrupted, the boundary conditions in `load_shared_memory()` could access invalid memory, causing GPU page fault → hang.

---

## Recommended Actions

### Action 1: Reboot system
Clear GPU state, reset VRAM, reload kernel module.

### Action 2: Check for shader compile errors
```bash
glslc --target-env=vulkan1.2 shaders/smft/kuramoto_step.comp -o /tmp/test.spv
spirv-val /tmp/test.spv
```

### Action 3: Add GPU timeout protection
In test code:
```cpp
// Set timeout for Vulkan fence wait
VkResult result = vkWaitForFences(device, 1, &fence, VK_TRUE, 5000000000); // 5 sec
if (result == VK_TIMEOUT) {
    std::cerr << "GPU computation timed out!" << std::endl;
    exit(1);
}
```

### Action 4: Simplify shader for testing
Create minimal shader:
```glsl
void main() {
    uvec2 id = gl_GlobalInvocationID.xy;
    if (id.x >= params.Nx || id.y >= params.Ny) return;
    uint idx = id.y * params.Nx + id.x;
    theta_out[idx] = theta[idx] + 0.01;  // Just add constant
}
```

Test if THIS hangs. If yes → dispatch problem. If no → shader logic problem.

---

## Critical Decision Point

**We cannot proceed with noise sweep until GPU hang is resolved.**

Two paths:

**Path A: System-level debugging**
1. Reboot
2. Check GPU temperature (`sensors`)
3. Update/downgrade AMD GPU driver
4. Test with different kernel

**Path B: Application-level debugging**
1. Simplify shader to absolute minimum
2. Add extensive logging to CPU side
3. Check push constants values before dispatch
4. Verify buffer sizes match dispatch

**Recommendation: Path A first (reboot), then Path B if problem persists.**

---

## Why This Wasn't Caught Earlier

The previous successful run (R=0.945) was **deterministic warmup only**.

Our stochastic runs may have triggered a latent bug by:
- Using different push constants structure
- Dispatching stochastic shader (new code path)
- Longer run time allowing thermal issues to manifest

But after reverting, even OLD code hangs → suggests **system state corruption**, not code bug.

---

## Next Steps

1. **User: Reboot system** (clear GPU state)
2. After reboot, test old code again
3. If still hangs → hardware/driver issue
4. If works → our changes revealed race condition or thermal problem

**Cannot continue until GPU hang resolved.**
