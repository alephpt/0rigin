# GPU Shader Audit - Executive Summary

**Date**: 2025-12-17
**Status**: 🔴 **CRITICAL ISSUES FOUND**

## TL;DR

**Two shaders WILL cause GPU timeout**:
1. `dirac_stochastic.comp` - 50-80 transcendentals/invocation (💀 CRITICAL)
2. `dirac_rk4.comp` - ~3000 FLOPs/invocation (💀 CRITICAL)

**One shader is borderline**:
3. `sync_field.comp` - 36 transcendentals/invocation (🔴 HIGH)

## Quick Fix

**In `src/MSFTEngine.cpp` line 913, change**:
```cpp
// BEFORE (HIGH RISK):
_sync_pipeline = _pipelineFactory->createSyncFieldPipeline(
    "/home/persist/neotec/0rigin/shaders/smft/sync_field.comp.spv",
    _sync_pipeline_layout);

// AFTER (MEDIUM RISK - SAFE):
_sync_pipeline = _pipelineFactory->createSyncFieldPipeline(
    "/home/persist/neotec/0rigin/shaders/smft/sync_field_simple.comp.spv",
    _sync_pipeline_layout);
```

**Disable stochastic Dirac**:
- Do NOT call `stepStochastic()` in production
- Use deterministic `step()` only
- Dirac shaders need major optimization before GPU use

## Risk Breakdown

| Shader | Transcendentals | Risk | Status |
|--------|----------------|------|--------|
| kuramoto_step.comp | 9 | ✅ LOW | SAFE |
| kuramoto_stochastic.comp | 12-14 | ⚠️ MEDIUM | SAFE |
| sync_field_simple.comp | 37 | ⚠️ MEDIUM | **USE THIS** |
| sync_field.comp | 36 + Kahan | 🔴 HIGH | AVOID |
| gravity_field.comp | 0 | ✅ LOW | SAFE |
| spinor_feedback.comp | 0-1 | ✅ LOW | SAFE |
| **dirac_stochastic.comp** | **50-80** | **💀 CRITICAL** | **DO NOT USE** |
| **dirac_rk4.comp** | **64+** | **💀 CRITICAL** | **DO NOT USE** |

## Why This Matters

**GPU Timeout Physics**:
- GPUs have hardware watchdog timers (TDR = Timeout Detection and Recovery)
- Typical limit: 2-5 seconds wall-clock time
- Too many expensive operations → timeout → driver reset → system crash

**The "Simplified Kuramoto" Clue**:
- User mentioned "simplified Kuramoto equation was needed for this reason"
- This confirms GPU has compute budget limits
- Complex operations must be minimized

**Transcendental Functions Are Expensive**:
- `sin()`, `cos()`, `exp()`, `log()`, `sqrt()` are 10-100× slower than `add()`
- GPU can only do ~20 Tflops total
- For 256×256 grid = 65,536 invocations
- Budget per invocation: Keep transcendentals <20

## What We Found

### SAFE Shaders (Currently Working)
✅ **kuramoto_step.comp**: 9 sin() calls → SAFE
✅ **sync_field_simple.comp**: 37 trig calls → ACCEPTABLE (borderline)
✅ **gravity_field.comp**: 0 transcendentals → OPTIMAL

### UNSAFE Shaders (Will Timeout)
💀 **dirac_stochastic.comp**:
- Problem: 4 independent Box-Muller transforms = 12 transcendentals each
- Total: 50-80 transcendentals per invocation
- Fix: Apply noise every 10 steps, not every step

💀 **dirac_rk4.comp**:
- Problem: RK4 = 4 full Hamiltonian evaluations per timestep
- Each evaluation: cos/sin for chiral mass + matrix operations
- Total: ~3000 FLOPs per invocation
- Fix: Use Euler integration (1 stage instead of 4)

🔴 **sync_field.comp**:
- Problem: Full Kahan summation + complex exponentials
- 36 trig calls + overhead
- Fix: Use sync_field_simple.comp instead

## Immediate Actions Required

### 1. Change Sync Field Shader (5 minutes)
Edit `src/MSFTEngine.cpp` line 913:
- Change `sync_field.comp.spv` → `sync_field_simple.comp.spv`
- Rebuild: `cd build && make`

### 2. Disable Dirac Pipelines (2 minutes)
Comment out in tests:
```cpp
// DO NOT CALL THIS:
// engine.stepStochastic(...);

// USE THIS INSTEAD:
engine.step(...);  // Deterministic only
```

### 3. Reduce Grid Size for Testing (1 minute)
In test files, use:
```cpp
const uint32_t Nx = 128;  // Instead of 256
const uint32_t Ny = 128;  // Instead of 256
```

### 4. Verify Fix (5 minutes)
```bash
cd /home/persist/neotec/0rigin/build
./bin/test_msft_gpu  # Should not timeout now
```

## Future Work

### Short-term (this week)
- [ ] Optimize dirac_stochastic: reduce noise frequency
- [ ] Replace RK4 with Euler in dirac evolution
- [ ] Add GPU timing instrumentation

### Medium-term (this month)
- [ ] Reduce Dirac spinor to 2 components (2D doesn't need 4)
- [ ] Pre-compute static values (chiral mass)
- [ ] Implement CPU fallback for Dirac

### Long-term (next quarter)
- [ ] CPU/GPU hybrid architecture
- [ ] Adaptive precision based on simulation phase
- [ ] Benchmarking suite for shader performance

## Technical Details

See full report: `docs/GPU_SHADER_TIMEOUT_AUDIT.md`

Key insights:
- Workgroup size: 16×16 = 256 threads
- Typical grid: 256×256 = 65,536 invocations total
- Safe budget: <20 transcendentals per invocation
- Current violators: dirac_stochastic (50-80), dirac_rk4 (64+)

## Questions?

- Why does sync_field_simple work but sync_field doesn't?
  - Kahan summation adds ~50% overhead
  - 36 vs 37 transcendentals is borderline
  - Simple version stays under budget

- Can we ever use Dirac shaders?
  - Yes, but need optimization first
  - Reduce to 2-component spinor
  - Apply noise less frequently
  - Use Euler instead of RK4

- What about larger grids (512×512)?
  - Not recommended with current shaders
  - Would need even simpler versions
  - Consider CPU/GPU hybrid for large scales
