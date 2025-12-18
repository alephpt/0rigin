# GPU Dirac Timesync Test Report

## Executive Summary

**FINDING: Timesync N is already implemented in the codebase!**

The operator splitting with timesync ratio is fully implemented in `MSFTEngine`:
- Default N = 10 (10 Kuramoto steps per 1 Dirac step)
- Can be set via `setSubstepRatio(N)` method
- Implementation at lines 328-360 of `MSFTEngine.cpp`
- Born-Oppenheimer approximation as described in `Feature-Not-Bug.md`

## 1. Timesync Implementation Found

### Code Evidence

**In MSFTEngine.h (lines 150-156)**:
```cpp
/**
 * Set the substep ratio N for operator splitting adiabatic approximation
 * N = ratio of fast (Kuramoto) to slow (Dirac) timescales
 * Typical values: 10 (testing), 100 (production)
 *
 * @param N Number of Kuramoto substeps per Dirac step
 */
void setSubstepRatio(int N);
```

**In MSFTEngine.cpp**:
- Line 120: `_substep_ratio(10),  // Default: 10 Kuramoto steps per Dirac step`
- Line 304: Check if operator splitting enabled with `_substep_ratio > 1`
- Line 328-360: Full operator splitting logic:
  - Accumulate theta and R fields over N substeps
  - Every N steps: compute time averages
  - Update Dirac with averaged fields at dt_slow = N * dt_fast
- Lines 1044-1055: `setSubstepRatio(int N)` implementation

### How It Works

1. **Fast subsystem (Kuramoto)**: Runs every timestep dt_fast = 0.01
2. **Accumulation**: Sum theta and R fields over N steps
3. **Time averaging**: After N steps, divide sums by N
4. **Slow subsystem (Dirac)**: Update with dt_slow = N * dt_fast using averaged fields
5. **Reset**: Clear accumulators and repeat

This matches the Born-Oppenheimer approximation exactly as described in Feature-Not-Bug.md.

## 2. GPU Shader Status

### Shaders Found
- ✅ `dirac_rk4.comp` - RK4 integration for Dirac
- ✅ `dirac_stochastic.comp` - Stochastic Dirac with noise
- ✅ `kuramoto_stochastic.comp` - Stochastic Kuramoto
- ✅ `accumulate.comp` - Accumulation for operator splitting
- ✅ `sync_field.comp` - Compute R field

### GPU Safety Status (from MSFTPipelineFactory.cpp)
- Kuramoto: ✅ SAFE (10 transcendentals)
- Sync field: ✅ SAFE
- Dirac RK4: ❌ DANGEROUS (~3000 FLOPs, 10× over budget)
- Dirac stochastic: ❌ DANGEROUS (50-80 transcendentals, 4× over budget)

### Current Implementation
- Kuramoto runs on GPU (safe)
- Dirac runs on CPU (via `stepWithDirac()` and `DiracEvolution` class)
- This is exactly the hybrid GPU-CPU approach needed!

## 3. Memory Mapping Issue

### Problem Identified
```
Failed to map memory at vkMapMemory
```

**Root Cause**: The accumulator buffers and some other buffers are allocated without proper memory flags for host access.

**Location**: `src/MSFTBufferManager.cpp` line 122

**Issue**: When trying to download accumulated sums for averaging, the memory cannot be mapped because it wasn't allocated with `VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT`.

## 4. Testing Results

### Test 1: test_dirac_gpu_timesync
- **Status**: ❌ Failed due to memory corruption (double-free in destructor)
- **Issue**: Resources destroyed twice - once in `destroyResources()`, again in component destructors

### Test 2: test_timesync_simple
- **Status**: ❌ Failed due to vkMapMemory error
- **Issue**: Cannot map GPU memory to download accumulated fields

### Valgrind Analysis
Found double-free issues:
1. Descriptor set layouts freed in pool destruction, then freed again in DescriptorManager
2. Buffers freed in destroyResources(), then freed again in BufferManager destructor

## 5. Recommendations

### Immediate Fixes Needed

1. **Fix Memory Allocation** (HIGH PRIORITY)
   - Ensure all buffers that need CPU access have `VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT`
   - Specifically: accumulator buffers need host visibility for downloading

2. **Fix Double-Free Issue** (HIGH PRIORITY)
   - Remove duplicate cleanup in either destroyResources() or component destructors
   - Use smart pointers or clear ownership model

3. **Complete Hybrid Implementation** (MEDIUM)
   - Verify initializeHybrid() properly sets up GPU-CPU split
   - Test with different N values (10, 100, 1000)

### Testing Strategy

Once fixes are applied:
1. Start with N=10 (minimal separation)
2. Run 100 timesteps, monitor for NaN/Inf
3. Increase to N=100 (production separation)
4. Verify norm conservation and energy stability
5. Document performance gains vs full GPU

## 6. Conclusion

**The timesync implementation already exists and follows the Born-Oppenheimer approximation correctly!**

Key findings:
- ✅ Operator splitting implemented with configurable N
- ✅ Default N=10, can be set to 100 or 1000
- ✅ GPU Kuramoto + CPU Dirac hybrid approach in place
- ❌ Memory allocation bugs prevent execution
- ❌ Double-free issues cause crashes

**Next Steps**:
1. Fix memory allocation flags for host-visible buffers
2. Fix double-free destructor issue
3. Test with N=100 to validate stability
4. Document performance improvement

The infrastructure is there - it just needs these two bugs fixed to work properly.