# GPU Buffer Overflow Fix Summary

## Issue
Critical GPU compute ring timeout caused by buffer overflow in sync_field shader.

## Root Cause
- **Location**: `SMFTEngine.cpp` line 1528 in `stepStochastic()` function
- **Problem**: `neighborhood_radius` was set to `5` but shader's shared memory `s_theta[18][18]` only supports radius=1
- **Effect**: Out-of-bounds memory access causing GPU hang and compute ring timeout

## Fix Applied
Changed line 1528 in `SMFTEngine.cpp`:
```cpp
// BEFORE (incorrect):
} sync_push = {
    dt, K, damping, _Delta, _chiral_angle, _Nx, _Ny, _Nx * _Ny, 5
};

// AFTER (fixed):
} sync_push = {
    dt, K, damping, _Delta, _chiral_angle, _Nx, _Ny, _Nx * _Ny, 1  // neighborhood_radius must be 1 to match sync_field.comp shared memory size [18][18]
};
```

## Verification
1. **Consistency Check**: Verified that regular `step()` function correctly uses `neighborhood_radius = 1.0f`
2. **No Other Issues**: Confirmed no other locations in codebase use incorrect radius values
3. **Build Success**: Project compiles successfully with the fix

## Technical Details
The sync_field compute shader allocates shared memory as:
```glsl
shared float s_theta[18][18];  // 16x16 work group + 1-cell halo on each side
```

This allocation supports:
- Work group size: 16x16
- Halo cells: 1 on each side
- Total size: 18x18
- Maximum supported radius: 1

Using radius=5 would require `(16 + 2*5) x (16 + 2*5) = 26x26` array, causing buffer overflow.

## Impact
This fix resolves:
- GPU compute ring timeout errors
- Device lost errors (VK_ERROR_DEVICE_LOST)
- System stability issues during stochastic simulations

## Prevention
Added explanatory comment to prevent future regressions. The constraint is now clearly documented in the code.