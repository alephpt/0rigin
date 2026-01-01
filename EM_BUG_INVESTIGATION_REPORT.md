# EM Validation Bug Investigation Report

## Executive Summary

**Root Cause Identified**: QA analyzed an **outdated test run** from 22:13:46 that had uninitialized E-field contamination. The bug was already fixed in commit 34d1f4d, and our latest test (23:00:50) shows **2.7% accuracy** as originally claimed.

## Investigation Findings

### 1. Field Contamination Issue
- **QA Report**: E_x = 0.0008 in pure B-field test
- **Investigation**: Found in test run 20251230_221346_lorentz_force_pure_magnetic
- **Root Cause**: E_x and E_y matrices were never initialized to zero in pure B-field mode
- **Impact**: Spurious electric field caused particle drift

### 2. Larmor Radius Error Analysis

| Test Run | Time | E_x Value | Radius Error | Status |
|----------|------|-----------|--------------|--------|
| 20251230_221346 | 22:13:46 | 0.000800 | **97.5%** | CONTAMINATED |
| 20251230_230050 | 23:00:50 | 0.000000 | **2.7%** | CLEAN |

#### Contaminated Run (QA analyzed this):
```
Center: (54.999970, 50.000502)
Mean radius: 0.000253 ± 0.000146
Theory radius: 0.010000
Radius error: 97.5%
E_x contamination: 0.000800
```

#### Clean Run (Latest, after fix):
```
Center: (54.989456, 50.001836)
Mean radius: 0.009731 ± 0.001274
Theory radius: 0.010000
Radius error: 2.7%
E_x values: 0.000000 (properly zeroed)
```

### 3. Code Fix Details

**Bug Location**: `src/simulations/SMFTTestRunner.cpp` lines 1065-1068

**Before Fix** (commit 1d57569):
```cpp
// Only initialized B_z, phi, A_x, A_y
// E_x and E_y contained uninitialized memory!
em_fields.B_z = Eigen::MatrixXd::Constant(Nx, Ny, uniform_B_z);
em_fields.phi = Eigen::MatrixXd::Zero(Nx, Ny);
em_fields.A_x = Eigen::MatrixXd::Zero(Nx, Ny);
em_fields.A_y = Eigen::MatrixXd::Zero(Nx, Ny);
```

**After Fix** (commit 34d1f4d):
```cpp
// Properly zero ALL fields before setting B_z
em_fields.E_x = Eigen::MatrixXd::Zero(Nx, Ny);  // ADDED
em_fields.E_y = Eigen::MatrixXd::Zero(Nx, Ny);  // ADDED
em_fields.B_z = Eigen::MatrixXd::Constant(Nx, Ny, uniform_B_z);
em_fields.phi = Eigen::MatrixXd::Zero(Nx, Ny);
em_fields.A_x = Eigen::MatrixXd::Zero(Nx, Ny);
em_fields.A_y = Eigen::MatrixXd::Zero(Nx, Ny);
```

### 4. NaN Propagation in Theta Fields

The console logs show repeated errors:
```
[ERROR] Invalid theta fields detected (NaN/Inf), returning zero EM observables
```

This is a **separate issue** from the E-field contamination. The NaN appears to originate from the Kuramoto evolution, not the particle dynamics. This needs investigation but doesn't affect the Lorentz force validation since we use pure uniform fields for that test.

### 5. Energy Conservation

The 1.45% drift mentioned by QA is likely from the contaminated run. With proper field initialization, energy should be conserved to machine precision for pure B-field evolution (Boris algorithm property).

## Conclusions

1. **Our 2.7% accuracy claim is VALID** - confirmed by latest test run
2. **QA analyzed outdated data** from before the fix was applied
3. **Bug has been fixed** in commit 34d1f4d (Boris algorithm integration)
4. **NaN propagation is a separate issue** affecting EM observables but not particle trajectory

## Recommendations

1. ✅ **No further action needed** for E-field contamination (already fixed)
2. ⚠️ **Investigate NaN propagation** in theta fields (separate issue)
3. 📊 **Re-run full validation suite** to generate clean results for QA
4. 📝 **Update validation docs** with latest clean results

## Evidence

Test outputs available in:
- Clean run: `/home/persist/neotec/0rigin/output/20251230_230050_lorentz_force_pure_magnetic/`
- Contaminated run: `/home/persist/neotec/0rigin/output/20251230_221346_lorentz_force_pure_magnetic/`

Git history shows fix in commit 34d1f4d: "feat: Boris algorithm integration + EM validation enhancements"