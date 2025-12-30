# SMFT Validation Framework - Next Steps

## Summary of Work Completed

### 1. Validation Framework Implementation ✅
Created C++ validation framework with full 6-criteria support:

**Files Created:**
- `src/validation/ValidationCriteria.h` - Data structures and enums
- `src/validation/ValidationCriteria.cpp` - Report serialization
- `src/validation/GlobalValidator.h/cpp` - 6 global requirements
- `src/validation/ScenarioValidator.h/cpp` - Test-specific validation (6 criteria)

**Files Modified:**
- `src/simulations/TestConfig.h` - Extended ValidationConfig with all parameters
- `src/simulations/TestConfig.cpp` - Updated parseValidation() to parse new fields
- `CMakeLists.txt` - Added validation sources to build

**Build Status:** ✅ Compiles successfully

### 2. Documentation Created ✅
- `docs/VALIDATION_FRAMEWORK.md` - Complete framework specification
- `docs/VALIDATION_IMPLEMENTATION.md` - Implementation guide with examples
- `docs/NEXT_STEPS.md` - This file

### 3. Python Verification Script ✅
- `verify_phase_2.3.Alpha.py` - Implements 6 criteria validation in Python

---

## Critical Issues Identified (Phase 2.2 Failures)

Running `verify_phase_2.3.Alpha.py` on Phase 2.2 data revealed **ALL 3 tests FAILED**:

### Issue 1: Missing Spatial Snapshots
**Problem**: No `theta_field_t0.00.csv` or `R_field_t0.00.csv` files saved
**Impact**: Cannot verify criteria 1-3 (vortex structure, winding number, R-field core)
**Status**: ❌ BLOCKING

### Issue 2: Initial Momentum ≈ 0
**Problem**: `⟨p⟩(t=0) = 0.0001` instead of expected `0.314`
**Expected**: For v=0.3c, γ=1.0483, m=1: p = γ·m·v = 0.3145
**Error**: 99.98%!
**Impact**: Boosted Gaussian initialization is broken
**Status**: ❌ BLOCKING - Criteria 5 fails

### Issue 3: γ_measured = 2.57 (145% error!)
**Problem**: Measured gamma factor way off theory
**Expected**: γ_theory = 1.0483
**Measured**: γ_measured ≈ 2.57
**Error**: 145%!
**Impact**: Entire relativistic mass validation fails
**Status**: ❌ BLOCKING - Criteria 6 fails

---

## Immediate Next Steps (Priority Order)

### Step 1: Fix Spatial Snapshot Saving ⏳
**Goal**: Save theta and R fields at t=0 for vortex verification

**Action**: Update `SMFTTestRunner` to call spatial snapshot saving
**Location**: `src/simulations/SMFTTestRunner.cpp`
**Method**: `saveSpatialFieldSnapshot()`

**Already exists in header but not called!**
```cpp
void saveSpatialFieldSnapshot(int N, int step) const;
```

**TODO**:
1. Implement `saveSpatialFieldSnapshot()` in `SMFTTestRunner.cpp`
2. Call at initialization: `saveSpatialFieldSnapshot(N, 0)`
3. Call at configured snapshots per `output.snapshot_steps`
4. Save format:
   ```
   theta_field_tXX.XX.csv  (columns: x, y, theta)
   R_field_tXX.XX.csv       (columns: x, y, R)
   ```

### Step 2: Fix Boosted Gaussian Initialization ⏳
**Goal**: Initialize spinor with correct momentum ⟨p⟩ = γmv

**Problem locations** (likely):
1. `src/DiracEvolution.cpp` - `initializeBoostedGaussian()` method
2. Initial condition setup in `SMFTTestRunner`

**Theory**: Boosted Gaussian should have phase factor:
```
ψ(x,y,t=0) = exp(i·p·r) * Gaussian(x-x0, y-y0, σ)
where p = γ·m·v
```

**TODO**:
1. Search for boosted Gaussian initialization code
2. Verify phase factor is being applied correctly
3. Test with simple case: v=0.3c, verify p=0.314

### Step 3: Debug γ_measured Calculation ⏳
**Goal**: Understand why γ factor is so far off

**Formula (from ObservableComputer.cpp)**:
```cpp
m_eff = √(E² - p²)
γ_measured = m_eff / (Δ · R_avg)
```

**Possible causes**:
- Energy calculation wrong
- Momentum calculation wrong (related to Issue 2)
- R_avg not being computed correctly
- Δ (delta) parameter mismatch

**TODO**:
1. Add debug output to print E, p, R_avg at t=0 and t=final
2. Verify energy-momentum relation manually
3. Check if Issue 2 fix resolves this

### Step 4: Integrate Validation into SMFTTestRunner ⏳
**Goal**: Hook validation framework into test execution loop

**Changes needed in `SMFTTestRunner.cpp`**:

```cpp
#include "validation/GlobalValidator.h"
#include "validation/ScenarioValidator.h"

class SMFTTestRunner {
private:
    Validation::ValidationReport m_initialValidation;
    Validation::ValidationReport m_finalValidation;

    void validateInitialState(int N);
    void validateRuntimeState(int N, int step);
    void validateFinalState(int N);
};

// In runSingleTest(int N):
if (_config.validation.validate_initial_state) {
    validateInitialState(N);
    if (!m_initialValidation.overall_pass && _config.validation.fail_on_critical) {
        std::cerr << "CRITICAL: Initial state validation failed!\n";
        std::cerr << m_initialValidation.toString() << "\n";
        return false;  // Abort test
    }
}

// In evolution loop:
if (_config.validation.validate_during_evolution &&
    step % _config.validation.validation_interval == 0) {
    validateRuntimeState(N, step);
}

// After evolution:
if (_config.validation.validate_final_state) {
    validateFinalState(N);
    m_finalValidation.saveToFile(output_dir + "/validation_report.txt");
}
```

### Step 5: Create Phase 2.3 Test Config with Validation ⏳
**Goal**: Create a properly configured YAML file for Phase 2.3 re-run

**File**: `config/phase_2.3_validated.yaml`

```yaml
test:
  name: "Phase 2.3 - Relativistic Mass Validation"
  description: "Validate m(v) = γ·Δ·R formula across velocities with 6 criteria"

# Validation configuration (NEW!)
validation:
  # Global requirements
  norm_tolerance: 0.005        # 0.5%
  energy_tolerance: 0.01       # 1%
  enforce_R_bounds: true
  enforce_causality: true
  check_numerical_stability: true

  # Scenario-specific (6 criteria for Phase 2.3 - includes all Phase 2.2 requirements)
  scenario: "relativistic_mass"
  require_vortex: true         # Criterion 1 & 2: W = ±1
  winding_tolerance: 0.2
  require_core: true           # Criterion 3: R_min < 0.5
  core_R_threshold: 0.5
  require_boost: true          # Criterion 5: ⟨p⟩ = γmv
  initial_momentum_tolerance: 0.05  # 5%
  validate_gamma_factor: true  # Criterion 6: γ_measured
  gamma_tolerance: 0.05        # 5%

  # Validation timing
  validate_initial_state: true
  validate_during_evolution: false
  validate_final_state: true
  fail_on_critical: true
  verbose: true

# Grid (same as before)
grid:
  size_x: 128
  size_y: 128
  L_domain: 100.0

# Physics (same as before)
physics:
  delta: 1.0
  coupling: 0.1
  dt: 0.01
  total_steps: 10000
  K: 1.0
  damping: 0.1

# Grid convergence (Phase 2.3 requirement)
grid:
  grid_sizes: [64, 128, 256]  # Test all three for convergence

# Initial conditions
initial_conditions:
  dirac:
    type: "boosted_gaussian"
    boost_velocities: [0.0, 0.3, 0.5, 0.7]  # Multiple velocities for Phase 2.3
    sigma_physical: 5.0        # Grid-independent
    x0_physical: 60.0          # Offset from vortex center
    y0_physical: 50.0
  kuramoto:
    type: "vortex"
    winding_number: 1
    vortex_core_radius: 3.0
    vortex_center_x: 50.0
    vortex_center_y: 50.0

# Operator splitting (Phase 2.3 convergence testing)
operator_splitting:
  enabled: true
  substep_ratios: [1, 10, 100]

# Output
output:
  directory: "output/phase_2.3_validated"
  save_every: 100
  save_spatial_snapshots: true
  snapshot_steps: [0, 100, 500, 1000, 5000, 9999]  # Include t=0!
  auto_visualize: true
```

---

## Testing Strategy

### Phase 1: Fix Critical Issues
1. Implement spatial snapshot saving
2. Fix boosted Gaussian initialization
3. Debug γ_measured calculation
4. Verify fixes with minimal test (1 configuration)

### Phase 2: Integrate Validation
1. Hook validation into SMFTTestRunner
2. Test with Phase 2.2 config
3. Verify all 6 criteria checked correctly

### Phase 3: Validate Minimal Test
**Verify fixes work on simple case first:**

```bash
cd /home/persist/neotec/0rigin
# Create minimal test: single grid (128), single N (10), single velocity (0.3c)
./build/bin/smft --test config/phase_2.2_minimal.yaml
python3 verify_phase_2.3.Alpha.py  # Should now PASS all 6 criteria
```

### Phase 4: Re-run Phase 2.3 with Full Validation
**Only after minimal test passes all 6 criteria:**

```bash
./build/bin/smft --test config/phase_2.3_validated.yaml
```

This will run:
- 3 grid sizes (64, 128, 256)
- 3 operator splitting ratios (N=1, 10, 100)
- 4 velocities (0.0, 0.3, 0.5, 0.7c)
- **Total: 36 individual configurations**

Each configuration will be validated against the 6 criteria and report individual results (not averaged!)

---

## Success Criteria

### Minimal Test (Phase 2.2 at v=0.3c) Must Pass ALL 6 Criteria:

1. ✅ θ(x,y,t=0) shows vortex structure (W ≈ ±1)
2. ✅ W = ±1 computed from phase winding (|W-1| < 0.2)
3. ✅ R(x,y,t=0) shows core (R_min < 0.5)
4. ✅ ψ(x,y,t=0) Gaussian at offset (localized wavepacket)
5. ✅ ⟨p⟩(t=0) = 0.314 ± 5% (boosted correctly)
6. ✅ γ_measured = 1.048 ± 5% (relativistic mass correct)

**If ANY criterion fails on minimal test: DEBUG until it passes. Do NOT proceed to full Phase 2.3.**

### Phase 2.3 Full Validation Success:
- Each of the 36 individual configurations must be validated separately
- Report pass/fail for each configuration
- Identify patterns:
  - Which velocities cause failures?
  - Which grid sizes are insufficient?
  - Which N values fail?
- Generate convergence analysis showing grid and N-dependence

---

## File Roadmap

### Core Framework (Completed ✅)
- `src/validation/ValidationCriteria.{h,cpp}`
- `src/validation/GlobalValidator.{h,cpp}`
- `src/validation/ScenarioValidator.{h,cpp}`
- `src/simulations/TestConfig.{h,cpp}` (extended)

### Next Files to Modify (TODO ⏳)
- `src/simulations/SMFTTestRunner.cpp` - Add validation hooks
- `src/DiracEvolution.cpp` - Fix boosted Gaussian init
- `config/phase_2.2_validated.yaml` - New test config

### Verification Scripts
- `verify_phase_2.3.Alpha.py` - Python validation (already works!)

---

## Current Status

✅ **Framework implemented and compiled**
✅ **YAML schema extended**
✅ **Documentation complete**
⏳ **Integration with test runner** - TODO
⏳ **Fix spatial snapshots** - TODO (BLOCKING)
⏳ **Fix boosted Gaussian** - TODO (BLOCKING)
⏳ **Debug γ calculation** - TODO (BLOCKING)

**Next Immediate Action**: Implement `saveSpatialFieldSnapshot()` in SMFTTestRunner.cpp

---

**Date**: 2025-12-20
**Author**: SMFT Validation Team
