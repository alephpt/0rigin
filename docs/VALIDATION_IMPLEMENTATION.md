# SMFT Validation Framework - Implementation Summary

## Overview

This document describes the implemented C++ validation framework for SMFT simulations, based on the requirements in `VALIDATION_FRAMEWORK.md`.

**Status**: ✅ Framework implemented and compiled successfully
**Date**: 2025-12-20
**Next Steps**: Integrate with SMFTTestRunner and add YAML configuration support

---

## Architecture

### File Structure

```
src/validation/
├── ValidationCriteria.h         # Data structures and enums
├── ValidationCriteria.cpp       # ValidationReport serialization
├── GlobalValidator.h            # 6 global requirements
├── GlobalValidator.cpp          # Global validation implementation
├── ScenarioValidator.h          # Test-specific requirements
└── ScenarioValidator.cpp        # Scenario validation implementation
```

### Build Integration

Added to `CMakeLists.txt` as part of the main SMFT executable:
```cmake
set(SMFT_SOURCES
    ...
    src/validation/ValidationCriteria.cpp
    src/validation/GlobalValidator.cpp
    src/validation/ScenarioValidator.cpp
    ...
)
```

**Status**: ✅ Compiles successfully with no errors

---

## I. Global Validation (6 Requirements)

### Implementation: `GlobalValidator` class

#### 1. Probability Conservation
```cpp
CriterionResult checkProbabilityConservation(
    double norm,
    double norm_initial,
    double tolerance);
```
- **Check**: `|norm - norm_initial| < tolerance`
- **Default tolerance**: 0.5% (0.005)
- **Critical**: YES - simulation invalid if failed

#### 2. Energy Conservation
```cpp
CriterionResult checkEnergyConservation(
    double E_current,
    double E_initial,
    double tolerance);
```
- **Check**: `|ΔE/E₀| < tolerance`
- **Default tolerance**: 1% (0.01)
- **Critical**: YES - simulation invalid if failed

#### 3. Order Parameter Bounds
```cpp
CriterionResult checkOrderParameterBounds(
    const std::vector<double>& R_field);
```
- **Check**: `0 ≤ R(x,y,t) ≤ 1` for all points
- **Critical**: YES - unphysical if violated

#### 4. Gauge Invariance
```cpp
CriterionResult checkGaugeInvariance(
    const ObservableComputer::Observables& obs);
```
- **Status**: Placeholder (advanced feature)
- **Critical**: NO - optional check

#### 5. Causality
```cpp
CriterionResult checkCausality(
    const ObservableComputer::Observables& obs,
    double c_light);
```
- **Check**: `v = |p|/E ≤ c`
- **Default**: c = 1.0 (Planck units)
- **Critical**: YES - relativistic framework broken if violated

#### 6. Numerical Stability
```cpp
CriterionResult checkNumericalStability(
    const ObservableComputer::Observables& obs,
    const std::vector<double>& R_field);
```
- **Check**: All fields finite (no NaN, no Inf)
- **Critical**: YES - numerical breakdown if failed

### Public Interface

```cpp
// Pre-flight check (before evolution)
ValidationReport GlobalValidator::validateInitialState(
    const DiracEvolution& dirac,
    const std::vector<double>& R_field,
    const std::vector<double>& theta_field,
    const GlobalCriteria& criteria);

// Runtime monitoring (during evolution)
ValidationReport GlobalValidator::validateRuntimeState(
    const ObservableComputer::Observables& current_obs,
    const ObservableComputer::Observables& initial_obs,
    const std::vector<double>& R_field,
    const GlobalCriteria& criteria,
    double time);

// Final validation (after evolution)
ValidationReport GlobalValidator::validateFinalState(
    const ObservableComputer::Observables& final_obs,
    const ObservableComputer::Observables& initial_obs,
    const std::vector<double>& R_field,
    const GlobalCriteria& criteria);
```

---

## II. Scenario-Specific Validation

### Implementation: `ScenarioValidator` class

#### Scenario Types

```cpp
enum class ScenarioType {
    NONE,
    DEFECT_LOCALIZATION,       // Scenario 2.1
    TRAVELING_WAVE,            // Scenario 2.2
    RELATIVISTIC_MASS,         // Scenario 2.3
    BREAKDOWN_INVESTIGATION    // Scenario 2.4
};
```

### Scenario 2.1: Defect Localization

```cpp
ValidationReport ScenarioValidator::validateDefectLocalization(
    const DiracEvolution& dirac,
    const std::vector<double>& R_field,
    const std::vector<double>& theta_field,
    const ObservableComputer::Observables& obs,
    const DefectLocalizationCriteria& criteria);
```

**Required checks**:
- ✓ Vortex present: W = ±1
- ✓ R-field core: R_min < 0.5

**NOT required**:
- ✗ Particle present (vacuum only)
- ✗ Boosted conditions

### Scenario 2.2: Traveling Wave Surfing

```cpp
ValidationReport ScenarioValidator::validateTravelingWave(
    const DiracEvolution& dirac,
    const std::vector<double>& R_field,
    const std::vector<double>& theta_field,
    const ObservableComputer::Observables& initial_obs,
    const ObservableComputer::Observables& final_obs,
    double gamma_theory,
    const TravelingWaveCriteria& criteria);
```

**6 Critical Criteria** (as requested by user):

1. **θ(x,y,t=0) shows vortex structure**
   ```cpp
   CriterionResult checkVortexStructure(
       const std::vector<double>& theta_field,
       int Nx, int Ny,
       double tolerance);
   ```
   - Computes winding number: W = (1/2π) ∮ ∇θ · dl
   - Passes if |W - 1| < 0.2

2. **W = ±1 computed from phase winding**
   - Included in criterion 1

3. **R(x,y,t=0) shows core (R_min < 0.5)**
   ```cpp
   CriterionResult checkRFieldCore(
       const std::vector<double>& R_field,
       double threshold);
   ```
   - Passes if R_min < 0.5

4. **ψ(x,y,t=0) is Gaussian at offset position**
   ```cpp
   CriterionResult checkGaussianWavepacket(
       const DiracEvolution& dirac,
       double tolerance);
   ```
   - Checks spatial localization via variance

5. **⟨p⟩(t=0) = γmv (within 5%)**
   ```cpp
   CriterionResult checkInitialMomentum(
       const ObservableComputer::Observables& obs,
       double p_expected,
       double tolerance);
   ```
   - Computes p_mag = √(p_x² + p_y²)
   - Passes if |p_mag - p_expected| / p_expected < 0.05

6. **γ_measured(t=final) within 5% of theory**
   ```cpp
   CriterionResult checkGammaFactor(
       const ObservableComputer::Observables& obs,
       double gamma_theory,
       double delta,
       double tolerance);
   ```
   - Computes m_eff = √(E² - p²)
   - γ_measured = m_eff / (Δ · R_avg)
   - Passes if |γ_measured - γ_theory| / γ_theory < 0.05

### Scenario 2.3: Relativistic Mass Validation

```cpp
ValidationReport ScenarioValidator::validateRelativisticMass(
    const DiracEvolution& dirac,
    const std::vector<double>& R_field,
    const std::vector<double>& theta_field,
    const ObservableComputer::Observables& initial_obs,
    const ObservableComputer::Observables& final_obs,
    double gamma_theory,
    const RelativisticMassCriteria& criteria);
```

**Includes all Scenario 2.2 requirements, plus**:
- Grid convergence (checked externally by comparing multiple runs)
- N-convergence (checked externally)

---

## III. Validation Configuration

### Data Structures

```cpp
struct GlobalCriteria {
    double norm_tolerance = 0.005;         // 0.5%
    double energy_tolerance = 0.01;        // 1%
    bool enforce_R_bounds = true;
    bool check_gauge_invariance = false;
    bool enforce_causality = true;
    double c_light = 1.0;
    bool check_numerical_stability = true;
};

struct TravelingWaveCriteria {
    bool require_vortex = true;
    double winding_tolerance = 0.2;
    bool require_core = true;
    double core_R_threshold = 0.5;
    bool require_gaussian = true;
    double gaussian_width_tolerance = 0.2;
    bool require_initial_boost = true;
    double initial_momentum_tolerance = 0.05;  // 5%
    bool check_particle_tracking = true;
    bool validate_gamma_factor = true;
    double gamma_tolerance = 0.05;             // 5%
};

struct RelativisticMassCriteria {
    TravelingWaveCriteria base_criteria;
    bool require_grid_convergence = true;
    double grid_convergence_tolerance = 0.02;
    bool require_N_convergence = true;
    double N_convergence_tolerance = 0.05;
};
```

### YAML Configuration (Planned)

```yaml
validation:
  # GLOBAL (always enforced)
  norm_tolerance: 0.005        # 0.5%
  energy_tolerance: 0.01       # 1%
  R_bounds_check: true
  causality_check: true

  # TEST-SPECIFIC
  scenario: "traveling_wave"   # Determines which checks to run

  # Traveling wave specific
  gamma_tolerance: 0.05        # 5% for γ factor
  require_vortex: true
  require_core: true
  require_boost: true

  # Validation timing
  validate_initial_state: true
  validate_during_evolution: true
  validation_interval: 100     # Steps between runtime checks
  validate_final_state: true
```

**Implementation**: `ValidationConfig::fromYAML()` method provided

---

## IV. Validation Reports

### Report Structure

```cpp
struct CriterionResult {
    std::string name;
    bool passed;
    double measured_value;
    double expected_value;
    double tolerance;
    std::string message;
    bool is_critical;
};

struct ValidationReport {
    std::vector<CriterionResult> global_results;
    bool global_pass;

    std::vector<CriterionResult> scenario_results;
    bool scenario_pass;

    bool overall_pass;
    std::string summary;

    // Serialization
    std::string toString() const;
    void saveToFile(const std::string& filepath) const;
};
```

### Example Output

```
================================================================================
VALIDATION REPORT
================================================================================

GLOBAL REQUIREMENTS (Apply to ALL simulations):
--------------------------------------------------------------------------------
✓ Probability Conservation
  ||ψ||² = 1.0002 (drift: 0.02% from initial 1.0000) ✓ PASS

✓ Energy Conservation
  E_total = 0.9987 (drift: 0.13% from initial 1.0000) ✓ PASS

✓ Order Parameter Bounds
  R ∈ [0.0124, 0.9998] ✓ PASS (within [0,1])

✓ Causality
  v = |p|/E = 0.314 (c = 1.0) ✓ PASS (v ≤ c)

✓ Numerical Stability
  All fields finite ✓ PASS

Global status: ✅ PASS

SCENARIO-SPECIFIC REQUIREMENTS:
--------------------------------------------------------------------------------
✓ 1. θ(x,y,t=0) shows vortex structure (W = ±1)
  Winding number W = 1.03 ✓ PASS (|W| ≈ 1)

✓ 3. R(x,y,t=0) shows core (R_min < 0.5)
  R_min = 0.0124 ✓ PASS (< 0.5)

✓ 4. ψ(x,y,t=0) is Gaussian at offset position
  Wavepacket variance = 0.0023 ✓ PASS (localized)

✓ 5. ⟨p⟩(t=0) = γmv (within 5%)
  ⟨p⟩ = 0.3145 (expect 0.3145, error = 0.01%) ✓ PASS

✓ 6. γ_measured(t=final) within 5% of theory
  γ_measured = 1.0521 (γ_theory = 1.0483, error = 0.36%) ✓ PASS

Scenario status: ✅ PASS

================================================================================
OVERALL STATUS: ✅ PASS
Traveling Wave validation: ✅ PASS (6 criteria checked)
================================================================================
```

---

## V. Integration Plan

### Step 1: Add to SMFTTestRunner (TODO)

```cpp
class SMFTTestRunner {
private:
    Validation::ValidationConfig m_validationConfig;
    Validation::ValidationReport m_initialValidation;
    Validation::ValidationReport m_finalValidation;

public:
    void loadValidationConfig(const YAML::Node& yaml);

    void validateInitialState();
    void validateRuntimeState(int step);
    void validateFinalState();

    void saveValidationReport(const std::string& filepath);
};
```

### Step 2: Add YAML Schema (TODO)

Extend `TestConfig` to parse `validation:` section from YAML.

### Step 3: Hook into Test Loop (TODO)

```cpp
// In SMFTTestRunner::runSingleConfiguration()

// Before evolution
m_initialValidation = GlobalValidator::validateInitialState(...);
if (!m_initialValidation.overall_pass && m_validationConfig.fail_on_critical) {
    std::cerr << "CRITICAL: Initial state validation failed!\n";
    std::cerr << m_initialValidation.toString() << "\n";
    return;  // Abort test
}

// During evolution (every N steps)
if (step % m_validationConfig.validation_interval == 0) {
    auto runtime_report = GlobalValidator::validateRuntimeState(...);
    if (!runtime_report.overall_pass) {
        std::cerr << "WARNING: Runtime validation failed at step " << step << "\n";
    }
}

// After evolution
m_finalValidation = ScenarioValidator::validateTravelingWave(...);
m_finalValidation.saveToFile(output_dir + "/validation_report.txt");

if (!m_finalValidation.overall_pass) {
    std::cerr << "FINAL: Validation failed - test results unreliable\n";
}
```

---

## VI. Testing the Framework

### Phase 2.3.Alpha Verification Script

The Python script `verify_phase_2.3.Alpha.py` demonstrates the expected validation checks.

Once integrated with SMFTTestRunner, the C++ framework will enforce these same checks automatically.

### Example Test with Validation

```yaml
# config/traveling_wave_validated.yaml

test:
  name: "Phase 2.2 Traveling Wave (Validated)"

validation:
  scenario: "traveling_wave"
  norm_tolerance: 0.005
  energy_tolerance: 0.01
  gamma_tolerance: 0.05
  require_vortex: true
  require_core: true
  require_boost: true
  validation_interval: 100
  verbose: true
  fail_on_critical: true

grid:
  size_x: 128
  size_y: 128

initial_conditions:
  dirac:
    type: "boosted_gaussian"
    boost_velocities: [0.3]
  kuramoto:
    type: "vortex"
    vortex_charge: 1
```

**Expected behavior**:
- ✅ Initial validation runs before any evolution
- ✅ Runtime validation runs every 100 steps
- ✅ Final validation runs after evolution completes
- ✅ Validation report saved to `output/.../validation_report.txt`
- ✅ Test aborts if critical validation fails

---

## VII. Summary

### Implemented ✅

1. **Global validation** (6 requirements) - `GlobalValidator`
2. **Scenario-specific validation** (2.1, 2.2, 2.3) - `ScenarioValidator`
3. **6 criteria checklist** for Phase 2.2/2.3 (as requested by user)
4. **Validation configuration** structures
5. **Validation report** generation and serialization
6. **YAML parsing** support (fromYAML method)
7. **CMake integration** and successful compilation

### TODO 🔨

1. **Integrate with SMFTTestRunner** - Add validation hooks to test loop
2. **Extend TestConfig** - Parse `validation:` section from YAML
3. **Save spatial snapshots** - Required for criteria 1-3 (θ, R fields)
4. **Fix boosted Gaussian initialization** - Required for criteria 4-5
5. **Test framework end-to-end** - Run Phase 2.2 with validation enabled
6. **Re-run Phase 2.2** - Verify all 6 criteria pass
7. **Re-run Phase 2.3** - Apply validation to all 50 configurations

---

**Status**: Validation framework foundation complete, ready for integration testing

**Next Action**: Integrate validation calls into SMFTTestRunner and extend YAML schema
