# Test Framework Consolidation - Critical Architecture Fix

## Problem

**Current state**: 43 separate test executables with duplicated code
- test_smft_gpu.cpp
- test_smft_validation.cpp
- test_descriptor_bindings.cpp
- test_stochastic_cpu.cpp
- ... (39 more)

**Issues**:
1. Code duplication (Kuramoto evolution, R computation, output helpers)
2. No shared testing infrastructure
3. Inconsistent parameter sets
4. Difficult to maintain
5. Build time bloat (43 separate link steps)

---

## Solution: Single Test Framework with Modular Design

### Architecture

```
bin/smft_test [scenario] [options]

Scenarios:
  phase1/split_operator       - Phase 1 validation (Dirac split-operator)
  phase1/ehrenfest            - Force direction validation
  phase1/norm_conservation    - Long-time stability

  phase2/defect_localization  - Static defect attraction
  phase2/traveling_wave       - Wave surfing
  phase2/defect_collision     - Merger dynamics

  validation/gpu_compute      - GPU pipeline validation
  validation/descriptors      - Descriptor binding tests

  stochastic/baseline         - MSR formalism baseline
  stochastic/critical_noise   - Noise sweep
```

### Core Structure

```cpp
// src/testing/TestFramework.h
class SMFTTestFramework {
public:
    // Shared infrastructure
    void initializeKuramoto(params);
    void initializeDirac(params);
    void evolveStep(dt);

    // Diagnostics
    void trackTrajectory();
    void computeForceAlignment();
    void measureCoreD ensity();
    void computeEnergy();

    // Output
    void saveResults(scenario_name);
    void generateVisualization();
};

// test/scenarios/Phase1SplitOperator.cpp
class Phase1SplitOperator : public Scenario {
    void run() override {
        framework.initializeDirac(...);
        framework.evolveStep(dt);
        framework.validateNormConservation();
    }
};

// test/main.cpp
int main(int argc, char** argv) {
    string scenario = argv[1];

    auto test = ScenarioFactory::create(scenario);
    test->run();
    test->report();
}
```

### Migration Plan

**Phase 1: Create Core Framework**
1. `src/testing/TestFramework.{h,cpp}` - Shared infrastructure
2. `src/testing/KuramotoSim.{h,cpp}` - CPU Kuramoto evolution
3. `src/testing/Diagnostics.{h,cpp}` - Force, energy, correlation analysis
4. `src/testing/OutputManager.{h,cpp}` - Standardized data output

**Phase 2: Extract Scenarios**
1. `test/scenarios/Phase1*.cpp` - Phase 1 validations
2. `test/scenarios/Phase2*.cpp` - Phase 2 SMFT coupling
3. `test/scenarios/Validation*.cpp` - GPU/infrastructure tests
4. `test/scenarios/Stochastic*.cpp` - MSR formalism tests

**Phase 3: Consolidate**
1. `test/main.cpp` - Single entry point
2. Update CMakeLists.txt - Single executable
3. Delete 43 old test files
4. Update documentation

---

## Implementation Priority

### CRITICAL (Do Now - Before Scenario 2)
- Extract common Kuramoto evolution → `KuramotoSim`
- Extract diagnostics → `Diagnostics` (force alignment, core density, energy)
- Create single entry point `smft_test`
- Migrate Phase 2 Scenario 1 to new framework

### HIGH (Before Phase 3)
- Migrate all Phase 1 tests
- Migrate all Phase 2 tests
- Standardize output format

### MEDIUM (Cleanup)
- Delete deprecated test files
- Update CI/CD scripts
- Consolidate documentation

---

## Benefits

1. **Code Reuse**: Kuramoto/Dirac evolution written once
2. **Consistency**: Same parameter sets across all tests
3. **Maintainability**: Single codebase to update
4. **Build Speed**: 1 executable vs 43
5. **Discoverability**: `smft_test --list` shows all scenarios
6. **Reproducibility**: Standardized output format

---

## Example Usage

```bash
# Run Phase 2 Scenario 1 with custom parameters
$ ./bin/smft_test phase2/defect_localization \
    --delta 0.5 \
    --defect_radius 15 \
    --steps 50000 \
    --output output/10/

# List all available scenarios
$ ./bin/smft_test --list

# Run validation suite
$ ./bin/smft_test phase1/all
$ ./bin/smft_test phase2/all
```

---

## Next Steps

1. Implement `src/testing/TestFramework` (2-3 hours)
2. Migrate Scenario 1 to use it (1 hour)
3. Delete old `test_scenario1_defect.cpp`
4. Verify results match
5. Proceed with Scenario 2 using new framework

---

**Status**: Design complete, awaiting implementation

**Blocking**: Should complete before Scenario 2/3 to avoid further proliferation
