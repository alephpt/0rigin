# SMFT Experimental Roadmap

**Status**: Implementation Phase - Ready for Verification
**Date**: 2025-12-16

---

## Overview

This document tracks the two major experimental campaigns for validating Synchronization Mass Field Theory (SMFT):

1. **Noise Sweep Experiment** - Path A vs Path B decision
2. **Dirac Coupling Experiment** - Particle generation validation

---

## Experiment 1: Noise Sweep (Path B Validation)

### Objective
Determine if SMFT can tolerate thermal noise at the Planck scale, validating the stochastic formalism (Path B) vs deterministic formalism (Path A).

### Critical Question
**What is σ_c, the critical noise amplitude where synchronization breaks down?**

### Decision Criterion
```
σ_c > 10^-5  →  Path B (stochastic) validated ✓
σ_c < 10^-5  →  Path A (deterministic) required ✗
```

### Implementation Status

#### ✓ Completed
- [x] PRNG implementation (PCG algorithm) in `kuramoto_stochastic.comp`
- [x] Box-Muller Gaussian transform
- [x] Euler-Maruyama stochastic integrator
- [x] Test driver: `test_noise_sweep.cpp`
- [x] Documentation: `docs/noise_sweep_experiment.md`
- [x] Added `stepStochastic()` method signature to SMFTEngine.h

#### ⏳ Pending
- [ ] Implement `stepStochastic()` in SMFTEngine.cpp
- [ ] Compile stochastic shader (kuramoto_stochastic.comp.spv)
- [ ] Run Verification Test 1: PRNG quality
- [ ] Run Verification Test 2: Noise scaling
- [ ] Run Verification Test 3: Deterministic limit
- [ ] Execute noise sweep (8 σ values)
- [ ] Fit L(σ) curve to extract σ_c
- [ ] Make Path A vs Path B decision

### Timeline
- **Week 1** (Dec 16-22): Implementation + Verification
- **Week 2** (Dec 23-29): Noise sweep execution
- **Week 3** (Dec 30 - Jan 5): Analysis + Decision

### Key Files
```
shaders/smft/kuramoto_stochastic.comp   - Stochastic integrator shader
test/test_noise_sweep.cpp               - Experiment driver
docs/noise_sweep_experiment.md          - Full protocol
build/output/noise_sweep/               - Results (to be generated)
```

---

## Experiment 2: Dirac Coupling (Particle Generation)

### Objective
Test if Dirac spinors localize in vacuum defects and produce discrete energy spectrum, validating the core SMFT mechanism for particle mass.

### Five Success Criteria
1. **Localization**: O > 0.7 (Dirac field concentrates at defects)
2. **Stabilization**: ΔS > 10% (Defect survival rate increases)
3. **Discrete Energies**: 2-5 peaks in energy histogram
4. **Particle Count**: N ∈ [10, 200] (selective binding)
5. **Stability**: τ > 500 steps (bound states don't decay)

### Implementation Status

#### ✓ Completed
- [x] Dirac RK4 shader exists: `dirac_rk4.comp`
- [x] Spinor feedback shader exists: `spinor_feedback.comp`
- [x] Test driver: `test_dirac_coupling.cpp`
- [x] Documentation: `docs/dirac_coupling_experiment.md`
- [x] Success criteria analysis functions
- [x] Energy spectrum computation
- [x] Defect counting algorithms

#### ⏳ Pending
- [ ] Implement `initializeDiracField()` in SMFTEngine
- [ ] Implement `stepWithDirac()` in SMFTEngine
- [ ] Implement `getDiracDensity()` in SMFTEngine
- [ ] Connect Dirac pipeline to test driver
- [ ] Run Phase 1: Vacuum equilibration
- [ ] Run Phase 2: Dirac initialization
- [ ] Run Phase 3: Coupled evolution (3 λ values)
- [ ] Analyze energy histogram for discrete peaks
- [ ] Validate all 5 criteria

### Timeline
- **Week 1** (Dec 16-22): Implementation
- **Week 2** (Dec 23-29): Initial tests (λ = 0.1, 1.0, 10.0)
- **Week 3** (Dec 30 - Jan 5): Validation + Final report

### Key Files
```
shaders/smft/dirac_rk4.comp             - Dirac evolution (exists)
shaders/smft/spinor_feedback.comp       - Feedback mechanism (exists)
test/test_dirac_coupling.cpp            - Experiment driver
docs/dirac_coupling_experiment.md       - Full protocol
build/output/dirac_coupling/            - Results (to be generated)
```

---

## Dependencies & Prerequisites

### For Both Experiments

**Hardware**:
- Vulkan 1.2+ capable GPU
- Compute shader support
- 8GB+ VRAM recommended (256² grid = 256KB per field)

**Software**:
- Nova engine with compute pipeline
- SMFTEngine with GPU buffer infrastructure
- glslc shader compiler
- CMake build system

### Experiment-Specific

**Noise Sweep**:
- Requires: `SMFTEngine::stepStochastic()` implementation
- Shader: kuramoto_stochastic.comp compiled
- ~8 hours compute time (8 σ values × 6000 steps × 256² grid)

**Dirac Coupling**:
- Requires: Dirac field management in SMFTEngine
- Methods: `initializeDiracField()`, `stepWithDirac()`, `getDiracDensity()`
- ~12 hours compute time (3 λ values × 2000 steps × 256² grid × RK4)

---

## Build Instructions

### Add to CMakeLists.txt

```cmake
# Noise sweep test
add_executable(test_noise_sweep test/test_noise_sweep.cpp)
target_link_libraries(test_noise_sweep Nova ${VULKAN_LIBRARIES} SDL2)

# Dirac coupling test
add_executable(test_dirac_coupling test/test_dirac_coupling.cpp)
target_link_libraries(test_dirac_coupling Nova ${VULKAN_LIBRARIES} SDL2)

# Compile stochastic shader
add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/shaders/smft/kuramoto_stochastic.comp.spv
    COMMAND glslc
        ${CMAKE_SOURCE_DIR}/shaders/smft/kuramoto_stochastic.comp
        -o ${CMAKE_BINARY_DIR}/shaders/smft/kuramoto_stochastic.comp.spv
    DEPENDS ${CMAKE_SOURCE_DIR}/shaders/smft/kuramoto_stochastic.comp
)
```

### Compile & Run

```bash
cd build
cmake ..
make test_noise_sweep
make test_dirac_coupling

# Run noise sweep
./bin/test_noise_sweep

# Run Dirac coupling
./bin/test_dirac_coupling
```

---

## Expected Outputs

### Noise Sweep Results

```
build/output/noise_sweep/
├── results.csv                         # σ, L_mean, R_mean, etc.
├── sigma_0.00e+00/                     # Deterministic baseline
│   ├── L_timeseries.dat
│   ├── R_field_step_0000.dat
│   └── ...
├── sigma_1.00e-05/                     # Critical threshold
└── sigma_1.00e+00/                     # High noise
```

**Key Result**: `σ_c = ??? ± ???` extracted from L(σ) fit

### Dirac Coupling Results

```
build/output/dirac_coupling/
├── results_summary.txt                 # All 5 criteria: PASS/FAIL
├── lambda_0.1/
│   ├── energy_spectrum.dat            # KEY: Discrete peaks?
│   ├── overlap_timeseries.dat
│   └── step_0000/
│       ├── R_field.dat
│       └── psi_density.dat
├── lambda_1.0/
└── lambda_10.0/
```

**Key Result**: Energy histogram showing 2-5 discrete peaks

---

## Success Metrics

### Scientific Success (Both Experiments)
- ✓ Hypothesis tested rigorously
- ✓ Results interpreted honestly
- ✓ Theory matches experiment
- ✓ Falsifiability demonstrated

**Either outcome (pass/fail) is publishable.**

### Technical Success
- ✓ All verification tests pass
- ✓ Clear measurements with <20% uncertainty
- ✓ Reproducible results
- ✓ No implementation bugs

### Failure Modes to Avoid
- ✗ Verification tests fail → broken implementation
- ✗ Ambiguous results → need better resolution
- ✗ Data contradicts theory but theory kept anyway → dishonest
- ✗ Claim validation without evidence → unscientific

---

## Integration with Existing Tests

### Current Test Suite
```
test/test_smft_headless.cpp            - Basic headless test
test/test_smft_compute_only.cpp        - Compute pipeline validation
test/test_smft_phase0.cpp              - Deterministic characterization
test/test_defect_detection.cpp         - Vortex detection
test/test_defect_evolution.cpp         - Defect annihilation
test/test_mass_quantization.cpp        - Mass distribution analysis
```

### New Tests (This Work)
```
test/test_noise_sweep.cpp              - Path B validation ← NEW
test/test_dirac_coupling.cpp           - Particle generation ← NEW
```

---

## Next Steps (Priority Order)

1. **Implement SMFTEngine::stepStochastic()** (1-2 hours)
   - Load kuramoto_stochastic pipeline
   - Pass sigma parameter via push constants
   - Dispatch compute shader

2. **Compile stochastic shader** (10 minutes)
   - `glslc kuramoto_stochastic.comp -o kuramoto_stochastic.comp.spv`
   - Add to CMake

3. **Run Verification Tests** (2 hours)
   - Test 1: PRNG quality
   - Test 2: Noise scaling
   - Test 3: Deterministic limit
   - **GATE**: Must pass before noise sweep

4. **Execute Noise Sweep** (8 hours compute)
   - 8 σ values
   - 1000 measurement steps each
   - Save all outputs

5. **Analyze σ_c** (1 day)
   - Fit L(σ) curve
   - Extract critical threshold
   - **DECISION**: Path A or Path B

6. **Implement Dirac coupling** (2-3 days)
   - initializeDiracField()
   - stepWithDirac()
   - getDiracDensity()

7. **Run Dirac experiment** (12 hours compute + 2 days analysis)
   - 3 λ values
   - Validate 5 criteria
   - Analyze energy spectrum

---

## Risk Assessment

### High Risk
- ⚠️ **Noise sweep shows σ_c ≈ 10^-5** (ambiguous region)
  - Mitigation: Finer σ sampling around threshold
  - Timeline impact: +1 week

- ⚠️ **Dirac field doesn't localize** (Criterion 1 fails)
  - Mitigation: Increase Δ, check potential depth
  - May indicate fundamental theory issue

### Medium Risk
- ⚠️ **Energy spectrum continuous** (Criterion 3 fails)
  - Mitigation: Higher grid resolution, longer runs
  - Alternative: Quantization may require 3D

- ⚠️ **Verification tests fail** (implementation bug)
  - Mitigation: Debug PRNG, check numerical stability
  - Timeline impact: +3 days

### Low Risk
- ⚠️ **Compute time exceeds estimate**
  - Mitigation: Reduce grid size (256² → 128²)
  - Or run overnight batches

---

## References

1. **Directive.md** - Experimental protocol authorization
2. **Determinism.md** - Effective stochasticity hypothesis
3. **Dirac-Anomaly.md** - Defect evolution + Dirac predictions
4. **noise_sweep_experiment.md** - Detailed protocol for Experiment 1
5. **dirac_coupling_experiment.md** - Detailed protocol for Experiment 2

---

## Contact

**Lead**: [Your Name]
**Platform**: Nova Vulkan Compute Engine
**GPU**: [Model]
**Start Date**: 2025-12-16
**Expected Completion**: 2026-01-05

**Status Dashboard**: See `docs/experiment_roadmap.md` (this file)

---

**Last Updated**: 2025-12-16 21:45 UTC
**Next Milestone**: Implement stepStochastic() - Due Dec 17
