- ./trd is the single executable we should have. No additional executables. All testing should be ran through the --test flag, and all core tests that are crucial to TRD should be part of the engine to ensure 100% correctness for our theory

      IMPORTANT: this context may or may not be relevant to your tasks. You should not respond to this context unless it is highly relevant to your task.

## TRD-Specific Standards

### 1. Single Unified Executable
- **ABSOLUTE RULE**: `./trd --test <config.yaml>` is the ONLY execution pattern
- NO standalone test binaries (no test_*, no separate executables per test)
- All physics tests integrated into single TRD executable via test harness
- Legacy smft binary deprecated - all work uses trd
- Verify with: `ls build/bin/` should only show `trd` executable

### 2. Core Framework Integration
- **MANDATORY**: All tests MUST use TRDCore3D/TRDEngine3D infrastructure
- NO isolated test implementations bypassing core engine
- Wave 4 architecture violations discovered: tests implementing own integrators had >1% energy drift
- Root cause: Tests bypassed proven symplectic TRDCore3D framework
- Prevention: All tests MUST inherit from unified engine classes

### 3. Symplectic Integration Standard
- **REQUIRED**: All field evolution MUST use symplectic integrators for conservative physics
- Approved methods:
  - RK2 Midpoint Method (TRDCore3D default implementation)
  - Velocity Verlet (particles, Sine-Gordon wave equations)
  - Half-Strang split-stepping (θ-R phase-magnitude decoupling)
- **FORBIDDEN**: Forward Euler (dissipative, becomes heat equation)
- **FORBIDDEN**: RK4 (rejected - 0.0002% drift exceeds our 0.01% standard)
- Energy conservation requirement: ΔE/E < 0.01% (GO/NO-GO criterion)
- Time reversibility validation required for all integrators

### 4. Conservative Field Dynamics
- Field evolution MUST preserve energy functionals
- Sine-Gordon for topological solitons: ∂²θ/∂t² = ∇²θ - sin(θ)
- Klein-Gordon for scalar fields (with documented limitations)
- NO diffusion equations (∂θ/∂t = K·∇²θ) for conservative physics
- Velocity Verlet kick-drift-kick pattern for wave equations
- Reference: D4 scattering fix - migrated from dissipative (85% loss) to Sine-Gordon (0.127% drift)

### 5. YAML-Based Configuration
- **ABSOLUTE RULE**: All tests configured via .yaml files in config/
- No hardcoded physics parameters in test implementations
- Config must specify:
  - test_file: implementation file to execute
  - golden_key: 246 GeV standard reference
  - physics_params: all simulation parameters
  - quality_gates: energy conservation thresholds
- Existing configs: josephson_junction.yaml, spin_magnetism.yaml, electroweak.yaml
- Results documented in config (see josephson_junction.yaml:171-195)

### 6. Quality Gates & Validation
- Energy conservation: ΔE/E < 0.01% (proven in benchmarks A2, A3, C1)
- Time reversibility: <1e-4 rad phase error after forward+backward evolution
- Symplectic structure: Verified via time-reversal test
- QA approval REQUIRED before any physics commit (@qa agent must validate)
- Test results documented in YAML config files with timestamps
- Reference implementations: TRDTestFramework.h, TRDValidationSuite.cpp

### 7. Architecture Consistency
- **CRITICAL**: Validated infrastructure MUST be reused, not reimplemented
- Wave 4 architecture violation example:
  - D4 scattering test used dissipative diffusion equation (85% energy loss)
  - Fix: Migrated to TRDCore3D + Sine-Gordon (0.127% drift - 670× improvement)
  - Root cause: Test bypassed proven symplectic TRDCore3D framework
- Prevention: All tests MUST use TRDEngine3D::runSimulation() entry point
- Documentation: B1_ARCHITECTURAL_FAILURE_ANALYSIS.md

### 8. Migration Protocol
- Legacy tests on smft binary MUST migrate to trd executable
- Verify zero standalone binaries: `find build/ -type f -executable | grep -v trd` should be empty
- Check CMakeLists.txt for test_* executables (should be removed)
- All new tests: implement via TRD test harness + YAML config
- Migration tracking in TODO.md under "Migration to Unified TRD"

### 9. Forbidden Patterns
- ❌ Standalone test binaries (./test_foo, ./test_bar)
- ❌ Isolated physics implementations (not using TRDCore3D base classes)
- ❌ Dissipative integrators for conservative physics (Forward Euler)
- ❌ Hardcoded physics parameters (use YAML configs)
- ❌ Manual output directory creation (use OutputManager class)
- ❌ RK4 integration (0.0002% drift exceeds 0.01% standard)
- ❌ Direct file I/O (use OutputManager for all data export)
- ❌ Custom wave equation solvers (use proven Sine-Gordon/Klein-Gordon)

### 10. Documentation Standards
- Test results MUST be documented in YAML configs with:
  - Timestamp of run
  - Energy conservation achieved
  - Time reversibility error
  - Pass/fail status
- Example: config/josephson_junction.yaml has complete results section
- Architecture failures documented: B1_ARCHITECTURAL_FAILURE_ANALYSIS.md
- Symplectic methods documented: SYMPLECTIC_INVESTIGATION_REPORT.md
- Update TODO.md validation status after test completion
- Maintain EM_INTEGRATION_REPORT.md for electromagnetic validations

### 11. Validation Hierarchy
- **Level 1**: Energy conservation < 0.01% (mandatory)
- **Level 2**: Time reversibility < 1e-4 rad (mandatory)
- **Level 3**: Symplectic structure preservation (mandatory for conservative systems)
- **Level 4**: Physical observable accuracy (domain-specific thresholds)
- All levels must pass for production approval

### 12. Code Review Requirements
- NO physics changes without:
  - Energy conservation benchmark
  - Time reversibility test
  - Comparison to analytical solution (when available)
  - QA agent approval (@qa must run validation suite)
- Pull requests must reference:
  - YAML config used for validation
  - Energy drift achieved
  - TODO.md item being addressed

**Rationale**: These standards emerged from systematic validation failures where tests bypassing the core framework showed 670× worse energy conservation. The unified TRD executable with mandatory framework integration ensures all physics simulations benefit from proven symplectic integrators and validated numerical methods. The 0.01% energy conservation threshold is our GO/NO-GO criterion for theoretical validity.