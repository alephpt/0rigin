# TRD Engine Developer Guide

## Overview

This guide covers the development workflow, build system, testing infrastructure, and contribution process for the TRD Engine. It is intended for developers who will be writing new physics tests, extending the core framework, or maintaining the codebase.

For physics theory, see [PHYSICS.md](./PHYSICS.md). For API details, see [API.md](./API.md).

---

## 1. Development Setup

### Prerequisites

```bash
# Ubuntu/Debian
sudo apt update
sudo apt install -y \
    build-essential \
    cmake \
    vulkan-tools \
    libvulkan-dev \
    libyaml-cpp-dev \
    libglm-dev \
    libsdl2-dev \
    libfftw3-dev \
    libomp-dev \
    glslang-tools \
    git

# Verify Vulkan
vulkaninfo | grep "Vulkan Instance Version"
```

### Building from Source

```bash
git clone https://github.com/alephpt/0rigin.git
cd 0rigin

mkdir build && cd build

# Release build (default, optimized)
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)

# Debug build (with symbols, assertions enabled)
cmake .. -DCMAKE_BUILD_TYPE=Debug
make -j$(nproc)
```

The single output binary is `build/bin/trd`.

### Verifying the Build

```bash
./bin/trd --help
./bin/trd --test ../config/weak_field_3d.yaml
# Expected: Energy drift < 0.01%
```

### Build Options

| Option | Default | Description |
|--------|---------|-------------|
| `CMAKE_BUILD_TYPE` | Release | Debug, Release, or RelWithDebInfo |
| OpenMP | Auto-detected | Multi-core CPU parallelism |

CMake finds Vulkan, SDL2, FFTW3, and yaml-cpp automatically. If any are missing, the configure step will report an error with instructions.

---

## 2. Project Architecture

### Directory Structure

```
0rigin/
  main.cpp                     # Entry point: routes to interactive or test mode
  src/
    TRDCore3D.cpp/.h           # 3D Kuramoto phase dynamics (vacuum solver)
    TRDEngine3D.cpp/.h         # GPU compute orchestration, dual-solver routing
    ConservativeSolver.cpp     # Conservative particle dynamics
    Dirac3D.cpp               # 3D+1 Dirac equation (FFT split-operator)
    Maxwell3D.cpp             # 3D Maxwell equations
    TRDCore.cpp               # 2D core with Stueckelberg EM coupling
    TRDEngine.cpp/.h          # 2D engine (test runner integration)
    DiracEvolution.cpp/.h     # 2D Dirac evolution
    GeodesicIntegrator.cpp    # Geodesic motion integrator
    TRDPipelineFactory.cpp/.h # Vulkan pipeline creation
    TRDBufferManager.cpp/.h   # Vulkan buffer management
    TRDCompute.cpp/.h         # Vulkan compute dispatch
    TRDDescriptorManager.cpp/.h # Vulkan descriptor management
    physics/
      StuckelbergEM.cpp       # Stueckelberg gauge-restored EM
    simulations/
      TRDTestRunner.cpp/.h    # Unified test execution framework
      TestConfig.cpp/.h       # YAML configuration loading
      ObservableComputer.cpp/.h # Centralized observable computation
    output/
      OutputManager.cpp/.h    # Data export and directory management
  include/
    TRDCore3D.h               # Public 3D core header
    ConservativeSolver.h      # Public conservative solver header
    Dirac3D.h                 # Public Dirac solver header
    Maxwell3D.h               # Public Maxwell solver header
    TRDCore.h                 # Public 2D core header
    physics/
      StuckelbergEM.h         # Stueckelberg EM header
      GaugeTheory.h           # Abstract gauge theory base
  test/
    test_*.cpp                # Physics test implementations (40+ files)
  config/
    *.yaml                    # YAML test configurations (40+ files)
  shaders/
    smft/
      *.comp                  # Vulkan GLSL compute shaders
  lib/
    Nova/                     # Vulkan graphics/compute engine
    imgui/                    # ImGui UI framework
  docs/
    API.md                    # API reference
    PHYSICS.md                # Physics theory and validation
    DEVELOPER_GUIDE.md        # This file
    reports/                  # Validation reports archive
```

### Dual-Solver Architecture

TRDEngine3D routes physics evolution through two solvers:

```
TRDEngine3D::runSimulation(dt)
    |
    +-- "vacuum_kuramoto" --> TRDCore3D
    |     Dissipative Kuramoto gradient flow
    |     Energy NOT conserved (thermodynamic)
    |
    +-- "particle_sine_gordon" --> ConservativeSolver
    |     Conservative Sine-Gordon solitons
    |     Energy conserved < 0.01%
    |
    +-- "particle_dirac" --> ConservativeSolver --> Dirac3D
    |     Conservative Dirac fermion
    |     Mass from uniform vacuum fields
    |
    +-- "coupled_vacuum_particle" --> TRDCore3D + ConservativeSolver
          1. Evolve vacuum (Kuramoto, dissipative)
          2. Extract R-field and theta-field
          3. Evolve particle (Dirac) with chiral mass from vacuum
```

### Build System

CMake manages the build with these targets:

- **Nova** (static library): Vulkan engine core
- **imgui** (static library): UI framework
- **Physics** (static library): Stueckelberg EM implementation
- **TRD** (executable): The unified `trd` binary, linking all libraries

All test source files in `test/` are compiled directly into the TRD executable, not as separate binaries. The CMakeLists.txt `TRD_SOURCES` list includes every test file.

Compute shaders in `shaders/smft/` are compiled to SPIR-V by `glslc` as a custom CMake target (`CompileShaders`).

---

## 3. Adding a New Physics Test

This is the most common development task. Follow these steps exactly.

### Step 1: Define the Physics

Before writing code, identify:
- What equation(s) are being solved
- Which solver to use (TRDCore3D for vacuum, ConservativeSolver for particles)
- What observables to measure
- What validation criteria apply (energy conservation threshold, expected values)

### Step 2: Create the YAML Configuration

Create `config/<test_name>.yaml`:

```yaml
# config/my_new_test.yaml
test_name: "My New Physics Test"
description: |
  Tests [specific physics phenomenon].
  Validates [what is being validated].

validation:
  framework: "Category - Test Name"
  test_file: "test/test_my_new_test.cpp"
  golden_key: "1 TRD unit = 246 GeV"

physics:
  grid_size: 64
  time_step: 0.01
  total_steps: 1000
  coupling_strength: 1.0

quality_gates:
  energy_conservation_threshold: 0.01   # 0.01% maximum drift
  time_reversibility_threshold: 1e-4    # radians
```

The YAML structure is flexible. Different tests use different parameter names depending on the physics. The key requirements are:
- `test_name` and `description` for documentation
- Physics parameters (no hardcoding allowed)
- Quality gates for validation

### Step 3: Write the Test Implementation

Create `test/test_my_new_test.cpp`. The test must use the framework classes:

```cpp
// test/test_my_new_test.cpp
#include "ConservativeSolver.h"
#include <iostream>
#include <cmath>
#include <cassert>

/**
 * Test: [Description of what this validates]
 *
 * Physics: [Equation being solved]
 * Method: [Integration method used]
 * Quality gate: Energy conservation < 0.01%
 */
void runMyNewTest() {
    // 1. Configure the solver
    ConservativeSolver solver;
    ConservativeSolver::Config config;
    config.nx = 64;
    config.ny = 64;
    config.nz = 64;
    config.dx = 1.0f;
    config.dt = 0.01f;
    config.method = ConservativeSolver::IntegrationMethod::STRANG_SPLITTING;
    config.spatial_order = ConservativeSolver::SpatialOrder::FOURTH_ORDER;
    solver.initialize(config);

    // 2. Set initial conditions
    solver.initializeGaussian(32.0f, 32.0f, 32.0f, 5.0f, 1.0f);

    // 3. Record initial state
    float E0 = solver.computeTotalEnergy();
    std::cout << "Initial energy: " << E0 << std::endl;

    // 4. Time evolution
    const int num_steps = 1000;
    for (int step = 0; step < num_steps; ++step) {
        solver.evolveSineGordon(config.dt);

        // Optional: periodic logging
        if (step % 100 == 0) {
            float drift = solver.measureEnergyDrift(E0);
            std::cout << "Step " << step
                      << " energy drift: " << drift * 100.0f << "%"
                      << std::endl;
        }
    }

    // 5. Validate
    float final_drift = solver.measureEnergyDrift(E0);
    std::cout << "Final energy drift: " << final_drift * 100.0f << "%"
              << std::endl;

    bool passed = final_drift < 0.0001f;  // < 0.01%
    std::cout << (passed ? "PASS" : "FAIL") << std::endl;
}
```

### Step 4: Register in CMakeLists.txt

Add the test file to the `TRD_SOURCES` list in `CMakeLists.txt`:

```cmake
set(TRD_SOURCES
    ...
    test/test_my_new_test.cpp
    ...
)
```

### Step 5: Wire into the Test Harness

In `main.cpp`, the test harness dispatches to tests based on the YAML config filename. Add your test to the dispatch table:

```cpp
// In the test dispatch logic of main.cpp
if (test_name == "my_new_test") {
    runMyNewTest();
}
```

### Step 6: Build, Run, and Validate

```bash
cd build
make -j$(nproc)
./bin/trd --test ../config/my_new_test.yaml
```

Check that:
- Energy drift < 0.01%
- No NaN values
- Observables match expected behavior
- Output files are generated in `output/`

### Step 7: Document Results

Update the YAML config with validation results:

```yaml
validation_results:
  timestamp: "2026-03-03T14:00:00"
  energy_drift_percent: 0.0042
  time_reversal_error_rad: 1.2e-9
  status: "PASS"
  notes: "Validated with 4th-order spatial, Strang splitting"
```

---

## 4. YAML Configuration Reference

### Full Configuration Structure

The `TestConfig` class supports the following YAML structure:

```yaml
test_name: "string"
description: "string"

grid:
  size_x: 64       # Grid width
  size_y: 64       # Grid height

physics:
  delta: 2.5       # Mass gap parameter
  coupling: 0.1    # Kuramoto-Dirac coupling
  dt: 0.01         # Timestep
  total_steps: 100 # Number of evolution steps
  K: 1.0           # Kuramoto coupling strength
  damping: 0.1     # Phase damping

  # EM configuration (optional)
  em_coupling_enabled: false
  em_coupling_strength: 1.0
  em_coupling_type: "stuckelberg"
  photon_mass: 0.0

initial_conditions:
  dirac:
    type: "gaussian"       # "gaussian", "plane_wave"
    x0: 32.0
    y0: 32.0
    sigma: 3.0
    amplitude: 1.0
  kuramoto:
    phase_distribution: "uniform"  # "uniform", "random", "vortex", "phase_gradient"
    omega_distribution: "zero"     # "zero", "gaussian", "random"
    omega_mean: 0.0
    omega_std: 0.1
    wave_vector_x: 0.0
    wave_vector_y: 0.0
  em_pulse:
    type: "none"                   # "none", "gaussian_pulse"
    component: "A_y"
    center_x: 32.0
    center_y: 32.0
    width_x: 8.0
    width_y: 8.0
    amplitude: 0.1
    wave_vector_x: 10.0
    wave_vector_y: 0.0

operator_splitting:
  enabled: false
  substep_ratios: [1, 10, 100]

validation:
  norm_tolerance: 1.0e-4
  energy_tolerance: 1.0e-2
  convergence_tolerance: 0.05
  expected_wave_speed: 1.0
  wave_speed_tolerance: 0.01

output:
  directory: "output/my_test"
  save_every: 10
  formats: ["csv"]
  auto_visualize: false
  track_wave_position: false
  measure_phase_velocity: false
```

Not all fields are required for every test. Many tests use only a subset plus custom physics parameters defined directly in their YAML.

---

## 5. The Test Framework

### TRDTestRunner

The `TRDTestRunner` class orchestrates test execution:

1. **Load**: Parse YAML configuration via `TestConfig`
2. **Initialize**: Create Nova instance, TRD engine, output directories
3. **Execute**: Run simulation(s) with configured parameters
4. **Validate**: Check norm conservation, energy conservation, convergence
5. **Report**: Generate text reports and CSV data

```cpp
TRDTestRunner runner("config/my_test.yaml");
runner.initialize();
runner.run();
runner.generateReport();
```

### ObservableComputer

Centralized observable computation ensures consistent measurements across all tests:

```cpp
ObservableComputer::Observables obs;
ObservableComputer::compute(&obs, dirac, R_field, delta, time, E0);

// Access results
double norm = obs.norm;           // Should be ~1.0
double energy = obs.energy_total; // Total system energy
double R_avg = obs.R_avg;        // Average synchronization
```

### OutputManager

All file I/O should go through `TRD::OutputManager`:

```cpp
TRD::OutputManager output;
std::string dir = output.createExperimentDirectory("my_test");

// Write timeseries
std::map<std::string, std::string> metadata;
metadata["grid_size"] = "64x64";
metadata["dt"] = "0.01";
output.writeTimeseries(dir + "/energy.dat", times, energies, metadata, "energy");

// Write CSV
output.writeCSV(dir + "/observables.csv", headers, data);
```

Output directories follow the naming convention: `output/YYYYMMDD_HHMMSS_<name>/`

---

## 6. Coding Standards

### C++ Style

```cpp
// File header
/**
 * @file MyNewSolver.cpp
 * @brief Description of what this file implements
 */

// Include guard (use #pragma once for headers)
#pragma once

// Class naming: PascalCase
class MyNewSolver {
public:
    // Method naming: camelCase
    void computeEnergy();

    // Constants: UPPER_SNAKE_CASE
    static constexpr float MAX_DRIFT = 0.0001f;

private:
    // Member variables: underscore prefix or suffix
    float _coupling_strength;
    std::vector<float> theta_;
};
```

### Size Limits

These are mandatory and enforced during review:

| Metric | Maximum |
|--------|---------|
| File length | 500 lines |
| Function length | 50 lines |
| Nesting depth | 3 levels |

When a file exceeds 500 lines, split it into logical modules. When a function exceeds 50 lines, extract subfunctions. When nesting exceeds 3 levels, restructure the logic.

### Quality Requirements

- Zero compiler warnings (build with `-Wall -Wextra -Wpedantic`)
- No hardcoded physics parameters (use YAML or constructor arguments)
- No `TODO` comments in production code
- Comprehensive error handling (no silent failures)
- Descriptive naming (no single-letter variables except loop indices)
- Doxygen comments on all public methods

### Energy Conservation Mandate

Every new physics implementation must include an energy conservation test demonstrating drift < 0.01%. This is not optional.

---

## 7. Debugging

### Common Issues

#### Energy Divergence

If energy grows exponentially:
1. Check that you are using a symplectic integrator (not Euler)
2. Reduce the timestep `dt`
3. Verify initial conditions (uninitialized velocity fields cause 164% drift; see `initializeVortexWithProperVelocity`)
4. Ensure the energy functional uses the same spatial order as the evolution operator

#### NaN Values

If fields become NaN:
1. Check for division by zero in coupling terms
2. Verify grid indices are within bounds (use `wrapX/Y/Z` for periodic BC)
3. Check that FFT plans are properly initialized (Dirac3D)
4. Reduce timestep -- CFL condition may be violated

#### Memory Issues

For large grids:
- 64^3 grid: ~1 MB per field, ~10-20 MB total
- 128^3 grid: ~8 MB per field, ~80-160 MB total
- Use Valgrind: `valgrind --leak-check=full ./bin/trd --test config/test.yaml`

### Debug Build

```bash
cmake .. -DCMAKE_BUILD_TYPE=Debug
make -j$(nproc)

# Run with gdb
gdb ./bin/trd
(gdb) run --test ../config/test.yaml
(gdb) bt    # Backtrace on crash
(gdb) print _theta_data[0]    # Inspect field values
```

### Address Sanitizer

```bash
cmake .. -DCMAKE_CXX_FLAGS="-fsanitize=address -fno-omit-frame-pointer"
make -j$(nproc)
./bin/trd --test ../config/test.yaml
```

This catches buffer overflows, use-after-free, and memory leaks at runtime.

---

## 8. Profiling

### CPU Profiling

```bash
# Using perf
perf record -g ./bin/trd --test ../config/test.yaml
perf report

# Using gprof
cmake .. -DCMAKE_CXX_FLAGS="-pg"
make -j$(nproc)
./bin/trd --test ../config/test.yaml
gprof ./bin/trd gmon.out > profile.txt
```

### GPU Profiling

```bash
# NVIDIA Nsight Systems
nsys profile ./bin/trd --test ../config/test.yaml
nsys-ui report.nsys-rep
```

### Optimization Tips

1. **Vectorization**: Use `#pragma omp simd` for inner loops over field arrays
2. **Cache locality**: Access arrays in row-major order (k outer, j middle, i inner)
3. **Shared memory**: GPU compute shaders should use shared memory tiles with halos for stencil operations
4. **Workgroup size**: Use 8x8x8 local workgroups for 3D compute shaders

---

## 9. Version Control

### Branch Strategy

```
main              # Stable, all tests pass
  feature/*       # New physics implementations
  fix/*           # Bug fixes
  validation/*    # Physics validation campaigns
```

### Commit Message Format

```
<type>(<scope>): <subject>

# Types: feat, fix, refactor, test, docs, perf
# Scope: the module affected

# Examples:
feat(dirac): implement chiral mass coupling with eigenvalue decomposition
fix(energy): resolve 4th-order discretization drift in Sine-Gordon solver
test(validation): add electroweak unification quality gates
docs(api): update ConservativeSolver method documentation
refactor(core): extract Laplacian computation to helper function
perf(fft): optimize FFT plan reuse in Dirac3D
```

### Pull Request Checklist

Before opening a PR for physics code:

- [ ] Builds without warnings (`-Wall -Wextra -Wpedantic`)
- [ ] All existing tests still pass
- [ ] New test(s) included with YAML config
- [ ] Energy conservation < 0.01% verified
- [ ] Time reversibility < 1e-4 rad verified
- [ ] Uses TRDCore3D/ConservativeSolver framework (no custom integrators)
- [ ] YAML configuration present (no hardcoded parameters)
- [ ] Files < 500 lines, functions < 50 lines, nesting < 3 levels
- [ ] Doxygen comments on all new public methods
- [ ] Validation results documented in YAML config

---

## 10. Troubleshooting

### Build Errors

| Error | Solution |
|-------|----------|
| `Could not find Vulkan` | `sudo apt install libvulkan-dev` |
| `yaml-cpp not found` | `sudo apt install libyaml-cpp-dev` |
| `FFTW3 not found` | `sudo apt install libfftw3-dev` |
| `SDL2 not found` | `sudo apt install libsdl2-dev` |
| `C++20 required` | Install GCC 9+ or Clang 10+: `sudo apt install g++-9` |
| `glslc not found` | `sudo apt install glslang-tools` (optional, for shader dev) |

### Runtime Errors

| Error | Solution |
|-------|----------|
| `No Vulkan devices available` | Check GPU drivers, or run in CPU fallback mode |
| `Energy conservation violated` | Reduce `dt` in YAML config |
| `NaN detected` | Check initial conditions and CFL stability |
| `Failed to allocate field arrays` | Reduce grid size (`nx`, `ny`, `nz`) |
| `Failed to load configuration` | Check YAML syntax; verify file path |

### Performance Issues

- **Slow evolution**: Check that you are building in Release mode, not Debug
- **High memory usage**: Reduce grid size or check for memory leaks with Valgrind
- **GPU not utilized**: Verify Vulkan device initialization; some tests run CPU-only

---

## 11. Key Design Decisions

### Why Symplectic Integrators?

Non-symplectic integrators (Euler, RK4) introduce systematic energy drift that accumulates over thousands of timesteps. A 2nd-order symplectic integrator has bounded energy error that oscillates rather than growing, making it superior for long-time physics simulations even though its local truncation error is larger than RK4.

Benchmarking confirmed: RK4 gave 0.0002% drift (accumulating), while Velocity Verlet gave 0.0038% drift (bounded). Over 100,000 steps, Velocity Verlet remains stable while RK4 diverges.

### Why a Single Executable?

Before consolidation, the project had 20+ standalone test binaries. This caused:
- Inconsistent framework usage (some tests bypassed TRDCore3D)
- 670x worse energy conservation in bypassing tests
- Build system complexity
- Difficulty running the full test suite

The unified `./trd --test` pattern ensures all tests use the validated framework.

### Why 4th-Order Spatial Discretization?

The 2nd-to-4th order upgrade provided 18x improvement in energy conservation at minimal computational cost (12 neighbors instead of 6 per dimension). The critical insight was that energy measurement must match evolution order -- using 4th-order evolution with 2nd-order energy creates artificial drift.

### Why Stueckelberg Over Proca?

The Proca approach to massive photons breaks gauge invariance, resulting in B-fields that are ~10^9 too weak. The Stueckelberg mechanism restores gauge invariance through a compensator scalar phi that couples directly to the TRD phase field theta. This produces physically realistic electromagnetic fields.

---

For physics theory: [PHYSICS.md](./PHYSICS.md)
For API reference: [API.md](./API.md)
For contribution guidelines: [CONTRIBUTING.md](../CONTRIBUTING.md)
