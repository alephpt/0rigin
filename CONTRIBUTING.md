# Contributing to TRD Engine

## Project Overview

TRD (Topological Resonance Dynamics) Engine is a C++20 physics simulation framework exploring mass generation through topological phase synchronization. The core equation is M = Delta R e^(i theta gamma^5), where synchronized vacuum fields generate fermion masses dynamically.

The project builds a single unified executable (`./trd`) that serves as both an interactive visualization tool and a comprehensive physics test harness. All physics tests run through YAML configuration files.

**Theory paper**: "Mass Synchronization Field Theory" by R. Christopher

## Prerequisites

### Required

- **C++20 compiler**: GCC 9+ or Clang 10+
- **CMake**: 3.20 or later
- **Vulkan SDK**: 1.3+ (GPU compute backend)
- **SDL2**: Window management and input
- **yaml-cpp**: YAML configuration parsing
- **FFTW3**: Fast Fourier transforms (Dirac solver, single-precision: `fftw3f`)
- **GLM**: OpenGL Mathematics library

### Optional

- **OpenMP**: Multi-core CPU parallelism (HPC scaling tests)
- **glslc**: GLSL to SPIR-V shader compiler (for compute shader development)
- **Valgrind**: Memory debugging

### Ubuntu/Debian Installation

```bash
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

# Verify Vulkan installation
vulkaninfo | grep "Vulkan Instance Version"
```

## Building

```bash
git clone https://github.com/alephpt/0rigin.git
cd 0rigin

mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

The build produces a single executable at `build/bin/trd`.

### Build Options

| Option | Default | Description |
|--------|---------|-------------|
| `CMAKE_BUILD_TYPE` | Release | Build type: Debug, Release, RelWithDebInfo |
| `ENABLE_OPENMP` | auto | OpenMP multi-core CPU support |

### Verifying the Build

```bash
./bin/trd --help
./bin/trd --test ../config/weak_field_3d.yaml
```

## Running Tests

All physics tests are invoked through the unified executable with a YAML configuration file.

```bash
# From the build directory:
./bin/trd --test ../config/<test_name>.yaml

# Examples:
./bin/trd --test ../config/weak_field_3d.yaml
./bin/trd --test ../config/electroweak.yaml
./bin/trd --test ../config/josephson_junction.yaml
./bin/trd --test ../config/causality.yaml
```

Output is written to `output/YYYYMMDD_HHMMSS_<test_name>/` containing:
- `console.log` -- full simulation log
- `observables.csv` -- time series of physical observables
- `energy_history.csv` -- energy conservation tracking
- `test_report.txt` -- pass/fail validation results

## Code Structure

```
0rigin/
  src/
    TRDCore3D.cpp          # 3D field dynamics (Kuramoto phase evolution)
    TRDEngine3D.cpp/.h     # GPU compute orchestration, dual-solver routing
    ConservativeSolver.cpp  # Conservative particle dynamics (Sine-Gordon, Dirac)
    Dirac3D.cpp            # 3D+1 Dirac equation solver (split-operator FFT)
    Maxwell3D.cpp          # 3D Maxwell electromagnetic solver
    TRDCore.cpp            # 2D TRD core with Stueckelberg EM coupling
    TRDEngine.cpp/.h       # 2D engine (test runner integration)
    DiracEvolution.cpp     # 2D Dirac evolution
    GeodesicIntegrator.cpp # Geodesic motion in curved spacetime
    physics/
      StuckelbergEM.cpp    # Stueckelberg gauge-restored EM fields
    simulations/
      TRDTestRunner.cpp/.h     # Unified test execution framework
      TestConfig.cpp/.h        # YAML configuration loading and validation
      ObservableComputer.cpp/.h # Centralized observable computation
    output/
      OutputManager.cpp/.h     # Systematic data export
  include/
    TRDCore3D.h            # 3D grid infrastructure, phase field storage
    ConservativeSolver.h   # Conservative solver API
    Dirac3D.h              # 3D Dirac equation API
    Maxwell3D.h            # 3D Maxwell equations API
    TRDCore.h              # 2D core with EM coupling
    physics/
      StuckelbergEM.h      # Stueckelberg EM interface
      GaugeTheory.h        # Abstract gauge theory base class
  test/
    test_*.cpp             # Physics test implementations
  config/
    *.yaml                 # YAML test configurations
  shaders/
    smft/
      kuramoto3d.comp      # 3D Kuramoto phase evolution (GPU)
      dirac_velocity_verlet.comp  # Dirac spinor dynamics (GPU)
      gravity_field.comp   # Gravitational field computation (GPU)
      ...                  # Additional compute shaders
  lib/
    Nova/                  # Vulkan graphics/compute engine
    imgui/                 # ImGui UI framework
  docs/
    API.md                 # API reference
    PHYSICS.md             # Physics theory and validation
    DEVELOPER_GUIDE.md     # Development guide
```

## Coding Standards

### Size Limits (Mandatory)

| Metric | Limit |
|--------|-------|
| File length | 500 lines maximum |
| Function length | 50 lines maximum |
| Nesting depth | 3 levels maximum |

If a file or function exceeds these limits, refactor into smaller modules or extracted subfunctions.

### C++ Style

- **Class names**: PascalCase (`ConservativeSolver`, `TRDCore3D`)
- **Method names**: camelCase (`computeEnergy`, `evolveSymplecticCPU`)
- **Constants**: UPPER_SNAKE_CASE (`ENERGY_THRESHOLD`)
- **Member variables**: underscore prefix or suffix (`_theta_data`, `nx_`)
- **Include guards**: `#pragma once`
- **Documentation**: Doxygen-style comments for all public APIs

```cpp
/**
 * @brief Compute 4th-order Laplacian using 12-neighbor stencil
 *
 * Error: O(dx^4) - provides 16x better accuracy than 2nd-order.
 *
 * @param field Input field array
 * @param i X-index
 * @param j Y-index
 * @param k Z-index
 * @return 4th-order accurate Laplacian value
 */
float computeLaplacian4thOrder(
    const std::vector<float>& field, int i, int j, int k) const;
```

### Zero Tolerance Policy

These are never acceptable in production code:

- Duplicate files or components
- Stub implementations or mock data in production paths
- Hardcoded credentials or secrets
- Orphaned code or commented-out blocks
- Unhandled edge cases
- Copy-paste code duplication
- Files or functions exceeding size limits
- Compiler warnings

## Physics Standards

### Framework Integration (Mandatory)

All physics simulations must use the TRDCore3D / TRDEngine3D infrastructure. Custom integrators that bypass the validated framework are prohibited. The framework provides proven symplectic integrators with validated energy conservation.

### Symplectic Integration (Required)

Conservative physics must use one of the approved symplectic integration methods:

| Method | Use Case | Implementation |
|--------|----------|----------------|
| **RK2 Midpoint** | First-order field evolution (Kuramoto) | `TRDCore3D::evolveSymplecticCPU()` |
| **Velocity Verlet** | Wave equations (Sine-Gordon, Klein-Gordon) | `ConservativeSolver::velocityVerletStep()` |
| **Strang Splitting** | Nonlinear PDEs (Sine-Gordon with strong nonlinearity) | `ConservativeSolver::strangSplittingStep()` |
| **Half-Strang** | Phase-magnitude decoupling | `ConservativeSolver` |
| **Split-Operator** | Dirac equation (kinetic-mass splitting) | `Dirac3D::step()` |

### Forbidden Integrators

- **Forward Euler**: Dissipative; destroys energy conservation (becomes heat equation)
- **RK4**: 0.0002% drift exceeds the 0.01% standard; rejected after benchmarking

### Energy Conservation (GO/NO-GO Criterion)

All conservative physics tests must demonstrate energy drift below 0.01%:

```
Delta_E / E < 0.0001  (0.01%)
```

This is a hard gate. Tests that fail this criterion do not pass.

### Time Reversibility

Symplectic integrators must demonstrate time reversibility with phase error below 1e-4 radians after forward-backward evolution.

### Validation Hierarchy

1. **Level 1**: Energy conservation < 0.01% (mandatory)
2. **Level 2**: Time reversibility < 1e-4 rad (mandatory)
3. **Level 3**: Symplectic structure preservation (mandatory for conservative systems)
4. **Level 4**: Physical observable accuracy (domain-specific thresholds)

## YAML Configuration Format

All tests are configured through YAML files in the `config/` directory. No physics parameters should be hardcoded in test implementations.

```yaml
# config/my_test.yaml
test_name: "My Physics Test"
description: "Description of what this test validates"

# Required
validation:
  framework: "Category - Test Name"
  test_file: "test/test_my_physics.cpp"
  golden_key: "1 TRD unit = 246 GeV"

# Grid configuration
grid:
  size_x: 64
  size_y: 64

# Physics parameters
physics:
  delta: 2.5
  coupling: 0.1
  dt: 0.01
  total_steps: 1000
  coupling_strength: 1.0

# Validation thresholds
validation:
  norm_tolerance: 1.0e-4
  energy_tolerance: 1.0e-2

# Output configuration
output:
  directory: "output/my_test"
  save_every: 10
```

## How to Add a New Physics Test

### Step 1: Create the Test Implementation

Write a test file in `test/` that uses the TRDCore3D/TRDEngine3D framework:

```cpp
// test/test_my_physics.cpp
#include "TRDEngine3D.h"
#include "ConservativeSolver.h"
#include <cassert>
#include <cmath>
#include <iostream>

void runMyPhysicsTest(/* parameters from YAML */) {
    // Initialize using framework classes
    ConservativeSolver solver;
    ConservativeSolver::Config config;
    config.nx = 64;
    config.ny = 64;
    config.nz = 64;
    config.method = ConservativeSolver::IntegrationMethod::VELOCITY_VERLET;
    config.spatial_order = ConservativeSolver::SpatialOrder::FOURTH_ORDER;
    solver.initialize(config);

    // Set initial conditions
    solver.initializeGaussian(32.0f, 32.0f, 32.0f, 5.0f, 1.0f);

    // Record initial energy
    float E0 = solver.computeTotalEnergy();

    // Evolve
    for (int step = 0; step < 1000; ++step) {
        solver.evolveSineGordon(0.01f);
    }

    // Validate energy conservation
    float drift = solver.measureEnergyDrift(E0);
    std::cout << "Energy drift: " << drift * 100.0f << "%" << std::endl;
    assert(drift < 0.0001f);  // <0.01% required
}
```

### Step 2: Create the YAML Configuration

Add a configuration file in `config/`:

```yaml
# config/my_physics.yaml
test_name: "My Physics Validation"
description: "Validates [specific physics] using [method]"

validation:
  framework: "Category - Test Name"
  test_file: "test/test_my_physics.cpp"
  golden_key: "1 TRD unit = 246 GeV"

physics:
  grid_size: 64
  time_step: 0.01
  coupling_strength: 1.0

quality_gates:
  energy_conservation_threshold: 0.01
  time_reversibility_threshold: 1e-4
```

### Step 3: Register in CMakeLists.txt

Add the test source file to the `TRD_SOURCES` list in `CMakeLists.txt`:

```cmake
set(TRD_SOURCES
    ...
    test/test_my_physics.cpp
    ...
)
```

### Step 4: Build and Validate

```bash
cd build
make -j$(nproc)
./bin/trd --test ../config/my_physics.yaml
```

Verify that:
- Energy drift is below 0.01%
- Time reversibility error is below 1e-4 rad
- All quality gates pass

### Step 5: Document Results

Add validation results to the YAML configuration file:

```yaml
validation_results:
  timestamp: "2026-03-03T12:00:00"
  energy_drift: 0.0038
  time_reversal_error: 1.0e-9
  status: "PASS"
```

## Forbidden Patterns

- **Standalone test binaries**: All tests run through `./trd --test`; no `test_*` executables
- **Forward Euler integration**: Dissipative; violates energy conservation
- **RK4 integration**: Drift exceeds the 0.01% threshold
- **Hardcoded physics parameters**: All parameters must come from YAML configs
- **Custom integrators bypassing the framework**: Use TRDCore3D/ConservativeSolver
- **Direct file I/O**: Use OutputManager for all data export
- **Custom wave equation solvers**: Use the proven Sine-Gordon/Klein-Gordon implementations

## Git Workflow

### Branch Strategy

```
main              # Stable, validated physics
  feature/*       # New features
  fix/*           # Bug fixes
  validation/*    # Physics validation runs
```

### Commit Messages

```
<type>(<scope>): <description>

# Examples:
feat(dirac): implement chiral mass coupling with eigenvalue decomposition
fix(energy): resolve 4th-order discretization drift in Sine-Gordon solver
test(validation): add electroweak unification quality gates
docs(api): update ConservativeSolver method documentation
refactor(core): extract Laplacian computation into separate method
```

Types: `feat`, `fix`, `refactor`, `test`, `docs`, `perf`

### Pull Request Requirements

Every pull request that modifies physics code must include:

1. **Energy conservation benchmark**: Measured drift with the relevant test configuration
2. **Time reversibility test**: Forward-backward evolution error measurement
3. **YAML config reference**: Which configuration file was used for validation
4. **Comparison to analytical solution** (when available)

Before submitting:

- [ ] All tests pass locally
- [ ] No compiler warnings
- [ ] Files < 500 lines, functions < 50 lines, nesting < 3 levels
- [ ] Energy conservation < 0.01% verified
- [ ] Framework integration confirmed (uses TRDCore3D/TRDEngine3D)
- [ ] YAML configuration present and valid
- [ ] Documentation updated where needed

## Getting Help

- **Questions**: Open a GitHub issue with the `question` label
- **Bugs**: Open an issue with the `bug` label and include reproduction steps
- **Feature requests**: Open an issue with the `enhancement` label

## License

Research use -- 0rigin Project. By contributing, you agree that your contributions will be licensed under the same terms.
