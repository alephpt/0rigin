# TRD Engine Developer Guide

## Overview

This guide covers development practices, build system, testing procedures, and contribution workflow for the TRD Engine project.

## Development Setup

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
    libomp-dev \
    git

# Verify Vulkan
vulkaninfo | grep "Vulkan Instance Version"
```

### Building from Source

```bash
# Clone repository
git clone https://github.com/alephpt/0rigin.git
cd 0rigin

# Create build directory
mkdir build && cd build

# Configure with CMake
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DENABLE_VULKAN=ON \
    -DENABLE_OPENMP=ON

# Build (parallel)
make -j$(nproc)

# Run tests
./bin/trd --test config/trdcore_symplectic.yaml
```

### Build Options

| Option | Default | Description |
|--------|---------|-------------|
| `CMAKE_BUILD_TYPE` | Release | Debug/Release/RelWithDebInfo |
| `ENABLE_VULKAN` | ON | GPU compute support |
| `ENABLE_OPENMP` | ON | Multi-core CPU support |
| `BUILD_TESTS` | ON | Build test suite |
| `ENABLE_PROFILING` | OFF | Performance profiling |

## Code Organization

### Directory Structure

```
0rigin/
├── src/                    # Core implementation
│   ├── TRDEngine3D.cpp    # Main orchestration
│   ├── TRDCore3D.cpp      # Base field dynamics
│   ├── ConservativeSolver.cpp  # Conservative physics
│   └── Dirac3D.cpp        # Dirac equation solver
├── include/                # Public headers
│   ├── TRDEngine3D.h
│   ├── TRDCore3D.h
│   └── ConservativeSolver.h
├── test/                   # Test implementations
│   ├── test_sine_gordon.cpp
│   └── test_dirac_vacuum.cpp
├── shaders/               # Vulkan compute shaders
│   ├── kuramoto.comp
│   └── dirac_velocity_verlet.comp
├── config/                # YAML test configurations
│   ├── weak_field_3d.yaml
│   └── electroweak.yaml
├── docs/                  # Documentation
│   ├── API.md
│   └── PHYSICS.md
└── build/                 # Build artifacts
    └── bin/
        └── trd           # Single executable
```

## Coding Standards

### C++ Style Guide

```cpp
// File header (mandatory)
/**
 * @file ConservativeSolver.cpp
 * @brief Conservative physics solver with symplectic integration
 * @author TRD Team
 * @date 2026-01-17
 */

// Include guard
#ifndef CONSERVATIVE_SOLVER_H
#define CONSERVATIVE_SOLVER_H

// Class naming: PascalCase
class ConservativeSolver : public TRDCore3D {
public:
    // Method naming: camelCase
    void stepSineGordon(float dt);

    // Constants: UPPER_SNAKE_CASE
    static constexpr float ENERGY_THRESHOLD = 0.0001f;

private:
    // Member variables: m_prefix
    float m_sineGordonCoupling;
    std::vector<float> m_thetaField;

    // Private methods: underscore prefix
    void _computeLaplacian();
};
```

### Quality Requirements

1. **File Length**: ≤500 lines per file
2. **Function Length**: ≤50 lines per function
3. **Nesting Depth**: ≤3 levels
4. **Cyclomatic Complexity**: ≤10
5. **No Magic Numbers**: Use named constants
6. **No Compiler Warnings**: Build must be warning-free

### Energy Conservation Mandate

**ABSOLUTE REQUIREMENT**: All physics implementations MUST demonstrate:

```cpp
// Energy conservation test
float E0 = solver.computeEnergy();
for (int i = 0; i < 10000; i++) {
    solver.step(0.01f);
}
float E1 = solver.computeEnergy();

float drift = std::abs(E1 - E0) / E0;
ASSERT(drift < 0.0001);  // <0.01% required
```

## Testing Framework

### Writing Tests

All tests MUST use the TRDCore3D framework:

```cpp
// test/test_my_physics.cpp
#include "TRDEngine3D.h"
#include <cassert>

void testMyPhysics() {
    // Setup
    TRDEngine3D engine;
    engine.loadConfig("config/my_physics.yaml");
    engine.initialize();

    // Record initial state
    float initialEnergy = engine.computeEnergy();

    // Run simulation
    for (int i = 0; i < 1000; i++) {
        engine.step(0.01f);
    }

    // Validate energy conservation
    float finalEnergy = engine.computeEnergy();
    float drift = std::abs(finalEnergy - initialEnergy) / initialEnergy;

    assert(drift < 0.0001);  // <0.01% mandatory

    std::cout << "Test passed. Energy drift: "
              << drift * 100 << "%\n";
}

// Register with test harness
int main(int argc, char** argv) {
    testMyPhysics();
    return 0;
}
```

### YAML Test Configuration

```yaml
# config/my_physics.yaml
test_name: "My Physics Test"
test_file: "test/test_my_physics.cpp"
golden_key: 246.0  # GeV

physics_params:
  # Grid setup
  nx: 64
  ny: 64
  nz: 64
  dx: 0.1
  dy: 0.1
  dz: 0.1

  # Time evolution
  dt: 0.01
  num_steps: 10000

  # Physics parameters
  coupling_constant: 1.0
  mass_scale: 1.0

quality_gates:
  energy_conservation_threshold: 0.0001  # 0.01%
  time_reversibility_threshold: 1e-4     # rad

validation_results:
  # Filled after test run
  timestamp: null
  energy_drift: null
  time_reversal_error: null
  status: null
```

### Running Tests

```bash
# Single test
./bin/trd --test config/my_physics.yaml

# All tests
./bin/trd --test-all

# Specific category
./bin/trd --test-category "Standard Model"

# With profiling
./bin/trd --test config/test.yaml --profile
```

## Debugging

### Common Issues

#### 1. Energy Divergence

```cpp
// Problem: Energy grows exponentially
// Solution: Reduce timestep
float dt = 0.01;
while (energyDrift > 0.0001) {
    dt *= 0.5;
    // Re-run with smaller dt
}
```

#### 2. NaN Detection

```cpp
// Add NaN checks
if (std::isnan(field[i])) {
    std::cerr << "NaN detected at index " << i << "\n";
    std::cerr << "Previous value: " << field_prev[i] << "\n";
    abort();
}
```

#### 3. Memory Leaks

```bash
# Use valgrind
valgrind --leak-check=full ./bin/trd --test config/test.yaml
```

### Debug Build

```bash
# Configure debug build
cmake .. -DCMAKE_BUILD_TYPE=Debug

# Enable sanitizers
cmake .. -DENABLE_SANITIZERS=ON

# Run with gdb
gdb ./bin/trd
(gdb) run --test config/test.yaml
(gdb) bt  # Backtrace on crash
```

## Performance Optimization

### Profiling

```bash
# CPU profiling with perf
perf record -g ./bin/trd --test config/test.yaml
perf report

# GPU profiling with nsight
nsys profile ./bin/trd --test config/test.yaml
```

### Optimization Techniques

1. **Vectorization**
```cpp
// Enable auto-vectorization
#pragma omp simd
for (int i = 0; i < n; i++) {
    result[i] = a[i] + b[i];
}
```

2. **Cache Optimization**
```cpp
// Access patterns matter
// Good: Sequential access
for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
        field[i][j] = ...

// Bad: Strided access
for (int j = 0; j < ny; j++)
    for (int i = 0; i < nx; i++)
        field[i][j] = ...
```

3. **GPU Optimization**
```glsl
// shaders/optimized.comp
layout(local_size_x = 8, local_size_y = 8, local_size_z = 8) in;

// Shared memory for tile
shared float tile[10][10][10];  // +1 for halo

void main() {
    // Load to shared memory
    tile[lx][ly][lz] = field[gx][gy][gz];
    barrier();

    // Compute using shared memory
    float laplacian = ...
}
```

## Version Control

### Branch Strategy

```bash
main              # Stable, validated physics
├── develop       # Integration branch
├── feature/*     # New features
├── fix/*         # Bug fixes
└── validation/*  # Physics validation
```

### Commit Messages

```bash
# Format
<type>(<scope>): <subject>

# Examples
feat(dirac): Implement chiral mass coupling
fix(energy): Resolve 4th-order discretization drift
test(validation): Add Sine-Gordon soliton stability test
docs(api): Update ConservativeSolver documentation
```

### Pull Request Process

1. **Create feature branch**
```bash
git checkout -b feature/my-feature
```

2. **Implement with tests**
```bash
# Write code
vim src/MyFeature.cpp

# Write test
vim test/test_my_feature.cpp

# Add YAML config
vim config/my_feature.yaml
```

3. **Validate locally**
```bash
# Build and test
make -j$(nproc)
./bin/trd --test config/my_feature.yaml

# Verify energy conservation
grep "Energy drift" output/*/console.log
```

4. **Submit PR**
```bash
git push origin feature/my-feature
# Create PR on GitHub with:
# - Description of physics
# - Energy conservation results
# - Test configuration used
```

## Continuous Integration

### GitHub Actions Workflow

```yaml
# .github/workflows/physics-validation.yml
name: Physics Validation

on: [push, pull_request]

jobs:
  validate:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2

      - name: Build
        run: |
          mkdir build && cd build
          cmake ..
          make -j$(nproc)

      - name: Test Energy Conservation
        run: |
          ./build/bin/trd --test-all
          grep "PASS" results.log

      - name: Validate <0.01% drift
        run: |
          python3 scripts/validate_energy.py results.log
```

## Documentation

### Code Documentation

```cpp
/**
 * @brief Compute 4th-order Laplacian
 *
 * Uses 18-neighbor stencil for O(dx^4) accuracy.
 * Energy conservation improved by 18x over 2nd-order.
 *
 * @param field Input field
 * @param laplacian Output Laplacian
 * @param dx Spatial step size
 *
 * @note Boundary conditions: Periodic
 * @see ConservativeSolver::stepSineGordon
 */
void computeLaplacian4thOrder(
    const std::vector<float>& field,
    std::vector<float>& laplacian,
    float dx
);
```

### Physics Documentation

Document in YAML configs:

```yaml
# config/test.yaml
test_description: |
  Tests Sine-Gordon equation for topological soliton stability.
  Validates energy conservation <0.01% over 10,000 timesteps.

  Physics:
  - Equation: ∂²θ/∂t² = ∇²θ - sin(θ)
  - Soliton solution: θ = 4·atan(exp((x-vt)/√(1-v²)))
  - Conserved: Energy, momentum, topological charge

expected_results:
  energy_drift: <0.0001
  soliton_velocity: 0.5 ± 0.01
  topological_charge: 1.0 ± 0.001
```

## Troubleshooting

### Build Errors

```bash
# Missing Vulkan
CMake Error: Could not find Vulkan
Solution: sudo apt install libvulkan-dev

# Missing yaml-cpp
CMake Error: yaml-cpp not found
Solution: sudo apt install libyaml-cpp-dev

# Compiler too old
Error: C++17 required
Solution: sudo apt install g++-9
export CXX=g++-9
```

### Runtime Errors

```bash
# Vulkan device not found
Error: No Vulkan devices available
Solution: Check GPU drivers, use CPU fallback

# Energy divergence
Error: Energy conservation violated: 0.15%
Solution: Reduce timestep in YAML config

# Out of memory
Error: Failed to allocate field arrays
Solution: Reduce grid size (nx, ny, nz)
```

## Best Practices Summary

1. **Always use TRDCore3D framework** - No custom integrators
2. **Validate energy conservation** - Must be <0.01%
3. **Use YAML configuration** - No hardcoded parameters
4. **Write comprehensive tests** - Include time reversibility
5. **Document physics** - Equations and expected behavior
6. **Follow coding standards** - Clean, readable, maintainable
7. **Profile before optimizing** - Measure, don't guess

---

For physics theory: [PHYSICS.md](./PHYSICS.md)
For API reference: [API.md](./API.md)
For contribution: [CONTRIBUTING.md](../CONTRIBUTING.md)