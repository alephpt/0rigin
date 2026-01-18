# TRD Engine API Reference

## Overview

The TRD Engine provides a C++ API for simulating topological field theories with validated energy conservation. All simulations must use the TRDEngine3D/TRDCore3D framework.

## Core Classes

### TRDEngine3D

Main orchestration class for running simulations.

```cpp
class TRDEngine3D {
public:
    // Lifecycle
    TRDEngine3D();
    ~TRDEngine3D();

    // Configuration
    bool loadConfig(const std::string& yamlPath);
    void initialize();

    // Simulation
    void runSimulation();
    void step(float dt);

    // Observables
    float computeEnergy() const;
    float computeHelicity() const;
    std::vector<float> getObservables() const;

    // Output
    void saveState(const std::string& filename);
    void exportCSV(const std::string& path);
};
```

#### Usage Example

```cpp
#include "TRDEngine3D.h"

int main() {
    TRDEngine3D engine;
    engine.loadConfig("config/my_test.yaml");
    engine.initialize();
    engine.runSimulation();
    return 0;
}
```

### TRDCore3D

Base class for 3D field dynamics with symplectic integration.

```cpp
class TRDCore3D {
protected:
    // Field storage
    std::vector<std::complex<float>> psi1, psi2;  // Complex fields
    std::vector<float> E_field, B_field;           // EM fields
    std::vector<float> theta_field, R_field;       // Phase-magnitude

    // Grid parameters
    int nx, ny, nz;
    float dx, dy, dz;

public:
    // Evolution methods (CPU)
    void evolveSymplecticCPU(float dt);
    void evolveSineGordonCPU(float dt);
    void evolveKuramotoCPU(float dt);
    void evolveMaxwell3DCPU(float dt);
    void evolveDirac3DCPU(float dt);

    // GPU compute
    void evolveGPU(float dt);

    // Energy computation
    virtual float computeEnergy() const;
    virtual float computeHelicity() const;
};
```

### ConservativeSolver

Specialized solver for conservative physics with 4th-order spatial discretization.

```cpp
class ConservativeSolver : public TRDCore3D {
public:
    ConservativeSolver(int nx, int ny, int nz);

    // Wave equation solvers
    void stepSineGordon(float dt);
    void stepKleinGordon(float dt);

    // Spatial operators
    void computeLaplacian4thOrder();
    void computeGradient4thOrder();

    // Energy functionals
    float computeSineGordonEnergy() const;
    float computeKleinGordonEnergy() const;

    // Configuration
    void setSineGordonCoupling(float g);
    void setKleinGordonMass(float m);
};
```

#### Key Methods

**stepSineGordon(dt)**
- Velocity Verlet integration for Sine-Gordon equation
- Energy conservation: <0.01%
- Supports topological solitons

**computeLaplacian4thOrder()**
- 4th-order accurate spatial discretization
- 18-neighbor stencil
- O(dx⁴) truncation error

### Dirac3D

Dirac equation solver with chiral mass coupling.

```cpp
class Dirac3D : public TRDCore3D {
public:
    Dirac3D(int nx, int ny, int nz);

    // Main evolution
    void stepWithChiralMass(float dt);

    // Mass operator
    void computeMassDerivative();
    void setChiralCoupling(float Delta);

    // Observables
    float computeProbability() const;
    std::complex<float> computeChiralCondensate() const;

private:
    // Eigenvalue decomposition for stability
    void eigenvalueDecomposition(
        const std::array<std::complex<float>, 16>& M,
        std::array<float, 4>& eigenvalues,
        std::array<std::complex<float>, 16>& eigenvectors
    );
};
```

#### Chiral Mass Coupling

The mass operator implements M = Δ·R·e^(iθγ⁵):

```cpp
// Set coupling strength
dirac.setChiralCoupling(0.1f);  // Δ parameter

// Evolution preserves unitarity
dirac.stepWithChiralMass(dt);

// Measure observables
float probability = dirac.computeProbability();  // Should be 1.0
```

## Configuration System

### YAML Configuration

All simulations configured via YAML files:

```yaml
# Required fields
test_name: "My Physics Test"
test_file: "test/my_test.cpp"
golden_key: 246.0  # GeV

# Physics parameters
physics_params:
  nx: 64
  ny: 64
  nz: 64
  dx: 0.1
  dy: 0.1
  dz: 0.1
  dt: 0.01
  num_steps: 10000

  # Field-specific
  sine_gordon_coupling: 1.0
  klein_gordon_mass: 1.0
  dirac_chiral_coupling: 0.1

# Quality gates
quality_gates:
  energy_conservation_threshold: 0.01  # Percentage
  time_reversibility_threshold: 1e-4   # Radians

# Output configuration
output:
  directory: "output/"
  save_interval: 100
  export_csv: true
  visualize: false
```

### Loading Configuration

```cpp
TRDEngine3D engine;
if (!engine.loadConfig("config/test.yaml")) {
    std::cerr << "Failed to load configuration\n";
    return 1;
}
```

## Output Management

### CSV Export

Automatic CSV export for observables:

```cpp
// Exported to output/YYYYMMDD_HHMMSS_test_name/
observables.csv:
  - time
  - energy
  - helicity
  - custom observables...

energy_history.csv:
  - step
  - total_energy
  - kinetic_energy
  - potential_energy
  - energy_drift_percent
```

### State Serialization

```cpp
// Save current state
engine.saveState("checkpoint.bin");

// Load state (in new session)
TRDEngine3D engine2;
engine2.loadState("checkpoint.bin");
engine2.runSimulation();  // Continue from checkpoint
```

## Integration Methods

### Symplectic Integrators

All conservative physics uses symplectic methods:

```cpp
// RK2 Midpoint (default)
void evolveRK2(float dt) {
    auto k1 = computeDerivative(state);
    auto k2 = computeDerivative(state + dt/2 * k1);
    state += dt * k2;
}

// Velocity Verlet (wave equations)
void evolveVelocityVerlet(float dt) {
    x += dt * v + 0.5 * dt * dt * a;
    float a_new = computeAcceleration(x);
    v += 0.5 * dt * (a + a_new);
    a = a_new;
}
```

### Energy Conservation

Monitor energy drift:

```cpp
float initial_energy = engine.computeEnergy();
engine.runSimulation();
float final_energy = engine.computeEnergy();

float drift = abs(final_energy - initial_energy) / initial_energy;
assert(drift < 0.0001);  // Must be <0.01%
```

## GPU Compute

### Vulkan Backend

GPU acceleration via Vulkan compute shaders:

```cpp
// Automatic GPU selection
engine.enableGPU(true);

// Manual device selection
VulkanCompute compute;
compute.selectDevice(0);  // Use first GPU
engine.setComputeBackend(&compute);
```

### Compute Shaders

Located in `shaders/`:
- `kuramoto.comp` - Phase synchronization
- `maxwell3d.comp` - Electromagnetic evolution
- `dirac_velocity_verlet.comp` - Spinor dynamics
- `sine_gordon.comp` - Topological solitons

## Best Practices

### 1. Always Use Framework

```cpp
// ✅ CORRECT: Use framework
TRDEngine3D engine;
engine.loadConfig("config.yaml");
engine.runSimulation();

// ❌ WRONG: Custom integrator
for (int i = 0; i < steps; i++) {
    field += dt * derivative;  // NO! Violates energy conservation
}
```

### 2. Validate Energy Conservation

```cpp
// After simulation
float drift = engine.getEnergyDrift();
if (drift > 0.0001) {
    std::cerr << "Energy conservation violated: "
              << drift * 100 << "%\n";
    return 1;
}
```

### 3. Use YAML Configuration

```cpp
// ✅ CORRECT: Configuration file
engine.loadConfig("config/test.yaml");

// ❌ WRONG: Hardcoded parameters
engine.nx = 64;  // NO! Use YAML
engine.dt = 0.01;  // NO! Use YAML
```

### 4. Check Quality Gates

```cpp
// Automatic validation
engine.setQualityGate("energy_conservation", 0.0001);
engine.setQualityGate("time_reversibility", 1e-4);

if (!engine.passesQualityGates()) {
    std::cerr << "Quality gates failed\n";
    return 1;
}
```

## Error Handling

### Common Errors

```cpp
try {
    engine.runSimulation();
} catch (const std::runtime_error& e) {
    if (strstr(e.what(), "Energy divergence")) {
        // Timestep too large
        config.dt *= 0.5;
    } else if (strstr(e.what(), "NaN detected")) {
        // Numerical instability
        config.use_adaptive_timestep = true;
    }
}
```

### Validation Errors

```cpp
// Energy conservation violation
if (engine.getEnergyDrift() > 0.0001) {
    throw ValidationError("Energy drift exceeds threshold");
}

// Time reversibility failure
if (engine.getPhaseError() > 1e-4) {
    throw ValidationError("Time reversal symmetry broken");
}
```

## Examples

### Complete Example: Sine-Gordon Soliton

```cpp
#include "TRDEngine3D.h"

int main() {
    // Load configuration
    TRDEngine3D engine;
    engine.loadConfig("config/sine_gordon_soliton.yaml");

    // Set quality gates
    engine.setQualityGate("energy_conservation", 0.0001);
    engine.setQualityGate("soliton_stability", 0.99);

    // Initialize
    engine.initialize();

    // Run simulation
    engine.runSimulation();

    // Validate results
    if (!engine.passesQualityGates()) {
        std::cerr << "Validation failed\n";
        return 1;
    }

    // Export results
    engine.exportCSV("output/soliton_results.csv");

    std::cout << "Simulation complete. Energy drift: "
              << engine.getEnergyDrift() * 100 << "%\n";

    return 0;
}
```

## Thread Safety

The API is NOT thread-safe by default. For parallel simulations:

```cpp
#pragma omp parallel for
for (int i = 0; i < num_simulations; i++) {
    TRDEngine3D engine;  // Each thread gets own instance
    engine.loadConfig(configs[i]);
    engine.runSimulation();
    results[i] = engine.getObservables();
}
```

## Performance

### Benchmarks (64³ grid, 10,000 steps)

| Method | Time (s) | Energy Drift |
|--------|----------|--------------|
| CPU (single core) | 45.2 | 0.0038% |
| CPU (8 cores OpenMP) | 6.8 | 0.0038% |
| GPU (RTX 3080) | 2.1 | 0.0041% |

### Optimization Tips

1. Use GPU for grids larger than 32³
2. Enable OpenMP for multi-core CPU: `export OMP_NUM_THREADS=8`
3. Adjust timestep for stability: smaller dt = better conservation
4. Use 4th-order spatial discretization for <0.01% drift

## Version History

- **v1.0**: Initial TRDCore3D framework
- **v2.0**: Dirac equation with chiral coupling
- **v3.0**: 4th-order spatial discretization
- **v4.0**: Unified executable architecture

---

For detailed physics theory, see [PHYSICS.md](./PHYSICS.md)
For development guidelines, see [CONTRIBUTING.md](../CONTRIBUTING.md)