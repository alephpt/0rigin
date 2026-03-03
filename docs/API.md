# TRD Engine API Reference

## Overview

The TRD Engine provides a C++ API for simulating topological field theories with validated energy conservation. The architecture follows a dual-solver pattern: dissipative vacuum dynamics (TRDCore3D) and conservative particle dynamics (ConservativeSolver), orchestrated by TRDEngine3D.

All simulations are configured via YAML files and executed through the unified `./trd` executable.

---

## Class Hierarchy

```
TRDEngine3D              # GPU compute orchestration, dual-solver routing
  |-- TRDCore3D          # 3D Kuramoto phase dynamics (vacuum/dissipative)
  |-- ConservativeSolver # Conservative particle dynamics (Sine-Gordon/Dirac)
  |     |-- Dirac3D     # 3D+1 Dirac equation solver
  |-- Maxwell3D          # 3D Maxwell electromagnetic solver
  |-- TRDPipelineFactory # Vulkan compute pipeline creation
  |-- TRDBufferManager   # Vulkan buffer allocation and transfer
  |-- TRDCompute         # Vulkan compute dispatch
  |-- TRDDescriptorManager # Vulkan descriptor set management

TRDCore                  # 2D TRD core with Stueckelberg EM coupling
  |-- StuckelbergEM      # Gauge-restored massive photon field

TRDTestRunner            # Unified test execution framework
  |-- TestConfig         # YAML configuration loading
  |-- ObservableComputer # Centralized observable computation

TRD::OutputManager       # Systematic data export
```

---

## TRDEngine3D

**Header**: `src/TRDEngine3D.h`
**Source**: `src/TRDEngine3D.cpp`

Main orchestration class for 3D GPU-accelerated simulations. Implements the dual-solver architecture that routes physics evolution to the appropriate solver based on the configured physics model.

### Constructor

```cpp
TRDEngine3D(Nova* nova);
```

**Parameters**:
- `nova` -- Pointer to an initialized Nova graphics/compute engine instance. Provides Vulkan device handles for GPU resource allocation.

Creates both the vacuum solver (`TRDCore3D`) and the particle solver (`ConservativeSolver`). Initializes Vulkan pipeline factory, buffer manager, compute dispatcher, and descriptor manager.

### Destructor

```cpp
~TRDEngine3D();
```

Cleans up all Vulkan resources: buffers, pipelines, pipeline layouts, descriptor sets, descriptor layouts, and descriptor pools.

### Initialization

```cpp
void initialize(uint32_t Nx, uint32_t Ny, uint32_t Nz, float Delta);
```

**Parameters**:
- `Nx`, `Ny`, `Nz` -- Grid dimensions in each spatial direction
- `Delta` -- Mass gap parameter (vacuum potential strength)

Allocates the 3D grid for both CPU and GPU computation. Initializes TRDCore3D with the specified grid configuration, allocates GPU buffers for all field data, and creates Vulkan compute pipelines.

**Memory**: Each scalar field requires `Nx * Ny * Nz * sizeof(float)` bytes. A 64^3 grid uses approximately 1 MB per field.

### Physics Model Selection

```cpp
void setPhysicsModel(const std::string& model);
```

**Parameters**:
- `model` -- One of the following model identifiers:

| Model | Solver | Physics |
|-------|--------|---------|
| `"vacuum_kuramoto"` | TRDCore3D | Dissipative Kuramoto gradient flow |
| `"particle_sine_gordon"` | ConservativeSolver | Conservative Sine-Gordon solitons |
| `"particle_dirac"` | ConservativeSolver | Conservative Dirac fermion |
| `"coupled_vacuum_particle"` | Both | Hybrid vacuum + particle coupling |

When a particle model is selected, automatically initializes the ConservativeSolver with the current grid dimensions and default Strang splitting integration.

```cpp
void setIntegrationMethod(const std::string& method);
```

**Parameters**:
- `method` -- One of `"velocity_verlet"`, `"rk2_symplectic"`, `"strang_splitting"`, `"half_strang"`

Re-initializes the ConservativeSolver with the specified integration method.

### Simulation Execution

```cpp
void runSimulation(float dt);
```

**Parameters**:
- `dt` -- Time step size

Routes the physics evolution to the appropriate solver based on the current physics model:

- **vacuum_kuramoto**: Calls `TRDCore3D::evolveKuramotoCPU(dt)` for dissipative vacuum synchronization.
- **particle_sine_gordon**: Calls `ConservativeSolver::evolveSineGordon(dt)` and validates energy conservation at each step.
- **particle_dirac**: Calls `ConservativeSolver::evolveDirac(dt, R_field, theta_field, Delta)` with uniform vacuum fields.
- **coupled_vacuum_particle**: First evolves vacuum (Kuramoto), then extracts R-field and theta-field, then evolves the Dirac particle with chiral mass coupling from the vacuum fields.

### Field Initialization

```cpp
void setInitialPhases(const std::vector<float>& theta);
void setNaturalFrequencies(const std::vector<float>& omega);
```

Set the phase field and natural frequency distribution for Kuramoto dynamics. Vectors must have size `Nx * Ny * Nz`. Data is uploaded to GPU if buffers are allocated.

### Kuramoto Dynamics

```cpp
void stepKuramoto3D(float dt, float K, float damping);
void computeSyncField3D();
```

Execute one GPU time step of 3D Kuramoto dynamics and compute the synchronization order parameter R(x,y,z). Currently falls back to CPU via TRDCore3D when GPU dispatch is not yet configured.

### Field Access

```cpp
std::vector<float> getPhaseField3D() const;
std::vector<float> getSyncField3D() const;
std::vector<float> getMassField3D() const;
```

Return copies of the current phase field theta(x,y,z), synchronization field R(x,y,z), and mass field m(x,y,z) = Delta * R(x,y,z). Downloads from GPU if GPU buffers are active; otherwise returns CPU-side data.

### Electromagnetic Fields

```cpp
void initializeEM3D();
void stepMaxwell3D(float dt);
void initializeStuckelberg3D();
void applyGaugeTransform3D();
```

Initialize and evolve 3D electromagnetic fields. The Stueckelberg gauge transformation restores gauge invariance: A'_mu = A_mu + d_mu phi / e.

### Dirac Spinor

```cpp
void initializeDirac3D();
void stepDirac3D(float dt);
```

Initialize and evolve the 4-component 3D+1 Dirac spinor field with electromagnetic coupling: i d_t psi = (-i gamma^i d_i + m) psi with D_mu = d_mu - ieA'_mu.

### Vortex Initialization

```cpp
void initializeVortexLine(float center_x, float center_y, float center_z,
                          float radius, int axis = 2);
```

Initialize a closed vortex line (ring) in the phase field with the specified center, radius, and orientation axis (0=x, 1=y, 2=z).

### Accessors

```cpp
uint32_t getNx() const;
uint32_t getNy() const;
uint32_t getNz() const;
uint32_t getTotalPoints() const;
bool isGPUReady() const;
std::string getPhysicsModel() const;
```

---

## TRDCore3D

**Header**: `include/TRDCore3D.h`
**Source**: `src/TRDCore3D.cpp`

Foundational 3D grid infrastructure for TRD simulation. Manages the phase field theta(x,y,z), natural frequency field omega(x,y,z), and synchronization order parameter R(x,y,z) on a periodic 3D lattice. Provides both dissipative (Euler) and symplectic (RK2 midpoint) integrators for Kuramoto phase dynamics.

### Integration Modes

```cpp
enum class IntegrationMode {
    EULER,      // Fast but dissipative (legacy)
    SYMPLECTIC  // Energy-conserving RK2 midpoint (recommended)
};
```

### Configuration

```cpp
struct Config {
    uint32_t Nx = 32;
    uint32_t Ny = 32;
    uint32_t Nz = 32;
    float dx = 1.0f;
    float dt = 0.01f;
    float coupling_strength = 1.0f;
    IntegrationMode mode = IntegrationMode::SYMPLECTIC;
};
```

### Constructors

```cpp
TRDCore3D(VkDevice device, VkPhysicalDevice physicalDevice);
TRDCore3D();  // CPU-only mode
```

The CPU-only constructor sets Vulkan handles to `VK_NULL_HANDLE` and operates entirely on the host.

### Initialization

```cpp
void initialize(const Config& config);
void initializeUniform(float phase = 0.0f);
void initializeRandom(uint32_t seed = 42);
```

- `initialize` -- Allocate storage for all 3D fields according to the configuration
- `initializeUniform` -- Set all phases to a constant value, omega to zero, R to 1.0
- `initializeRandom` -- Set phases to uniform random in [-pi, pi], omega to Gaussian(0, 0.1)

### Evolution Methods

```cpp
void evolveKuramotoCPU(float dt);
void evolveSymplecticCPU(float dt);
void evolveEulerCPU(float dt);
```

- `evolveKuramotoCPU` -- Dispatches to either `evolveSymplecticCPU` or `evolveEulerCPU` based on `Config::mode`
- `evolveSymplecticCPU` -- RK2 midpoint method for the Kuramoto model. Computes k1 = f(theta), theta_mid = theta + k1 * dt/2, k2 = f(theta_mid), theta += k2 * dt. Maintains excellent time reversibility (phase error < 1e-5 rad). Note: Kuramoto is gradient flow, not Hamiltonian, so energy is not conserved by design.
- `evolveEulerCPU` -- Legacy forward Euler. Kept for backward compatibility only.

### Grid Utilities

```cpp
uint32_t index3D(uint32_t i, uint32_t j, uint32_t k) const;
void coords3D(uint32_t idx, uint32_t& i, uint32_t& j, uint32_t& k) const;
Neighbors3D getNeighbors(uint32_t i, uint32_t j, uint32_t k) const;
uint32_t wrapX(int32_t x) const;
uint32_t wrapY(int32_t y) const;
uint32_t wrapZ(int32_t z) const;
```

- `index3D` -- Convert (i,j,k) coordinates to row-major linear index: `k * (Nx * Ny) + j * Nx + i`
- `coords3D` -- Convert linear index back to (i,j,k) coordinates
- `getNeighbors` -- Return the 6-neighbor stencil indices with periodic boundary wrapping
- `wrapX/Y/Z` -- Apply periodic boundary conditions to a single coordinate

### Observables

```cpp
float computeEnergy() const;
void computeRField();
float getAverageR() const;
```

- `computeEnergy` -- Total system energy: E = sum[ 0.5 * omega_i^2 + 0.5 * K * coupling_i^2 ]
- `computeRField` -- Compute local synchronization order parameter: R = |<e^(i theta)>| using nearest-neighbor averaging
- `getAverageR` -- Return the spatially averaged R value over the entire grid

### Field Access

```cpp
std::vector<float>& getTheta();
const std::vector<float>& getTheta() const;
std::vector<float>& getOmega();
const std::vector<float>& getOmega() const;
std::vector<float>& getRField();
const std::vector<float>& getRField() const;

void setPhaseField(const float* data);
void getPhaseField(float* data) const;
void getRField(float* data) const;
```

Direct mutable and const access to internal field arrays for GPU upload/download or external manipulation.

### GPU Control

```cpp
void enableGPU();
bool isGPUEnabled() const;
```

---

## ConservativeSolver

**Header**: `include/ConservativeSolver.h`
**Source**: `src/ConservativeSolver.cpp`

Conservative/unitary particle dynamics solver. Implements symplectic integration for wave equations (Sine-Gordon) and spinor dynamics (Dirac). Provides validated energy conservation below 0.01%.

### Integration Methods

```cpp
enum class IntegrationMethod {
    VELOCITY_VERLET,   // Kick-drift-kick for wave equations
    RK2_SYMPLECTIC,    // Conservative field theories
    STRANG_SPLITTING,  // T-V-T splitting for nonlinear PDEs
    HALF_STRANG        // Phase-magnitude decoupling
};
```

### Spatial Discretization

```cpp
enum class SpatialOrder {
    SECOND_ORDER,   // 6-neighbor stencil: O(dx^2)
    FOURTH_ORDER    // 12-neighbor stencil: O(dx^4)
};
```

Fourth-order is the default. It uses a 12-neighbor stencil per dimension (two neighbors in each direction) and achieves 16x better accuracy than the 2nd-order stencil.

### Configuration

```cpp
struct Config {
    uint32_t nx = 64, ny = 64, nz = 64;
    float dx = 1.0f;
    float dt = 0.005f;
    IntegrationMethod method = IntegrationMethod::VELOCITY_VERLET;
    SpatialOrder spatial_order = SpatialOrder::FOURTH_ORDER;
    float sine_gordon_mass = 1.0f;
    float coupling_K = 1.0f;
};
```

### Initialization

```cpp
ConservativeSolver();
~ConservativeSolver();
void initialize(const Config& config);
```

### Field Evolution

```cpp
void evolveSineGordon(float dt);
```

Evolve the Sine-Gordon equation: d^2 theta / dt^2 = nabla^2 theta - sin(theta).

Uses Velocity Verlet (kick-drift-kick):
1. Half-step velocity: `v += 0.5 * dt * F(x)`
2. Full-step position: `x += dt * v`
3. Recompute force: `F(x_new)`
4. Half-step velocity: `v += 0.5 * dt * F(x_new)`

When configured for Strang splitting, uses the T-V-T pattern instead:
1. Half-step kinetic: `theta += (dt/2) * pi`
2. Full-step potential: `pi += dt * (nabla^2 theta - sin(theta))`
3. Half-step kinetic: `theta += (dt/2) * pi`

```cpp
void evolveDirac(float dt, const std::vector<float>& R_field,
                 const std::vector<float>& theta_field, float Delta);
```

Evolve the Dirac equation with chiral mass coupling from vacuum fields. Uses the internal Dirac3D solver with split-step integration: kinetic half-step, chiral mass step, kinetic half-step.

**Parameters**:
- `dt` -- Time step
- `R_field` -- Vacuum R(x) synchronization field from TRDCore3D
- `theta_field` -- Vacuum theta(x) phase field from TRDCore3D
- `Delta` -- Coupling strength

### Initial Conditions

```cpp
void initializeVortexWithProperVelocity(
    float x0, float y0, float z0, int charge);
```

Initialize a topological vortex with proper velocity field. The phase field is set to theta(r, phi) = n * phi (topological winding number `charge`), and the velocity field is set to d_theta/dt = 0 for a stationary vortex. This initialization was critical for resolving a 164% to <0.01% energy drift improvement.

```cpp
void initializeCollisionScenario(
    float x1, float y1, float x2, float y2, float velocity);
```

Initialize a two-vortex collision scenario with proper velocity fields for moving vortices using linear superposition.

```cpp
void initializeGaussian(
    float x0, float y0, float z0, float sigma, float amplitude);
```

Initialize a Gaussian wave packet. Validated: 0.127% energy drift over 1000 steps.

### Energy Computation

```cpp
float computeTotalEnergy() const;
```

Compute the total energy functional:

```
E = integral[ (d_theta/dt)^2 + (nabla theta)^2 + (1 - cos(theta)) ] dV
```

Components: kinetic energy (d_theta/dt)^2, gradient energy (nabla theta)^2, potential energy 1 - cos(theta) (Sine-Gordon).

```cpp
float measureEnergyDrift(float E_initial) const;
```

Return the relative energy drift: |E_current - E_initial| / E_initial.

### Validation

```cpp
bool validateEnergyConservation(float threshold = 0.0001f);
bool validateTimeReversibility(float threshold = 1e-4f);
```

- `validateEnergyConservation` -- Returns true if the measured drift is below the threshold (default 0.01%)
- `validateTimeReversibility` -- Runs forward evolution then backward evolution and measures the phase error. Returns true if error is below the threshold.

### Field Access

```cpp
const std::vector<float>& getTheta() const;
const std::vector<float>& getThetaDot() const;
const std::vector<float>& getPsi() const;  // Real part of Dirac spinor

uint32_t getNx() const;
uint32_t getNy() const;
uint32_t getNz() const;
uint32_t getTotalPoints() const;

void setMass(float mass);
```

---

## Dirac3D

**Header**: `include/Dirac3D.h`
**Source**: `src/Dirac3D.cpp`

3D+1 Dirac equation solver implementing 4-component spinor evolution using the split-operator method with FFT.

**Equation**: i hbar d_psi/dt = H psi, where H = -i hbar c alpha . nabla + beta m c^2

The alpha matrices (alpha_x, alpha_y, alpha_z) and beta matrix are stored as static 4x4 complex arrays in the Dirac representation. The gamma5 matrix is used for chiral mass coupling.

### Constructor

```cpp
Dirac3D(uint32_t Nx, uint32_t Ny, uint32_t Nz);
~Dirac3D();
```

Allocates 4-component spinor fields in both position and momentum space, sets up FFTW3 plans for forward and inverse transforms, and computes the momentum grid.

### Initialization

```cpp
void initialize(const std::vector<std::complex<float>>& psi_init);
void initializeGaussian(float x0, float y0, float z0, float sigma);
```

- `initialize` -- Set the spinor field from a flat vector of 4 * N_total complex values
- `initializeGaussian` -- Create a Gaussian wavepacket centered at (x0, y0, z0) with width sigma. Populates all four spinor components.

### Evolution

```cpp
void step(const std::vector<float>& mass_field, float dt);
```

Split-operator evolution step with scalar mass:
1. Kinetic half-step: apply exp(-i alpha . p dt/2) in momentum space via FFT
2. Mass step: apply exp(-i beta m dt) in position space
3. Kinetic half-step: apply exp(-i alpha . p dt/2) in momentum space via FFT

```cpp
void stepWithChiralMass(const std::vector<float>& R_field,
                        const std::vector<float>& theta_field,
                        float Delta, float dt);
```

Split-operator step with the full chiral mass operator M = Delta R e^(i theta gamma5):
1. Kinetic half-step (FFT-based, momentum space)
2. Chiral mass step using Velocity Verlet with eigenvalue decomposition for numerical stability
3. Kinetic half-step (FFT-based, momentum space)

The eigenvalue decomposition of the mass matrix ensures stable time evolution even for large coupling strengths. The eigenvalues are:
- lambda_+ = Delta R (1 + cos(theta)) [positive chirality]
- lambda_- = Delta R (1 - cos(theta)) [negative chirality]

### Observables

```cpp
std::vector<float> getDensity() const;
std::vector<float> getCurrent(int component) const;
const std::vector<std::complex<float>>& getComponent(int component) const;
float getNorm() const;
```

- `getDensity` -- Spinor probability density rho = sum_alpha |psi_alpha|^2 at each grid point
- `getCurrent` -- Probability current j^i = psi_dagger alpha^i psi for component 0=x, 1=y, 2=z
- `getComponent` -- Direct access to one of the 4 spinor components (0-3)
- `getNorm` -- Total norm integral(psi_dagger psi d^3x), which should be conserved

### Static Data Members

```cpp
static const std::array<std::complex<float>, 16> alpha_x;  // Dirac alpha_x (4x4)
static const std::array<std::complex<float>, 16> alpha_y;  // Dirac alpha_y (4x4)
static const std::array<std::complex<float>, 16> alpha_z;  // Dirac alpha_z (4x4)
static const std::array<std::complex<float>, 16> beta;     // Dirac beta (4x4)
static const std::array<std::complex<float>, 16> gamma5;   // Chirality operator (4x4)
```

All stored in row-major order: M[i,j] = data[i*4 + j].

---

## Maxwell3D

**Header**: `include/Maxwell3D.h`
**Source**: `src/Maxwell3D.cpp`

3D Maxwell electromagnetic field solver implementing the full set of vacuum Maxwell equations with periodic boundary conditions.

**Equations**:
- Faraday: dE/dt = curl(B)
- Ampere: dB/dt = -curl(E)

### Constructor

```cpp
Maxwell3D(uint32_t Nx, uint32_t Ny, uint32_t Nz);
```

Allocates storage for all 6 field components: Ex, Ey, Ez, Bx, By, Bz.

### Initialization

```cpp
void initialize(const std::vector<float>& Ex_init,
                const std::vector<float>& Ey_init,
                const std::vector<float>& Ez_init,
                const std::vector<float>& Bx_init,
                const std::vector<float>& By_init,
                const std::vector<float>& Bz_init);
void initializeSphericalWave(float wavelength, float amplitude);
```

- `initialize` -- Set all 6 field components from external data
- `initializeSphericalWave` -- Create a spherical electromagnetic wave for testing

### Evolution

```cpp
void step(float dt);
void evolveElectricField(float dt);
void evolveMagneticField(float dt);
```

- `step` -- Full Maxwell step using second-order accurate Strang splitting
- `evolveElectricField` -- Advance E-field: dE/dt = curl(B)
- `evolveMagneticField` -- Advance B-field: dB/dt = -curl(E)

### Curl Operator

```cpp
std::vector<float> curl(const std::vector<float>& Fx,
                        const std::vector<float>& Fy,
                        const std::vector<float>& Fz,
                        int component) const;
```

Compute a single component (0=x, 1=y, 2=z) of the curl of a vector field using central differences with periodic boundary conditions.

### Field Access

```cpp
const std::vector<float>& getEx() const;
const std::vector<float>& getEy() const;
const std::vector<float>& getEz() const;
const std::vector<float>& getBx() const;
const std::vector<float>& getBy() const;
const std::vector<float>& getBz() const;
```

### Energy

```cpp
std::vector<float> getEnergyDensity() const;
float getTotalEnergy() const;
```

- `getEnergyDensity` -- Local energy density u = (E^2 + B^2) / 2 at each grid point
- `getTotalEnergy` -- Integrated total electromagnetic energy over the grid

---

## TRDCore

**Header**: `include/TRDCore.h`
**Source**: `src/TRDCore.cpp`

Minimal 2D demonstration core for TRD + Stueckelberg EM coupled evolution. Provides the integration pattern for coupling Stueckelberg EM fields to TRD phase dynamics.

### Configuration

```cpp
struct Config {
    uint32_t nx = 64;
    uint32_t ny = 64;
    float dx = 1.0f;
    float dt = 0.01f;
};
```

### Methods

```cpp
TRDCore(VkDevice device, VkPhysicalDevice physicalDevice);
void initialize(const Config& config);
void enableEM(float photon_mass_coupling);
void evolveFields(float dt);
void evolveEM(float dt);
physics::FieldTensor getEMFieldAt(int i, int j) const;
const std::vector<float>& getPhaseField() const;
const std::vector<float>& getSyncField() const;
```

- `enableEM` -- Initialize the Stueckelberg EM system with the specified photon mass coupling
- `evolveFields` -- Evolve TRD phase fields one time step
- `evolveEM` -- Evolve EM fields using Stueckelberg dynamics
- `getEMFieldAt` -- Return the electromagnetic field tensor at grid point (i,j)

---

## StuckelbergEM

**Header**: `include/physics/StuckelbergEM.h`
**Source**: `src/physics/StuckelbergEM.cpp`
**Namespace**: `physics`

Stueckelberg gauge-restored electromagnetic field implementation. Unlike the Proca approach, this formulation maintains gauge invariance through a scalar compensator field.

**Gauge transformation**: Under theta -> theta + e alpha:
- A_mu -> A_mu + d_mu alpha
- phi -> phi - e alpha
- A'_mu = A_mu + d_mu phi / e (gauge invariant)

**Evolution equations**:
- Box A_mu = 0 (massless Maxwell for the gauge potential)
- Box phi + m^2 phi = 0 (Klein-Gordon for the Stueckelberg scalar)
- B = curl(A + d phi) (field from both A and phi)

### Constructor

```cpp
StuckelbergEM(int nx, int ny, float dx, float photon_mass);
```

### Core Methods (GaugeTheory Interface)

```cpp
void computePotentials(const float* theta_field, const float* R_field,
                       int nx, int ny, float dx, float dt) override;
void computeFieldStrengths() override;
FieldTensor getFieldAt(int i, int j) const override;
bool isGaugeInvariant() const override;  // Returns true
float computeFieldEnergy() const override;
```

### Stueckelberg-Specific Methods

```cpp
void evolveStuckelbergField(float dt);
void evolveMaxwell(float dt);
void initializeGaussianPulse(float center_x, float center_y,
                             float width_x, float width_y,
                             float amplitude, float k_x, float k_y,
                             const std::string& component);
float getPhiAt(int i, int j) const;
float getAprimeX(int i, int j) const;
float getAprimeY(int i, int j) const;
```

- `evolveStuckelbergField` -- Evolve the scalar field phi using Klein-Gordon dynamics
- `evolveMaxwell` -- Evolve the gauge potential A_mu using Maxwell dynamics
- `initializeGaussianPulse` -- Set up a Gaussian EM pulse for wave propagation testing

---

## TRDTestRunner

**Header**: `src/simulations/TRDTestRunner.h`
**Source**: `src/simulations/TRDTestRunner.cpp`

Unified test execution framework. Loads YAML configurations, initializes the engine, runs simulations with multiple substep ratios, validates results, and generates reports.

### Constructors

```cpp
TRDTestRunner(const std::string& config_path);
TRDTestRunner(const TestConfig& config);
```

### Lifecycle

```cpp
bool initialize();
bool run();
void generateReport(const std::string& output_path = "") const;
```

- `initialize` -- Load configuration, create Nova instance, initialize TRD engine, create output directories
- `run` -- Execute all tests. If operator splitting is enabled, runs once per substep ratio (e.g., N=1, 10, 100) and validates convergence between them.
- `generateReport` -- Write a detailed text report with pass/fail status, quantitative metrics, and EM field observables

### Validation

```cpp
bool allTestsPassed() const;
const std::vector<ObservableComputer::Observables>& getResults(int N) const;
```

Internal validation covers:
- **Norm conservation**: |psi|^2 should remain 1.0 within the configured tolerance
- **Energy conservation**: |Delta E / E0| below the configured threshold
- **Convergence**: When testing multiple substep ratios, observables at different N should converge

---

## TestConfig

**Header**: `src/simulations/TestConfig.h`
**Source**: `src/simulations/TestConfig.cpp`

YAML-driven test configuration structure with nested sub-configurations for grid, physics, initial conditions, operator splitting, validation, and output.

### Sub-Configurations

| Struct | Key Fields |
|--------|------------|
| `GridConfig` | `size_x`, `size_y` |
| `PhysicsConfig` | `delta`, `coupling`, `dt`, `total_steps`, `K`, `damping`, `em_coupling_enabled`, `photon_mass` |
| `DiracInitialCondition` | `type`, `x0`, `y0`, `sigma`, `amplitude` |
| `KuramotoInitialCondition` | `phase_distribution`, `omega_distribution`, `omega_mean`, `omega_std`, `wave_vector_x`, `wave_vector_y` |
| `EMPulseConfig` | `type`, `component`, `center_x`, `center_y`, `width_x`, `width_y`, `amplitude`, `wave_vector_x`, `wave_vector_y` |
| `OperatorSplittingConfig` | `enabled`, `substep_ratios` |
| `ValidationConfig` | `norm_tolerance`, `energy_tolerance`, `convergence_tolerance`, `expected_wave_speed` |
| `OutputConfig` | `directory`, `save_every`, `formats`, `auto_visualize`, `track_wave_position` |

### Methods

```cpp
bool loadFromYAML(const std::string& yaml_path);
bool saveToYAML(const std::string& yaml_path) const;
static TestConfig createDefault(const std::string& test_name = "default_test");
bool validate() const;
void print() const;
```

---

## ObservableComputer

**Header**: `src/simulations/ObservableComputer.h`
**Source**: `src/simulations/ObservableComputer.cpp`

Centralized computation of all physical observables. All methods are static.

### Observables Structure

```cpp
struct Observables {
    double time;
    double norm, norm_error;
    double energy_total, energy_kinetic, energy_potential;
    std::complex<double> position_x, position_y;
    std::complex<double> momentum_x, momentum_y;
    double R_avg, R_max, R_min, R_variance;
    double EM_B_max, EM_B_rms, EM_energy;
    bool norm_valid, energy_valid;
};
```

### Methods

```cpp
static void compute(Observables* result,
                    const DiracEvolution& dirac,
                    const std::vector<double>& R_field,
                    double delta, double time,
                    double E0 = 0.0,
                    double norm_tolerance = 1e-4,
                    double energy_tolerance = 1e-2,
                    const TRDEngine* engine = nullptr);
```

Compute all observables for the current simulation state and populate the output struct. Uses an output parameter rather than return-by-value to avoid a compiler optimization issue where struct fields were corrupted during return.

```cpp
static double computeNorm(const DiracEvolution& dirac);
static double computeEnergy(const DiracEvolution& dirac,
                            const std::vector<double>& R_field, double delta);
static double computeKineticEnergy(const DiracEvolution& dirac);
static double computePotentialEnergy(const DiracEvolution& dirac,
                                     const std::vector<double>& R_field, double delta);
static std::complex<double> computePositionExpectation(
    const DiracEvolution& dirac, int component);
static std::complex<double> computeMomentumExpectation(
    const DiracEvolution& dirac, int component);
static std::tuple<double, double, double, double>
    computeSyncFieldStats(const std::vector<double>& R_field);
static std::string toCSVLine(const Observables& obs);
static std::string getCSVHeader();
```

---

## TRD::OutputManager

**Header**: `src/output/OutputManager.h`
**Source**: `src/output/OutputManager.cpp`
**Namespace**: `TRD`

Systematic output generation for experiments. Creates timestamped directories, writes metadata-annotated data files, and manages CSV export.

### Constructor

```cpp
explicit OutputManager(
    const std::string& base_path = "/home/persist/neotec/0rigin/output");
```

### Directory Management

```cpp
std::string createExperimentDirectory(const std::string& experiment_name);
std::string createSubdirectory(const std::string& parent_dir,
                                const std::string& subdir_name);
static std::string getTimestamp();
```

Creates directories of the form: `output/YYYYMMDD_HHMMSS_<experiment_name>/`

### Data Writing

```cpp
void writeTimeseries(const std::string& filepath,
                     const std::vector<float>& time_data,
                     const std::vector<float>& observable_data,
                     const std::map<std::string, std::string>& metadata,
                     const std::string& observable_name = "observable");

void writeSnapshot(const std::string& filepath,
                   const std::vector<float>& field_data,
                   int nx, int ny,
                   const std::map<std::string, std::string>& metadata);

void writeMultipleTimeseries(const std::string& filepath,
                              const std::vector<float>& time_data,
                              const std::map<std::string, std::vector<float>>& observables,
                              const std::map<std::string, std::string>& metadata);

void writeMetadata(const std::string& filepath,
                   const std::map<std::string, std::string>& parameters);

void writeCSV(const std::string& filepath,
              const std::vector<std::string>& headers,
              const std::vector<std::vector<float>>& data);

void appendCSV(const std::string& filepath,
               const std::vector<float>& row,
               bool create_if_missing = true,
               const std::vector<std::string>& headers = {});
```

All `.dat` files include metadata headers (lines prefixed with `#`) containing simulation parameters, making them self-documenting and compatible with Python visualization scripts.

---

## Vulkan Infrastructure

These classes manage the GPU compute backend and are internal to TRDEngine3D.

### TRDPipelineFactory

**Header/Source**: `src/TRDPipelineFactory.h/.cpp`

Creates Vulkan compute pipelines from SPIR-V shader modules. Handles pipeline layout creation with push constant ranges and descriptor set layouts.

### TRDBufferManager

**Header/Source**: `src/TRDBufferManager.h/.cpp`

Manages Vulkan buffer allocation, memory binding, and host-device data transfers. Provides `createBuffer()`, `uploadData()`, and `downloadData()` for field data movement between CPU and GPU.

### TRDCompute

**Header/Source**: `src/TRDCompute.h/.cpp`

Dispatches Vulkan compute shaders. Manages the compute command pool and command buffers, submits workgroups, and synchronizes execution.

### TRDDescriptorManager

**Header/Source**: `src/TRDDescriptorManager.h/.cpp`

Manages Vulkan descriptor pools, descriptor set layouts, and descriptor set allocation/updates for binding storage buffers to compute shaders.

---

## GPU Compute Shaders

Located in `shaders/smft/`:

| Shader | Purpose |
|--------|---------|
| `kuramoto3d.comp` | 3D Kuramoto phase evolution |
| `kuramoto_step.comp` | Single-step Kuramoto dynamics |
| `kuramoto_stochastic.comp` | Stochastic Kuramoto model |
| `r_field_evolution.comp` | R-field (synchronization) evolution |
| `dirac_velocity_verlet.comp` | Dirac spinor dynamics (Velocity Verlet) |
| `gravity_field.comp` | Gravitational field computation |
| `em_stress_energy.comp` | EM stress-energy tensor |
| `spinor_feedback.comp` | Spinor field feedback |
| `accumulate.comp` | Field accumulation operations |

Shaders are compiled from GLSL to SPIR-V by CMake using `glslc`.

---

## Performance Benchmarks

### CPU Performance (64^3 grid, 10,000 steps)

| Configuration | Wall Time | Energy Drift |
|---------------|-----------|--------------|
| Single core | ~45 s | 0.0038% |
| 8 cores (OpenMP) | ~7 s | 0.0038% |

### Optimization Guidance

- Use GPU for grids larger than 32^3
- Set `OMP_NUM_THREADS` for CPU parallelism
- Smaller timesteps improve energy conservation
- 4th-order spatial discretization provides 16x better accuracy than 2nd-order at minimal cost

---

## Thread Safety

TRDEngine3D instances are not thread-safe. For parallel parameter sweeps, create independent instances per thread:

```cpp
#pragma omp parallel for
for (int i = 0; i < num_simulations; i++) {
    ConservativeSolver solver;
    solver.initialize(configs[i]);
    // ... run independently
}
```

---

For physics theory, see [PHYSICS.md](./PHYSICS.md).
For development workflow, see [DEVELOPER_GUIDE.md](./DEVELOPER_GUIDE.md).
For contribution guidelines, see [CONTRIBUTING.md](../CONTRIBUTING.md).
