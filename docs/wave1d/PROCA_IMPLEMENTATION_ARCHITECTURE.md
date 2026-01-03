# Proca EM Implementation Architecture
## C++ Class Design for Gauge-Covariant Electromagnetic Theory

**Author**: Claude Code (Operations Tier 1)
**Date**: 2025-12-31
**Status**: DESIGN PHASE
**Dependency**: GAUGE_COVARIANT_EM_THEORY.md (mathematical foundation)

---

## 1. Architecture Overview

### 1.1 Design Principles

**Modularity**: EM theory as plugin to TRD, not core dependency
**Abstraction**: GaugeTheory interface supports multiple implementations (Proca, Stückelberg, Maxwell)
**GPU-first**: All field evolution on GPU via Vulkan compute shaders
**Type safety**: Strong typing for physical quantities (Vec3, Vec4, FieldTensor)
**Zero-cost abstractions**: Interface overhead eliminated at compile time

### 1.2 Component Diagram

```
┌─────────────────────────────────────────────────────────────┐
│                        TRDCore                              │
│  - Kuramoto evolution (θ, R fields)                         │
│  - Manages simulation state                                  │
│  - Optional EM coupling                                      │
└────────────────────┬────────────────────────────────────────┘
                     │ has-a (optional)
                     ▼
┌─────────────────────────────────────────────────────────────┐
│            GaugeTheory (abstract interface)                  │
│  + computePotentials(FieldState&) = 0                       │
│  + computeFieldStrengths(FieldState&) = 0                   │
│  + evolveFields(dt) = 0                                      │
│  + isGaugeInvariant() const = 0                             │
│  + getCouplingStrength() const = 0                          │
└────────────────────┬────────────────────────────────────────┘
                     │ implements
                     ▼
┌─────────────────────────────────────────────────────────────┐
│                      ProcaEM                                 │
│  - m_photon_mass_coupling: float (g)                        │
│  - m_current_coupling: float (α)                            │
│  - m_gauge_parameter: float (ξ for gauge fixing)           │
│  - m_A_field: Vec4Field (A_μ at t, t-dt, t-2dt)            │
│  - m_j_current: Vec4Field (j_μ source)                      │
│  - m_EM_buffers: VulkanBufferSet                            │
│  + computePotentials() override                             │
│  + computeFieldStrengths() override                         │
│  + evolveFields(dt) override                                │
│  + isGaugeInvariant() const override { return false; }     │
│  + quantifyGaugeViolation() const                           │
└─────────────────────────────────────────────────────────────┘
                     │ uses
                     ▼
┌─────────────────────────────────────────────────────────────┐
│                  GPU Compute Shaders                         │
│  - computeEMCurrent.comp (j_μ from θ, R)                    │
│  - evolveProcaField.comp (A_μ evolution)                    │
│  - computeFieldStrengths.comp (E, B from A_μ)               │
│  - enforceGaugeCondition.comp (Lorenz gauge ∂_μA^μ=0)      │
└─────────────────────────────────────────────────────────────┘
```

---

## 2. Core Class Definitions

### 2.1 GaugeTheory (Abstract Interface)

**File**: `include/physics/GaugeTheory.h`

```cpp
#pragma once

#include "FieldState.h"
#include "PhysicsTypes.h"
#include <memory>

namespace physics {

/**
 * @brief Abstract interface for electromagnetic gauge theories
 *
 * Supports multiple implementations:
 * - Maxwell (massless photon, U(1) gauge invariant)
 * - Proca (massive photon, breaks U(1))
 * - Stückelberg (massive photon, restores U(1))
 */
class GaugeTheory {
public:
    virtual ~GaugeTheory() = default;

    /**
     * @brief Compute gauge potentials A_μ from field state
     * @param state Current field configuration (θ, R, etc.)
     */
    virtual void computePotentials(const FieldState& state) = 0;

    /**
     * @brief Compute field strengths E, B from potentials A_μ
     * @param state Field state (updated with E, B)
     */
    virtual void computeFieldStrengths(FieldState& state) = 0;

    /**
     * @brief Evolve EM fields forward by timestep dt
     * @param dt Timestep
     */
    virtual void evolveFields(float dt) = 0;

    /**
     * @brief Check if theory is U(1) gauge invariant
     * @return true if A_μ → A_μ + ∂_μα leaves physics unchanged
     */
    virtual bool isGaugeInvariant() const = 0;

    /**
     * @brief Get EM-TRD coupling strength
     * @return Dimensionful coupling constant α
     */
    virtual float getCouplingStrength() const = 0;

    /**
     * @brief Get photon mass (0 for Maxwell, m_γ for Proca)
     * @param R Local order parameter value
     * @return Effective photon mass
     */
    virtual float getPhotonMass(float R) const = 0;
};

} // namespace physics
```

### 2.2 ProcaEM (Concrete Implementation)

**File**: `include/physics/ProcaEM.h`

```cpp
#pragma once

#include "GaugeTheory.h"
#include "vulkan/VulkanBufferSet.h"
#include "vulkan/VulkanComputePipeline.h"
#include <array>

namespace physics {

/**
 * @brief Proca electromagnetic theory with TRD coupling
 *
 * Lagrangian:
 *   ℒ = -1/4 F_μν F^μν + 1/2 m_γ²(R) A_μ A^μ + j_μ(θ,R) A^μ
 *
 * Key features:
 * - Massive photon: m_γ(R) = g·(1-R)
 * - Breaks U(1) gauge invariance
 * - Allows B ≠ 0 from scalar source
 * - Lorenz gauge: ∂_μ A^μ = 0
 */
class ProcaEM : public GaugeTheory {
public:
    /**
     * @brief Construct Proca EM module
     * @param device Vulkan device
     * @param gridSize (Nx, Ny, Nz)
     * @param spatialStep (dx, dy, dz)
     * @param photonMassCoupling g (mass scale)
     * @param currentCoupling α (EM-TRD coupling)
     */
    ProcaEM(
        VkDevice device,
        const std::array<uint32_t, 3>& gridSize,
        const std::array<float, 3>& spatialStep,
        float photonMassCoupling,
        float currentCoupling
    );

    ~ProcaEM() override;

    // GaugeTheory interface implementation
    void computePotentials(const FieldState& state) override;
    void computeFieldStrengths(FieldState& state) override;
    void evolveFields(float dt) override;
    bool isGaugeInvariant() const override { return false; }
    float getCouplingStrength() const override { return m_current_coupling; }
    float getPhotonMass(float R) const override;

    /**
     * @brief Quantify gauge violation magnitude
     * @return δB/B when θ → θ + α (dimensionless)
     */
    float quantifyGaugeViolation() const;

    /**
     * @brief Compute magnetic flux through domain
     * @return Φ_B = ∫ B·dA
     */
    float computeMagneticFlux() const;

    /**
     * @brief Check Lorenz gauge condition ∂_μ A^μ = 0
     * @return RMS violation |∂_μ A^μ|
     */
    float checkLorenzGauge() const;

private:
    // Physical parameters
    float m_photon_mass_coupling;  // g
    float m_current_coupling;      // α
    float m_gauge_parameter;       // ξ for R_ξ gauge fixing

    // Grid parameters
    std::array<uint32_t, 3> m_grid_size;  // (Nx, Ny, Nz)
    std::array<float, 3> m_spatial_step;  // (dx, dy, dz)
    float m_dt;                            // Current timestep

    // Field storage (GPU buffers)
    std::unique_ptr<vulkan::VulkanBufferSet> m_buffers;

    // A_μ field at multiple time levels (for leap-frog)
    VkBuffer m_A_field_current;   // A^μ(t)
    VkBuffer m_A_field_previous;  // A^μ(t-dt)
    VkBuffer m_A_field_next;      // A^μ(t+dt)

    // Current j_μ and field strengths E, B
    VkBuffer m_j_current;         // j^μ
    VkBuffer m_E_field;           // E^i
    VkBuffer m_B_field;           // B^i

    // Compute pipelines
    std::unique_ptr<vulkan::VulkanComputePipeline> m_compute_current_pipeline;
    std::unique_ptr<vulkan::VulkanComputePipeline> m_evolve_field_pipeline;
    std::unique_ptr<vulkan::VulkanComputePipeline> m_compute_strengths_pipeline;
    std::unique_ptr<vulkan::VulkanComputePipeline> m_enforce_gauge_pipeline;

    // Internal methods
    void initializeBuffers();
    void createComputePipelines();
    void computeEMCurrent(const FieldState& state);
    void enforceLorenzGauge();
};

} // namespace physics
```

### 2.3 PhysicsTypes (Type Definitions)

**File**: `include/physics/PhysicsTypes.h`

```cpp
#pragma once

#include <array>
#include <cstdint>

namespace physics {

/**
 * @brief 3-vector (spatial)
 */
struct Vec3 {
    float x, y, z;

    Vec3() : x(0), y(0), z(0) {}
    Vec3(float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}

    float norm() const;
    float normSquared() const;
    Vec3 normalized() const;
};

/**
 * @brief 4-vector (spacetime)
 */
struct Vec4 {
    float t, x, y, z;

    Vec4() : t(0), x(0), y(0), z(0) {}
    Vec4(float t_, float x_, float y_, float z_) : t(t_), x(x_), y(y_), z(z_) {}

    // Minkowski inner product: η_μν = diag(+1,-1,-1,-1)
    float minkowskiDot(const Vec4& other) const;
};

/**
 * @brief Antisymmetric field strength tensor F_μν
 */
struct FieldTensor {
    // F_μν = -F_νμ, so store 6 independent components
    // F_01 = E_x, F_02 = E_y, F_03 = E_z (electric field)
    // F_12 = B_z, F_23 = B_x, F_31 = B_y (magnetic field)

    float E_x, E_y, E_z;
    float B_x, B_y, B_z;

    FieldTensor() : E_x(0), E_y(0), E_z(0), B_x(0), B_y(0), B_z(0) {}

    // Compute invariants
    float electricInvariant() const;  // E·E
    float magneticInvariant() const;  // B·B
    float pseudoscalarInvariant() const;  // E·B
};

/**
 * @brief Field state containing TRD + EM fields
 */
struct FieldState {
    // TRD fields (host-side metadata, GPU buffers)
    VkBuffer theta_field;   // θ(x,y,z)
    VkBuffer R_field;       // R(x,y,z)

    // EM fields (managed by GaugeTheory implementation)
    VkBuffer A_mu;          // A_μ(x,y,z)
    VkBuffer E_field;       // E(x,y,z)
    VkBuffer B_field;       // B(x,y,z)

    // Grid metadata
    uint32_t Nx, Ny, Nz;
    float dx, dy, dz;
    float dt;
};

} // namespace physics
```

---

## 3. GPU Shader Architecture

### 3.1 Shader: computeEMCurrent.comp

**Purpose**: Compute j_μ = α ∂_μθ f(R) from TRD fields

**File**: `build/shaders/proca/computeEMCurrent.comp`

```glsl
#version 450
#extension GL_GOOGLE_include_directive : enable

#include "../include/precision.glsl"

/**
 * Compute electromagnetic current j_μ from TRD fields
 *
 * Prescription:
 *   j_μ = α · ∂_μθ · f(R)
 *
 * where:
 *   α = current coupling (push constant)
 *   f(R) = R² (form factor, only synchronized regions source)
 *   ∂_μθ = (∂_tθ, ∂_xθ, ∂_yθ, ∂_zθ)
 */

layout(local_size_x = 8, local_size_y = 8, local_size_z = 8) in;

layout(set = 0, binding = 0) readonly buffer ThetaCurrent {
    float theta_current[];
};

layout(set = 0, binding = 1) readonly buffer ThetaPrevious {
    float theta_previous[];
};

layout(set = 0, binding = 2) readonly buffer RField {
    float R[];
};

layout(set = 0, binding = 3) writeonly buffer JCurrent {
    vec4 j_mu[];  // (j^0, j^1, j^2, j^3)
};

layout(push_constant) uniform Params {
    uint Nx, Ny, Nz;
    float dx, dy, dz, dt;
    float alpha;  // Current coupling
} params;

const float PI = 3.14159265359;

void main() {
    uint i = gl_GlobalInvocationID.x;
    uint j = gl_GlobalInvocationID.y;
    uint k = gl_GlobalInvocationID.z;

    if (i >= params.Nx || j >= params.Ny || k >= params.Nz) return;

    uint idx = k * (params.Nx * params.Ny) + j * params.Nx + i;

    // Get R and form factor f(R) = R²
    float R_val = R[idx];
    float f_R = R_val * R_val;

    // Compute ∂_tθ (temporal derivative)
    float theta_diff = theta_current[idx] - theta_previous[idx];

    // Wrap to [-π, π]
    while (theta_diff > PI) theta_diff -= 2.0 * PI;
    while (theta_diff < -PI) theta_diff += 2.0 * PI;

    float d_t_theta = theta_diff / params.dt;

    // Compute ∂_xθ (spatial derivative, centered difference)
    uint i_left = (i == 0) ? (params.Nx - 1) : (i - 1);
    uint i_right = (i == params.Nx - 1) ? 0 : (i + 1);
    uint idx_left = k * (params.Nx * params.Ny) + j * params.Nx + i_left;
    uint idx_right = k * (params.Nx * params.Ny) + j * params.Nx + i_right;

    float theta_left = theta_current[idx_left];
    float theta_right = theta_current[idx_right];
    float theta_center = theta_current[idx];

    // Phase-aware gradient (handle 2π wraps)
    float dtheta_x = theta_right - theta_left;
    while (dtheta_x > PI) dtheta_x -= 2.0 * PI;
    while (dtheta_x < -PI) dtheta_x += 2.0 * PI;
    float d_x_theta = dtheta_x / (2.0 * params.dx);

    // Similar for ∂_yθ, ∂_zθ (code omitted for brevity)
    // ...

    // Compute j_μ = α · ∂_μθ · f(R)
    vec4 j_mu_val;
    j_mu_val.x = params.alpha * d_t_theta * f_R;  // j^0
    j_mu_val.y = params.alpha * d_x_theta * f_R;  // j^1
    // j_mu_val.z = ... (j^2)
    // j_mu_val.w = ... (j^3)

    j_mu[idx] = j_mu_val;
}
```

### 3.2 Shader: evolveProcaField.comp

**Purpose**: Evolve A_μ via Proca equation (□ + m_γ²) A^μ = j^μ

**File**: `build/shaders/proca/evolveProcaField.comp`

```glsl
#version 450
#extension GL_GOOGLE_include_directive : enable

/**
 * Evolve Proca field A_μ forward in time
 *
 * Equation (Lorenz gauge):
 *   (∂_t² - ∇² + m_γ²) A^μ = j^μ
 *
 * Discretization (leap-frog):
 *   A^μ(t+dt) = 2A^μ(t) - A^μ(t-dt) + dt² [∇²A^μ - m_γ² A^μ + j^μ]
 */

layout(local_size_x = 8, local_size_y = 8, local_size_z = 8) in;

layout(set = 0, binding = 0) readonly buffer ACurrent {
    vec4 A_current[];
};

layout(set = 0, binding = 1) readonly buffer APrevious {
    vec4 A_previous[];
};

layout(set = 0, binding = 2) readonly buffer RField {
    float R[];
};

layout(set = 0, binding = 3) readonly buffer JCurrent {
    vec4 j_mu[];
};

layout(set = 0, binding = 4) writeonly buffer ANext {
    vec4 A_next[];
};

layout(push_constant) uniform Params {
    uint Nx, Ny, Nz;
    float dx, dy, dz, dt;
    float g;  // Photon mass coupling
} params;

void main() {
    uint i = gl_GlobalInvocationID.x;
    uint j = gl_GlobalInvocationID.y;
    uint k = gl_GlobalInvocationID.z;

    if (i >= params.Nx || j >= params.Ny || k >= params.Nz) return;

    uint idx = k * (params.Nx * params.Ny) + j * params.Nx + i;

    // Get photon mass m_γ = g·(1-R)
    float R_val = R[idx];
    float m_gamma = params.g * (1.0 - R_val);
    float m_gamma_sq = m_gamma * m_gamma;

    // Get A_μ at current and previous times
    vec4 A_curr = A_current[idx];
    vec4 A_prev = A_previous[idx];

    // Compute Laplacian ∇²A_μ (3D, periodic BC)
    // Need neighbors in all 6 directions
    uint i_left = (i == 0) ? (params.Nx - 1) : (i - 1);
    uint i_right = (i == params.Nx - 1) ? 0 : (i + 1);
    // ... (similar for j, k)

    uint idx_left = k * (params.Nx * params.Ny) + j * params.Nx + i_left;
    uint idx_right = k * (params.Nx * params.Ny) + j * params.Nx + i_right;
    // ... (6 total neighbor indices)

    vec4 A_left = A_current[idx_left];
    vec4 A_right = A_current[idx_right];
    // ... (6 total neighbor values)

    // Laplacian (7-point stencil):
    //   ∇²A ≈ (A_left + A_right - 2A_center)/dx²
    //       + (A_down + A_up - 2A_center)/dy²
    //       + (A_back + A_front - 2A_center)/dz²

    vec4 laplacian_A =
        (A_left + A_right - 2.0*A_curr) / (params.dx * params.dx) +
        (A_down + A_up - 2.0*A_curr) / (params.dy * params.dy) +
        (A_back + A_front - 2.0*A_curr) / (params.dz * params.dz);

    // Get source j_μ
    vec4 j = j_mu[idx];

    // Leap-frog update:
    //   A(t+dt) = 2A(t) - A(t-dt) + dt²·[∇²A - m²A + j]
    vec4 A_next_val = 2.0*A_curr - A_prev
                    + params.dt*params.dt * (laplacian_A - m_gamma_sq*A_curr + j);

    A_next[idx] = A_next_val;
}
```

### 3.3 Shader: enforceGaugeCondition.comp

**Purpose**: Project A_μ to satisfy Lorenz gauge ∂_μ A^μ = 0

**Strategy**: Solve Poisson equation for gauge function χ, then A_μ → A_μ - ∂_μχ

```glsl
/**
 * Enforce Lorenz gauge condition ∂_μ A^μ = 0
 *
 * Method:
 * 1. Compute divergence div_A = ∂_μ A^μ
 * 2. Solve Poisson equation ∇²χ = div_A for gauge function χ
 * 3. Project: A_μ → A_μ - ∂_μχ
 *
 * After projection: ∂_μ(A^μ - ∂^μχ) = div_A - ∇²χ = 0
 */

// Implementation details: Use Jacobi iteration for Poisson solve
// (Could upgrade to multigrid if needed for performance)
```

---

## 4. Integration with TRDCore

### 4.1 Modified TRDCore Class

**File**: `include/TRDCore.h` (modifications)

```cpp
class TRDCore {
public:
    // ... existing TRD methods ...

    /**
     * @brief Enable electromagnetic coupling
     * @param emTheory Gauge theory implementation (Proca, Maxwell, etc.)
     */
    void enableEMCoupling(std::unique_ptr<physics::GaugeTheory> emTheory);

    /**
     * @brief Disable EM coupling (pure TRD mode)
     */
    void disableEMCoupling();

    /**
     * @brief Check if EM coupling is active
     */
    bool isEMCoupled() const { return m_em_theory != nullptr; }

    /**
     * @brief Get EM theory (const access)
     */
    const physics::GaugeTheory* getEMTheory() const { return m_em_theory.get(); }

private:
    // EM coupling (optional)
    std::unique_ptr<physics::GaugeTheory> m_em_theory;

    // Modified evolution step (operator splitting)
    void evolveStepWithEM(float dt);
};
```

### 4.2 Operator Splitting Evolution

**File**: `src/TRDCore.cpp`

```cpp
void TRDCore::evolveStepWithEM(float dt) {
    if (!isEMCoupled()) {
        // Pure TRD evolution (existing code path)
        evolveStepPure(dt);
        return;
    }

    // Strang splitting: S(dt/2) → EM(dt) → S(dt/2)

    // 1. TRD half-step (dt/2) without EM backreaction
    evolveKuramotoHalfStep(dt / 2.0f);

    // 2. EM full step (dt) with current j_μ from θ, R
    FieldState state;
    state.theta_field = m_theta_buffer;
    state.R_field = m_R_buffer;
    state.Nx = m_grid_size[0];
    state.Ny = m_grid_size[1];
    state.Nz = m_grid_size[2];
    state.dx = m_spatial_step[0];
    state.dy = m_spatial_step[1];
    state.dz = m_spatial_step[2];
    state.dt = dt;

    m_em_theory->computePotentials(state);
    m_em_theory->evolveFields(dt);
    m_em_theory->computeFieldStrengths(state);

    // 3. TRD half-step (dt/2) WITH EM backreaction (if α ≠ 0)
    if (m_em_backreaction_enabled) {
        applyEMBackreaction(state);
    }
    evolveKuramotoHalfStep(dt / 2.0f);
}
```

---

## 5. Build System Integration

### 5.1 CMakeLists.txt Additions

```cmake
# Proca EM module
set(PROCA_SOURCES
    src/physics/GaugeTheory.cpp
    src/physics/ProcaEM.cpp
    src/physics/PhysicsTypes.cpp
)

# Add to TRD target
target_sources(TRD PRIVATE ${PROCA_SOURCES})

# Include directories
target_include_directories(TRD PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/include/physics
)

# Compile Proca shaders
set(PROCA_SHADERS
    shaders/proca/computeEMCurrent.comp
    shaders/proca/evolveProcaField.comp
    shaders/proca/computeFieldStrengths.comp
    shaders/proca/enforceGaugeCondition.comp
)

foreach(SHADER ${PROCA_SHADERS})
    get_filename_component(SHADER_NAME ${SHADER} NAME_WE)
    add_custom_command(
        OUTPUT ${CMAKE_BINARY_DIR}/shaders/proca/${SHADER_NAME}.spv
        COMMAND glslangValidator -V ${CMAKE_SOURCE_DIR}/${SHADER}
                -o ${CMAKE_BINARY_DIR}/shaders/proca/${SHADER_NAME}.spv
        DEPENDS ${CMAKE_SOURCE_DIR}/${SHADER}
    )
    list(APPEND PROCA_SHADER_SPIRV ${CMAKE_BINARY_DIR}/shaders/proca/${SHADER_NAME}.spv)
endforeach()

add_custom_target(proca_shaders DEPENDS ${PROCA_SHADER_SPIRV})
add_dependencies(TRD proca_shaders)
```

---

## 6. Testing Infrastructure

### 6.1 Unit Test: Vortex B-field

**File**: `test/test_proca_vortex_bfield.cpp`

```cpp
#include "physics/ProcaEM.h"
#include "FieldState.h"
#include <gtest/gtest.h>

/**
 * Critical test: Can θ vortex generate B ≠ 0?
 *
 * Setup:
 * 1. Create θ(x,y) = atan2(y-y_c, x-x_c) (vortex winding n=1)
 * 2. Uniform R = 0.9 (synchronized)
 * 3. Evolve Proca equation
 * 4. Measure B_z in neighborhood of vortex
 *
 * Success: |B_z|_max > 0.01, Φ_B within 20% of theory
 */
TEST(ProcaEM, VortexGeneratesBField) {
    // Grid setup
    const uint32_t Nx = 128, Ny = 128, Nz = 1;
    const float dx = 0.1f, dy = 0.1f, dz = 0.1f;

    // Proca parameters
    const float g = 0.1f;      // Photon mass coupling
    const float alpha = 0.01f; // Current coupling

    // Create Proca module
    ProcaEM proca(device, {Nx, Ny, Nz}, {dx, dy, dz}, g, alpha);

    // Initialize θ vortex
    const float x_core = Nx * dx / 2.0f;
    const float y_core = Ny * dy / 2.0f;

    FieldState state;
    // ... (initialize theta_field with vortex, R_field uniform) ...

    // Evolve EM fields
    proca.computePotentials(state);
    for (int step = 0; step < 100; ++step) {
        proca.evolveFields(0.01f);
    }
    proca.computeFieldStrengths(state);

    // Measure B field
    float B_z_max = measureMaxBField(state.B_field, Nx, Ny);
    float Phi_B = proca.computeMagneticFlux();

    // Theoretical predictions
    float m_gamma = g * (1.0f - 0.9f);  // R = 0.9
    float Phi_B_theory = 2.0f * M_PI / m_gamma;

    // Assertions
    EXPECT_GT(B_z_max, 0.01f) << "B field must be non-zero!";
    EXPECT_NEAR(Phi_B, Phi_B_theory, 0.2f * Phi_B_theory)
        << "Magnetic flux within 20% of Proca prediction";
}
```

### 6.2 Integration Test: Boris with Emergent B

**File**: `test/test_boris_emergent_bfield.cpp`

```cpp
/**
 * DEFINITIVE TEST: Boris algorithm with emergent B field
 *
 * This is THE test that determines GO/NO-GO for EM emergence.
 *
 * Setup:
 * 1. Create θ vortex → emergent B field
 * 2. Place test particle in region with B ≠ 0
 * 3. Evolve with Boris pusher using emergent E, B (not external!)
 * 4. Verify cyclotron motion
 *
 * Success criteria:
 * - Circular trajectory (r_std/r_mean < 0.05)
 * - Larmor radius r_L = mv/(qB) within 10%
 * - Energy conservation <0.1%
 * - Stable over 10+ orbits
 */
TEST(ProcaEM, BorisTestEmergentBField) {
    // ... (implementation) ...

    // CRITICAL ASSERTION
    ASSERT_TRUE(trajectory_is_circular)
        << "FAIL: Cannot produce cyclotron motion from emergent B";
    ASSERT_NEAR(r_measured, r_L_theory, 0.1f * r_L_theory)
        << "FAIL: Larmor radius incorrect";
}
```

---

## 7. Performance Considerations

### 7.1 Computational Cost Estimate

**Proca evolution cost per timestep**:
- Compute j_μ: ~1 GFLOP (gradient computation)
- Evolve A_μ: ~10 GFLOP (Laplacian + leap-frog)
- Enforce gauge: ~5 GFLOP (Jacobi iteration)
- Compute E, B: ~2 GFLOP (gradients)
- **Total: ~18 GFLOP/step**

**TRD evolution cost**:
- Kuramoto step: ~5 GFLOP
- Order parameter: ~2 GFLOP
- **Total: ~7 GFLOP/step**

**Performance ratio**: Proca/TRD ≈ 2.5×

**Mitigation strategies**:
1. Adaptive stepping: Only compute EM every N TRD steps
2. GPU optimization: Fuse kernels, optimize memory access
3. Conditional EM: Only enable in regions with R > threshold

### 7.2 Memory Footprint

**Per grid point**:
- TRD: θ (4B), R (4B) = 8 bytes
- EM: A_μ×3 (48B), j_μ (16B), E (12B), B (12B) = 88 bytes
- **Total with EM: 96 bytes/point**

**For 128³ grid**:
- Without EM: 16 MB
- With EM: 192 MB
- **Increase: 12×**

**Impact**: Still fits comfortably in GPU memory (modern GPUs have 8-24 GB)

---

## 8. Implementation Timeline

### Week 2 (Current):
- ✅ Design complete (this document)
- Day 1-2: Implement GaugeTheory interface + ProcaEM class skeleton
- Day 3-4: Write GPU shaders (computeEMCurrent, evolveProcaField)
- Day 5: Integration with TRDCore (operator splitting)
- Day 6-7: **CRITICAL TEST**: Vortex B-field (Test 2)

### Week 3:
- Day 8-9: Boris test with emergent B (Test 3) **← GO/NO-GO decision**
- Day 10-11: Gauge invariance test (Test 4)
- Day 12-14: Full regression testing, bug fixes

### Week 4:
- Day 15-16: Performance profiling and optimization
- Day 17-19: QA review, documentation
- Day 20-21: Final report, roadmap update

---

## 9. Risk Mitigation

### High-Risk Item: B = 0 even with Proca

**Early detection** (Week 2, Day 6):
- If vortex test gives B = 0 → IMMEDIATE pivot to Stückelberg
- Do NOT wait until Week 3

**Fallback plan**:
- Implement Stückelberg mechanism (adds 1 week)
- If that also fails → Document failure honestly, pivot to Option A

### Medium-Risk Item: Performance unacceptable

**Threshold**: >3× slowdown unacceptable

**Optimization plan**:
1. Profile GPU shaders (identify bottleneck)
2. Fuse kernels (reduce memory bandwidth)
3. Adaptive EM (only compute where needed)
4. If still >3×: Make EM optional feature, not default

---

## 10. Conclusion

**Architecture is modular, extensible, and testable.**

**Key design decisions**:
1. GaugeTheory interface supports multiple theories (Proca, Stückelberg, Maxwell)
2. GPU-first implementation (all field evolution on GPU)
3. Operator splitting (2nd-order Strang) for TRD-EM coupling
4. Early testing (vortex B-field Week 2) to detect failures quickly

**Next action**: Begin implementation of GaugeTheory.h and ProcaEM.h

**Status**: DESIGN COMPLETE - Ready for development.
