# Stochastic Shader Design for GPU Implementation

**Date:** 2025-12-17
**Status:** Implementation Architecture
**Context:** Building on validated noise parameters (σ_c ≈ 0.65-0.80)

---

## Executive Summary

This document outlines the GPU shader architecture for implementing stochastic MSFT with Dirac coupling. The design leverages the existing validated `kuramoto_stochastic.comp` shader and extends it to include stochastic Dirac evolution while maintaining the robust synchronization demonstrated in our noise sweep experiments.

---

## 1. Shader Architecture Overview

### 1.1 Shader Pipeline Structure

```
┌─────────────────────────────────────────────┐
│              Initialization                 │
│         (Random phases, spinors)            │
└─────────────────┬───────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────┐
│        kuramoto_stochastic.comp             │
│   (Phase evolution with noise σ_θ)          │
└─────────────────┬───────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────┐
│           sync_field.comp                   │
│   (Compute R(θ) order parameter)            │
└─────────────────┬───────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────┐
│         dirac_stochastic.comp               │
│   (Spinor evolution with noise σ_Ψ)         │
└─────────────────┬───────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────┐
│          gravity_field.comp                 │
│   (Compute gravitational corrections)        │
└─────────────────┬───────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────┐
│            Measurement                      │
│   (R_global, particle properties)           │
└─────────────────────────────────────────────┘
```

### 1.2 Buffer Layout

```glsl
// Phase field buffers
layout(set = 0, binding = 0) buffer ThetaBuffer {
    float theta[];      // Current phases [-π, π]
};

layout(set = 0, binding = 1) buffer ThetaOutBuffer {
    float theta_out[];  // Updated phases
};

// Synchronization field
layout(set = 0, binding = 2) buffer SyncBuffer {
    vec2 R_field[];     // Complex order parameter R = ⟨e^{iθ}⟩
};

// Spinor field buffers
layout(set = 0, binding = 3) buffer SpinorBuffer {
    vec4 psi[];         // Current spinor (4 complex components = 8 floats)
};

layout(set = 0, binding = 4) buffer SpinorOutBuffer {
    vec4 psi_out[];     // Updated spinor
};

// Coupling field
layout(set = 0, binding = 5) buffer CouplingBuffer {
    float lambda[];     // |Ψ|² feedback to phases
};

// Gravitational field (optional)
layout(set = 0, binding = 6) buffer GravityBuffer {
    vec4 g_field[];     // Gravitational 4-potential
};
```

---

## 2. Enhanced Kuramoto Stochastic Shader

### 2.1 Modifications to Existing `kuramoto_stochastic.comp`

The current implementation is already well-structured. Key enhancements needed:

```glsl
// Add spinor feedback coupling
layout(set = 0, binding = 5) readonly buffer CouplingBuffer {
    float lambda[];
};

// In main evolution:
void main() {
    // ... existing code ...

    // Get spinor feedback (NEW)
    float spinor_density = lambda[idx];

    // Modified drift term
    float drift = omega_i
                + coupling_force           // Kuramoto coupling
                + damping_force            // Phase damping
                - params.lambda_coupling * spinor_density;  // Spinor feedback

    // Stochastic term (unchanged)
    float noise = params.sigma * sqrt(params.dt) * randn();

    // Update
    float theta_new = theta_i + drift * params.dt + noise;

    // ... rest unchanged ...
}
```

### 2.2 Push Constants Structure

```glsl
layout(push_constant) uniform PushConstants {
    float dt;              // Time step (0.01)
    float K;               // Coupling strength (1.0)
    float sigma_theta;     // Phase noise amplitude (0.01-0.05)
    float sigma_psi;       // Spinor noise amplitude
    float damping;         // Phase damping γ (0.1)
    float lambda_coupling; // Spinor-phase coupling λ
    float Delta;           // Mass gap parameter (2.5)
    float omega_mean;      // Mean frequency
    uint Nx;               // Grid width
    uint Ny;               // Grid height
    uint time_step;        // Current timestep for PRNG
} params;
```

---

## 3. New Dirac Stochastic Shader

### 3.1 `dirac_stochastic.comp` Structure

```glsl
#version 450
#extension GL_ARB_separate_shader_objects : enable

#include "../include/complex.glsl"
#include "../include/precision.glsl"

/**
 * Stochastic Dirac evolution with MSR noise
 * i∂_t Ψ = H_Dirac Ψ + σ_Ψ ξ(t)
 */

layout(local_size_x = 16, local_size_y = 16, local_size_z = 1) in;

// Dirac matrices (Pauli representation)
const mat2 sigma_x = mat2(0, 1, 1, 0);
const mat2 sigma_y = mat2(0, -1, 1, 0);  // Without i factor
const mat2 sigma_z = mat2(1, 0, 0, -1);
const mat2 I_2 = mat2(1, 0, 0, 1);

// PCG PRNG (same as kuramoto_stochastic.comp)
uint pcg_state;

uint pcg_hash(uint seed) {
    uint state = seed * 747796405u + 2891336453u;
    uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
    return (word >> 22u) ^ word;
}

void pcg_init(uint seed) {
    pcg_state = pcg_hash(seed);
}

float pcg_random() {
    pcg_state = pcg_state * 747796405u + 2891336453u;
    uint word = ((pcg_state >> ((pcg_state >> 28u) + 4u)) ^ pcg_state) * 277803737u;
    word = (word >> 22u) ^ word;
    return float(word) / 4294967296.0;
}

// Complex Gaussian noise
vec2 complex_randn() {
    float u1 = max(pcg_random(), 1e-10);
    float u2 = pcg_random();

    const float PI = 3.14159265359;
    float r = sqrt(-2.0 * log(u1));
    float theta = 2.0 * PI * u2;

    return vec2(r * cos(theta), r * sin(theta));
}

// Compute spatial derivatives using finite differences
vec4 compute_gradient(uint x, uint y) {
    // Load neighboring spinor values
    uint idx_center = y * params.Nx + x;
    uint idx_left   = y * params.Nx + ((x + params.Nx - 1) % params.Nx);
    uint idx_right  = y * params.Nx + ((x + 1) % params.Nx);
    uint idx_up     = ((y + params.Ny - 1) % params.Ny) * params.Nx + x;
    uint idx_down   = ((y + 1) % params.Ny) * params.Nx + x;

    vec4 psi_center = psi[idx_center];
    vec4 psi_left   = psi[idx_left];
    vec4 psi_right  = psi[idx_right];
    vec4 psi_up     = psi[idx_up];
    vec4 psi_down   = psi[idx_down];

    // Centered differences
    vec4 grad_x = (psi_right - psi_left) / 2.0;
    vec4 grad_y = (psi_down - psi_up) / 2.0;

    // Apply Dirac matrices: -i γ^0 γ^i ∂_i
    // In 2D: γ^0 γ^1 = σ_x ⊗ σ_z, γ^0 γ^2 = σ_y ⊗ σ_z

    vec4 result = vec4(0.0);

    // x-derivative contribution
    result.xy += mat2(sigma_x) * grad_x.xy;  // Upper components
    result.zw -= mat2(sigma_x) * grad_x.zw;  // Lower components (σ_z factor)

    // y-derivative contribution
    result.xy += mat2(0, -1, 1, 0) * grad_y.xy;  // σ_y without i
    result.zw -= mat2(0, -1, 1, 0) * grad_y.zw;  // With σ_z

    return result;
}

void main() {
    uvec2 global_id = gl_GlobalInvocationID.xy;

    if (global_id.x >= params.Nx || global_id.y >= params.Ny) {
        return;
    }

    uint idx = global_id.y * params.Nx + global_id.x;

    // Initialize PRNG
    uint seed = global_id.x + global_id.y * 65536u + params.time_step * 16777216u;
    pcg_init(seed);

    // Load current spinor
    vec4 psi_current = psi[idx];

    // Load synchronization field R(θ)
    vec2 R = R_field[idx];
    float R_magnitude = length(R);

    // Compute Hamiltonian terms
    vec4 kinetic = compute_gradient(global_id.x, global_id.y);

    // Mass term: Δ·R(θ)·Ψ
    vec4 mass = params.Delta * R_magnitude * psi_current;

    // Total Hamiltonian evolution
    vec4 H_psi = kinetic + mass;

    // Stochastic term: σ_Ψ·√(dt)·ξ(t)
    vec4 noise;
    noise.xy = complex_randn() * params.sigma_psi * sqrt(params.dt);
    noise.zw = complex_randn() * params.sigma_psi * sqrt(params.dt);

    // Euler-Maruyama update: Ψ(t+dt) = Ψ(t) - i·H·Ψ·dt + noise
    vec4 psi_new = psi_current;

    // Deterministic evolution (complex multiplication by -i)
    psi_new.x -= H_psi.y * params.dt;  // Real part
    psi_new.y += H_psi.x * params.dt;  // Imaginary part
    psi_new.z -= H_psi.w * params.dt;
    psi_new.w += H_psi.z * params.dt;

    // Add noise
    psi_new += noise;

    // Normalize to conserve probability (optional, for stability)
    float norm = length(psi_new);
    if (norm > 0.001) {
        psi_new /= sqrt(norm);  // Soft normalization
    }

    // Write output
    psi_out[idx] = psi_new;

    // Compute |Ψ|² for feedback to phase field
    lambda[idx] = dot(psi_new, psi_new);
}
```

---

## 4. Noise Implementation Details

### 4.1 PRNG Strategy

**PCG (Permuted Congruential Generator):**
- Period: 2³²
- Quality: Passes all statistical tests
- Performance: Single multiplication + XOR per number
- Seeding: Spatial position + timestep ensures decorrelation

### 4.2 Box-Muller Transform

For Gaussian noise generation:
```glsl
float randn() {
    float u1 = max(pcg_random(), 1e-10);  // Avoid log(0)
    float u2 = pcg_random();

    float r = sqrt(-2.0 * log(u1));
    float theta = 2.0 * PI * u2;

    return r * cos(theta);  // N(0,1)
}
```

### 4.3 Euler-Maruyama Scaling

**Critical:** Noise must be scaled by √(dt) for proper Brownian motion:
```glsl
// CORRECT (Euler-Maruyama)
update = drift * dt + sigma * sqrt(dt) * randn();

// INCORRECT (wrong scaling)
update = drift * dt + sigma * randn();  // Missing sqrt(dt)!
```

---

## 5. Synchronization Field Computation

### 5.1 Existing `sync_field.comp` (No Changes Needed)

The current implementation correctly computes:
```glsl
// Local order parameter
vec2 R_local = vec2(0.0);
for (int dy = -radius; dy <= radius; dy++) {
    for (int dx = -radius; dx <= radius; dx++) {
        uint nx = (x + dx + Nx) % Nx;
        uint ny = (y + dy + Ny) % Ny;
        float theta_j = theta[ny * Nx + nx];
        R_local += vec2(cos(theta_j), sin(theta_j));
    }
}
R_local /= float((2*radius+1) * (2*radius+1));
R_field[idx] = R_local;
```

This provides the R(θ) field needed for Dirac mass generation.

---

## 6. Gravity Field Integration

### 6.1 Existing `gravity_field.comp` (Minor Updates)

Add stochastic corrections:
```glsl
// Compute Ricci scalar (existing)
float R = compute_ricci_scalar();

// Add noise correction to account for fluctuations
float noise_correction = params.sigma_theta * params.sigma_theta;
R *= (1.0 - noise_correction);  // Reduced effective curvature

// Compute gravitational field (existing)
vec4 g = compute_gravity_from_ricci(R);
```

---

## 7. Integration Scheme

### 7.1 Time-Stepping Order

For each timestep:
1. **Phase evolution** (`kuramoto_stochastic.comp`)
   - Input: θ(t), λ(t)
   - Output: θ(t+dt)

2. **Synchronization** (`sync_field.comp`)
   - Input: θ(t+dt)
   - Output: R(t+dt)

3. **Spinor evolution** (`dirac_stochastic.comp`)
   - Input: Ψ(t), R(t+dt)
   - Output: Ψ(t+dt), λ(t+dt)

4. **Gravity** (`gravity_field.comp`) - Optional
   - Input: R(t+dt), |Ψ(t+dt)|²
   - Output: g_μν corrections

### 7.2 Synchronization Requirements

```cpp
// C++ dispatch code
void MSFTEngine::stepStochastic(float dt, float sigma_theta, float sigma_psi) {
    // Step 1: Phase evolution
    dispatchKuramatoStochastic();
    vkCmdPipelineBarrier(...);  // Ensure completion

    // Step 2: Compute sync field
    dispatchSyncField();
    vkCmdPipelineBarrier(...);

    // Step 3: Spinor evolution
    dispatchDiracStochastic();
    vkCmdPipelineBarrier(...);

    // Step 4: Gravity (if enabled)
    if (gravity_enabled) {
        dispatchGravityField();
        vkCmdPipelineBarrier(...);
    }

    // Swap buffers for next iteration
    swapBuffers();
}
```

---

## 8. Performance Optimization

### 8.1 Shared Memory Usage

Leverage shared memory for local computations:
```glsl
shared float s_theta[18][18];   // 16+2 border for phases
shared vec4 s_psi[18][18];      // For spinor neighbors
```

### 8.2 Workgroup Sizing

Optimal for modern GPUs:
- Workgroup: 16×16 = 256 threads
- Occupancy: 4-8 workgroups per SM
- Memory: ~32KB shared per workgroup

### 8.3 Memory Access Patterns

```glsl
// Coalesced reads
float theta_block[16] = theta[base_idx : base_idx + 16];

// Bank conflict avoidance
s_theta[ly][lx] = theta_block[lx];  // No conflicts
```

---

## 9. Validation and Testing

### 9.1 Noise Amplitude Tests

Validate implementation with known results:
```cpp
// Test 1: Zero noise -> deterministic evolution
engine.stepStochastic(dt=0.01, sigma_theta=0.0, sigma_psi=0.0);
assert(R_final == 1.000);  // Perfect sync

// Test 2: Small noise -> maintained sync
engine.stepStochastic(dt=0.01, sigma_theta=0.05, sigma_psi=0.05);
assert(R_final > 0.95);  // Still synchronized

// Test 3: Critical noise -> transition
engine.stepStochastic(dt=0.01, sigma_theta=0.65, sigma_psi=0.65);
assert(R_final ≈ 0.7);  // Near critical point

// Test 4: Large noise -> desync
engine.stepStochastic(dt=0.01, sigma_theta=1.0, sigma_psi=1.0);
assert(R_final < 0.4);  // Thermal state
```

### 9.2 Conservation Checks

```glsl
// Spinor norm conservation (approximately)
float norm_initial = sum(|Ψ|²);
stepStochastic(...);
float norm_final = sum(|Ψ|²);
assert(abs(norm_final - norm_initial) / norm_initial < 0.01);
```

### 9.3 Statistical Tests

Verify noise properties:
```cpp
// Collect noise samples over many steps
vector<float> noise_samples;
for (int i = 0; i < 10000; i++) {
    noise_samples.push_back(generate_noise());
}

// Check Gaussian distribution
assert(mean(noise_samples) ≈ 0.0);
assert(std_dev(noise_samples) ≈ 1.0);
assert(passes_kolmogorov_smirnov_test(noise_samples));
```

---

## 10. Debugging Support

### 10.1 Diagnostic Buffers

```glsl
// Debug output buffer
layout(set = 0, binding = 7) buffer DebugBuffer {
    float debug[];
};

// In shader:
debug[idx] = pcg_state;        // Check PRNG state
debug[idx] = drift;            // Verify deterministic terms
debug[idx] = noise;            // Check noise amplitude
debug[idx] = length(psi_new);  // Monitor spinor norm
```

### 10.2 Visualization Modes

```cpp
enum VisualizationMode {
    PHASE_FIELD,        // θ(x,y)
    SYNC_FIELD,         // |R(x,y)|
    SPINOR_DENSITY,     // |Ψ(x,y)|²
    NOISE_AMPLITUDE,    // Local σ_eff
    ENERGY_DENSITY,     // H(x,y)
    GRADIENT_FLOW       // ∇θ, ∇|Ψ|²
};
```

---

## 11. Summary

This shader architecture provides:

1. **Robust noise implementation** using PCG + Box-Muller
2. **Proper Euler-Maruyama** integration with √(dt) scaling
3. **Efficient GPU parallelization** with shared memory
4. **Flexible parameter control** via push constants
5. **Conservation properties** through careful numerics

The design leverages the validated σ_c ≈ 0.65 threshold to operate safely at σ = 0.05 (13× margin), ensuring stable particle formation while incorporating realistic vacuum fluctuations.

**Implementation priority:**
1. Enhance existing `kuramoto_stochastic.comp` with spinor feedback
2. Create new `dirac_stochastic.comp` shader
3. Update C++ dispatch code for proper synchronization
4. Validate against CPU reference implementation
5. Measure particle formation under stochastic dynamics

---

**End of Document**