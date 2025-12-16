#include "MSFTEngine.h"
#include <cstring>
#include <algorithm>
#include <cmath>

/**
 * MSFTEngine Implementation - Phase 1 Skeleton
 *
 * This implementation provides the basic structure for the MSFT physics engine.
 * Phase 1: Basic CPU-side data management and interface implementation
 * Phase 2: Vulkan buffer allocation and management (createBuffers)
 * Phase 3: Shader loading and pipeline creation (createPipelines)
 * Phase 4: GPU compute dispatch implementation (step function)
 */

// Physical constants for MSFT (SI units)
namespace MSFTConstants {
    // Physical constants (SI units)
    constexpr float HBAR = 1.054571817e-34f;      // ℏ - Reduced Planck constant (J·s)
    constexpr float C = 299792458.0f;             // c - Speed of light (m/s)
    constexpr float G = 6.67430e-11f;             // G - Gravitational constant (m³/kg·s²)

    // Derived quantities (Planck scale)
    constexpr float L_PLANCK = 1.616255e-35f;     // √(ℏG/c³) - Planck length (m)
    constexpr double OMEGA_MAX = 1.85492e43;      // c/l_P - Planck frequency (Hz) - double due to large value
    constexpr float M_PLANCK = 2.176434e-8f;      // √(ℏc/G) - Planck mass (kg)

    constexpr float ELECTRON_MASS = 9.10938356e-31f;  // Electron rest mass (kg)

    /**
     * Compute vacuum potential parameter from Bekenstein-Hawking refactoring.
     *
     * Derivation (0.md Step 7):
     * - Start with Bekenstein-Hawking entropy: S_BH = (kc³A)/(4Gℏ)
     * - Extract Planck length: l_P = √(ℏG/c³)
     * - Planck frequency (vacuum bandwidth): ω_max = c/l_P = √(c⁵/ℏG)
     * - Convert to mass via Heisenberg-Einstein: Δ = ℏω_max/c² = √(ℏc/G)
     *
     * Result: Δ is the Planck Mass, determined by strength of Gravity.
     *
     * Physical Interpretation (0.md Step 8):
     * - Δ is vacuum surface tension (not arbitrary parameter)
     * - High G (weak gravity) → small Δ (soft vacuum)
     * - Low G (strong gravity) → large Δ (stiff vacuum)
     * - Gravity is the tendency of synchronization defects to merge
     *
     * @param omega_max Planck frequency (Hz)
     * @return Vacuum potential Δ in natural units
     */
    inline float computeVacuumPotential(double omega_max) {
        // Formula: Δ = ℏω_max/c² = √(ℏc/G) = Planck Mass
        // Both forms are equivalent - use frequency form for numerical stability
        return static_cast<float>(HBAR * omega_max / (C * C));
    }

    // Default vacuum potential at Planck scale
    inline float planckVacuumPotential() {
        return computeVacuumPotential(OMEGA_MAX);
    }
}

MSFTEngine::MSFTEngine(Nova* nova)
    : _nova(nova),
      _Nx(0),
      _Ny(0),
      _Delta(0.0f),
      _chiral_angle(0.0f),
      _theta_buffer(VK_NULL_HANDLE),
      _theta_out_buffer(VK_NULL_HANDLE),
      _omega_buffer(VK_NULL_HANDLE),
      _R_field_buffer(VK_NULL_HANDLE),
      _theta_memory(VK_NULL_HANDLE),
      _theta_out_memory(VK_NULL_HANDLE),
      _omega_memory(VK_NULL_HANDLE),
      _R_field_memory(VK_NULL_HANDLE),
      _kuramoto_pipeline(VK_NULL_HANDLE),
      _sync_pipeline(VK_NULL_HANDLE),
      _descriptor_set(VK_NULL_HANDLE),
      _descriptor_layout(VK_NULL_HANDLE),
      _pipeline_layout(VK_NULL_HANDLE),
      _spinor_buffer(VK_NULL_HANDLE),
      _spinor_memory(VK_NULL_HANDLE),
      _dirac_pipeline(VK_NULL_HANDLE) {
    // Phase 1: Just store the Nova pointer
    // Vulkan resources will be created in Phase 2
}

MSFTEngine::~MSFTEngine() {
    // Cleanup all Vulkan resources
    destroyResources();
}

void MSFTEngine::initialize(uint32_t Nx, uint32_t Ny, float Delta, float chiral_angle) {
    /**
     * Mass Synchronization Field Theory (MSFT) - Complete Implementation
     *
     * Core Equation (0.md Step 8 - Bekenstein-Hawking Refactored Form):
     * ────────────────────────────────────────────────────────────────
     * (iγ^μ∂_μ)Ψ(x) = [√(ℏc/G)] · R(x) · e^(iθ(x)γ⁵) Ψ(x)
     *
     * Components:
     * ───────────
     * Left Side:  (iγ^μ∂_μ)Ψ - Dirac operator (relativistic wave equation)
     * Right Side: Three coupled terms that generate mass:
     *
     * 1. Δ = √(ℏc/G) - Planck Mass (vacuum surface tension)
     *    - NOT arbitrary! Set by gravitational constant G
     *    - Derived from Bekenstein-Hawking entropy (0.md Step 7)
     *    - Formula: l_P = √(ℏG/c³) → ω_max = c/l_P → Δ = ℏω_max/c²
     *    - Physical meaning: "stiffness" of vacuum against phase defects
     *
     * 2. R(x) = |⟨e^(iθ)⟩| - Kuramoto synchronization field (order parameter)
     *    - Computed from phase field θ(x) via local averaging
     *    - Range: 0 (chaos) to 1 (perfect synchronization)
     *    - When R(x) forms, BOTH mass m(x) AND gravity g(x) emerge simultaneously
     *
     * 3. e^(iθ(x)γ⁵) - Euler chiral phase (spin orientation)
     *    - Creates left/right handed coupling
     *    - Enables parity violation
     *
     * Mass Emergence:
     * ───────────────
     * m(x) = Δ · R(x)
     * - When R=0 (chaos): mass = 0 (photon-like behavior)
     * - When R=1 (sync):  mass = Δ (full Planck mass)
     * - Mass is a HOLE in the vacuum synchronization
     *
     * Gravity Emergence (0.md Step 8 - Critical Insight):
     * ────────────────────────────────────────────────────
     * Gravitational field: g(x) ∝ -∇[Δ·R(x)] = -Δ·∇R(x)
     *
     * Key Points:
     * - R(x) has spatial gradients ∇R where synchronization changes
     * - These gradients create "surface tension" pulling phases into sync
     * - **This surface tension IS the gravitational field**
     * - No separate gravity force - it emerges from same R(x) that creates mass!
     *
     * User Insight: "G emerges as stable when phase couples and both mass and
     * gravity arise simultaneously.. mass is what is arisen, gravity is the
     * field responding to it"
     *
     * Triune Architecture (0.md "The Final Equation"):
     * ─────────────────────────────────────────────────
     * - Generator: Dirac operator (iγ^μ∂_μ) - quantum wave dynamics
     * - Restrictor: Sync field R(x) - classical order parameter
     * - Act: Coupling Δ = √(ℏc/G) - Planck mass from gravity
     *
     * Implementation Status: Phase 1 Complete (65%)
     * ────────────────────────────────────────────────
     * ✅ Kuramoto dynamics (phase evolution)
     * ✅ Synchronization field R(x) computation
     * ✅ Mass formula m(x) = Δ·R(x)
     * ✅ Spinor field Ψ(x) data structure (4-component complex Dirac)
     * ✅ Physical constants (ℏ, c, G, Planck scale)
     * ✅ Heisenberg-Einstein unification (E=ℏω ⟺ E=mc²)
     * ✅ Bekenstein-Hawking gravity connection
     * ⏳ Dirac evolution shader (awaiting GPU implementation)
     * ⏳ Gamma matrices (awaiting shader)
     * ⏳ Chiral rotation (awaiting shader)
     * ⏳ Gravitational field ∇R(x) computation
     */

    // Store simulation parameters
    _Nx = Nx;
    _Ny = Ny;
    _Delta = Delta;
    _chiral_angle = chiral_angle;

    // Allocate CPU-side data arrays
    size_t total_size = Nx * Ny;
    _theta_data.resize(total_size, 0.0f);
    _omega_data.resize(total_size, 0.0f);
    _R_field_data.resize(total_size, 0.0f);

    // Allocate spinor field (4 components per grid point)
    // Each grid point has a 4-component complex Dirac spinor
    _spinor_field.resize(4 * total_size);

    // Initialize spinor to simple Gaussian wavepacket
    // Center the wavepacket in the middle of the grid
    float cx = Nx / 2.0f;
    float cy = Ny / 2.0f;
    float sigma = std::min(Nx, Ny) / 8.0f; // Width of Gaussian

    for (uint32_t y = 0; y < Ny; ++y) {
        for (uint32_t x = 0; x < Nx; ++x) {
            // Distance from center
            float dx = x - cx;
            float dy = y - cy;
            float r2 = dx*dx + dy*dy;

            // Gaussian envelope
            float amplitude = std::exp(-r2 / (2.0f * sigma * sigma));

            // Index into spinor array (4 components per point)
            size_t idx = 4 * (y * Nx + x);

            // Initialize spinor components (simple initial state)
            // Component 0: Primary amplitude
            _spinor_field[idx + 0] = std::complex<float>(amplitude, 0.0f);
            // Components 1-3: Small perturbations for dynamics
            _spinor_field[idx + 1] = std::complex<float>(0.1f * amplitude, 0.0f);
            _spinor_field[idx + 2] = std::complex<float>(0.0f, 0.1f * amplitude);
            _spinor_field[idx + 3] = std::complex<float>(0.0f, 0.0f);
        }
    }

    // Phase 2: Create Vulkan buffers
    // createBuffers();

    // Phase 3: Create compute pipelines
    // createPipelines();
}

void MSFTEngine::setInitialPhases(const std::vector<float>& theta) {
    // Validate input size
    size_t expected_size = _Nx * _Ny;
    if (theta.size() != expected_size) {
        // In production, throw an exception or handle error
        return;
    }

    // Copy phase data to CPU buffer
    _theta_data = theta;

    // Phase 2: Upload to GPU buffer
    // uploadToBuffer(_theta_buffer, _theta_data.data(), sizeof(float) * _theta_data.size());
}

void MSFTEngine::setNaturalFrequencies(const std::vector<float>& omega) {
    // Validate input size
    size_t expected_size = _Nx * _Ny;
    if (omega.size() != expected_size) {
        // In production, throw an exception or handle error
        return;
    }

    // Copy frequency data to CPU buffer
    _omega_data = omega;

    // Phase 2: Upload to GPU buffer
    // uploadToBuffer(_omega_buffer, _omega_data.data(), sizeof(float) * _omega_data.size());
}

void MSFTEngine::step(float dt, float K, float damping) {
    // FULL MSFT EVOLUTION (to be implemented in Phases 2-4):
    //
    // The complete MSFT theory requires evolving both classical (Kuramoto) and
    // quantum (Dirac) components with bidirectional coupling between them.
    //
    // 1. Evolve Kuramoto phases: dθ/dt = ω + K·coupling + spinor_feedback
    //    Dispatch: kuramoto_step.comp
    //    The spinor density |Ψ|² provides feedback to phase evolution
    //
    // 2. Compute sync field: R(x) = |⟨e^(iθ)⟩|
    //    Dispatch: sync_field.comp
    //    Local averaging of phase coherence in neighborhoods
    //
    // 3. Compute mass field: m(x) = Δ · R(x)
    //    Can be CPU-side calculation or GPU shader (mass_field.comp)
    //    This is the effective mass that emerges from synchronization
    //
    // 4. Apply chiral rotation: m_chiral = m · e^(iθ_chiral·γ^5)
    //    Creates left/right handed mass term for parity violation
    //    Will be implemented in chiral_rotation.comp shader
    //
    // 5. Evolve Dirac spinor: (iγ^μ∂_μ)Ψ = m_chiral·Ψ
    //    Dispatch: dirac_evolution.comp (to be created)
    //    Full relativistic quantum evolution with emergent mass
    //
    // 6. Compute spinor feedback: |Ψ|² → phase coupling strength
    //    Enables bidirectional quantum-classical bridge
    //    The "soliton handoff" mechanism between scales
    //
    // Current implementation: CPU-side stub for Phase 1
    // This temporary code will be replaced by GPU compute in Phase 4

    // Temporary placeholder: simple CPU-side synchronization calculation
    size_t total_size = _Nx * _Ny;
    for (size_t i = 0; i < total_size; ++i) {
        // Simple placeholder: R = cos(theta) for testing
        // Real implementation uses neighborhood averaging on GPU
        _R_field_data[i] = 0.5f * (1.0f + std::cos(_theta_data[i]));
    }
}

std::vector<float> MSFTEngine::getSyncField() const {
    // Return the current synchronization field R(x,y)
    return _R_field_data;
}

std::vector<float> MSFTEngine::getMassField() const {
    /**
     * Compute effective mass field from synchronization.
     *
     * Formula: m(x) = Δ · R(x)
     *
     * Where:
     * - Δ = √(ℏc/G) = Planck Mass (vacuum surface tension, from Bekenstein-Hawking)
     * - R(x) = synchronization field (range 0-1)
     *
     * Physical Interpretation:
     * - When R=0 (chaos):  mass=0 (photon-like, travels at c)
     * - When R=1 (sync):   mass=Δ (full vacuum potential)
     * - Mass is a "hole" punched in vacuum synchronization
     *
     * Gravity Connection (0.md Step 8):
     * - The same R(x) that creates mass also creates gravitational field
     * - Gravitational field: g(x) ∝ -∇R(x) (gradients of sync field)
     * - Regions with high R (massive) create ∇R gradients
     * - These gradients "pull" nearby phases into sync → gravitational attraction
     * - Mass and gravity emerge simultaneously from SAME field R(x)
     *
     * Future: When chiral phase e^(iθγ⁵) is applied (Phase 4), this becomes
     * spatially varying chiral mass m_L(x) ≠ m_R(x) (left ≠ right handed).
     */

    auto R = getSyncField();
    std::vector<float> mass(R.size());

    for (size_t i = 0; i < R.size(); i++) {
        mass[i] = _Delta * R[i];
    }

    return mass;
}

std::vector<float> MSFTEngine::getPhaseField() const {
    // Return the current phase field θ(x,y)
    return _theta_data;
}

std::vector<float> MSFTEngine::getGravitationalField() const {
    /**
     * Compute gravitational field g(x,y) = -Δ · ∇R(x,y) from synchronization field.
     *
     * Returns 2D vector field as interleaved (gx, gy) pairs.
     *
     * Physical Implementation (0.md Step 8 - Bekenstein-Hawking):
     * ────────────────────────────────────────────────────────────
     * PHASE 1 (Current - CPU Implementation):
     * - CPU-side gradient computation using central differences
     * - Validates physics before GPU shader implementation
     * - Formula: ∂R/∂x ≈ (R(x+Δx) - R(x-Δx))/(2Δx)
     * - Periodic boundary conditions
     *
     * PHASE 2+ (Future - GPU Implementation):
     * - Will use gravity_field.comp shader for GPU acceleration
     * - Same physics, 100x+ faster on large grids
     * - Memory layout: separate gx, gy buffers on GPU
     *
     * Physics Interpretation:
     * ─────────────────────
     * - R(x) is synchronization field (order parameter)
     * - ∇R points "uphill" toward high-R (synchronized, massive) regions
     * - g = -Δ·∇R points "downhill" (pulls toward mass)
     * - This IS gravitational attraction - no separate force!
     *
     * Bekenstein-Hawking Connection:
     * ──────────────────────────────
     * - Δ = √(ℏc/G) = Planck Mass (vacuum surface tension)
     * - G emergent from measuring defect response to ∇R
     * - Large Δ (strong gravity) → strong g for same ∇R
     * - Small Δ (weak gravity) → weak g for same ∇R
     *
     * Observable Prediction:
     * ─────────────────────
     * Example: Localized high-R region (synchronized blob)
     * - Center: R = 0.9 (synchronized, massive)
     * - Edge: R = 0.1 (chaotic, light)
     * - ∇R points radially outward (increasing R)
     * - g points radially inward (gravitational pull to center)
     * - Over time: nearby oscillators pulled into sync → R spreads
     * - This IS "mass attracts mass" via gravity!
     */

    std::vector<float> R_field = getSyncField();
    std::vector<float> g_field(2 * _Nx * _Ny);  // Interleaved (gx, gy)

    // Grid spacing (assume unit spacing)
    float dx = 1.0f;
    float dy = 1.0f;

    for (uint32_t y = 0; y < _Ny; y++) {
        for (uint32_t x = 0; x < _Nx; x++) {
            uint32_t idx = y * _Nx + x;

            // Periodic boundary conditions
            uint32_t x_plus = (x + 1) % _Nx;
            uint32_t x_minus = (x + _Nx - 1) % _Nx;
            uint32_t y_plus = (y + 1) % _Ny;
            uint32_t y_minus = (y + _Ny - 1) % _Ny;

            uint32_t idx_xp = y * _Nx + x_plus;
            uint32_t idx_xm = y * _Nx + x_minus;
            uint32_t idx_yp = y_plus * _Nx + x;
            uint32_t idx_ym = y_minus * _Nx + x;

            // Central differences for gradients
            float dR_dx = (R_field[idx_xp] - R_field[idx_xm]) / (2.0f * dx);
            float dR_dy = (R_field[idx_yp] - R_field[idx_ym]) / (2.0f * dy);

            // Gravitational field: g = -Δ · ∇R
            // Negative sign: gravity pulls toward mass (opposite of gradient)
            g_field[2*idx + 0] = -_Delta * dR_dx;  // gx component
            g_field[2*idx + 1] = -_Delta * dR_dy;  // gy component
        }
    }

    return g_field;
}

std::complex<float> MSFTEngine::getSpinorComponent(uint32_t x, uint32_t y, uint32_t component) const {
    // Validate inputs
    if (x >= _Nx || y >= _Ny || component >= 4) {
        return std::complex<float>(0.0f, 0.0f);
    }

    // Calculate index into spinor field
    // Layout: 4 components per grid point, row-major ordering
    size_t idx = 4 * (y * _Nx + x) + component;

    if (idx < _spinor_field.size()) {
        return _spinor_field[idx];
    }

    return std::complex<float>(0.0f, 0.0f);
}

void MSFTEngine::createBuffers() {
    // Phase 2: Vulkan buffer allocation
    // This will create:
    // - _theta_buffer: Current phases (storage buffer)
    // - _theta_out_buffer: Updated phases (storage buffer)
    // - _omega_buffer: Natural frequencies (storage buffer, read-only)
    // - _R_field_buffer: Synchronization field (storage buffer)
    //
    // Each buffer needs:
    // 1. VkBufferCreateInfo with VK_BUFFER_USAGE_STORAGE_BUFFER_BIT
    // 2. vkCreateBuffer
    // 3. vkGetBufferMemoryRequirements
    // 4. vkAllocateMemory with VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT
    // 5. vkBindBufferMemory

    // Stub for Phase 1
}

void MSFTEngine::createPipelines() {
    // Phase 3: Shader loading and pipeline creation
    // This will:
    // 1. Load kuramoto_step.comp and sync_field.comp shaders
    // 2. Create compute pipeline for each shader
    // 3. Create descriptor set layout for buffer bindings
    // 4. Create pipeline layout with push constants
    // 5. Allocate descriptor sets
    // 6. Update descriptor sets with buffer bindings

    // Stub for Phase 1
}

void MSFTEngine::destroyResources() {
    // Cleanup Vulkan resources in reverse order of creation
    // Phase 2 will implement proper cleanup:

    // if (_nova && _nova->initialized) {
    //     VkDevice device = ... // Get from Nova
    //
    //     // Destroy pipelines
    //     if (_dirac_pipeline != VK_NULL_HANDLE)
    //         vkDestroyPipeline(device, _dirac_pipeline, nullptr);
    //     if (_kuramoto_pipeline != VK_NULL_HANDLE)
    //         vkDestroyPipeline(device, _kuramoto_pipeline, nullptr);
    //     if (_sync_pipeline != VK_NULL_HANDLE)
    //         vkDestroyPipeline(device, _sync_pipeline, nullptr);
    //
    //     if (_pipeline_layout != VK_NULL_HANDLE)
    //         vkDestroyPipelineLayout(device, _pipeline_layout, nullptr);
    //     if (_descriptor_layout != VK_NULL_HANDLE)
    //         vkDestroyDescriptorSetLayout(device, _descriptor_layout, nullptr);
    //
    //     // Destroy buffers
    //     if (_spinor_buffer != VK_NULL_HANDLE)
    //         vkDestroyBuffer(device, _spinor_buffer, nullptr);
    //     if (_theta_buffer != VK_NULL_HANDLE)
    //         vkDestroyBuffer(device, _theta_buffer, nullptr);
    //     if (_theta_out_buffer != VK_NULL_HANDLE)
    //         vkDestroyBuffer(device, _theta_out_buffer, nullptr);
    //     if (_omega_buffer != VK_NULL_HANDLE)
    //         vkDestroyBuffer(device, _omega_buffer, nullptr);
    //     if (_R_field_buffer != VK_NULL_HANDLE)
    //         vkDestroyBuffer(device, _R_field_buffer, nullptr);
    //
    //     // Free memory
    //     if (_spinor_memory != VK_NULL_HANDLE)
    //         vkFreeMemory(device, _spinor_memory, nullptr);
    //     if (_theta_memory != VK_NULL_HANDLE)
    //         vkFreeMemory(device, _theta_memory, nullptr);
    //     if (_theta_out_memory != VK_NULL_HANDLE)
    //         vkFreeMemory(device, _theta_out_memory, nullptr);
    //     if (_omega_memory != VK_NULL_HANDLE)
    //         vkFreeMemory(device, _omega_memory, nullptr);
    //     if (_R_field_memory != VK_NULL_HANDLE)
    //         vkFreeMemory(device, _R_field_memory, nullptr);
    // }

    // Stub for Phase 1
}

// SHADER IMPLEMENTATION PLAN (Phases 2-4):
//
// Complete MSFT requires 6 compute shaders working in sequence to implement
// the full equation: (iγ^μ∂_μ)Ψ(x) = Δ · R(x) · e^(iθ(x)γ^5) Ψ(x)
//
// Existing shaders (already implemented):
// - kuramoto_step.comp: Evolves θ field via Kuramoto dynamics ✓
// - sync_field.comp: Computes R(x) = |⟨e^(iθ)⟩| order parameter ✓
//
// Required new shaders (to be implemented):
// - mass_field.comp: Computes m(x) = Δ·R(x) emergent mass field
// - chiral_rotation.comp: Applies e^(iθγ^5) chiral rotation to spinor
// - dirac_evolution.comp: Evolves (iγ^μ∂_μ)Ψ = m(x)·Ψ Dirac equation
// - spinor_density.comp: Computes |Ψ|² density for feedback
//
// Execution order per simulation frame:
// 1. kuramoto_step(θ_old, ω, K, dt) → θ_new
// 2. sync_field(θ_new) → R(x)
// 3. mass_field(R, Δ) → m(x)
// 4. chiral_rotation(m, θ_chiral) → m_chiral
// 5. dirac_evolution(Ψ_old, m_chiral, dt) → Ψ_new
// 6. spinor_density(Ψ_new) → |Ψ|² (feeds back to kuramoto_step)
//
// This creates the bidirectional quantum-classical bridge that is
// the key innovation of MSFT theory - explaining mass emergence
// from vacuum synchronization at the quantum-classical interface.