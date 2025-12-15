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
    constexpr float HBAR = 1.054571817e-34f;          // Reduced Planck constant (J·s)
    constexpr float C = 299792458.0f;                 // Speed of light (m/s)
    constexpr float ELECTRON_MASS = 9.10938356e-31f;  // Electron rest mass (kg)
    constexpr double PLANCK_FREQUENCY = 1.85492e43;   // Planck frequency (Hz) = sqrt(c^5/(ℏG))

    // Compute vacuum potential (Heisenberg-Einstein unification)
    // Δ = ℏω_max/c² - represents the maximum mass capacity of vacuum fluctuations
    inline float computeVacuumPotential(float omega_max) {
        return HBAR * omega_max / (C * C);
    }

    // Default vacuum potential at Planck scale
    inline float planckVacuumPotential() {
        return computeVacuumPotential(PLANCK_FREQUENCY);
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
    // MSFT Theory: Mass Synchronization Field Theory
    // Formula: (iγ^μ∂_μ)Ψ(x) = Δ · R(x) · e^(iθ(x)γ^5) Ψ(x)
    //
    // Components:
    // - Δ (Delta): Vacuum potential limit = ℏω_max/c² (Heisenberg-Einstein unification)
    //              Represents the maximum mass capacity of vacuum fluctuations
    // - R(x): Kuramoto synchronization field = |⟨e^(iθ)⟩| (local order parameter)
    //         Measures coherence of phase oscillators in local neighborhood
    // - θ(x): Phase field from oscillator array (Kuramoto dynamics)
    //         Evolves via dθ/dt = ω + K·Σsin(θⱼ - θᵢ)
    // - e^(iθγ^5): Euler chiral phase rotation (determines spin/handedness)
    //              Creates left/right handed fermions through γ^5 chirality operator
    // - Ψ(x): 4-component Dirac spinor (quantum wavefunction)
    //         Describes relativistic fermion with spin-1/2
    //
    // Mass emerges as: m(x) = Δ · R(x)
    // When R→0 (chaos): m→0 (massless, photon-like)
    // When R→1 (sync): m→Δ (mass emerges from vacuum structure)
    //
    // This bridges quantum (Dirac) and classical (Kuramoto) physics,
    // explaining mass as emergent from vacuum synchronization

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
    // Compute mass field: m(x,y) = Δ · R(x,y)
    // This implements the core MSFT equation for mass emergence
    auto R = getSyncField();
    std::vector<float> mass(R.size());

    // Basic mass from synchronization: m = Δ · R
    // When R=0 (chaos): mass=0 (particle is massless like photon)
    // When R=1 (sync): mass=Δ (particle gains full vacuum potential)
    for (size_t i = 0; i < R.size(); i++) {
        mass[i] = _Delta * R[i];
    }

    // TODO Phase 4: Apply chiral modulation for parity violation
    // The full mass operator includes chiral rotation: m · e^(iθ·γ^5)
    // For scalar mass magnitude: m_scalar = Δ · R · cos(θ_chiral)
    // For vector form: m_vector = Δ · R · e^(iθ·γ^5) (requires spinor context)
    // This enables left/right handedness in weak interactions

    return mass;
}

std::vector<float> MSFTEngine::getPhaseField() const {
    // Return the current phase field θ(x,y)
    return _theta_data;
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