#include "MSFTEngine.h"
#include <cstring>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <vector>

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
      _time_step(0),
      _theta_buffer(VK_NULL_HANDLE),
      _theta_out_buffer(VK_NULL_HANDLE),
      _omega_buffer(VK_NULL_HANDLE),
      _R_field_buffer(VK_NULL_HANDLE),
      _gravity_x_buffer(VK_NULL_HANDLE),
      _gravity_y_buffer(VK_NULL_HANDLE),
      _spinor_density_buffer(VK_NULL_HANDLE),
      _theta_memory(VK_NULL_HANDLE),
      _theta_out_memory(VK_NULL_HANDLE),
      _omega_memory(VK_NULL_HANDLE),
      _R_field_memory(VK_NULL_HANDLE),
      _gravity_x_memory(VK_NULL_HANDLE),
      _gravity_y_memory(VK_NULL_HANDLE),
      _spinor_density_memory(VK_NULL_HANDLE),
      _kuramoto_pipeline(VK_NULL_HANDLE),
      _sync_pipeline(VK_NULL_HANDLE),
      _gravity_pipeline(VK_NULL_HANDLE),
      _kuramoto_descriptor_set(VK_NULL_HANDLE),
      _sync_descriptor_set(VK_NULL_HANDLE),
      _gravity_descriptor_set(VK_NULL_HANDLE),
      _kuramoto_descriptor_layout(VK_NULL_HANDLE),
      _sync_descriptor_layout(VK_NULL_HANDLE),
      _gravity_descriptor_layout(VK_NULL_HANDLE),
      _kuramoto_pipeline_layout(VK_NULL_HANDLE),
      _sync_pipeline_layout(VK_NULL_HANDLE),
      _gravity_pipeline_layout(VK_NULL_HANDLE),
      _descriptor_pool(VK_NULL_HANDLE),
      _spinor_buffer(VK_NULL_HANDLE),
      _spinor_memory(VK_NULL_HANDLE),
      _dirac_pipeline(VK_NULL_HANDLE),
      _kuramoto_stochastic_pipeline(VK_NULL_HANDLE),
      _dirac_stochastic_pipeline(VK_NULL_HANDLE),
      _dirac_descriptor_set(VK_NULL_HANDLE),
      _dirac_descriptor_layout(VK_NULL_HANDLE),
      _dirac_pipeline_layout(VK_NULL_HANDLE) {
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
    _gravity_x_data.resize(total_size, 0.0f);
    _gravity_y_data.resize(total_size, 0.0f);

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
    createBuffers();

    // Upload initial data to GPU
    uploadToGPU();

    // Phase 3: Create compute pipelines
    createPipelines();
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

    // Upload to GPU buffer if it exists
    if (_theta_memory != VK_NULL_HANDLE) {
        VkDevice device = _nova->_architect->logical_device;
        size_t gridSizeBytes = sizeof(float) * _Nx * _Ny;
        void* mappedMemory = nullptr;
        vkMapMemory(device, _theta_memory, 0, gridSizeBytes, 0, &mappedMemory);
        memcpy(mappedMemory, _theta_data.data(), gridSizeBytes);
        vkUnmapMemory(device, _theta_memory);
    }
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

    // Upload to GPU buffer if it exists
    if (_omega_memory != VK_NULL_HANDLE) {
        VkDevice device = _nova->_architect->logical_device;
        size_t gridSizeBytes = sizeof(float) * _Nx * _Ny;
        void* mappedMemory = nullptr;
        vkMapMemory(device, _omega_memory, 0, gridSizeBytes, 0, &mappedMemory);
        memcpy(mappedMemory, _omega_data.data(), gridSizeBytes);
        vkUnmapMemory(device, _omega_memory);
    }
}

void MSFTEngine::step(float dt, float K, float damping) {
    /**
     * Phase 4: GPU Compute Dispatch Implementation
     *
     * Executes the MSFT simulation step using GPU compute shaders:
     * 1. kuramoto_step: Evolve phases θ(t) → θ(t+dt)
     * 2. sync_field: Compute synchronization field R(x)
     * 3. gravity_field: Compute gravitational field g(x) = -Δ·∇R(x)
     */

    // Verify pipelines are created
    if (_kuramoto_pipeline == VK_NULL_HANDLE ||
        _sync_pipeline == VK_NULL_HANDLE ||
        _gravity_pipeline == VK_NULL_HANDLE) {
        // Pipelines not ready, fall back to CPU placeholder
        size_t total_size = _Nx * _Ny;
        for (size_t i = 0; i < total_size; ++i) {
            _R_field_data[i] = 0.5f * (1.0f + std::cos(_theta_data[i]));
        }
        return;
    }

    VkDevice device = _nova->_architect->logical_device;

    // 1. Upload current data to GPU
    uploadToGPU();

    // 2. Create command pool for compute operations
    VkCommandPool commandPool;
    VkCommandPoolCreateInfo poolInfo{};
    poolInfo.sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
    poolInfo.queueFamilyIndex = _nova->_architect->queues.indices.compute_family.value_or(
        _nova->_architect->queues.indices.graphics_family.value());  // Use compute or graphics queue
    poolInfo.flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;

    if (vkCreateCommandPool(device, &poolInfo, nullptr, &commandPool) != VK_SUCCESS) {
        // In production, throw exception
        return;
    }

    // 3. Create command buffer
    VkCommandBufferAllocateInfo allocInfo{};
    allocInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;
    allocInfo.commandPool = commandPool;
    allocInfo.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
    allocInfo.commandBufferCount = 1;

    VkCommandBuffer commandBuffer;
    if (vkAllocateCommandBuffers(device, &allocInfo, &commandBuffer) != VK_SUCCESS) {
        vkDestroyCommandPool(device, commandPool, nullptr);
        return;
    }

    // 4. Begin recording commands
    VkCommandBufferBeginInfo beginInfo{};
    beginInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    beginInfo.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;

    VkResult result = vkBeginCommandBuffer(commandBuffer, &beginInfo);
    if (result != VK_SUCCESS) {
        vkFreeCommandBuffers(device, commandPool, 1, &commandBuffer);
        vkDestroyCommandPool(device, commandPool, nullptr);
        return;
    }

    // 5. Prepare push constants
    struct PushConstants {
        float dt;
        float K;
        float damping;
        float Delta;
        float chiral_angle;
        float Nx;
        float Ny;
        float N_total;
        float neighborhood_radius;
    } pushConstants;

    pushConstants.dt = dt;
    pushConstants.K = K;
    pushConstants.damping = damping;
    pushConstants.Delta = _Delta;
    pushConstants.chiral_angle = _chiral_angle;
    pushConstants.Nx = static_cast<float>(_Nx);
    pushConstants.Ny = static_cast<float>(_Ny);
    pushConstants.N_total = static_cast<float>(_Nx * _Ny);
    pushConstants.neighborhood_radius = 1.0f;  // Default to immediate neighbors

    // Calculate workgroup counts (shaders use local_size_x = 16, local_size_y = 16)
    uint32_t workgroupsX = (_Nx + 15) / 16;
    uint32_t workgroupsY = (_Ny + 15) / 16;

    // 6. Dispatch kuramoto_step shader with its specific descriptor set and pipeline layout
    vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, _kuramoto_pipeline);
    vkCmdBindDescriptorSets(commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE,
                           _kuramoto_pipeline_layout, 0, 1, &_kuramoto_descriptor_set, 0, nullptr);
    vkCmdPushConstants(commandBuffer, _kuramoto_pipeline_layout, VK_SHADER_STAGE_COMPUTE_BIT,
                      0, sizeof(pushConstants), &pushConstants);
    vkCmdDispatch(commandBuffer, workgroupsX, workgroupsY, 1);

    // 7. Memory barrier between shaders
    VkMemoryBarrier memBarrier{};
    memBarrier.sType = VK_STRUCTURE_TYPE_MEMORY_BARRIER;
    memBarrier.srcAccessMask = VK_ACCESS_SHADER_WRITE_BIT;
    memBarrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT;

    vkCmdPipelineBarrier(commandBuffer,
                        VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                        VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                        0, 1, &memBarrier, 0, nullptr, 0, nullptr);

    // 8. Copy theta_out back to theta for next iteration
    VkBufferCopy copyRegion{};
    copyRegion.size = sizeof(float) * _Nx * _Ny;
    vkCmdCopyBuffer(commandBuffer, _theta_out_buffer, _theta_buffer, 1, &copyRegion);

    // Another barrier after copy
    vkCmdPipelineBarrier(commandBuffer,
                        VK_PIPELINE_STAGE_TRANSFER_BIT,
                        VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                        0, 1, &memBarrier, 0, nullptr, 0, nullptr);

    // 9. Dispatch sync_field shader with its specific descriptor set and pipeline layout
    vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, _sync_pipeline);
    vkCmdBindDescriptorSets(commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE,
                           _sync_pipeline_layout, 0, 1, &_sync_descriptor_set, 0, nullptr);
    vkCmdPushConstants(commandBuffer, _sync_pipeline_layout, VK_SHADER_STAGE_COMPUTE_BIT,
                      0, sizeof(pushConstants), &pushConstants);
    vkCmdDispatch(commandBuffer, workgroupsX, workgroupsY, 1);

    // 10. Memory barrier
    vkCmdPipelineBarrier(commandBuffer,
                        VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                        VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                        0, 1, &memBarrier, 0, nullptr, 0, nullptr);

    // 11. Dispatch gravity_field shader with its specific descriptor set and pipeline layout
    vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, _gravity_pipeline);
    vkCmdBindDescriptorSets(commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE,
                           _gravity_pipeline_layout, 0, 1, &_gravity_descriptor_set, 0, nullptr);
    vkCmdPushConstants(commandBuffer, _gravity_pipeline_layout, VK_SHADER_STAGE_COMPUTE_BIT,
                      0, sizeof(pushConstants), &pushConstants);
    vkCmdDispatch(commandBuffer, workgroupsX, workgroupsY, 1);

    // 12. End command buffer recording
    result = vkEndCommandBuffer(commandBuffer);
    if (result != VK_SUCCESS) {
        vkFreeCommandBuffers(device, commandPool, 1, &commandBuffer);
        vkDestroyCommandPool(device, commandPool, nullptr);
        return;
    }

    // 13. Submit command buffer to compute queue
    VkSubmitInfo submitInfo{};
    submitInfo.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;
    submitInfo.commandBufferCount = 1;
    submitInfo.pCommandBuffers = &commandBuffer;

    VkQueue computeQueue = _nova->_architect->queues.compute ?
        _nova->_architect->queues.compute : _nova->_architect->queues.graphics;

    // Create fence for synchronization
    VkFenceCreateInfo fenceInfo{};
    fenceInfo.sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO;

    VkFence fence;
    result = vkCreateFence(device, &fenceInfo, nullptr, &fence);
    if (result != VK_SUCCESS) {
        vkFreeCommandBuffers(device, commandPool, 1, &commandBuffer);
        vkDestroyCommandPool(device, commandPool, nullptr);
        return;
    }

    // Submit work to GPU
    result = vkQueueSubmit(computeQueue, 1, &submitInfo, fence);
    if (result != VK_SUCCESS) {
        vkDestroyFence(device, fence, nullptr);
        vkFreeCommandBuffers(device, commandPool, 1, &commandBuffer);
        vkDestroyCommandPool(device, commandPool, nullptr);
        return;
    }

    // 14. Wait for GPU to complete
    vkWaitForFences(device, 1, &fence, VK_TRUE, UINT64_MAX);

    // 15. Download results from GPU
    downloadFromGPU();

    // 16. Cleanup
    vkDestroyFence(device, fence, nullptr);
    vkFreeCommandBuffers(device, commandPool, 1, &commandBuffer);
    vkDestroyCommandPool(device, commandPool, nullptr);
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
     * Return gravitational field g(x,y) = -Δ · ∇R(x,y) from GPU computation.
     *
     * Returns 2D vector field as interleaved (gx, gy) pairs.
     *
     * Phase 4 Implementation (GPU-accelerated):
     * ─────────────────────────────────────────
     * - GPU computes gradients via gravity_field.comp shader
     * - Results stored in _gravity_x_data and _gravity_y_data
     * - This method combines them into interleaved format
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

    std::vector<float> g_field(2 * _Nx * _Ny);  // Interleaved (gx, gy)

    // Check if GPU computation has been performed
    if (_gravity_x_data.empty() || _gravity_y_data.empty()) {
        // Fallback: CPU computation if GPU not ready
        std::vector<float> R_field = getSyncField();

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
    } else {
        // Use GPU-computed gravity field
        for (uint32_t i = 0; i < _Nx * _Ny; i++) {
            g_field[2*i + 0] = _gravity_x_data[i];  // gx component
            g_field[2*i + 1] = _gravity_y_data[i];  // gy component
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
    // Get Vulkan device from Nova
    VkDevice device = _nova->_architect->logical_device;
    VkPhysicalDevice physicalDevice = _nova->_architect->physical_device;

    // Calculate buffer sizes
    VkDeviceSize gridSize = sizeof(float) * _Nx * _Ny;
    VkDeviceSize spinorSize = sizeof(float) * 8 * _Nx * _Ny;  // 4 complex components = 8 floats

    // Helper lambda to find memory type
    auto findMemoryType = [physicalDevice](uint32_t typeFilter, VkMemoryPropertyFlags properties) -> uint32_t {
        VkPhysicalDeviceMemoryProperties memProperties;
        vkGetPhysicalDeviceMemoryProperties(physicalDevice, &memProperties);

        for (uint32_t i = 0; i < memProperties.memoryTypeCount; i++) {
            if ((typeFilter & (1 << i)) &&
                (memProperties.memoryTypes[i].propertyFlags & properties) == properties) {
                return i;
            }
        }

        // Should handle error properly in production
        return 0;
    };

    // Helper lambda to create buffer
    auto createBuffer = [device, findMemoryType](
        VkDeviceSize size,
        VkBufferUsageFlags usage,
        VkMemoryPropertyFlags properties,
        VkBuffer& buffer,
        VkDeviceMemory& bufferMemory) {

        // Create buffer
        VkBufferCreateInfo bufferInfo{};
        bufferInfo.sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;
        bufferInfo.size = size;
        bufferInfo.usage = usage;
        bufferInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;

        if (vkCreateBuffer(device, &bufferInfo, nullptr, &buffer) != VK_SUCCESS) {
            // In production, throw exception
            return;
        }

        // Get memory requirements
        VkMemoryRequirements memRequirements;
        vkGetBufferMemoryRequirements(device, buffer, &memRequirements);

        // Allocate memory
        VkMemoryAllocateInfo allocInfo{};
        allocInfo.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
        allocInfo.allocationSize = memRequirements.size;
        allocInfo.memoryTypeIndex = findMemoryType(memRequirements.memoryTypeBits, properties);

        if (vkAllocateMemory(device, &allocInfo, nullptr, &bufferMemory) != VK_SUCCESS) {
            // In production, throw exception
            vkDestroyBuffer(device, buffer, nullptr);
            buffer = VK_NULL_HANDLE;
            return;
        }

        // Bind memory to buffer
        vkBindBufferMemory(device, buffer, bufferMemory, 0);
    };

    // Create buffers for physics fields
    // Using STORAGE_BUFFER for compute shader access
    // Adding TRANSFER flags for buffer copies
    // Using HOST_VISIBLE | HOST_COHERENT for CPU upload/download
    VkBufferUsageFlags storageUsage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT |
                                      VK_BUFFER_USAGE_TRANSFER_SRC_BIT |
                                      VK_BUFFER_USAGE_TRANSFER_DST_BIT;
    VkMemoryPropertyFlags hostProperties = VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT |
                                          VK_MEMORY_PROPERTY_HOST_COHERENT_BIT;

    // Phase field buffers
    createBuffer(gridSize, storageUsage, hostProperties,
                _theta_buffer, _theta_memory);

    createBuffer(gridSize, storageUsage, hostProperties,
                _theta_out_buffer, _theta_out_memory);

    // Natural frequency buffer
    createBuffer(gridSize, storageUsage, hostProperties,
                _omega_buffer, _omega_memory);

    // Synchronization field buffer
    createBuffer(gridSize, storageUsage, hostProperties,
                _R_field_buffer, _R_field_memory);

    // Gravity field buffers (separate x and y components)
    createBuffer(gridSize, storageUsage, hostProperties,
                _gravity_x_buffer, _gravity_x_memory);

    createBuffer(gridSize, storageUsage, hostProperties,
                _gravity_y_buffer, _gravity_y_memory);

    // Spinor field buffer (4 complex components per grid point)
    createBuffer(spinorSize, storageUsage, hostProperties,
                _spinor_buffer, _spinor_memory);

    // Spinor density buffer (|Ψ|² for quantum-classical feedback)
    createBuffer(gridSize, storageUsage, hostProperties,
                _spinor_density_buffer, _spinor_density_memory);
}

void MSFTEngine::createPipelines() {
    /**
     * Phase 3: Compute Pipeline Creation
     *
     * Creates 3 compute pipelines for MSFT GPU execution with SEPARATE descriptor sets:
     * 1. kuramoto_step.comp - Phase evolution using Kuramoto dynamics
     *    Bindings: theta, theta_out, omega, spinor_density
     * 2. sync_field.comp - Synchronization field R(x) calculation
     *    Bindings: theta, R_field
     * 3. gravity_field.comp - Gravitational field g(x) = -Δ·∇R(x)
     *    Bindings: R_field, gravity_x, gravity_y
     */

    VkDevice device = _nova->_architect->logical_device;

    // Helper lambda to load SPIR-V shader file
    auto loadShaderFile = [](const std::string& path) -> std::vector<uint32_t> {
        std::ifstream file(path, std::ios::ate | std::ios::binary);

        if (!file.is_open()) {
            // In production, throw exception
            return {};
        }

        size_t fileSize = static_cast<size_t>(file.tellg());
        std::vector<uint32_t> buffer(fileSize / sizeof(uint32_t));

        file.seekg(0);
        file.read(reinterpret_cast<char*>(buffer.data()), fileSize);
        file.close();

        return buffer;
    };

    // Helper lambda to create shader module from SPIR-V code
    auto createShaderModule = [device](const std::vector<uint32_t>& code) -> VkShaderModule {
        VkShaderModuleCreateInfo createInfo{};
        createInfo.sType = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO;
        createInfo.codeSize = code.size() * sizeof(uint32_t);
        createInfo.pCode = code.data();

        VkShaderModule shaderModule;
        if (vkCreateShaderModule(device, &createInfo, nullptr, &shaderModule) != VK_SUCCESS) {
            return VK_NULL_HANDLE;
        }

        return shaderModule;
    };

    // 1. Create descriptor set layouts for each shader (each has different bindings)

    // Kuramoto shader layout: theta, theta_out, omega, spinor_density (4 bindings)
    std::vector<VkDescriptorSetLayoutBinding> kuramoto_bindings = {
        // Binding 0: theta_buffer (input phases)
        {
            .binding = 0,
            .descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
            .descriptorCount = 1,
            .stageFlags = VK_SHADER_STAGE_COMPUTE_BIT,
            .pImmutableSamplers = nullptr
        },
        // Binding 1: theta_out_buffer (output phases)
        {
            .binding = 1,
            .descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
            .descriptorCount = 1,
            .stageFlags = VK_SHADER_STAGE_COMPUTE_BIT,
            .pImmutableSamplers = nullptr
        },
        // Binding 2: omega_buffer (natural frequencies)
        {
            .binding = 2,
            .descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
            .descriptorCount = 1,
            .stageFlags = VK_SHADER_STAGE_COMPUTE_BIT,
            .pImmutableSamplers = nullptr
        },
        // Binding 3: spinor_density_buffer (quantum feedback)
        {
            .binding = 3,
            .descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
            .descriptorCount = 1,
            .stageFlags = VK_SHADER_STAGE_COMPUTE_BIT,
            .pImmutableSamplers = nullptr
        }
    };

    VkDescriptorSetLayoutCreateInfo kuramoto_layout_info{};
    kuramoto_layout_info.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
    kuramoto_layout_info.bindingCount = static_cast<uint32_t>(kuramoto_bindings.size());
    kuramoto_layout_info.pBindings = kuramoto_bindings.data();

    if (vkCreateDescriptorSetLayout(device, &kuramoto_layout_info, nullptr, &_kuramoto_descriptor_layout) != VK_SUCCESS) {
        // In production, throw exception
        return;
    }

    // Sync shader layout: theta, R_field (2 bindings)
    std::vector<VkDescriptorSetLayoutBinding> sync_bindings = {
        // Binding 0: theta_buffer (input phases)
        {
            .binding = 0,
            .descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
            .descriptorCount = 1,
            .stageFlags = VK_SHADER_STAGE_COMPUTE_BIT,
            .pImmutableSamplers = nullptr
        },
        // Binding 1: R_field_buffer (output sync field)
        {
            .binding = 1,
            .descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
            .descriptorCount = 1,
            .stageFlags = VK_SHADER_STAGE_COMPUTE_BIT,
            .pImmutableSamplers = nullptr
        }
    };

    VkDescriptorSetLayoutCreateInfo sync_layout_info{};
    sync_layout_info.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
    sync_layout_info.bindingCount = static_cast<uint32_t>(sync_bindings.size());
    sync_layout_info.pBindings = sync_bindings.data();

    if (vkCreateDescriptorSetLayout(device, &sync_layout_info, nullptr, &_sync_descriptor_layout) != VK_SUCCESS) {
        // In production, throw exception
        return;
    }

    // Gravity shader layout: R_field, gravity_x, gravity_y (3 bindings)
    std::vector<VkDescriptorSetLayoutBinding> gravity_bindings = {
        // Binding 0: R_field_buffer (input sync field)
        {
            .binding = 0,
            .descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
            .descriptorCount = 1,
            .stageFlags = VK_SHADER_STAGE_COMPUTE_BIT,
            .pImmutableSamplers = nullptr
        },
        // Binding 1: gravity_x_buffer (output x-component)
        {
            .binding = 1,
            .descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
            .descriptorCount = 1,
            .stageFlags = VK_SHADER_STAGE_COMPUTE_BIT,
            .pImmutableSamplers = nullptr
        },
        // Binding 2: gravity_y_buffer (output y-component)
        {
            .binding = 2,
            .descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
            .descriptorCount = 1,
            .stageFlags = VK_SHADER_STAGE_COMPUTE_BIT,
            .pImmutableSamplers = nullptr
        }
    };

    VkDescriptorSetLayoutCreateInfo gravity_layout_info{};
    gravity_layout_info.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
    gravity_layout_info.bindingCount = static_cast<uint32_t>(gravity_bindings.size());
    gravity_layout_info.pBindings = gravity_bindings.data();

    if (vkCreateDescriptorSetLayout(device, &gravity_layout_info, nullptr, &_gravity_descriptor_layout) != VK_SUCCESS) {
        // In production, throw exception
        return;
    }

    // 2. Define push constants structure for shader parameters
    VkPushConstantRange pushConstantRange{};
    pushConstantRange.stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
    pushConstantRange.offset = 0;
    pushConstantRange.size = sizeof(float) * 9;  // 9 floats: dt, K, damping, Delta, chiral_angle, Nx, Ny, N_total, neighborhood_radius

    // 3. Create pipeline layouts for each shader

    // Kuramoto pipeline layout
    VkPipelineLayoutCreateInfo kuramoto_pipeline_layout_info{};
    kuramoto_pipeline_layout_info.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
    kuramoto_pipeline_layout_info.setLayoutCount = 1;
    kuramoto_pipeline_layout_info.pSetLayouts = &_kuramoto_descriptor_layout;
    kuramoto_pipeline_layout_info.pushConstantRangeCount = 1;
    kuramoto_pipeline_layout_info.pPushConstantRanges = &pushConstantRange;

    if (vkCreatePipelineLayout(device, &kuramoto_pipeline_layout_info, nullptr, &_kuramoto_pipeline_layout) != VK_SUCCESS) {
        // In production, throw exception
        return;
    }

    // Sync pipeline layout
    VkPipelineLayoutCreateInfo sync_pipeline_layout_info{};
    sync_pipeline_layout_info.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
    sync_pipeline_layout_info.setLayoutCount = 1;
    sync_pipeline_layout_info.pSetLayouts = &_sync_descriptor_layout;
    sync_pipeline_layout_info.pushConstantRangeCount = 1;
    sync_pipeline_layout_info.pPushConstantRanges = &pushConstantRange;

    if (vkCreatePipelineLayout(device, &sync_pipeline_layout_info, nullptr, &_sync_pipeline_layout) != VK_SUCCESS) {
        // In production, throw exception
        return;
    }

    // Gravity pipeline layout
    VkPipelineLayoutCreateInfo gravity_pipeline_layout_info{};
    gravity_pipeline_layout_info.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
    gravity_pipeline_layout_info.setLayoutCount = 1;
    gravity_pipeline_layout_info.pSetLayouts = &_gravity_descriptor_layout;
    gravity_pipeline_layout_info.pushConstantRangeCount = 1;
    gravity_pipeline_layout_info.pPushConstantRanges = &pushConstantRange;

    if (vkCreatePipelineLayout(device, &gravity_pipeline_layout_info, nullptr, &_gravity_pipeline_layout) != VK_SUCCESS) {
        // In production, throw exception
        return;
    }

    // 4. Create descriptor pool for all descriptor sets
    // Total buffers needed: kuramoto(4) + sync(2) + gravity(3) = 9
    std::vector<VkDescriptorPoolSize> poolSizes = {
        { VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 9 }  // 9 total storage buffers
    };

    VkDescriptorPoolCreateInfo poolInfo{};
    poolInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
    poolInfo.poolSizeCount = static_cast<uint32_t>(poolSizes.size());
    poolInfo.pPoolSizes = poolSizes.data();
    poolInfo.maxSets = 3;  // 3 descriptor sets (kuramoto, sync, gravity)

    if (vkCreateDescriptorPool(device, &poolInfo, nullptr, &_descriptor_pool) != VK_SUCCESS) {
        // In production, throw exception
        return;
    }

    // 5. Allocate descriptor sets for each pipeline

    // Allocate kuramoto descriptor set
    VkDescriptorSetAllocateInfo kuramoto_alloc_info{};
    kuramoto_alloc_info.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO;
    kuramoto_alloc_info.descriptorPool = _descriptor_pool;
    kuramoto_alloc_info.descriptorSetCount = 1;
    kuramoto_alloc_info.pSetLayouts = &_kuramoto_descriptor_layout;

    if (vkAllocateDescriptorSets(device, &kuramoto_alloc_info, &_kuramoto_descriptor_set) != VK_SUCCESS) {
        // In production, throw exception
        vkDestroyDescriptorPool(device, _descriptor_pool, nullptr);
        _descriptor_pool = VK_NULL_HANDLE;
        return;
    }

    // Allocate sync descriptor set
    VkDescriptorSetAllocateInfo sync_alloc_info{};
    sync_alloc_info.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO;
    sync_alloc_info.descriptorPool = _descriptor_pool;
    sync_alloc_info.descriptorSetCount = 1;
    sync_alloc_info.pSetLayouts = &_sync_descriptor_layout;

    if (vkAllocateDescriptorSets(device, &sync_alloc_info, &_sync_descriptor_set) != VK_SUCCESS) {
        // In production, throw exception
        vkDestroyDescriptorPool(device, _descriptor_pool, nullptr);
        _descriptor_pool = VK_NULL_HANDLE;
        return;
    }

    // Allocate gravity descriptor set
    VkDescriptorSetAllocateInfo gravity_alloc_info{};
    gravity_alloc_info.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO;
    gravity_alloc_info.descriptorPool = _descriptor_pool;
    gravity_alloc_info.descriptorSetCount = 1;
    gravity_alloc_info.pSetLayouts = &_gravity_descriptor_layout;

    if (vkAllocateDescriptorSets(device, &gravity_alloc_info, &_gravity_descriptor_set) != VK_SUCCESS) {
        // In production, throw exception
        vkDestroyDescriptorPool(device, _descriptor_pool, nullptr);
        _descriptor_pool = VK_NULL_HANDLE;
        return;
    }

    // 6. Update descriptor sets with buffer bindings

    // Update kuramoto descriptor set (theta, theta_out, omega, spinor_density)
    {
        std::vector<VkDescriptorBufferInfo> kuramoto_buffer_infos = {
            { _theta_buffer, 0, VK_WHOLE_SIZE },
            { _theta_out_buffer, 0, VK_WHOLE_SIZE },
            { _omega_buffer, 0, VK_WHOLE_SIZE },
            { _spinor_density_buffer, 0, VK_WHOLE_SIZE }
        };

        std::vector<VkWriteDescriptorSet> kuramoto_writes;
        for (size_t i = 0; i < kuramoto_buffer_infos.size(); ++i) {
            VkWriteDescriptorSet write{};
            write.sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
            write.dstSet = _kuramoto_descriptor_set;
            write.dstBinding = static_cast<uint32_t>(i);
            write.dstArrayElement = 0;
            write.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
            write.descriptorCount = 1;
            write.pBufferInfo = &kuramoto_buffer_infos[i];
            kuramoto_writes.push_back(write);
        }

        vkUpdateDescriptorSets(device, static_cast<uint32_t>(kuramoto_writes.size()),
                              kuramoto_writes.data(), 0, nullptr);
    }

    // Update sync descriptor set (theta, R_field)
    {
        std::vector<VkDescriptorBufferInfo> sync_buffer_infos = {
            { _theta_buffer, 0, VK_WHOLE_SIZE },
            { _R_field_buffer, 0, VK_WHOLE_SIZE }
        };

        std::vector<VkWriteDescriptorSet> sync_writes;
        for (size_t i = 0; i < sync_buffer_infos.size(); ++i) {
            VkWriteDescriptorSet write{};
            write.sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
            write.dstSet = _sync_descriptor_set;
            write.dstBinding = static_cast<uint32_t>(i);
            write.dstArrayElement = 0;
            write.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
            write.descriptorCount = 1;
            write.pBufferInfo = &sync_buffer_infos[i];
            sync_writes.push_back(write);
        }

        vkUpdateDescriptorSets(device, static_cast<uint32_t>(sync_writes.size()),
                              sync_writes.data(), 0, nullptr);
    }

    // Update gravity descriptor set (R_field, gravity_x, gravity_y)
    {
        std::vector<VkDescriptorBufferInfo> gravity_buffer_infos = {
            { _R_field_buffer, 0, VK_WHOLE_SIZE },
            { _gravity_x_buffer, 0, VK_WHOLE_SIZE },
            { _gravity_y_buffer, 0, VK_WHOLE_SIZE }
        };

        std::vector<VkWriteDescriptorSet> gravity_writes;
        for (size_t i = 0; i < gravity_buffer_infos.size(); ++i) {
            VkWriteDescriptorSet write{};
            write.sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
            write.dstSet = _gravity_descriptor_set;
            write.dstBinding = static_cast<uint32_t>(i);
            write.dstArrayElement = 0;
            write.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
            write.descriptorCount = 1;
            write.pBufferInfo = &gravity_buffer_infos[i];
            gravity_writes.push_back(write);
        }

        vkUpdateDescriptorSets(device, static_cast<uint32_t>(gravity_writes.size()),
                              gravity_writes.data(), 0, nullptr);
    }

    // 7. Load and create compute pipelines for each shader

    // Helper lambda to create compute pipeline with specific layout
    auto createComputePipeline = [device, &loadShaderFile, &createShaderModule]
        (const std::string& shaderPath, VkPipelineLayout pipelineLayout) -> VkPipeline {
        // Load SPIR-V shader
        auto shaderCode = loadShaderFile(shaderPath);
        if (shaderCode.empty()) {
            return VK_NULL_HANDLE;
        }

        // Create shader module
        VkShaderModule shaderModule = createShaderModule(shaderCode);
        if (shaderModule == VK_NULL_HANDLE) {
            return VK_NULL_HANDLE;
        }

        // Create compute pipeline
        VkPipelineShaderStageCreateInfo shaderStageInfo{};
        shaderStageInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
        shaderStageInfo.stage = VK_SHADER_STAGE_COMPUTE_BIT;
        shaderStageInfo.module = shaderModule;
        shaderStageInfo.pName = "main";

        VkComputePipelineCreateInfo pipelineInfo{};
        pipelineInfo.sType = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO;
        pipelineInfo.stage = shaderStageInfo;
        pipelineInfo.layout = pipelineLayout;  // Use specific layout for this pipeline

        VkPipeline pipeline;
        if (vkCreateComputePipelines(device, VK_NULL_HANDLE, 1, &pipelineInfo,
                                     nullptr, &pipeline) != VK_SUCCESS) {
            vkDestroyShaderModule(device, shaderModule, nullptr);
            return VK_NULL_HANDLE;
        }

        // Cleanup shader module (pipeline keeps reference internally)
        vkDestroyShaderModule(device, shaderModule, nullptr);

        return pipeline;
    };

    // Create the three compute pipelines with their specific layouts
    _kuramoto_pipeline = createComputePipeline("/home/persist/neotec/0rigin/shaders/smft/kuramoto_step.comp.spv",
                                               _kuramoto_pipeline_layout);
    _sync_pipeline = createComputePipeline("/home/persist/neotec/0rigin/shaders/smft/sync_field.comp.spv",
                                           _sync_pipeline_layout);
    _gravity_pipeline = createComputePipeline("/home/persist/neotec/0rigin/shaders/smft/gravity_field.comp.spv",
                                              _gravity_pipeline_layout);

    // Note: We keep descriptor pool alive for lifetime of engine (cleanup in destructor)
    // Store it as member variable if needed for cleanup
}

void MSFTEngine::uploadToGPU() {
    // Upload CPU-side data to GPU buffers
    VkDevice device = _nova->_architect->logical_device;

    size_t gridSizeBytes = sizeof(float) * _Nx * _Ny;

    // Upload theta data
    if (_theta_memory != VK_NULL_HANDLE && !_theta_data.empty()) {
        void* mappedMemory = nullptr;
        vkMapMemory(device, _theta_memory, 0, gridSizeBytes, 0, &mappedMemory);
        memcpy(mappedMemory, _theta_data.data(), gridSizeBytes);
        vkUnmapMemory(device, _theta_memory);
    }

    // Upload omega data
    if (_omega_memory != VK_NULL_HANDLE && !_omega_data.empty()) {
        void* mappedMemory = nullptr;
        vkMapMemory(device, _omega_memory, 0, gridSizeBytes, 0, &mappedMemory);
        memcpy(mappedMemory, _omega_data.data(), gridSizeBytes);
        vkUnmapMemory(device, _omega_memory);
    }

    // Initialize R_field buffer to zero (will be computed by shader)
    if (_R_field_memory != VK_NULL_HANDLE) {
        void* mappedMemory = nullptr;
        vkMapMemory(device, _R_field_memory, 0, gridSizeBytes, 0, &mappedMemory);
        memset(mappedMemory, 0, gridSizeBytes);
        vkUnmapMemory(device, _R_field_memory);
    }

    // Initialize gravity buffers to zero
    if (_gravity_x_memory != VK_NULL_HANDLE) {
        void* mappedMemory = nullptr;
        vkMapMemory(device, _gravity_x_memory, 0, gridSizeBytes, 0, &mappedMemory);
        memset(mappedMemory, 0, gridSizeBytes);
        vkUnmapMemory(device, _gravity_x_memory);
    }

    if (_gravity_y_memory != VK_NULL_HANDLE) {
        void* mappedMemory = nullptr;
        vkMapMemory(device, _gravity_y_memory, 0, gridSizeBytes, 0, &mappedMemory);
        memset(mappedMemory, 0, gridSizeBytes);
        vkUnmapMemory(device, _gravity_y_memory);
    }

    // Upload spinor field data (4 complex components = 8 floats per grid point)
    if (_spinor_memory != VK_NULL_HANDLE && !_spinor_field.empty()) {
        size_t spinorSizeBytes = sizeof(std::complex<float>) * _spinor_field.size();
        void* mappedMemory = nullptr;
        vkMapMemory(device, _spinor_memory, 0, spinorSizeBytes, 0, &mappedMemory);
        memcpy(mappedMemory, _spinor_field.data(), spinorSizeBytes);
        vkUnmapMemory(device, _spinor_memory);
    }

    // Initialize spinor density buffer to zero (will be computed from spinor field)
    if (_spinor_density_memory != VK_NULL_HANDLE) {
        void* mappedMemory = nullptr;
        vkMapMemory(device, _spinor_density_memory, 0, gridSizeBytes, 0, &mappedMemory);
        memset(mappedMemory, 0, gridSizeBytes);
        vkUnmapMemory(device, _spinor_density_memory);
    }
}

void MSFTEngine::downloadFromGPU() {
    // Download GPU results back to CPU-side arrays
    VkDevice device = _nova->_architect->logical_device;

    size_t gridSizeBytes = sizeof(float) * _Nx * _Ny;

    // Download updated theta field (from theta_out buffer)
    if (_theta_out_memory != VK_NULL_HANDLE) {
        void* mappedMemory = nullptr;
        vkMapMemory(device, _theta_out_memory, 0, gridSizeBytes, 0, &mappedMemory);
        memcpy(_theta_data.data(), mappedMemory, gridSizeBytes);
        vkUnmapMemory(device, _theta_out_memory);
    }

    // Download synchronization field
    if (_R_field_memory != VK_NULL_HANDLE) {
        void* mappedMemory = nullptr;
        vkMapMemory(device, _R_field_memory, 0, gridSizeBytes, 0, &mappedMemory);
        memcpy(_R_field_data.data(), mappedMemory, gridSizeBytes);
        vkUnmapMemory(device, _R_field_memory);
    }

    // Download gravity field components
    if (_gravity_x_memory != VK_NULL_HANDLE && !_gravity_x_data.empty()) {
        void* mappedMemory = nullptr;
        vkMapMemory(device, _gravity_x_memory, 0, gridSizeBytes, 0, &mappedMemory);
        memcpy(_gravity_x_data.data(), mappedMemory, gridSizeBytes);
        vkUnmapMemory(device, _gravity_x_memory);
    }

    if (_gravity_y_memory != VK_NULL_HANDLE && !_gravity_y_data.empty()) {
        void* mappedMemory = nullptr;
        vkMapMemory(device, _gravity_y_memory, 0, gridSizeBytes, 0, &mappedMemory);
        memcpy(_gravity_y_data.data(), mappedMemory, gridSizeBytes);
        vkUnmapMemory(device, _gravity_y_memory);
    }

    // Download updated spinor field
    if (_spinor_memory != VK_NULL_HANDLE && !_spinor_field.empty()) {
        size_t spinorSizeBytes = sizeof(std::complex<float>) * _spinor_field.size();
        void* mappedMemory = nullptr;
        vkMapMemory(device, _spinor_memory, 0, spinorSizeBytes, 0, &mappedMemory);
        memcpy(_spinor_field.data(), mappedMemory, spinorSizeBytes);
        vkUnmapMemory(device, _spinor_memory);
    }
}

void MSFTEngine::destroyResources() {
    // Cleanup Vulkan resources in reverse order of creation
    if (_nova && _nova->_architect) {
        VkDevice device = _nova->_architect->logical_device;

        // Destroy pipelines (Phase 3 - will be created later)
        if (_dirac_pipeline != VK_NULL_HANDLE) {
            vkDestroyPipeline(device, _dirac_pipeline, nullptr);
            _dirac_pipeline = VK_NULL_HANDLE;
        }
        if (_gravity_pipeline != VK_NULL_HANDLE) {
            vkDestroyPipeline(device, _gravity_pipeline, nullptr);
            _gravity_pipeline = VK_NULL_HANDLE;
        }
        if (_kuramoto_pipeline != VK_NULL_HANDLE) {
            vkDestroyPipeline(device, _kuramoto_pipeline, nullptr);
            _kuramoto_pipeline = VK_NULL_HANDLE;
        }
        if (_sync_pipeline != VK_NULL_HANDLE) {
            vkDestroyPipeline(device, _sync_pipeline, nullptr);
            _sync_pipeline = VK_NULL_HANDLE;
        }

        // Destroy pipeline layouts
        if (_kuramoto_pipeline_layout != VK_NULL_HANDLE) {
            vkDestroyPipelineLayout(device, _kuramoto_pipeline_layout, nullptr);
            _kuramoto_pipeline_layout = VK_NULL_HANDLE;
        }
        if (_sync_pipeline_layout != VK_NULL_HANDLE) {
            vkDestroyPipelineLayout(device, _sync_pipeline_layout, nullptr);
            _sync_pipeline_layout = VK_NULL_HANDLE;
        }
        if (_gravity_pipeline_layout != VK_NULL_HANDLE) {
            vkDestroyPipelineLayout(device, _gravity_pipeline_layout, nullptr);
            _gravity_pipeline_layout = VK_NULL_HANDLE;
        }

        // Destroy descriptor set layouts
        if (_kuramoto_descriptor_layout != VK_NULL_HANDLE) {
            vkDestroyDescriptorSetLayout(device, _kuramoto_descriptor_layout, nullptr);
            _kuramoto_descriptor_layout = VK_NULL_HANDLE;
        }
        if (_sync_descriptor_layout != VK_NULL_HANDLE) {
            vkDestroyDescriptorSetLayout(device, _sync_descriptor_layout, nullptr);
            _sync_descriptor_layout = VK_NULL_HANDLE;
        }
        if (_gravity_descriptor_layout != VK_NULL_HANDLE) {
            vkDestroyDescriptorSetLayout(device, _gravity_descriptor_layout, nullptr);
            _gravity_descriptor_layout = VK_NULL_HANDLE;
        }
        if (_descriptor_pool != VK_NULL_HANDLE) {
            vkDestroyDescriptorPool(device, _descriptor_pool, nullptr);
            _descriptor_pool = VK_NULL_HANDLE;
        }

        // Destroy buffers
        if (_spinor_buffer != VK_NULL_HANDLE) {
            vkDestroyBuffer(device, _spinor_buffer, nullptr);
            _spinor_buffer = VK_NULL_HANDLE;
        }
        if (_theta_buffer != VK_NULL_HANDLE) {
            vkDestroyBuffer(device, _theta_buffer, nullptr);
            _theta_buffer = VK_NULL_HANDLE;
        }
        if (_theta_out_buffer != VK_NULL_HANDLE) {
            vkDestroyBuffer(device, _theta_out_buffer, nullptr);
            _theta_out_buffer = VK_NULL_HANDLE;
        }
        if (_omega_buffer != VK_NULL_HANDLE) {
            vkDestroyBuffer(device, _omega_buffer, nullptr);
            _omega_buffer = VK_NULL_HANDLE;
        }
        if (_R_field_buffer != VK_NULL_HANDLE) {
            vkDestroyBuffer(device, _R_field_buffer, nullptr);
            _R_field_buffer = VK_NULL_HANDLE;
        }
        if (_gravity_x_buffer != VK_NULL_HANDLE) {
            vkDestroyBuffer(device, _gravity_x_buffer, nullptr);
            _gravity_x_buffer = VK_NULL_HANDLE;
        }
        if (_gravity_y_buffer != VK_NULL_HANDLE) {
            vkDestroyBuffer(device, _gravity_y_buffer, nullptr);
            _gravity_y_buffer = VK_NULL_HANDLE;
        }
        if (_spinor_density_buffer != VK_NULL_HANDLE) {
            vkDestroyBuffer(device, _spinor_density_buffer, nullptr);
            _spinor_density_buffer = VK_NULL_HANDLE;
        }

        // Free memory
        if (_spinor_memory != VK_NULL_HANDLE) {
            vkFreeMemory(device, _spinor_memory, nullptr);
            _spinor_memory = VK_NULL_HANDLE;
        }
        if (_theta_memory != VK_NULL_HANDLE) {
            vkFreeMemory(device, _theta_memory, nullptr);
            _theta_memory = VK_NULL_HANDLE;
        }
        if (_theta_out_memory != VK_NULL_HANDLE) {
            vkFreeMemory(device, _theta_out_memory, nullptr);
            _theta_out_memory = VK_NULL_HANDLE;
        }
        if (_omega_memory != VK_NULL_HANDLE) {
            vkFreeMemory(device, _omega_memory, nullptr);
            _omega_memory = VK_NULL_HANDLE;
        }
        if (_R_field_memory != VK_NULL_HANDLE) {
            vkFreeMemory(device, _R_field_memory, nullptr);
            _R_field_memory = VK_NULL_HANDLE;
        }
        if (_gravity_x_memory != VK_NULL_HANDLE) {
            vkFreeMemory(device, _gravity_x_memory, nullptr);
            _gravity_x_memory = VK_NULL_HANDLE;
        }
        if (_gravity_y_memory != VK_NULL_HANDLE) {
            vkFreeMemory(device, _gravity_y_memory, nullptr);
            _gravity_y_memory = VK_NULL_HANDLE;
        }
        if (_spinor_density_memory != VK_NULL_HANDLE) {
            vkFreeMemory(device, _spinor_density_memory, nullptr);
            _spinor_density_memory = VK_NULL_HANDLE;
        }
    }
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

void MSFTEngine::stepStochastic(float dt, float K, float damping,
                                float sigma_theta, float sigma_psi) {
    /**
     * Stochastic MSFT Evolution with MSR Noise Formalism
     *
     * Implements coupled stochastic dynamics:
     * 1. Kuramoto with noise: dθ = (ω + K·Σsin(θⱼ-θᵢ))dt + σ_θ·dW
     * 2. Dirac with noise: i∂_tΨ = H_DiracΨ + σ_Ψ·ξ(t)
     *
     * Pipeline sequence:
     * 1. kuramoto_stochastic: Evolve phases with noise
     * 2. sync_field: Compute R(x) order parameter
     * 3. gravity_field: Compute gravitational corrections
     * 4. dirac_stochastic: Evolve spinor with noise
     */

    // Increment timestep for PRNG seeding
    _time_step++;

    // Verify stochastic pipelines are created
    if (_kuramoto_stochastic_pipeline == VK_NULL_HANDLE ||
        _dirac_stochastic_pipeline == VK_NULL_HANDLE) {
        // Fall back to deterministic step if pipelines not ready
        step(dt, K, damping);
        return;
    }

    VkDevice device = _nova->_architect->logical_device;

    // 1. Upload current data to GPU
    uploadToGPU();

    // 2. Create command pool for compute operations
    VkCommandPool commandPool;
    VkCommandPoolCreateInfo poolInfo{};
    poolInfo.sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
    poolInfo.queueFamilyIndex = _nova->_architect->queues.indices.compute_family.value_or(
        _nova->_architect->queues.indices.graphics_family.value());
    poolInfo.flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;

    if (vkCreateCommandPool(device, &poolInfo, nullptr, &commandPool) != VK_SUCCESS) {
        return;
    }

    // 3. Create command buffer
    VkCommandBufferAllocateInfo allocInfo{};
    allocInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;
    allocInfo.commandPool = commandPool;
    allocInfo.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
    allocInfo.commandBufferCount = 1;

    VkCommandBuffer commandBuffer;
    if (vkAllocateCommandBuffers(device, &allocInfo, &commandBuffer) != VK_SUCCESS) {
        vkDestroyCommandPool(device, commandPool, nullptr);
        return;
    }

    // 4. Begin command buffer recording
    VkCommandBufferBeginInfo beginInfo{};
    beginInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    beginInfo.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;

    if (vkBeginCommandBuffer(commandBuffer, &beginInfo) != VK_SUCCESS) {
        vkFreeCommandBuffers(device, commandPool, 1, &commandBuffer);
        vkDestroyCommandPool(device, commandPool, nullptr);
        return;
    }

    // 5. Dispatch stochastic Kuramoto shader
    vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, _kuramoto_stochastic_pipeline);
    vkCmdBindDescriptorSets(commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE,
                           _kuramoto_pipeline_layout, 0, 1, &_kuramoto_descriptor_set, 0, nullptr);

    // Push constants for stochastic Kuramoto
    struct {
        float dt;
        float K;
        float sigma;       // sigma_theta
        float damping;
        float omega_mean;
        uint32_t Nx;
        uint32_t Ny;
        uint32_t time_step;
    } kuramoto_push = {
        dt, K, sigma_theta, damping, 0.0f, _Nx, _Ny, _time_step
    };

    vkCmdPushConstants(commandBuffer, _kuramoto_pipeline_layout, VK_SHADER_STAGE_COMPUTE_BIT,
                       0, sizeof(kuramoto_push), &kuramoto_push);

    // Dispatch with 16x16 workgroups
    uint32_t workgroup_x = (_Nx + 15) / 16;
    uint32_t workgroup_y = (_Ny + 15) / 16;
    vkCmdDispatch(commandBuffer, workgroup_x, workgroup_y, 1);

    // Memory barrier before sync field
    VkMemoryBarrier barrier{};
    barrier.sType = VK_STRUCTURE_TYPE_MEMORY_BARRIER;
    barrier.srcAccessMask = VK_ACCESS_SHADER_WRITE_BIT;
    barrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT;

    vkCmdPipelineBarrier(commandBuffer,
                        VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                        VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                        0, 1, &barrier, 0, nullptr, 0, nullptr);

    // Swap theta buffers
    std::swap(_theta_buffer, _theta_out_buffer);
    std::swap(_theta_memory, _theta_out_memory);

    // 6. Dispatch sync field shader
    vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, _sync_pipeline);
    vkCmdBindDescriptorSets(commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE,
                           _sync_pipeline_layout, 0, 1, &_sync_descriptor_set, 0, nullptr);

    // Push constants for sync field
    struct {
        float dt;
        float K;
        float damping;
        float Delta;
        float chiral_angle;
        uint32_t Nx;
        uint32_t Ny;
        uint32_t N_total;
        uint32_t neighborhood_radius;
    } sync_push = {
        dt, K, damping, _Delta, _chiral_angle, _Nx, _Ny, _Nx * _Ny, 5
    };

    vkCmdPushConstants(commandBuffer, _sync_pipeline_layout, VK_SHADER_STAGE_COMPUTE_BIT,
                       0, sizeof(sync_push), &sync_push);

    vkCmdDispatch(commandBuffer, workgroup_x, workgroup_y, 1);

    // Barrier before gravity field
    vkCmdPipelineBarrier(commandBuffer,
                        VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                        VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                        0, 1, &barrier, 0, nullptr, 0, nullptr);

    // 7. Dispatch gravity field shader
    vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, _gravity_pipeline);
    vkCmdBindDescriptorSets(commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE,
                           _gravity_pipeline_layout, 0, 1, &_gravity_descriptor_set, 0, nullptr);

    vkCmdPushConstants(commandBuffer, _gravity_pipeline_layout, VK_SHADER_STAGE_COMPUTE_BIT,
                       0, sizeof(sync_push), &sync_push);

    vkCmdDispatch(commandBuffer, workgroup_x, workgroup_y, 1);

    // Barrier before Dirac stochastic
    vkCmdPipelineBarrier(commandBuffer,
                        VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                        VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                        0, 1, &barrier, 0, nullptr, 0, nullptr);

    // 8. Dispatch stochastic Dirac shader
    if (_dirac_stochastic_pipeline != VK_NULL_HANDLE && _dirac_descriptor_set != VK_NULL_HANDLE) {
        vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, _dirac_stochastic_pipeline);
        vkCmdBindDescriptorSets(commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE,
                               _dirac_pipeline_layout, 0, 1, &_dirac_descriptor_set, 0, nullptr);

        // Push constants for stochastic Dirac
        struct {
            float dt;
            float K;
            float sigma_psi;
            float damping;
            float Delta;
            uint32_t Nx;
            uint32_t Ny;
            uint32_t time_step;
        } dirac_push = {
            dt, K, sigma_psi, damping, _Delta, _Nx, _Ny, _time_step
        };

        vkCmdPushConstants(commandBuffer, _dirac_pipeline_layout, VK_SHADER_STAGE_COMPUTE_BIT,
                          0, sizeof(dirac_push), &dirac_push);

        vkCmdDispatch(commandBuffer, workgroup_x, workgroup_y, 1);
    }

    // 9. End command buffer recording
    if (vkEndCommandBuffer(commandBuffer) != VK_SUCCESS) {
        vkFreeCommandBuffers(device, commandPool, 1, &commandBuffer);
        vkDestroyCommandPool(device, commandPool, nullptr);
        return;
    }

    // 10. Submit command buffer to compute queue
    VkSubmitInfo submitInfo{};
    submitInfo.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;
    submitInfo.commandBufferCount = 1;
    submitInfo.pCommandBuffers = &commandBuffer;

    VkQueue computeQueue = _nova->_architect->queues.compute ?
        _nova->_architect->queues.compute : _nova->_architect->queues.graphics;

    vkQueueSubmit(computeQueue, 1, &submitInfo, VK_NULL_HANDLE);
    vkQueueWaitIdle(computeQueue);

    // 11. Download results from GPU
    downloadFromGPU();

    // 12. Cleanup
    vkFreeCommandBuffers(device, commandPool, 1, &commandBuffer);
    vkDestroyCommandPool(device, commandPool, nullptr);
}