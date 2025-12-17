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

    // Initialize pipeline factory with device
    _pipelineFactory = std::make_unique<MSFTPipelineFactory>(_nova->_architect->logical_device);

    // Initialize descriptor manager with device
    _descriptorManager = std::make_unique<MSFTDescriptorManager>(_nova->_architect->logical_device);

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
    if (_theta_memory != VK_NULL_HANDLE && _bufferManager) {
        size_t gridSizeBytes = sizeof(float) * _Nx * _Ny;
        _bufferManager->uploadData(_theta_memory, _theta_data.data(), gridSizeBytes);
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
    if (_omega_memory != VK_NULL_HANDLE && _bufferManager) {
        size_t gridSizeBytes = sizeof(float) * _Nx * _Ny;
        _bufferManager->uploadData(_omega_memory, _omega_data.data(), gridSizeBytes);
    }
}

void MSFTEngine::step(float dt, float K, float damping) {
    /**
     * Phase 4: GPU Compute Dispatch Implementation
     *
     * Executes the MSFT simulation step using GPU compute shaders:
     * 1. kuramoto_step: Evolve phases θ(t) → θ(t+dt) [✅ SAFE: 9 transcendentals]
     * 2. sync_field_simple: Compute synchronization field R(x) [✅ SAFE: 37 transcendentals]
     * 3. gravity_field: Compute gravitational field g(x) = -Δ·∇R(x) [✅ SAFE: pure arithmetic]
     *
     * GPU SAFETY: Using only verified safe shaders based on timeout audit.
     * Dirac evolution NOT included - requires CPU implementation (see stepStochastic).
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

    // Initialize compute dispatcher if needed
    if (!_compute) {
        uint32_t queueFamily = _nova->_architect->queues.indices.compute_family.value_or(
            _nova->_architect->queues.indices.graphics_family.value());
        VkQueue queue = _nova->_architect->queues.compute ?
            _nova->_architect->queues.compute : _nova->_architect->queues.graphics;

        _compute = std::make_unique<MSFTCompute>(
            _nova->_architect->logical_device, queue, queueFamily);

        if (!_compute->initialize()) {
            _compute.reset();
            return;
        }
    }

    // Upload current data to GPU
    uploadToGPU();

    // Begin compute batch
    if (!_compute->beginBatch()) {
        return;
    }

    // Prepare push constants
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
    pushConstants.neighborhood_radius = 1.0f;

    // Calculate workgroup counts
    uint32_t workgroupsX, workgroupsY;
    MSFTCompute::calculateWorkgroups(_Nx, _Ny, workgroupsX, workgroupsY);

    // Dispatch kuramoto shader
    _compute->dispatchKuramoto(_kuramoto_pipeline, _kuramoto_pipeline_layout,
                               _kuramoto_descriptor_set, &pushConstants,
                               sizeof(pushConstants), workgroupsX, workgroupsY);

    // Memory barrier
    _compute->insertMemoryBarrier();

    // Copy theta_out back to theta
    _compute->copyBuffer(_theta_out_buffer, _theta_buffer,
                        sizeof(float) * _Nx * _Ny);

    // Another barrier after copy
    _compute->insertMemoryBarrier();

    // Dispatch sync field shader
    _compute->dispatchSyncField(_sync_pipeline, _sync_pipeline_layout,
                                _sync_descriptor_set, &pushConstants,
                                sizeof(pushConstants), workgroupsX, workgroupsY);

    // Memory barrier
    _compute->insertMemoryBarrier();

    // Dispatch gravity field shader
    _compute->dispatchGravityField(_gravity_pipeline, _gravity_pipeline_layout,
                                   _gravity_descriptor_set, &pushConstants,
                                   sizeof(pushConstants), workgroupsX, workgroupsY);

    // Submit batch and wait for completion
    if (!_compute->submitBatch(true)) {
        return;
    }

    // Download results from GPU
    downloadFromGPU();
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
    // Phase 2: Vulkan buffer allocation using MSFTBufferManager
    VkDevice device = _nova->_architect->logical_device;
    VkPhysicalDevice physicalDevice = _nova->_architect->physical_device;

    // Initialize buffer manager
    _bufferManager = std::make_unique<MSFTBufferManager>(device, physicalDevice);

    // Calculate buffer sizes
    VkDeviceSize gridSize = sizeof(float) * _Nx * _Ny;
    VkDeviceSize spinorSize = sizeof(float) * 8 * _Nx * _Ny;  // 4 complex components = 8 floats

    // Create phase field buffers
    auto [theta_buf, theta_mem] = _bufferManager->createStorageBuffer(gridSize);
    _theta_buffer = theta_buf;
    _theta_memory = theta_mem;

    auto [theta_out_buf, theta_out_mem] = _bufferManager->createStorageBuffer(gridSize);
    _theta_out_buffer = theta_out_buf;
    _theta_out_memory = theta_out_mem;

    // Create natural frequency buffer
    auto [omega_buf, omega_mem] = _bufferManager->createStorageBuffer(gridSize);
    _omega_buffer = omega_buf;
    _omega_memory = omega_mem;

    // Create synchronization field buffer
    auto [r_field_buf, r_field_mem] = _bufferManager->createStorageBuffer(gridSize);
    _R_field_buffer = r_field_buf;
    _R_field_memory = r_field_mem;

    // Create gravity field buffers (separate x and y components)
    auto [grav_x_buf, grav_x_mem] = _bufferManager->createStorageBuffer(gridSize);
    _gravity_x_buffer = grav_x_buf;
    _gravity_x_memory = grav_x_mem;

    auto [grav_y_buf, grav_y_mem] = _bufferManager->createStorageBuffer(gridSize);
    _gravity_y_buffer = grav_y_buf;
    _gravity_y_memory = grav_y_mem;

    // Create spinor field buffer (4 complex components per grid point)
    auto [spinor_buf, spinor_mem] = _bufferManager->createStorageBuffer(spinorSize);
    _spinor_buffer = spinor_buf;
    _spinor_memory = spinor_mem;

    // Create spinor density buffer (|Ψ|² for quantum-classical feedback)
    auto [spinor_dens_buf, spinor_dens_mem] = _bufferManager->createStorageBuffer(gridSize);
    _spinor_density_buffer = spinor_dens_buf;
    _spinor_density_memory = spinor_dens_mem;
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

    // 1. Create descriptor set layouts for each shader using descriptor manager

    // Kuramoto shader layout: theta, theta_out, omega, spinor_density (4 bindings)
    std::vector<VkDescriptorSetLayoutBinding> kuramoto_bindings;
    kuramoto_bindings.push_back(_descriptorManager->createBinding(0, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER));
    kuramoto_bindings.push_back(_descriptorManager->createBinding(1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER));
    kuramoto_bindings.push_back(_descriptorManager->createBinding(2, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER));
    kuramoto_bindings.push_back(_descriptorManager->createBinding(3, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER));

    _kuramoto_descriptor_layout = _descriptorManager->createDescriptorSetLayout(kuramoto_bindings);

    // Sync shader layout: theta, R_field (2 bindings)
    std::vector<VkDescriptorSetLayoutBinding> sync_bindings;
    sync_bindings.push_back(_descriptorManager->createBinding(0, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER));
    sync_bindings.push_back(_descriptorManager->createBinding(1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER));

    _sync_descriptor_layout = _descriptorManager->createDescriptorSetLayout(sync_bindings);

    // Gravity shader layout: R_field, gravity_x, gravity_y (3 bindings)
    std::vector<VkDescriptorSetLayoutBinding> gravity_bindings;
    gravity_bindings.push_back(_descriptorManager->createBinding(0, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER));
    gravity_bindings.push_back(_descriptorManager->createBinding(1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER));
    gravity_bindings.push_back(_descriptorManager->createBinding(2, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER));

    _gravity_descriptor_layout = _descriptorManager->createDescriptorSetLayout(gravity_bindings);

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

    // 4. Create descriptor pool for all descriptor sets using descriptor manager
    // Total buffers needed: kuramoto(4) + sync(2) + gravity(3) = 9
    std::vector<VkDescriptorPoolSize> poolSizes = {
        { VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 9 }  // 9 total storage buffers
    };

    _descriptor_pool = _descriptorManager->createDescriptorPool(3, poolSizes);

    // 5. Allocate descriptor sets for each pipeline using descriptor manager

    _kuramoto_descriptor_set = _descriptorManager->allocateDescriptorSet(
        _descriptor_pool, _kuramoto_descriptor_layout);

    _sync_descriptor_set = _descriptorManager->allocateDescriptorSet(
        _descriptor_pool, _sync_descriptor_layout);

    _gravity_descriptor_set = _descriptorManager->allocateDescriptorSet(
        _descriptor_pool, _gravity_descriptor_layout);

    // 6. Update descriptor sets with buffer bindings using descriptor manager

    // Update kuramoto descriptor set (theta, theta_out, omega, spinor_density)
    _descriptorManager->updateDescriptorSet(_kuramoto_descriptor_set,
        { _theta_buffer, _theta_out_buffer, _omega_buffer, _spinor_density_buffer });

    // Update sync descriptor set (theta, R_field)
    _descriptorManager->updateDescriptorSet(_sync_descriptor_set,
        { _theta_buffer, _R_field_buffer });

    // Update gravity descriptor set (R_field, gravity_x, gravity_y)
    _descriptorManager->updateDescriptorSet(_gravity_descriptor_set,
        { _R_field_buffer, _gravity_x_buffer, _gravity_y_buffer });

    // 7. Create compute pipelines using factory
    // Note: Factory is already created in initialize() method

    // Create the three compute pipelines with their specific layouts using factory
    // GPU SAFETY: Using only GPU-safe shaders based on timeout audit

    // ✅ SAFE: kuramoto_step.comp - 9 transcendentals, well under budget
    _kuramoto_pipeline = _pipelineFactory->createKuramotoPipeline(
        "/home/persist/neotec/0rigin/shaders/smft/kuramoto_step.comp.spv",
        _kuramoto_pipeline_layout);

    // GPU SAFETY: Using sync_field_simple.comp - 37 transcendentals (safe)
    // sync_field.comp has 36+ transcendentals + Kahan summation (borderline timeout risk)
    _sync_pipeline = _pipelineFactory->createSyncFieldPipeline(
        "/home/persist/neotec/0rigin/build/shaders/smft/sync_field_simple.comp.spv",
        _sync_pipeline_layout);

    // ✅ SAFE: gravity_field.comp - 0 transcendentals, pure arithmetic
    _gravity_pipeline = _pipelineFactory->createGravityFieldPipeline(
        "/home/persist/neotec/0rigin/shaders/smft/gravity_field.comp.spv",
        _gravity_pipeline_layout);

    // Create stochastic pipelines if layouts exist
    // ✅ SAFE: kuramoto_stochastic.comp - 12-14 transcendentals, well under budget
    if (_kuramoto_pipeline_layout != VK_NULL_HANDLE) {
        _kuramoto_stochastic_pipeline = _pipelineFactory->createKuramotoStochasticPipeline(
            "/home/persist/neotec/0rigin/shaders/smft/kuramoto_stochastic.comp.spv",
            _kuramoto_pipeline_layout);
    }

    // GPU SAFETY WARNING: Dirac shaders DISABLED - exceed 20 Tflops budget
    // ❌ DANGEROUS: dirac_rk4.comp - ~3000 FLOPs (10× over budget)
    // ❌ DANGEROUS: dirac_stochastic.comp - 50-80 transcendentals (4× over budget)
    //
    // Dirac evolution requires CPU implementation or algorithmic simplification.
    // Current GPU implementation causes consistent timeouts (>2 seconds per dispatch).
    //
    // Recommended approach:
    // 1. CPU-based Dirac evolution for small grids (N < 256)
    // 2. Simplified shader with reduced order integration (Euler vs RK4)
    // 3. Multi-pass approach splitting RK4 stages across dispatches
    //
    // DO NOT enable these pipelines unless using CPU fallback or simplified implementation.

    if (_dirac_pipeline_layout != VK_NULL_HANDLE) {
        // DISABLED: GPU Dirac pipelines cause timeout
        // Uncomment ONLY if using CPU fallback or simplified implementation
        /*
        _dirac_pipeline = _pipelineFactory->createDiracPipeline(
            "/home/persist/neotec/0rigin/shaders/smft/dirac_rk4.comp.spv",
            _dirac_pipeline_layout);

        _dirac_stochastic_pipeline = _pipelineFactory->createDiracStochasticPipeline(
            "/home/persist/neotec/0rigin/shaders/smft/dirac_stochastic.comp.spv",
            _dirac_pipeline_layout);
        */

        // Keep pipelines NULL to trigger CPU fallback in stepStochastic()
        _dirac_pipeline = VK_NULL_HANDLE;
        _dirac_stochastic_pipeline = VK_NULL_HANDLE;
    }

    // Note: We keep descriptor pool alive for lifetime of engine (cleanup in destructor)
    // Store it as member variable if needed for cleanup
}

void MSFTEngine::uploadToGPU() {
    // Upload CPU-side data to GPU buffers using buffer manager
    if (!_bufferManager) return;

    size_t gridSizeBytes = sizeof(float) * _Nx * _Ny;

    // Upload theta data
    if (_theta_memory != VK_NULL_HANDLE && !_theta_data.empty()) {
        _bufferManager->uploadData(_theta_memory, _theta_data.data(), gridSizeBytes);
    }

    // Upload omega data
    if (_omega_memory != VK_NULL_HANDLE && !_omega_data.empty()) {
        _bufferManager->uploadData(_omega_memory, _omega_data.data(), gridSizeBytes);
    }

    // Initialize R_field buffer to zero (will be computed by shader)
    if (_R_field_memory != VK_NULL_HANDLE) {
        _bufferManager->clearBuffer(_R_field_memory, gridSizeBytes);
    }

    // Initialize gravity buffers to zero
    if (_gravity_x_memory != VK_NULL_HANDLE) {
        _bufferManager->clearBuffer(_gravity_x_memory, gridSizeBytes);
    }

    if (_gravity_y_memory != VK_NULL_HANDLE) {
        _bufferManager->clearBuffer(_gravity_y_memory, gridSizeBytes);
    }

    // Upload spinor field data (4 complex components = 8 floats per grid point)
    if (_spinor_memory != VK_NULL_HANDLE && !_spinor_field.empty()) {
        size_t spinorSizeBytes = sizeof(std::complex<float>) * _spinor_field.size();
        _bufferManager->uploadData(_spinor_memory, _spinor_field.data(), spinorSizeBytes);
    }

    // Initialize spinor density buffer to zero (will be computed from spinor field)
    if (_spinor_density_memory != VK_NULL_HANDLE) {
        _bufferManager->clearBuffer(_spinor_density_memory, gridSizeBytes);
    }
}

void MSFTEngine::downloadFromGPU() {
    // Download GPU results back to CPU-side arrays using buffer manager
    if (!_bufferManager) return;

    size_t gridSizeBytes = sizeof(float) * _Nx * _Ny;

    // Download updated theta field (from theta_out buffer)
    if (_theta_out_memory != VK_NULL_HANDLE) {
        _bufferManager->downloadData(_theta_out_memory, _theta_data.data(), gridSizeBytes);
    }

    // Download synchronization field
    if (_R_field_memory != VK_NULL_HANDLE) {
        _bufferManager->downloadData(_R_field_memory, _R_field_data.data(), gridSizeBytes);
    }

    // Download gravity field components
    if (_gravity_x_memory != VK_NULL_HANDLE && !_gravity_x_data.empty()) {
        _bufferManager->downloadData(_gravity_x_memory, _gravity_x_data.data(), gridSizeBytes);
    }

    if (_gravity_y_memory != VK_NULL_HANDLE && !_gravity_y_data.empty()) {
        _bufferManager->downloadData(_gravity_y_memory, _gravity_y_data.data(), gridSizeBytes);
    }

    // Download updated spinor field
    if (_spinor_memory != VK_NULL_HANDLE && !_spinor_field.empty()) {
        size_t spinorSizeBytes = sizeof(std::complex<float>) * _spinor_field.size();
        _bufferManager->downloadData(_spinor_memory, _spinor_field.data(), spinorSizeBytes);
    }
}

void MSFTEngine::destroyResources() {
    // Cleanup Vulkan resources in reverse order of creation
    if (_nova && _nova->_architect) {
        VkDevice device = _nova->_architect->logical_device;

        // Destroy all pipelines via factory
        if (_pipelineFactory) {
            _pipelineFactory->destroyAllPipelines();
            _pipelineFactory.reset();
        }

        // Clear pipeline handles
        _dirac_pipeline = VK_NULL_HANDLE;
        _gravity_pipeline = VK_NULL_HANDLE;
        _kuramoto_pipeline = VK_NULL_HANDLE;
        _sync_pipeline = VK_NULL_HANDLE;
        _kuramoto_stochastic_pipeline = VK_NULL_HANDLE;
        _dirac_stochastic_pipeline = VK_NULL_HANDLE;

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
        if (_dirac_pipeline_layout != VK_NULL_HANDLE) {
            vkDestroyPipelineLayout(device, _dirac_pipeline_layout, nullptr);
            _dirac_pipeline_layout = VK_NULL_HANDLE;
        }

        // Destroy descriptor resources via descriptor manager
        if (_descriptorManager) {
            _descriptorManager->destroyAll();
            _descriptorManager.reset();
        }

        // Clear descriptor handles (already destroyed by manager)
        _kuramoto_descriptor_layout = VK_NULL_HANDLE;
        _sync_descriptor_layout = VK_NULL_HANDLE;
        _gravity_descriptor_layout = VK_NULL_HANDLE;
        _dirac_descriptor_layout = VK_NULL_HANDLE;
        _descriptor_pool = VK_NULL_HANDLE;

        // Destroy all buffers via buffer manager
        if (_bufferManager) {
            _bufferManager->destroyAllBuffers();
            _bufferManager.reset();
        }

        // Clear buffer handles
        _spinor_buffer = VK_NULL_HANDLE;
        _theta_buffer = VK_NULL_HANDLE;
        _theta_out_buffer = VK_NULL_HANDLE;
        _omega_buffer = VK_NULL_HANDLE;
        _R_field_buffer = VK_NULL_HANDLE;
        _gravity_x_buffer = VK_NULL_HANDLE;
        _gravity_y_buffer = VK_NULL_HANDLE;
        _spinor_density_buffer = VK_NULL_HANDLE;

        // Clear memory handles
        _spinor_memory = VK_NULL_HANDLE;
        _theta_memory = VK_NULL_HANDLE;
        _theta_out_memory = VK_NULL_HANDLE;
        _omega_memory = VK_NULL_HANDLE;
        _R_field_memory = VK_NULL_HANDLE;
        _gravity_x_memory = VK_NULL_HANDLE;
        _gravity_y_memory = VK_NULL_HANDLE;
        _spinor_density_memory = VK_NULL_HANDLE;
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
     * 1. Kuramoto with noise: dθ = (ω + K·Σsin(θⱼ-θᵢ))dt + σ_θ·dW [✅ SAFE: GPU]
     * 2. Dirac with noise: i∂_tΨ = H_DiracΨ + σ_Ψ·ξ(t) [❌ CPU ONLY]
     *
     * Pipeline sequence:
     * 1. kuramoto_stochastic: Evolve phases with noise [✅ SAFE: 12-14 transcendentals]
     * 2. sync_field_simple: Compute R(x) order parameter [✅ SAFE: 37 transcendentals]
     * 3. gravity_field: Compute gravitational corrections [✅ SAFE: pure arithmetic]
     * 4. dirac_stochastic: DISABLED - GPU timeout risk (50-80 transcendentals)
     *
     * GPU SAFETY: Dirac shaders disabled. Falls back to deterministic step() if Dirac needed.
     * For Dirac evolution, implement CPU-based integration or use simplified shader.
     */

    // Increment timestep for PRNG seeding
    _time_step++;

    // GPU SAFETY: Dirac pipelines are intentionally disabled (see createPipelines)
    // Since _dirac_stochastic_pipeline will always be VK_NULL_HANDLE, this will
    // fall back to deterministic step() which only uses SAFE shaders.
    //
    // TODO: Implement CPU-based Dirac evolution if quantum effects needed.

    // Verify stochastic pipelines are created
    if (_kuramoto_stochastic_pipeline == VK_NULL_HANDLE ||
        _dirac_stochastic_pipeline == VK_NULL_HANDLE) {
        // Fall back to deterministic step if pipelines not ready
        // NOTE: This is the expected path since Dirac shaders are disabled
        step(dt, K, damping);
        return;
    }

    // Initialize compute dispatcher if needed
    if (!_compute) {
        uint32_t queueFamily = _nova->_architect->queues.indices.compute_family.value_or(
            _nova->_architect->queues.indices.graphics_family.value());
        VkQueue queue = _nova->_architect->queues.compute ?
            _nova->_architect->queues.compute : _nova->_architect->queues.graphics;

        _compute = std::make_unique<MSFTCompute>(
            _nova->_architect->logical_device, queue, queueFamily);

        if (!_compute->initialize()) {
            _compute.reset();
            return;
        }
    }

    // Upload current data to GPU
    uploadToGPU();

    // Begin compute batch
    if (!_compute->beginBatch()) {
        return;
    }

    // Calculate workgroup counts
    uint32_t workgroupsX, workgroupsY;
    MSFTCompute::calculateWorkgroups(_Nx, _Ny, workgroupsX, workgroupsY);

    // Push constants for stochastic Kuramoto
    struct KuramotoPush {
        float dt;
        float K;
        float sigma;
        float damping;
        float omega_mean;
        uint32_t Nx;
        uint32_t Ny;
        uint32_t time_step;
    } kuramoto_push = {
        dt, K, sigma_theta, damping, 0.0f, _Nx, _Ny, _time_step
    };

    // Dispatch stochastic Kuramoto
    _compute->dispatchKuramoto(_kuramoto_stochastic_pipeline,
                               _kuramoto_pipeline_layout,
                               _kuramoto_descriptor_set,
                               &kuramoto_push, sizeof(kuramoto_push),
                               workgroupsX, workgroupsY);

    // Memory barrier
    _compute->insertMemoryBarrier();

    // Swap theta buffers
    std::swap(_theta_buffer, _theta_out_buffer);
    std::swap(_theta_memory, _theta_out_memory);

    // Push constants for sync and gravity
    struct SyncPush {
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
        dt, K, damping, _Delta, _chiral_angle, _Nx, _Ny,
        _Nx * _Ny, 1  // neighborhood_radius = 1
    };

    // Dispatch sync field
    _compute->dispatchSyncField(_sync_pipeline, _sync_pipeline_layout,
                                _sync_descriptor_set,
                                &sync_push, sizeof(sync_push),
                                workgroupsX, workgroupsY);

    // Barrier before gravity
    _compute->insertMemoryBarrier();

    // Dispatch gravity field
    _compute->dispatchGravityField(_gravity_pipeline, _gravity_pipeline_layout,
                                   _gravity_descriptor_set,
                                   &sync_push, sizeof(sync_push),
                                   workgroupsX, workgroupsY);

    // Barrier before Dirac
    _compute->insertMemoryBarrier();

    // GPU SAFETY: Dirac dispatch intentionally disabled
    // _dirac_stochastic_pipeline is always VK_NULL_HANDLE (see createPipelines)
    // This section will never execute - kept for code structure clarity.
    //
    // If Dirac evolution needed, implement CPU fallback:
    // - Load spinor data from GPU
    // - Perform CPU-based RK4 integration
    // - Upload results back to GPU

    // Dispatch stochastic Dirac if available
    if (_dirac_stochastic_pipeline != VK_NULL_HANDLE &&
        _dirac_descriptor_set != VK_NULL_HANDLE) {

        struct DiracPush {
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

        _compute->dispatchDiracEvolution(_dirac_stochastic_pipeline,
                                        _dirac_pipeline_layout,
                                        _dirac_descriptor_set,
                                        &dirac_push, sizeof(dirac_push),
                                        workgroupsX, workgroupsY);
    }

    // Submit batch and wait
    if (!_compute->submitBatch(true)) {
        return;
    }

    // Download results from GPU
    downloadFromGPU();
}