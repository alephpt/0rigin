#include "SMFTEngine.h"
#include "DiracEvolution.h"
#include "KleinGordonEvolution.h"
#include "SMFTCommon.h"
#include "physics/EMFieldComputer.h"
#include <cstring>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <random>
#include <chrono>

/**
 * SMFTEngine Implementation - Phase 1 Skeleton
 *
 * This implementation provides the basic structure for the SMFT physics engine.
 * Phase 1: Basic CPU-side data management and interface implementation
 * Phase 2: Vulkan buffer allocation and management (createBuffers)
 * Phase 3: Shader loading and pipeline creation (createPipelines)
 * Phase 4: GPU compute dispatch implementation (step function)
 */

// Physical constants for SMFT (SI units)
namespace SMFTConstants {
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
    static float computeVacuumPotential(double omega_max) {
        // Δ = ℏω_max/c² = Planck mass (in natural units where ℏ=c=1, this is just omega_max)
        return static_cast<float>(omega_max * HBAR / (C * C));
    }

    /**
     * Compute relativistic Doppler factor (0.md Step 10).
     * Synchronization follows Lorentz dynamics near c = R.
     *
     * @param v Velocity (m/s)
     * @param R Synchronization (0 to 1, acts as effective speed of light)
     * @return Doppler factor √((1 + v/cR)/(1 - v/cR))
     */
    static float computeDopplerFactor(float v, float R) {
        float beta = v / (C * std::max(R, 0.01f)); // v/c_effective, avoid div-by-zero
        if (std::abs(beta) >= 1.0f) {
            return 1e6f; // Return large value for v >= c_effective (causally disconnected)
        }
        return std::sqrt((1.0f + beta) / (1.0f - beta));
    }
}

SMFTEngine::SMFTEngine(Nova* nova)
    : _nova(nova),
      _Nx(0), _Ny(0),
      _Delta(0.0f),
      _chiral_angle(0.0f),
      _time_step(0),
      _cpu_rng(std::random_device{}()),  // Initialize RNG with random seed
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
      _dirac_pipeline_layout(VK_NULL_HANDLE),
      _R_history_index(0),
      _last_dt(0.0f),
      _dirac_evolution(nullptr),
      _dirac_antiparticle(nullptr),
      _kg_evolution(nullptr),
      _dirac_initialized(false),
      _kg_initialized(false),
      _two_particle_mode(false),
      _substep_count(0),
      _substep_ratio(10),  // Default: 10 Kuramoto steps per Dirac step
      _lambda_coupling(0.0f),
      _theta_sum_buffer(VK_NULL_HANDLE),
      _R_sum_buffer(VK_NULL_HANDLE),
      _theta_sum_memory(VK_NULL_HANDLE),
      _R_sum_memory(VK_NULL_HANDLE),
      _accumulation_pipeline(VK_NULL_HANDLE),
      _accumulation_descriptor_set(VK_NULL_HANDLE),
      _accumulation_descriptor_layout(VK_NULL_HANDLE),
      _accumulation_pipeline_layout(VK_NULL_HANDLE),
      _em_coupling_enabled(false),
      _em_coupling_strength(1.0f)
{
    // Ensure Nova is initialized before creating managers
    if (!nova || !nova->initialized) {
        std::cerr << "[SMFTEngine] Error: Nova not initialized" << std::endl;
        return;
    }

    // Initialize pipeline factory
    _pipelineFactory = std::make_unique<SMFTPipelineFactory>(nova->_architect->logical_device);

    // Initialize compute dispatcher
    _compute = std::make_unique<SMFTCompute>(nova->_architect->logical_device,
                                              nova->_architect->queues.compute,
                                              nova->_architect->queues.indices.compute_family.value_or(0));
    if (!_compute->initialize()) {
        std::cerr << "[SMFTEngine] Failed to initialize compute dispatcher" << std::endl;
        return;
    }

    // Initialize buffer manager
    _bufferManager = std::make_unique<SMFTBufferManager>(nova->_architect->logical_device,
                                                          nova->_architect->physical_device);

    // Initialize descriptor manager
    _descriptorManager = std::make_unique<SMFTDescriptorManager>(nova->_architect->logical_device);
}

SMFTEngine::~SMFTEngine() {
    // Clean up Dirac evolution
    if (_dirac_evolution) {
        delete _dirac_evolution;
        _dirac_evolution = nullptr;
    }

    // Clean up antiparticle evolution (Test 3.4)
    if (_dirac_antiparticle) {
        delete _dirac_antiparticle;
        _dirac_antiparticle = nullptr;
    }

    // Clean up Klein-Gordon evolution
    if (_kg_evolution) {
        delete _kg_evolution;
        _kg_evolution = nullptr;
    }

    // Clean up all Vulkan resources
    destroyResources();
}

void SMFTEngine::initialize(uint32_t Nx, uint32_t Ny, float Delta, float chiral_angle) {
    _Nx = Nx;
    _Ny = Ny;
    _Delta = Delta;
    _chiral_angle = chiral_angle;
    _time_step = 0;

    // Log physics interpretation
    if (_nova) {
        std::cout << "\n[SMFTEngine] Physical Interpretation:" << std::endl;
        std::cout << "  Vacuum potential Δ = " << Delta << " (Planck mass in natural units)" << std::endl;
        std::cout << "  Synchronization field R(x) → Emergent spacetime metric" << std::endl;
        std::cout << "  Mass field m(x) = Δ·R(x) → Local energy density" << std::endl;
        std::cout << "  ∇R(x) → Gravitational field (not separate force)" << std::endl;
        std::cout << "  Grid: " << Nx << " x " << Ny << " oscillators" << std::endl;
        std::cout << std::endl;
    }

    // Resize CPU-side data arrays
    size_t total = _Nx * _Ny;
    _theta_data.resize(total);
    _omega_data.resize(total);
    _R_field_data.resize(total);
    _gravity_x_data.resize(total);
    _gravity_y_data.resize(total);

    // Initialize R-field history ring buffer (Phase 4 Test 4.2)
    for (int i = 0; i < 3; i++) {
        _R_history[i].resize(total);
        std::fill(_R_history[i].begin(), _R_history[i].end(), 0.0f);
    }
    _R_history_index = 0;
    _last_dt = 0.0f;

    // Initialize spinor field for Phase 3
    // Each point has 4 complex components (8 floats total)
    _spinor_field.resize(4 * total);

    // Initialize all fields to zero
    std::fill(_theta_data.begin(), _theta_data.end(), 0.0f);
    std::fill(_omega_data.begin(), _omega_data.end(), 0.0f);
    std::fill(_R_field_data.begin(), _R_field_data.end(), 0.0f);
    std::fill(_gravity_x_data.begin(), _gravity_x_data.end(), 0.0f);
    std::fill(_gravity_y_data.begin(), _gravity_y_data.end(), 0.0f);
    std::fill(_spinor_field.begin(), _spinor_field.end(), std::complex<float>(0.0f, 0.0f));

    // Create GPU resources
    createBuffers();
    createPipelines();

    // Initialize EM field GPU resources if EM coupling is enabled
    if (_em_coupling_enabled && _nova && _nova->initialized) {
        initEMBuffers();
        createEMPipelines();
    }

    // Upload initial data to GPU
    uploadToGPU();
}

void SMFTEngine::setInitialPhases(const std::vector<float>& theta) {
    if (theta.size() != _Nx * _Ny) {
        std::cerr << "[SMFTEngine] Error: Invalid theta array size" << std::endl;
        return;
    }

    _theta_data = theta;
    uploadToGPU();
}

void SMFTEngine::setNaturalFrequencies(const std::vector<float>& omega) {
    if (omega.size() != _Nx * _Ny) {
        std::cerr << "[SMFTEngine] Error: Invalid omega array size" << std::endl;
        return;
    }

    _omega_data = omega;
    uploadToGPU();
}

void SMFTEngine::setInitialRField(const std::vector<float>& R_field) {
    if (R_field.size() != _Nx * _Ny) {
        std::cerr << "[SMFTEngine] Error: Invalid R_field array size ("
                  << R_field.size() << " != " << (_Nx * _Ny) << ")" << std::endl;
        return;
    }

    // Validate R values are in [0,1]
    for (size_t i = 0; i < R_field.size(); ++i) {
        if (R_field[i] < 0.0f || R_field[i] > 1.0f) {
            std::cerr << "[SMFTEngine] Warning: R_field[" << i << "] = " << R_field[i]
                      << " out of bounds [0,1], clamping" << std::endl;
        }
    }

    // Update host mirror
    _R_field_data = R_field;

    // Upload to GPU if buffers exist
    if (_bufferManager && _R_field_memory) {
        _bufferManager->uploadData(_R_field_memory, _R_field_data.data(),
                                   _R_field_data.size() * sizeof(float));
    }

    std::cout << "[SMFTEngine] R-field initialized with custom values" << std::endl;
    std::cout << "  R_min = " << *std::min_element(_R_field_data.begin(), _R_field_data.end()) << std::endl;
    std::cout << "  R_max = " << *std::max_element(_R_field_data.begin(), _R_field_data.end()) << std::endl;
}

void SMFTEngine::step(float dt, float K, float damping) {
    // Store timestep for derivative computation (Phase 4 Test 4.2)
    _last_dt = dt;

    // Ensure we have GPU resources
    if (!_kuramoto_pipeline || !_sync_pipeline || !_gravity_pipeline) {
        std::cerr << "[SMFTEngine] Pipelines not created, falling back to CPU simulation" << std::endl;
        return;
    }

    // Calculate workgroup counts (local size = 16x16 in shader)
    uint32_t workgroupsX = (_Nx + 15) / 16;
    uint32_t workgroupsY = (_Ny + 15) / 16;

    // Begin compute batch
    if (!_compute->beginBatch()) {
        std::cerr << "[SMFTEngine] Failed to begin compute batch" << std::endl;
        return;
    }

    // Step 1: Kuramoto phase evolution
    // Push constants for kuramoto shader
    struct KuramotoPush {
        float dt;
        float K;
        float damping;
        float Delta;
        float chiral_angle;
        uint32_t Nx;
        uint32_t Ny;
        uint32_t N_total;
        uint32_t enable_feedback;
    } kuramoto_push = {
        dt, K, damping, _Delta, 0.0f, _Nx, _Ny, _Nx * _Ny, 0
    };

    _compute->dispatchKuramoto(_kuramoto_pipeline,
                               _kuramoto_pipeline_layout,
                               _kuramoto_descriptor_set,
                               &kuramoto_push, sizeof(kuramoto_push),
                               workgroupsX, workgroupsY);

    // Memory barrier between kuramoto and sync
    _compute->insertMemoryBarrier();

    // Step 2: Calculate synchronization field
    // Push constants for sync shader (MUST match shader layout exactly)
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
        dt,              // dt
        K,               // K
        damping,         // damping
        _Delta,          // Delta
        0.0f,            // chiral_angle (unused for now)
        _Nx,             // Nx
        _Ny,             // Ny
        _Nx * _Ny,       // N_total
        1                // neighborhood_radius (1 = 3x3 Moore neighborhood)
    };

    _compute->dispatchSyncField(_sync_pipeline,
                                _sync_pipeline_layout,
                                _sync_descriptor_set,
                                &sync_push, sizeof(sync_push),
                                workgroupsX, workgroupsY);

    // Memory barrier between sync and gravity
    _compute->insertMemoryBarrier();

    // Step 3: Calculate gravitational field
    // Push constants for gravity shader
    struct GravityPush {
        float dt;
        float K;
        float damping;
        float Delta;
        float chiral_angle;
        uint32_t Nx;
        uint32_t Ny;
        uint32_t N_total;
        uint32_t neighborhood_radius;
    } gravity_push = {
        dt, K, damping, _Delta, 0.0f, _Nx, _Ny, _Nx * _Ny, 0
    };

    _compute->dispatchGravityField(_gravity_pipeline,
                                   _gravity_pipeline_layout,
                                   _gravity_descriptor_set,
                                   &gravity_push, sizeof(gravity_push),
                                   workgroupsX, workgroupsY);

    // Step 4: Accumulate for operator splitting (if enabled)
    if (_accumulation_pipeline && _substep_ratio > 1) {
        _compute->insertMemoryBarrier();

        struct AccumulatePush {
            uint32_t grid_size;
        } accumulate_push = { _Nx };

        _compute->dispatchAccumulation(_accumulation_pipeline,
                                      _accumulation_pipeline_layout,
                                      _accumulation_descriptor_set,
                                      &accumulate_push, sizeof(accumulate_push),
                                      workgroupsX, workgroupsY);
    }

    // Submit batch and wait
    if (!_compute->submitBatch(true)) {
        std::cerr << "[SMFTEngine] Failed to submit compute batch" << std::endl;
        return;
    }

    // Download results from GPU
    downloadFromGPU();

    // NOTE: Old internal operator splitting logic removed (lines 358-395)
    // Operator splitting is now handled exclusively by stepWithDirac() external loop
    // This eliminates the recursion bug and simplifies the codebase

    // Swap theta buffers for next iteration
    std::swap(_theta_buffer, _theta_out_buffer);
    std::swap(_theta_memory, _theta_out_memory);
}

void SMFTEngine::stepStochastic(float dt, float K, float damping,
                                float sigma_theta, float sigma_psi) {
    // Ensure we have GPU resources
    if (!_kuramoto_stochastic_pipeline) {
        std::cerr << "[SMFTEngine] Stochastic pipelines not created" << std::endl;
        return;
    }

    // Calculate workgroup counts (local size = 16x16 in shader)
    uint32_t workgroupsX = (_Nx + 15) / 16;
    uint32_t workgroupsY = (_Ny + 15) / 16;

    // Begin compute batch
    if (!_compute->beginBatch()) {
        return;
    }

    // Step 1: Stochastic Kuramoto evolution
    {
        struct KuramotoPush {
            float dt;
            float K;
            float sigma_theta;
            float damping;
            float Delta;
            uint32_t Nx;
            uint32_t Ny;
            uint32_t time_step;
        } kuramoto_push = {
            dt, K, sigma_theta, damping, _Delta, _Nx, _Ny, _time_step
        };

        _compute->dispatchKuramoto(_kuramoto_stochastic_pipeline,
                                   _kuramoto_pipeline_layout,
                                   _kuramoto_descriptor_set,
                                   &kuramoto_push, sizeof(kuramoto_push),
                                   workgroupsX, workgroupsY);
    }

    // Memory barrier
    _compute->insertMemoryBarrier();

    // Step 2: Calculate synchronization field (same as deterministic)
    {
        struct SyncPush {
            uint32_t Nx;
            uint32_t Ny;
        } sync_push = {
            _Nx, _Ny
        };

        _compute->dispatchSyncField(_sync_pipeline,
                                    _sync_pipeline_layout,
                                    _sync_descriptor_set,
                                    &sync_push, sizeof(sync_push),
                                    workgroupsX, workgroupsY);
    }

    // Memory barrier
    _compute->insertMemoryBarrier();

    // Step 3: Stochastic Dirac evolution (if available)
    if (_dirac_stochastic_pipeline) {
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

    // Swap theta buffers for next iteration
    std::swap(_theta_buffer, _theta_out_buffer);
    std::swap(_theta_memory, _theta_out_memory);

    // Increment timestep for PRNG seeding
    _time_step++;
}

void SMFTEngine::stepKuramotoCPUStochastic(float dt, float K, float damping, float sigma) {
    /**
     * CPU-only stochastic Kuramoto evolution step for Phase Transition tests
     *
     * Uses Euler-Maruyama integration:
     * dθ/dt = ω + K·∇²θ - γ·sin(θ) + σ·ξ(t)
     *
     * where ξ(t) is white noise with ⟨ξ⟩ = 0, ⟨ξ²⟩ = 1
     *
     * Implementation: Calls SMFT::stepKuramotoWithNoise from SMFTCommon
     * Then computes synchronization field R(x,y)
     */

    // Step phase field with noise
    SMFT::stepKuramotoWithNoise(_theta_data, _omega_data, dt, K, damping,
                                sigma, _Nx, _Ny, _cpu_rng);

    // Compute synchronization field R(x,y)
    _R_field_data = SMFT::computeLocalR(_theta_data, _Nx, _Ny);

    // Update R-field history ring buffer (Phase 4 Test 4.2)
    _R_history[_R_history_index] = _R_field_data;
    _R_history_index = (_R_history_index + 1) % 3;
    _last_dt = dt;

    // Compute gravitational field (if needed)
    // For phase transition tests, we primarily care about R_field
    // Gravity computation skipped to save CPU cycles

    _time_step++;
}

std::vector<float> SMFTEngine::getSyncField() const {
    return _R_field_data;
}

std::vector<float> SMFTEngine::getMassField() const {
    std::vector<float> mass_field(_R_field_data.size());
    for (size_t i = 0; i < _R_field_data.size(); ++i) {
        mass_field[i] = _Delta * _R_field_data[i];
    }
    return mass_field;
}

std::vector<float> SMFTEngine::getPhaseField() const {
    return _theta_data;
}

std::vector<float> SMFTEngine::getGravitationalField() const {
    std::vector<float> gravity_field(2 * _Nx * _Ny);
    for (size_t i = 0; i < _Nx * _Ny; ++i) {
        gravity_field[2*i + 0] = _gravity_x_data[i];
        gravity_field[2*i + 1] = _gravity_y_data[i];
    }
    return gravity_field;
}

std::complex<float> SMFTEngine::getSpinorComponent(uint32_t x, uint32_t y, uint32_t component) const {
    if (x >= _Nx || y >= _Ny || component >= 4) {
        return std::complex<float>(0.0f, 0.0f);
    }
    uint32_t idx = 4 * (y * _Nx + x) + component;
    return _spinor_field[idx];
}

void SMFTEngine::createBuffers() {
    if (!_nova || !_nova->initialized) {
        std::cerr << "[SMFTEngine] Nova not initialized, cannot create buffers" << std::endl;
        return;
    }

    size_t float_size = sizeof(float) * _Nx * _Ny;
    size_t spinor_size = sizeof(std::complex<float>) * 4 * _Nx * _Ny;

    // Create buffers using buffer manager - all need HOST_VISIBLE for CPU access
    auto theta_pair = _bufferManager->createStorageBuffer(float_size);
    _theta_buffer = theta_pair.first;
    _theta_memory = theta_pair.second;

    auto theta_out_pair = _bufferManager->createStorageBuffer(float_size);
    _theta_out_buffer = theta_out_pair.first;
    _theta_out_memory = theta_out_pair.second;

    auto omega_pair = _bufferManager->createStorageBuffer(float_size);
    _omega_buffer = omega_pair.first;
    _omega_memory = omega_pair.second;

    auto R_field_pair = _bufferManager->createStorageBuffer(float_size);
    _R_field_buffer = R_field_pair.first;
    _R_field_memory = R_field_pair.second;

    auto gravity_x_pair = _bufferManager->createStorageBuffer(float_size);
    _gravity_x_buffer = gravity_x_pair.first;
    _gravity_x_memory = gravity_x_pair.second;

    auto gravity_y_pair = _bufferManager->createStorageBuffer(float_size);
    _gravity_y_buffer = gravity_y_pair.first;
    _gravity_y_memory = gravity_y_pair.second;

    auto spinor_density_pair = _bufferManager->createStorageBuffer(float_size);
    _spinor_density_buffer = spinor_density_pair.first;
    _spinor_density_memory = spinor_density_pair.second;

    auto spinor_pair = _bufferManager->createStorageBuffer(spinor_size);
    _spinor_buffer = spinor_pair.first;
    _spinor_memory = spinor_pair.second;

    // Create accumulator buffers for operator splitting (GPU-side time averaging)
    auto accumulators = _bufferManager->createAccumulatorBuffers(float_size);
    if (accumulators.size() == 2) {
        _theta_sum_buffer = accumulators[0].first;
        _theta_sum_memory = accumulators[0].second;
        _R_sum_buffer = accumulators[1].first;
        _R_sum_memory = accumulators[1].second;
    }

    // Verify all buffers were created
    if (_theta_buffer == VK_NULL_HANDLE || _theta_out_buffer == VK_NULL_HANDLE ||
        _omega_buffer == VK_NULL_HANDLE || _R_field_buffer == VK_NULL_HANDLE ||
        _gravity_x_buffer == VK_NULL_HANDLE || _gravity_y_buffer == VK_NULL_HANDLE ||
        _spinor_density_buffer == VK_NULL_HANDLE || _spinor_buffer == VK_NULL_HANDLE ||
        _theta_sum_buffer == VK_NULL_HANDLE || _R_sum_buffer == VK_NULL_HANDLE) {
        std::cerr << "[SMFTEngine] Failed to create one or more buffers" << std::endl;
        destroyResources();
        return;
    }

    if (_nova) {
        std::cout << "[SMFTEngine] Created GPU buffers for " << _Nx << "x" << _Ny << " grid" << std::endl;
        std::cout << "  Phase buffers: " << float_size << " bytes each" << std::endl;
        std::cout << "  Spinor buffer: " << spinor_size << " bytes" << std::endl;
        std::cout << "  Accumulator buffers: " << float_size << " bytes each (operator splitting)" << std::endl;
    }
}

void SMFTEngine::createPipelines() {
    if (!_nova || !_nova->initialized) {
        std::cerr << "[SMFTEngine] Nova not initialized, cannot create pipelines" << std::endl;
        return;
    }

    // Create descriptor pool for all pipelines
    VkDescriptorPoolSize poolSizes[] = {
        { VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 32 }  // Enough for all pipelines
    };

    VkDescriptorPoolCreateInfo poolInfo{};
    poolInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
    poolInfo.flags = VK_DESCRIPTOR_POOL_CREATE_FREE_DESCRIPTOR_SET_BIT;
    poolInfo.maxSets = 10;  // kuramoto, sync, gravity, dirac, and stochastic variants
    poolInfo.poolSizeCount = 1;
    poolInfo.pPoolSizes = poolSizes;

    if (vkCreateDescriptorPool(_nova->_architect->logical_device, &poolInfo, nullptr, &_descriptor_pool) != VK_SUCCESS) {
        std::cerr << "[SMFTEngine] Failed to create descriptor pool" << std::endl;
        return;
    }

    // Create Kuramoto pipeline
    {
        // Define bindings for kuramoto shader
        std::vector<VkDescriptorSetLayoutBinding> bindings = {
            {0, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr}, // theta_in
            {1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr}, // theta_out
            {2, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr}, // omega
            {3, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr}  // spinor_density
        };

        _kuramoto_descriptor_layout = _descriptorManager->createDescriptorSetLayout(bindings);

        // Create pipeline layout manually
        VkPushConstantRange pushRange{};
        pushRange.stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
        pushRange.offset = 0;
        pushRange.size = sizeof(float) * 5 + sizeof(uint32_t) * 4;  // 36 bytes: dt,K,damping,Delta,chiral_angle,Nx,Ny,N_total,enable_feedback

        VkPipelineLayoutCreateInfo layoutInfo{};
        layoutInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
        layoutInfo.setLayoutCount = 1;
        layoutInfo.pSetLayouts = &_kuramoto_descriptor_layout;
        layoutInfo.pushConstantRangeCount = 1;
        layoutInfo.pPushConstantRanges = &pushRange;

        vkCreatePipelineLayout(_nova->_architect->logical_device, &layoutInfo, nullptr, &_kuramoto_pipeline_layout);

        _kuramoto_pipeline = _pipelineFactory->createKuramotoPipeline("shaders/smft/kuramoto.comp.spv", _kuramoto_pipeline_layout);

        // Allocate and update descriptor set
        _kuramoto_descriptor_set = _descriptorManager->allocateDescriptorSet(_descriptor_pool, _kuramoto_descriptor_layout);

        VkDescriptorBufferInfo bufferInfos[] = {
            {_theta_buffer, 0, VK_WHOLE_SIZE},
            {_theta_out_buffer, 0, VK_WHOLE_SIZE},
            {_omega_buffer, 0, VK_WHOLE_SIZE},
            {_spinor_density_buffer, 0, VK_WHOLE_SIZE}
        };

        VkWriteDescriptorSet descriptorWrites[] = {
            {VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET, nullptr, _kuramoto_descriptor_set, 0, 0, 1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, nullptr, &bufferInfos[0], nullptr},
            {VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET, nullptr, _kuramoto_descriptor_set, 1, 0, 1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, nullptr, &bufferInfos[1], nullptr},
            {VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET, nullptr, _kuramoto_descriptor_set, 2, 0, 1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, nullptr, &bufferInfos[2], nullptr},
            {VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET, nullptr, _kuramoto_descriptor_set, 3, 0, 1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, nullptr, &bufferInfos[3], nullptr}
        };

        vkUpdateDescriptorSets(_nova->_architect->logical_device, 4, descriptorWrites, 0, nullptr);
    }

    // Create Sync field pipeline
    {
        // Define bindings for sync shader
        std::vector<VkDescriptorSetLayoutBinding> bindings = {
            {0, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr}, // theta
            {1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr}  // R_field
        };

        _sync_descriptor_layout = _descriptorManager->createDescriptorSetLayout(bindings);

        // Create pipeline layout manually
        VkPushConstantRange pushRange{};
        pushRange.stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
        pushRange.offset = 0;
        pushRange.size = sizeof(float) * 5 + sizeof(uint32_t) * 4;  // 36 bytes: dt,K,damping,Delta,chiral_angle,Nx,Ny,N_total,neighborhood_radius

        VkPipelineLayoutCreateInfo layoutInfo{};
        layoutInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
        layoutInfo.setLayoutCount = 1;
        layoutInfo.pSetLayouts = &_sync_descriptor_layout;
        layoutInfo.pushConstantRangeCount = 1;
        layoutInfo.pPushConstantRanges = &pushRange;

        vkCreatePipelineLayout(_nova->_architect->logical_device, &layoutInfo, nullptr, &_sync_pipeline_layout);

        _sync_pipeline = _pipelineFactory->createSyncFieldPipeline("shaders/smft/sync_field.comp.spv", _sync_pipeline_layout);

        // Allocate and update descriptor set
        _sync_descriptor_set = _descriptorManager->allocateDescriptorSet(_descriptor_pool, _sync_descriptor_layout);

        VkDescriptorBufferInfo bufferInfos[] = {
            {_theta_buffer, 0, VK_WHOLE_SIZE},
            {_R_field_buffer, 0, VK_WHOLE_SIZE}
        };

        VkWriteDescriptorSet descriptorWrites[] = {
            {VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET, nullptr, _sync_descriptor_set, 0, 0, 1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, nullptr, &bufferInfos[0], nullptr},
            {VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET, nullptr, _sync_descriptor_set, 1, 0, 1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, nullptr, &bufferInfos[1], nullptr}
        };

        vkUpdateDescriptorSets(_nova->_architect->logical_device, 2, descriptorWrites, 0, nullptr);
    }

    // Create Gravity field pipeline
    {
        // Define bindings for gravity shader
        std::vector<VkDescriptorSetLayoutBinding> bindings = {
            {0, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr}, // R_field
            {1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr}, // gravity_x
            {2, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr}  // gravity_y
        };

        _gravity_descriptor_layout = _descriptorManager->createDescriptorSetLayout(bindings);

        // Create pipeline layout manually
        VkPushConstantRange pushRange{};
        pushRange.stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
        pushRange.offset = 0;
        pushRange.size = sizeof(float) * 5 + sizeof(uint32_t) * 4;  // 36 bytes: dt,K,damping,Delta,chiral_angle,Nx,Ny,N_total,neighborhood_radius

        VkPipelineLayoutCreateInfo layoutInfo{};
        layoutInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
        layoutInfo.setLayoutCount = 1;
        layoutInfo.pSetLayouts = &_gravity_descriptor_layout;
        layoutInfo.pushConstantRangeCount = 1;
        layoutInfo.pPushConstantRanges = &pushRange;

        vkCreatePipelineLayout(_nova->_architect->logical_device, &layoutInfo, nullptr, &_gravity_pipeline_layout);

        _gravity_pipeline = _pipelineFactory->createGravityFieldPipeline("shaders/smft/gravity.comp.spv", _gravity_pipeline_layout);

        // Allocate and update descriptor set
        _gravity_descriptor_set = _descriptorManager->allocateDescriptorSet(_descriptor_pool, _gravity_descriptor_layout);

        VkDescriptorBufferInfo bufferInfos[] = {
            {_R_field_buffer, 0, VK_WHOLE_SIZE},
            {_gravity_x_buffer, 0, VK_WHOLE_SIZE},
            {_gravity_y_buffer, 0, VK_WHOLE_SIZE}
        };

        VkWriteDescriptorSet descriptorWrites[] = {
            {VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET, nullptr, _gravity_descriptor_set, 0, 0, 1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, nullptr, &bufferInfos[0], nullptr},
            {VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET, nullptr, _gravity_descriptor_set, 1, 0, 1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, nullptr, &bufferInfos[1], nullptr},
            {VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET, nullptr, _gravity_descriptor_set, 2, 0, 1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, nullptr, &bufferInfos[2], nullptr}
        };

        vkUpdateDescriptorSets(_nova->_architect->logical_device, 3, descriptorWrites, 0, nullptr);
    }

    // Create stochastic Kuramoto pipeline
    {
        _kuramoto_stochastic_pipeline = _pipelineFactory->createKuramotoStochasticPipeline("shaders/smft/kuramoto_stochastic.comp.spv", _kuramoto_pipeline_layout);
        if (!_kuramoto_stochastic_pipeline) {
            std::cerr << "[SMFTEngine] Warning: Failed to load stochastic Kuramoto shader" << std::endl;
        }
    }

    // Create stochastic Dirac pipeline
    {
        // Define bindings for dirac shader
        std::vector<VkDescriptorSetLayoutBinding> bindings = {
            {0, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr}, // spinor field
            {1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr}, // R_field
            {2, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr}  // spinor_density
        };

        _dirac_descriptor_layout = _descriptorManager->createDescriptorSetLayout(bindings);

        // Create pipeline layout manually
        VkPushConstantRange pushRange{};
        pushRange.stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
        pushRange.offset = 0;
        pushRange.size = sizeof(float) * 8;

        VkPipelineLayoutCreateInfo layoutInfo{};
        layoutInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
        layoutInfo.setLayoutCount = 1;
        layoutInfo.pSetLayouts = &_dirac_descriptor_layout;
        layoutInfo.pushConstantRangeCount = 1;
        layoutInfo.pPushConstantRanges = &pushRange;

        vkCreatePipelineLayout(_nova->_architect->logical_device, &layoutInfo, nullptr, &_dirac_pipeline_layout);

        _dirac_stochastic_pipeline = _pipelineFactory->createDiracStochasticPipeline("shaders/smft/dirac_stochastic.comp.spv", _dirac_pipeline_layout);

        if (_dirac_stochastic_pipeline) {
            // Allocate and update descriptor set
            _dirac_descriptor_set = _descriptorManager->allocateDescriptorSet(_descriptor_pool, _dirac_descriptor_layout);

            VkDescriptorBufferInfo bufferInfos[] = {
                {_spinor_buffer, 0, VK_WHOLE_SIZE},
                {_R_field_buffer, 0, VK_WHOLE_SIZE},
                {_spinor_density_buffer, 0, VK_WHOLE_SIZE}
            };

            VkWriteDescriptorSet descriptorWrites[] = {
                {VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET, nullptr, _dirac_descriptor_set, 0, 0, 1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, nullptr, &bufferInfos[0], nullptr},
                {VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET, nullptr, _dirac_descriptor_set, 1, 0, 1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, nullptr, &bufferInfos[1], nullptr},
                {VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET, nullptr, _dirac_descriptor_set, 2, 0, 1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, nullptr, &bufferInfos[2], nullptr}
            };

            vkUpdateDescriptorSets(_nova->_architect->logical_device, 3, descriptorWrites, 0, nullptr);
        } else {
            std::cerr << "[SMFTEngine] Warning: Failed to load stochastic Dirac shader" << std::endl;
        }
    }

    // Create Accumulation pipeline for operator splitting
    {
        std::vector<VkDescriptorSetLayoutBinding> bindings = {
            {0, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr}, // theta (readonly)
            {1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr}, // R (readonly)
            {2, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr}, // theta_sum (read-write)
            {3, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr}  // R_sum (read-write)
        };

        _accumulation_descriptor_layout = _descriptorManager->createDescriptorSetLayout(bindings);

        VkPushConstantRange pushRange{};
        pushRange.stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
        pushRange.offset = 0;
        pushRange.size = sizeof(uint32_t);  // grid_size

        VkPipelineLayoutCreateInfo layoutInfo{};
        layoutInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
        layoutInfo.setLayoutCount = 1;
        layoutInfo.pSetLayouts = &_accumulation_descriptor_layout;
        layoutInfo.pushConstantRangeCount = 1;
        layoutInfo.pPushConstantRanges = &pushRange;

        vkCreatePipelineLayout(_nova->_architect->logical_device, &layoutInfo, nullptr, &_accumulation_pipeline_layout);

        _accumulation_pipeline = _pipelineFactory->createAccumulationPipeline(
            "shaders/smft/accumulate.comp.spv", _accumulation_pipeline_layout);

        _accumulation_descriptor_set = _descriptorManager->allocateDescriptorSet(
            _descriptor_pool, _accumulation_descriptor_layout);

        VkDescriptorBufferInfo bufferInfos[] = {
            {_theta_buffer, 0, VK_WHOLE_SIZE},
            {_R_field_buffer, 0, VK_WHOLE_SIZE},
            {_theta_sum_buffer, 0, VK_WHOLE_SIZE},
            {_R_sum_buffer, 0, VK_WHOLE_SIZE}
        };

        VkWriteDescriptorSet descriptorWrites[] = {
            {VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET, nullptr, _accumulation_descriptor_set, 0, 0, 1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, nullptr, &bufferInfos[0], nullptr},
            {VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET, nullptr, _accumulation_descriptor_set, 1, 0, 1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, nullptr, &bufferInfos[1], nullptr},
            {VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET, nullptr, _accumulation_descriptor_set, 2, 0, 1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, nullptr, &bufferInfos[2], nullptr},
            {VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET, nullptr, _accumulation_descriptor_set, 3, 0, 1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, nullptr, &bufferInfos[3], nullptr}
        };

        vkUpdateDescriptorSets(_nova->_architect->logical_device, 4, descriptorWrites, 0, nullptr);
    }

    if (_nova) {
        std::cout << "[SMFTEngine] Created compute pipelines" << std::endl;
        if (_kuramoto_pipeline) std::cout << "  ✓ Kuramoto pipeline" << std::endl;
        if (_sync_pipeline) std::cout << "  ✓ Sync field pipeline" << std::endl;
        if (_gravity_pipeline) std::cout << "  ✓ Gravity field pipeline" << std::endl;
        if (_kuramoto_stochastic_pipeline) std::cout << "  ✓ Stochastic Kuramoto pipeline" << std::endl;
        if (_dirac_stochastic_pipeline) std::cout << "  ✓ Stochastic Dirac pipeline" << std::endl;
        if (_accumulation_pipeline) std::cout << "  ✓ Accumulation pipeline (operator splitting)" << std::endl;
    }
}

void SMFTEngine::uploadToGPU() {
    if (!_bufferManager) return;

    // Upload phase field data
    _bufferManager->uploadData(_theta_memory, _theta_data.data(),
                              sizeof(float) * _theta_data.size());

    // Upload natural frequencies
    _bufferManager->uploadData(_omega_memory, _omega_data.data(),
                              sizeof(float) * _omega_data.size());

    // Upload spinor field data
    _bufferManager->uploadData(_spinor_memory, _spinor_field.data(),
                              sizeof(std::complex<float>) * _spinor_field.size());

    // Disabled verbose logging - too noisy during operator splitting
    // if (_nova) {
    //     std::cout << "[SMFTEngine] Uploaded data to GPU" << std::endl;
    // }
}

void SMFTEngine::downloadFromGPU() {
    if (!_bufferManager) return;

    // Download phase field (from theta_out_buffer)
    _bufferManager->downloadData(_theta_out_memory, _theta_data.data(),
                                sizeof(float) * _theta_data.size());

    // Download synchronization field
    _bufferManager->downloadData(_R_field_memory, _R_field_data.data(),
                                sizeof(float) * _R_field_data.size());

    // Update R-field history ring buffer (Phase 4 Test 4.2)
    // Store current R-field before updating to new timestep
    _R_history[_R_history_index] = _R_field_data;
    _R_history_index = (_R_history_index + 1) % 3;

    // Download gravitational field components
    _bufferManager->downloadData(_gravity_x_memory, _gravity_x_data.data(),
                                sizeof(float) * _gravity_x_data.size());
    _bufferManager->downloadData(_gravity_y_memory, _gravity_y_data.data(),
                                sizeof(float) * _gravity_y_data.size());

    // Disabled verbose logging - too noisy during operator splitting
    // if (_nova) {
    //     std::cout << "[SMFTEngine] Downloaded results from GPU" << std::endl;
    // }
}

void SMFTEngine::destroyResources() {
    // DOUBLE-FREE FIX: Let component managers handle resource destruction
    // The unique_ptr members (_pipelineFactory, _bufferManager, _descriptorManager)
    // already have proper cleanup logic in their destructors. Manual cleanup here
    // was causing resources to be freed twice, leading to crashes.

    if (_nova && _nova->_architect && _nova->_architect->logical_device) {
        vkDeviceWaitIdle(_nova->_architect->logical_device);

        // Clean up EM-specific resources before component managers
        cleanupEMResources();

        // Component managers will handle cleanup automatically via their destructors:
        // - SMFTPipelineFactory destroys pipelines and layouts
        // - SMFTBufferManager destroys buffers and frees memory
        // - SMFTDescriptorManager destroys descriptor layouts and pool
    }

    if (_nova) {
        std::cout << "[SMFTEngine] Resources cleanup delegated to component managers" << std::endl;
    }
}

// ============================================================================
// DIRAC FIELD METHODS (CPU-only implementation)
// ============================================================================

void SMFTEngine::initializeDiracField(float x0, float y0, float sigma, float amplitude) {
    /**
     * Initialize Dirac spinor field with Gaussian wavepacket
     * Uses split-operator DiracEvolution class
     */

    // Create DiracEvolution instance if not exists
    if (!_dirac_evolution) {
        _dirac_evolution = new DiracEvolution(_Nx, _Ny);
    }

    // Initialize with Gaussian wavepacket
    _dirac_evolution->initialize(x0, y0, sigma);

    _dirac_initialized = true;

    // Log initialization
    if (_nova) {
        float norm = _dirac_evolution->getNorm();
        std::cout << "[SMFTEngine] Dirac field initialized: "
                  << "center=(" << x0 << "," << y0 << "), "
                  << "sigma=" << sigma << ", "
                  << "norm=" << norm << std::endl;
    }
}

void SMFTEngine::initializeBoostedDiracField(float x0, float y0, float sigma,
                                             float vx, float vy, float R_bg) {
    /**
     * Initialize Dirac spinor field with boosted Gaussian wavepacket
     * Uses SMFT::initializeBoostedGaussian helper function
     */

    // Create DiracEvolution instance if not exists
    if (!_dirac_evolution) {
        _dirac_evolution = new DiracEvolution(_Nx, _Ny);
    }

    // Initialize with boosted Gaussian wavepacket using SMFT helper
    SMFT::initializeBoostedGaussian(*_dirac_evolution, x0, y0, sigma,
                                   vx, vy, _Delta, R_bg);

    _dirac_initialized = true;

    // Log initialization
    if (_nova) {
        float norm = _dirac_evolution->getNorm();
        float x_mean, y_mean;
        _dirac_evolution->getCenterOfMass(x_mean, y_mean);
        std::cout << "[SMFTEngine] Boosted Dirac field initialized: "
                  << "center=(" << x0 << "," << y0 << "), "
                  << "sigma=" << sigma << ", "
                  << "boost=(" << vx << "," << vy << ")c, "
                  << "norm=" << norm << ", "
                  << "<r>=(" << x_mean << "," << y_mean << ")" << std::endl;
    }
}

void SMFTEngine::stepWithDirac(float dt, float lambda_coupling, int substep_ratio, float K, float damping) {
    /**
     * Step coupled Kuramoto-Dirac evolution with mass coupling
     * CPU-only implementation - uses Strang splitting for symplectic evolution
     *
     * Physics:
     * - Dirac evolution: i·dΨ/dt = [-iα·∇ + β·m(x)]Ψ
     * - Klein-Gordon evolution: ∂²_tφ = ∇²φ - m²φ
     * - Mass coupling: m(x,y) = Δ·R(x,y) from synchronization field
     * - Feedback: field density influences phase field via coupling λ
     * - Operator splitting: N Kuramoto steps per 1 field step (Born-Oppenheimer)
     *
     * Strang Splitting (2nd order symplectic):
     * 1. Kuramoto half-step (dt/2): N/2 substeps with dt/(2·substep_ratio)
     * 2. Field full-step (dt): Full timestep for Dirac/Klein-Gordon evolution
     * 3. Kuramoto half-step (dt/2): N/2 substeps with dt/(2·substep_ratio)
     *
     * This preserves Hamiltonian structure and eliminates systematic energy drift.
     */

    // Check which solver is initialized
    bool has_dirac = (_dirac_initialized && _dirac_evolution != nullptr);
    bool has_kg = (_kg_initialized && _kg_evolution != nullptr);

    if (!has_dirac && !has_kg) {
        std::cerr << "[SMFTEngine::stepWithDirac] ERROR: No field initialized (neither Dirac nor Klein-Gordon)" << std::endl;
        return;
    }

    // Strang splitting parameters
    float substep_dt = dt / static_cast<float>(substep_ratio);
    int half_substeps = substep_ratio / 2;  // For even N, this is exact
    int remaining_substeps = substep_ratio - half_substeps;  // Handle odd N

    static int call_count = 0;
    if (call_count == 0) {
        std::cout << "[stepWithDirac] Strang splitting enabled:" << std::endl;
        std::cout << "  substep_ratio=" << substep_ratio << std::endl;
        std::cout << "  substep_dt=" << substep_dt << std::endl;
        std::cout << "  half_substeps=" << half_substeps << std::endl;
        std::cout << "  remaining_substeps=" << remaining_substeps << std::endl;
    }
    call_count++;

    // ========================================================================
    // STEP 1: Kuramoto half-step (dt/2)
    // ========================================================================
    for (int n = 0; n < half_substeps; ++n) {
        step(substep_dt, K, damping);
    }

    // ========================================================================
    // STEP 2: Field full-step (dt) - Dirac or Klein-Gordon
    // ========================================================================

    // Get current synchronization field for mass coupling
    std::vector<float> R_field = getSyncField();

    // Compute mass field: m(x,y) = Δ·R(x,y)
    std::vector<float> mass_field(_Nx * _Ny);

    // CRITICAL: For Klein-Gordon, use constant mass (R=1) to ensure norm conservation
    // Klein-Gordon evolution assumes constant mass - time-dependent mass breaks unitarity
    if (has_kg) {
        // Use R=1.0 everywhere for Klein-Gordon (constant background mass)
        for (uint32_t i = 0; i < _Nx * _Ny; i++) {
            mass_field[i] = _Delta * 1.0f;
        }
    } else {
        // Dirac can handle time-dependent mass (R-field coupling)
        for (uint32_t i = 0; i < _Nx * _Ny; i++) {
            mass_field[i] = _Delta * R_field[i];
        }
    }

    // ========================================================================
    // ELECTROMAGNETIC COUPLING (Phase 5 - Scenario 2.6B)
    // ========================================================================
    // Compute EM fields from Kuramoto phase gradients: A_μ = ∂_μ θ
    // Apply minimal coupling to Dirac evolution if enabled

    if (_em_coupling_enabled && has_dirac && _dirac_evolution) {
        // Initialize theta_previous on first call
        if (_theta_previous.empty()) {
            _theta_previous = _theta_data;
        }

        // Try GPU path first, fallback to CPU if GPU unavailable
        bool gpu_em_success = false;

        if (_em_potentials_pipeline != VK_NULL_HANDLE && _bufferManager && _compute) {
            // GPU PATH: Upload theta data and compute EM fields on GPU
            try {
                // Upload current theta to GPU
                _bufferManager->uploadData(_theta_memory, _theta_data.data(),
                                          _theta_data.size() * sizeof(float));

                // Upload previous theta for time derivatives
                _bufferManager->uploadData(_theta_previous_memory, _theta_previous.data(),
                                          _theta_previous.size() * sizeof(float));

                // Execute GPU EM computation
                computeEMFieldsGPU();

                // Download results from GPU
                _phi_data.resize(_Nx * _Ny);
                _A_x_data.resize(_Nx * _Ny);
                _A_y_data.resize(_Nx * _Ny);
                _E_x_data.resize(_Nx * _Ny);
                _E_y_data.resize(_Nx * _Ny);
                _B_z_data.resize(_Nx * _Ny);

                _bufferManager->downloadData(_phi_memory, _phi_data.data(), _phi_data.size() * sizeof(float));
                _bufferManager->downloadData(_A_x_memory, _A_x_data.data(), _A_x_data.size() * sizeof(float));
                _bufferManager->downloadData(_A_y_memory, _A_y_data.data(), _A_y_data.size() * sizeof(float));
                _bufferManager->downloadData(_E_x_memory, _E_x_data.data(), _E_x_data.size() * sizeof(float));
                _bufferManager->downloadData(_E_y_memory, _E_y_data.data(), _E_y_data.size() * sizeof(float));
                _bufferManager->downloadData(_B_z_memory, _B_z_data.data(), _B_z_data.size() * sizeof(float));

                gpu_em_success = true;

            } catch (const std::exception& e) {
                std::cerr << "[SMFTEngine::stepWithDirac] GPU EM computation failed: " << e.what() << std::endl;
                std::cerr << "  Falling back to CPU EM computation" << std::endl;
                gpu_em_success = false;
            }
        }

        // CPU FALLBACK: Compute EM fields on CPU if GPU unavailable or failed
        if (!gpu_em_success) {
            double dx = 1.0;  // Unit grid spacing in natural units
            double dy = 1.0;
            auto em_fields = EMFieldComputer::computeFromPhase(
                _theta_data,
                _theta_previous,
                _R_field_data,  // Pass R_field for conjugate product method
                static_cast<int>(_Nx),
                static_cast<int>(_Ny),
                dx, dy, static_cast<double>(dt)
            );

            // Extract fields and convert to std::vector<float>
            _phi_data.resize(_Nx * _Ny);
            _A_x_data.resize(_Nx * _Ny);
            _A_y_data.resize(_Nx * _Ny);
            _E_x_data.resize(_Nx * _Ny);
            _E_y_data.resize(_Nx * _Ny);
            _B_z_data.resize(_Nx * _Ny);

            for (uint32_t i = 0; i < _Nx; i++) {
                for (uint32_t j = 0; j < _Ny; j++) {
                    int idx = j * _Nx + i;  // Row-major indexing
                    _phi_data[idx] = static_cast<float>(em_fields.phi(i, j));
                    _A_x_data[idx] = static_cast<float>(em_fields.A_x(i, j));
                    _A_y_data[idx] = static_cast<float>(em_fields.A_y(i, j));
                    _E_x_data[idx] = static_cast<float>(em_fields.E_x(i, j));
                    _E_y_data[idx] = static_cast<float>(em_fields.E_y(i, j));
                    _B_z_data[idx] = static_cast<float>(em_fields.B_z(i, j));
                }
            }
        }

        // Apply EM coupling to Dirac evolution
        // Physics: Minimal coupling ∇ → ∇ - iqA and scalar potential φ
        _dirac_evolution->applyEMPotentialStep(_phi_data, dt);
        _dirac_evolution->applyMinimalCoupling(_A_x_data, _A_y_data);

        // Apply to antiparticle if in two-particle mode
        if (_two_particle_mode && _dirac_antiparticle) {
            _dirac_antiparticle->applyEMPotentialStep(_phi_data, dt);
            _dirac_antiparticle->applyMinimalCoupling(_A_x_data, _A_y_data);
        }

        // Update theta history for next timestep
        _theta_previous = _theta_data;
    }

    // Split-operator evolution step (unitary, preserves norm exactly)
    if (has_dirac) {
        _dirac_evolution->step(mass_field, dt);

        // If two-particle mode, also evolve antiparticle
        if (_two_particle_mode && _dirac_antiparticle) {
            _dirac_antiparticle->step(mass_field, dt);
        }
    } else if (has_kg) {
        _kg_evolution->step(mass_field, dt);
    }

    // ========================================================================
    // STEP 3: Kuramoto half-step (dt/2)
    // ========================================================================
    for (int n = 0; n < remaining_substeps; ++n) {
        step(substep_dt, K, damping);
    }

    // ========================================================================
    // Optional: Feedback to phase field (via lambda_coupling)
    // ========================================================================
    if (lambda_coupling > 0.0f && _theta_data.size() == _Nx * _Ny) {
        std::vector<float> density;
        if (has_dirac) {
            density = _dirac_evolution->getDensity();
        } else if (has_kg) {
            density = _kg_evolution->getDensity();
        }

        for (uint32_t i = 0; i < _Nx * _Ny; i++) {
            _theta_data[i] += lambda_coupling * dt * density[i];
        }

        // Re-upload phase data to GPU if modified
        if (_theta_buffer != VK_NULL_HANDLE) {
            uploadToGPU();
        }
    }
}

std::vector<float> SMFTEngine::getDiracDensity() const {
    /**
     * Get Dirac spinor density |Ψ|² for analysis
     * CPU-only implementation - uses split-operator DiracEvolution class
     *
     * Returns density at each grid point.
     * Full spinor: |Ψ|² = |ψ₀|² + |ψ₁|² + |ψ₂|² + |ψ₃|²
     */

    if (!_dirac_initialized || !_dirac_evolution) {
        return std::vector<float>(_Nx * _Ny, 0.0f);  // Return zeros if not initialized
    }

    return _dirac_evolution->getDensity();
}

const DiracEvolution* SMFTEngine::getDiracEvolution() const {
    /**
     * Get internal DiracEvolution object for observable computation
     * Returns nullptr if not initialized
     */
    return _dirac_evolution;
}

DiracEvolution* SMFTEngine::getDiracEvolutionNonConst() {
    /**
     * Get internal DiracEvolution object for modification (e.g., for specific analysis)
     * Returns nullptr if not initialized
     */
    return _dirac_evolution;
}

// ============================================================================
// OPERATOR SPLITTING METHODS (GPU-CPU Hybrid with Adiabatic Approximation)
// ============================================================================

void SMFTEngine::setSubstepRatio(int N) {
    /**
     * Set the substep ratio for operator splitting
     * N = number of fast Kuramoto steps per slow Dirac step
     *
     * Typical values based on timescale separation:
     * - N = 10 for testing/debugging
     * - N = 100 for production (matches ~100× physical timescale ratio)
     */
    _substep_ratio = N;
    _substep_count = 0;  // Reset counter
}

void SMFTEngine::updateAveragedFields(const std::vector<float>& theta_avg,
                                      const std::vector<float>& R_avg) {
    /**
     * Update time-averaged fields for Dirac evolution
     * Called internally after GPU accumulation and averaging completes
     *
     * These averaged fields represent the "slow" variables seen by the
     * fast-timescale Kuramoto subsystem during operator splitting.
     */
    _theta_avg = theta_avg;
    _R_avg = R_avg;
}

void SMFTEngine::initializeHybrid(float x0, float y0, float sigma) {
    /**
     * Initialize hybrid GPU-CPU system with operator splitting
     *
     * Sets up:
     * 1. Kuramoto (GPU) - already initialized via initialize()
     * 2. Dirac (CPU) - Gaussian wavepacket at defect location
     * 3. Accumulator buffers (GPU) - for time averaging
     * 4. Averaged field storage (CPU)
     */

    // Initialize Dirac field at defect location
    initializeDiracField(x0, y0, sigma, 1.0f);

    // Initialize averaged field storage
    size_t total = _Nx * _Ny;
    _theta_avg.resize(total);
    _R_avg.resize(total);

    // Initialize accumulators to current state
    std::vector<float> theta_init = getPhaseField();
    std::vector<float> R_init = getSyncField();
    _theta_avg = theta_init;
    _R_avg = R_init;

    // Upload initial Dirac density to GPU (frozen for N steps)
    auto psi_density = getDiracDensity();
    if (_spinor_density_buffer != VK_NULL_HANDLE && _spinor_density_memory != VK_NULL_HANDLE) {
        _bufferManager->uploadData(_spinor_density_memory, psi_density.data(),
                                    psi_density.size() * sizeof(float));
    }

    // Reset operator splitting state
    _substep_count = 0;
    if (_substep_ratio == 0) {
        _substep_ratio = 10;  // Default value
    }

    std::cout << "[SMFTEngine] Hybrid system initialized with N=" << _substep_ratio
              << " substeps" << std::endl;
}

void SMFTEngine::initializeDiracPlaneWave(float kx, float ky) {
    if (!_dirac_evolution) {
        _dirac_evolution = new DiracEvolution(_Nx, _Ny);
    }
    _dirac_evolution->initializePlaneWave(kx, ky);
    _dirac_initialized = true;
}

std::vector<std::complex<double>> SMFTEngine::getSpinorField() const {
    if (!_dirac_initialized || !_dirac_evolution) {
        return std::vector<std::complex<double>>(4 * _Nx * _Ny, std::complex<double>(0,0));
    }
    return _dirac_evolution->getSpinorField();
}

// ============================================================================
// Klein-Gordon Evolution Methods (Phase 2.5A)
// ============================================================================

void SMFTEngine::initializeKleinGordonField(float x0, float y0, float sigma, float amplitude) {
    /**
     * Initialize Klein-Gordon scalar field with Gaussian wavepacket
     * Uses split-operator KleinGordonEvolution class
     */

    // Create KleinGordonEvolution instance if not exists
    if (!_kg_evolution) {
        _kg_evolution = new KleinGordonEvolution(_Nx, _Ny);
    }

    // Initialize with Gaussian wavepacket
    _kg_evolution->initialize(x0, y0, sigma);

    _kg_initialized = true;

    // Log initialization
    if (_nova) {
        float norm = _kg_evolution->getNorm();
        std::cout << "[SMFTEngine] Klein-Gordon field initialized: "
                  << "center=(" << x0 << "," << y0 << "), "
                  << "sigma=" << sigma << ", "
                  << "norm=" << norm << std::endl;
    }
}

void SMFTEngine::initializeBoostedKleinGordonField(float x0, float y0, float sigma,
                                                   float vx, float vy, float R_bg) {
    /**
     * Initialize Klein-Gordon scalar field with boosted Gaussian wavepacket
     * Uses SMFT::initializeBoostedGaussian helper function (Phase 2.5A)
     */

    // Create KleinGordonEvolution instance if not exists
    if (!_kg_evolution) {
        _kg_evolution = new KleinGordonEvolution(_Nx, _Ny);
    }

    // Initialize with boosted Gaussian wavepacket using SMFT helper
    SMFT::initializeBoostedGaussian(*_kg_evolution, x0, y0, sigma,
                                   vx, vy, _Delta, R_bg);

    _kg_initialized = true;

    // Log initialization
    if (_nova) {
        float norm = _kg_evolution->getNorm();
        float x_mean, y_mean;
        _kg_evolution->getCenterOfMass(x_mean, y_mean);
        std::cout << "[SMFTEngine] Boosted Klein-Gordon field initialized: "
                  << "center=(" << x0 << "," << y0 << "), "
                  << "sigma=" << sigma << ", "
                  << "boost=(" << vx << "," << vy << ")c, "
                  << "norm=" << norm << ", "
                  << "<r>=(" << x_mean << "," << y_mean << ")" << std::endl;
    }
}

void SMFTEngine::initializeKleinGordonPlaneWave(float kx, float ky) {
    /**
     * Initialize Klein-Gordon scalar field with plane wave
     * Used for dispersion relation analysis
     */
    if (!_kg_evolution) {
        _kg_evolution = new KleinGordonEvolution(_Nx, _Ny);
    }
    _kg_evolution->initializePlaneWave(kx, ky);
    _kg_initialized = true;

    std::cout << "[SMFTEngine] Klein-Gordon plane wave initialized: k=("
              << kx << ", " << ky << ")" << std::endl;
}

const KleinGordonEvolution* SMFTEngine::getKleinGordonEvolution() const {
    return _kg_evolution;
}

KleinGordonEvolution* SMFTEngine::getKleinGordonEvolutionNonConst() {
    return _kg_evolution;
}

// ============================================================================
// TWO-PARTICLE SYSTEM METHODS (Test 3.4: Antiparticle Separation)
// ============================================================================

void SMFTEngine::initializeTwoParticleSystem(float x1, float y1, float x2, float y2, float sigma) {
    /**
     * Initialize particle + antiparticle system for Test 3.4
     *
     * Physics:
     * - Particle: DiracEvolution with β_sign = +1
     * - Antiparticle: DiracEvolution with β_sign = -1
     * - Both evolve in same mass field m(x) = Δ·R(x)
     * - Opposite force signs create spatial separation
     */

    // Create particle DiracEvolution (β = +1)
    if (!_dirac_evolution) {
        _dirac_evolution = new DiracEvolution(_Nx, _Ny, +1.0f);  // particle
    }
    _dirac_evolution->initialize(x1, y1, sigma);

    // Create antiparticle DiracEvolution (β = -1)
    if (!_dirac_antiparticle) {
        _dirac_antiparticle = new DiracEvolution(_Nx, _Ny, -1.0f);  // antiparticle
    }
    _dirac_antiparticle->initialize(x2, y2, sigma);

    _dirac_initialized = true;
    _two_particle_mode = true;

    // Log initialization
    if (_nova) {
        float norm1 = _dirac_evolution->getNorm();
        float norm2 = _dirac_antiparticle->getNorm();
        std::cout << "[SMFTEngine] Two-particle system initialized:" << std::endl;
        std::cout << "  Particle (β=+1): center=(" << x1 << "," << y1 << "), "
                  << "sigma=" << sigma << ", norm=" << norm1 << std::endl;
        std::cout << "  Antiparticle (β=-1): center=(" << x2 << "," << y2 << "), "
                  << "sigma=" << sigma << ", norm=" << norm2 << std::endl;
    }
}

const DiracEvolution* SMFTEngine::getAntiparticleEvolution() const {
    return _dirac_antiparticle;
}

DiracEvolution* SMFTEngine::getAntiparticleEvolutionNonConst() {
    return _dirac_antiparticle;
}

// ==============================================================================
// Phase 4 Test 4.2: R-Field Time Derivative
// ==============================================================================

std::vector<float> SMFTEngine::getRFieldDerivative() const {
    /**
     * Compute ∂R/∂t using centered finite differences from R-field history
     *
     * Ring buffer layout:
     * - idx_minus: R(t - dt)
     * - idx_center: R(t)
     * - idx_plus: R(t + dt)
     *
     * Centered difference: ∂R/∂t ≈ (R(t+dt) - R(t-dt)) / (2·dt)
     *
     * For first 2 timesteps (history incomplete), use forward difference:
     * ∂R/∂t ≈ (R(t) - R(t-dt)) / dt
     */
    std::vector<float> dR_dt(_Nx * _Ny, 0.0f);

    if (_last_dt <= 0.0f) {
        // No timestep information, return zeros
        return dR_dt;
    }

    // Determine which history slots are valid
    // Ring buffer: [0] = t-dt, [1] = t, [2] = t+dt (wraps around)
    int idx_center = (_R_history_index + 2) % 3;  // Most recent complete step
    int idx_minus = (_R_history_index + 1) % 3;   // Previous step (t-dt)
    int idx_plus = _R_history_index;              // Next step (not yet computed)

    // Check if we have enough history for centered difference
    bool has_full_history = true;
    for (int i = 0; i < 3; i++) {
        if (_R_history[i].empty()) {
            has_full_history = false;
            break;
        }
    }

    if (!has_full_history) {
        // Use forward difference: ∂R/∂t ≈ (R_center - R_minus) / dt
        const auto& R_minus = _R_history[idx_minus];
        const auto& R_center = _R_history[idx_center];

        if (!R_minus.empty() && !R_center.empty()) {
            for (size_t i = 0; i < dR_dt.size(); i++) {
                dR_dt[i] = (R_center[i] - R_minus[i]) / _last_dt;
            }
        }
    } else {
        // Use centered difference: ∂R/∂t ≈ (R_plus - R_minus) / (2·dt)
        const auto& R_minus = _R_history[idx_minus];
        const auto& R_plus = _R_history[idx_plus];

        for (size_t i = 0; i < dR_dt.size(); i++) {
            dR_dt[i] = (R_plus[i] - R_minus[i]) / (2.0f * _last_dt);
        }
    }

    return dR_dt;
}

// ============================================================================
// EM FIELD GPU COMPUTATION (Phase 5 - Sprint 3 Step 4)
// ============================================================================

void SMFTEngine::initEMBuffers() {
    /**
     * Initialize GPU buffers for EM field computation
     * Allocates device-local storage for field data and host-visible buffers for readback
     */
    if (!_nova || !_bufferManager) {
        std::cerr << "[SMFTEngine::initEMBuffers] ERROR: Nova or BufferManager not initialized" << std::endl;
        return;
    }

    size_t field_size = _Nx * _Ny * sizeof(float);

    try {
        // Allocate field buffers (DEVICE_LOCAL for GPU computation)
        auto phi_pair = _bufferManager->createStorageBuffer(field_size);
        _phi_buffer = phi_pair.first;
        _phi_memory = phi_pair.second;

        auto A_x_pair = _bufferManager->createStorageBuffer(field_size);
        _A_x_buffer = A_x_pair.first;
        _A_x_memory = A_x_pair.second;

        auto A_y_pair = _bufferManager->createStorageBuffer(field_size);
        _A_y_buffer = A_y_pair.first;
        _A_y_memory = A_y_pair.second;

        auto E_x_pair = _bufferManager->createStorageBuffer(field_size);
        _E_x_buffer = E_x_pair.first;
        _E_x_memory = E_x_pair.second;

        auto E_y_pair = _bufferManager->createStorageBuffer(field_size);
        _E_y_buffer = E_y_pair.first;
        _E_y_memory = E_y_pair.second;

        auto B_z_pair = _bufferManager->createStorageBuffer(field_size);
        _B_z_buffer = B_z_pair.first;
        _B_z_memory = B_z_pair.second;

        // Allocate theta_previous buffer for time derivatives
        auto theta_prev_pair = _bufferManager->createStorageBuffer(field_size);
        _theta_previous_buffer = theta_prev_pair.first;
        _theta_previous_memory = theta_prev_pair.second;

        // Allocate energy buffer (single float, HOST_VISIBLE for readback)
        // Note: createStorageBuffer creates DEVICE_LOCAL by default
        // For readback, we'll use downloadData which handles staging
        auto energy_pair = _bufferManager->createStorageBuffer(sizeof(float));
        _em_energy_buffer = energy_pair.first;
        _em_energy_memory = energy_pair.second;

        // Initialize energy buffer to 0
        float zero = 0.0f;
        _bufferManager->uploadData(_em_energy_memory, &zero, sizeof(float));

        // Allocate params buffer (64 bytes uniform, HOST_VISIBLE for updates)
        auto params_pair = _bufferManager->createStorageBuffer(64);
        _em_params_buffer = params_pair.first;
        _em_params_memory = params_pair.second;

        std::cout << "[SMFTEngine] EM field GPU buffers initialized ("
                  << (field_size / 1024) << " KB per field)" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "[SMFTEngine::initEMBuffers] Failed to create EM buffers: " << e.what() << std::endl;
        std::cerr << "  Falling back to CPU EM computation" << std::endl;
        // Mark EM GPU as unavailable (field computation will use CPU fallback)
        _em_potentials_pipeline = VK_NULL_HANDLE;
    }
}

void SMFTEngine::createEMPipelines() {
    /**
     * Create GPU compute pipelines for EM field computation
     * Pipelines: 1) computeEMPotentials, 2) computeFieldStrengths, 3) reduceFieldEnergy
     */
    if (!_nova || !_nova->initialized) {
        std::cerr << "[SMFTEngine::createEMPipelines] ERROR: Nova not initialized" << std::endl;
        return;
    }

    if (_phi_buffer == VK_NULL_HANDLE || _theta_buffer == VK_NULL_HANDLE) {
        std::cerr << "[SMFTEngine::createEMPipelines] ERROR: EM buffers not initialized, skipping pipeline creation" << std::endl;
        return;
    }

    try {
        VkDevice device = _nova->_architect->logical_device;

        // ========================================================================
        // Pipeline 1: computeEMPotentials (θ → A_μ) with conjugate product method
        // ========================================================================
        {
            // Create descriptor set layout
            std::vector<VkDescriptorSetLayoutBinding> bindings(7);  // Updated from 6 to 7 for R_field
            bindings[0] = {0, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr};  // theta_current
            bindings[1] = {1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr};  // theta_previous
            bindings[2] = {2, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr};  // R_field (NEW)
            bindings[3] = {3, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr};  // phi
            bindings[4] = {4, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr};  // A_x
            bindings[5] = {5, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr};  // A_y
            bindings[6] = {6, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr};  // params

            _em_potentials_desc_layout = _descriptorManager->createDescriptorSetLayout(bindings);

            // Create pipeline layout
            VkPipelineLayoutCreateInfo layoutInfo{};
            layoutInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
            layoutInfo.setLayoutCount = 1;
            layoutInfo.pSetLayouts = &_em_potentials_desc_layout;
            layoutInfo.pushConstantRangeCount = 0;

            if (vkCreatePipelineLayout(device, &layoutInfo, nullptr, &_em_potentials_layout) != VK_SUCCESS) {
                throw std::runtime_error("Failed to create EM potentials pipeline layout");
            }

            // Load shader and create pipeline
            _em_potentials_pipeline = _pipelineFactory->createEMPotentialsPipeline(
                "shaders/smft/computeEMPotentials.comp.spv",
                _em_potentials_layout
            );

            if (_em_potentials_pipeline == VK_NULL_HANDLE) {
                throw std::runtime_error("Failed to create EM potentials compute pipeline");
            }

            // Allocate and update descriptor set
            _em_potentials_desc_set = _descriptorManager->allocateDescriptorSet(_descriptor_pool, _em_potentials_desc_layout);

            VkDescriptorBufferInfo bufferInfos[7] = {  // Updated from 6 to 7
                {_theta_buffer, 0, VK_WHOLE_SIZE},
                {_theta_previous_buffer, 0, VK_WHOLE_SIZE},
                {_R_field_buffer, 0, VK_WHOLE_SIZE},  // NEW: R_field for conjugate product method
                {_phi_buffer, 0, VK_WHOLE_SIZE},
                {_A_x_buffer, 0, VK_WHOLE_SIZE},
                {_A_y_buffer, 0, VK_WHOLE_SIZE},
                {_em_params_buffer, 0, VK_WHOLE_SIZE}
            };

            std::vector<VkWriteDescriptorSet> writes(7);  // Updated from 6 to 7
            for (int i = 0; i < 7; i++) {
                writes[i].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
                writes[i].dstSet = _em_potentials_desc_set;
                writes[i].dstBinding = i;
                writes[i].dstArrayElement = 0;
                writes[i].descriptorCount = 1;
                writes[i].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
                writes[i].pBufferInfo = &bufferInfos[i];
            }

            vkUpdateDescriptorSets(device, writes.size(), writes.data(), 0, nullptr);
        }

        // ========================================================================
        // Pipeline 2: computeFieldStrengths (A_μ → E, B)
        // ========================================================================
        {
            std::vector<VkDescriptorSetLayoutBinding> bindings(6);
            bindings[0] = {0, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr};  // phi
            bindings[1] = {1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr};  // A_x
            bindings[2] = {2, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr};  // A_y
            bindings[3] = {3, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr};  // E_x
            bindings[4] = {4, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr};  // E_y
            bindings[5] = {5, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr};  // B_z

            _em_field_strengths_desc_layout = _descriptorManager->createDescriptorSetLayout(bindings);

            VkPipelineLayoutCreateInfo layoutInfo{};
            layoutInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
            layoutInfo.setLayoutCount = 1;
            layoutInfo.pSetLayouts = &_em_field_strengths_desc_layout;

            if (vkCreatePipelineLayout(device, &layoutInfo, nullptr, &_em_field_strengths_layout) != VK_SUCCESS) {
                throw std::runtime_error("Failed to create field strengths pipeline layout");
            }

            _em_field_strengths_pipeline = _pipelineFactory->createEMFieldStrengthsPipeline(
                "shaders/smft/computeFieldStrengths.comp.spv",
                _em_field_strengths_layout
            );

            if (_em_field_strengths_pipeline == VK_NULL_HANDLE) {
                throw std::runtime_error("Failed to create field strengths pipeline");
            }

            _em_field_strengths_desc_set = _descriptorManager->allocateDescriptorSet(_descriptor_pool, _em_field_strengths_desc_layout);

            VkDescriptorBufferInfo bufferInfos[6] = {
                {_phi_buffer, 0, VK_WHOLE_SIZE},
                {_A_x_buffer, 0, VK_WHOLE_SIZE},
                {_A_y_buffer, 0, VK_WHOLE_SIZE},
                {_E_x_buffer, 0, VK_WHOLE_SIZE},
                {_E_y_buffer, 0, VK_WHOLE_SIZE},
                {_B_z_buffer, 0, VK_WHOLE_SIZE}
            };

            std::vector<VkWriteDescriptorSet> writes(6);
            for (int i = 0; i < 6; i++) {
                writes[i].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
                writes[i].dstSet = _em_field_strengths_desc_set;
                writes[i].dstBinding = i;
                writes[i].dstArrayElement = 0;
                writes[i].descriptorCount = 1;
                writes[i].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
                writes[i].pBufferInfo = &bufferInfos[i];
            }

            vkUpdateDescriptorSets(device, writes.size(), writes.data(), 0, nullptr);
        }

        // ========================================================================
        // Pipeline 3: reduceFieldEnergy (E, B → scalar energy)
        // ========================================================================
        {
            std::vector<VkDescriptorSetLayoutBinding> bindings(4);
            bindings[0] = {0, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr};  // E_x
            bindings[1] = {1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr};  // E_y
            bindings[2] = {2, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr};  // B_z
            bindings[3] = {3, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr};  // energy_out

            _em_reduce_energy_desc_layout = _descriptorManager->createDescriptorSetLayout(bindings);

            VkPipelineLayoutCreateInfo layoutInfo{};
            layoutInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
            layoutInfo.setLayoutCount = 1;
            layoutInfo.pSetLayouts = &_em_reduce_energy_desc_layout;

            if (vkCreatePipelineLayout(device, &layoutInfo, nullptr, &_em_reduce_energy_layout) != VK_SUCCESS) {
                throw std::runtime_error("Failed to create reduce energy pipeline layout");
            }

            _em_reduce_energy_pipeline = _pipelineFactory->createEMReduceEnergyPipeline(
                "shaders/smft/reduceFieldEnergy.comp.spv",
                _em_reduce_energy_layout
            );

            if (_em_reduce_energy_pipeline == VK_NULL_HANDLE) {
                throw std::runtime_error("Failed to create reduce energy pipeline");
            }

            _em_reduce_energy_desc_set = _descriptorManager->allocateDescriptorSet(_descriptor_pool, _em_reduce_energy_desc_layout);

            VkDescriptorBufferInfo bufferInfos[4] = {
                {_E_x_buffer, 0, VK_WHOLE_SIZE},
                {_E_y_buffer, 0, VK_WHOLE_SIZE},
                {_B_z_buffer, 0, VK_WHOLE_SIZE},
                {_em_energy_buffer, 0, VK_WHOLE_SIZE}
            };

            std::vector<VkWriteDescriptorSet> writes(4);
            for (int i = 0; i < 4; i++) {
                writes[i].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
                writes[i].dstSet = _em_reduce_energy_desc_set;
                writes[i].dstBinding = i;
                writes[i].dstArrayElement = 0;
                writes[i].descriptorCount = 1;
                writes[i].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
                writes[i].pBufferInfo = &bufferInfos[i];
            }

            vkUpdateDescriptorSets(device, writes.size(), writes.data(), 0, nullptr);
        }

        std::cout << "[SMFTEngine] EM field GPU pipelines created successfully" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "[SMFTEngine::createEMPipelines] Failed to create EM pipelines: " << e.what() << std::endl;
        std::cerr << "  Falling back to CPU EM computation" << std::endl;
        _em_potentials_pipeline = VK_NULL_HANDLE;
    }
}

void SMFTEngine::computeEMFieldsGPU() {
    /**
     * Execute GPU EM field computation pipeline
     * Dispatch sequence: potentials → field strengths → energy reduction
     */
    if (_em_potentials_pipeline == VK_NULL_HANDLE || !_compute) {
        return;  // GPU EM not available, caller will use CPU fallback
    }

    // Performance timing (Step 5 Test 2)
    auto start_time = std::chrono::high_resolution_clock::now();

    // Update uniform params (Nx, Ny, dx, dy, dt)
    struct EMParams {
        uint32_t Nx, Ny;
        float dx, dy, dt;
        float padding[11];  // Pad to 64 bytes for alignment
    } params;

    params.Nx = _Nx;
    params.Ny = _Ny;
    params.dx = 1.0f;  // Grid spacing in natural units
    params.dy = 1.0f;
    params.dt = 0.01f;  // Will be overwritten by actual dt during stepWithDirac
    std::memset(params.padding, 0, sizeof(params.padding));

    _bufferManager->uploadData(_em_params_memory, &params, sizeof(EMParams));

    // Reset energy buffer to 0 before accumulation
    float zero = 0.0f;
    _bufferManager->uploadData(_em_energy_memory, &zero, sizeof(float));

    // Calculate workgroup counts (16x16 local size)
    uint32_t num_groups_x = (_Nx + 15) / 16;
    uint32_t num_groups_y = (_Ny + 15) / 16;

    // Begin compute batch
    if (!_compute->beginBatch()) {
        std::cerr << "[SMFTEngine::computeEMFieldsGPU] Failed to begin compute batch" << std::endl;
        return;
    }

    // Kernel 1: Compute A_μ = ∂_μ θ (potentials from phase gradients)
    _compute->dispatchEMPotentials(
        _em_potentials_pipeline,
        _em_potentials_layout,
        _em_potentials_desc_set,
        nullptr, 0,  // No push constants
        num_groups_x, num_groups_y
    );

    _compute->insertMemoryBarrier();

    // Kernel 2: Compute E, B from A_μ (field strengths from potentials)
    _compute->dispatchEMFieldStrengths(
        _em_field_strengths_pipeline,
        _em_field_strengths_layout,
        _em_field_strengths_desc_set,
        nullptr, 0,
        num_groups_x, num_groups_y
    );

    _compute->insertMemoryBarrier();

    // Kernel 3: Reduce to energy scalar E_field = ∫(E² + B²)/2 dA
    uint32_t total_elements = _Nx * _Ny;
    uint32_t num_groups_1d = (total_elements + 255) / 256;

    _compute->dispatchEMReduceEnergy(
        _em_reduce_energy_pipeline,
        _em_reduce_energy_layout,
        _em_reduce_energy_desc_set,
        nullptr, 0,
        num_groups_1d, 1
    );

    // Submit and wait for completion
    if (!_compute->submitBatch(true)) {
        std::cerr << "[SMFTEngine::computeEMFieldsGPU] Failed to submit compute batch" << std::endl;
    }

    // Performance timing report (Step 5 Test 2)
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration_us = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    static int call_count = 0;
    static long total_us = 0;
    call_count++;
    total_us += duration_us;
    if (call_count == 1 || call_count % 50 == 0) {
        std::cout << "[GPU EM Compute] " << call_count << " calls | Last: " << duration_us << " μs | Avg: "
                  << (total_us / call_count) << " μs | Grid: " << _Nx << "×" << _Ny << std::endl;
    }
}

float SMFTEngine::downloadEMEnergy() {
    /**
     * Download EM field energy from GPU
     * Returns total field energy: E_field = ∫(E² + B²)/2 dA
     */
    if (_em_energy_buffer == VK_NULL_HANDLE || !_bufferManager) {
        return 0.0f;
    }

    float energy = 0.0f;
    _bufferManager->downloadData(_em_energy_memory, &energy, sizeof(float));
    return energy;
}

void SMFTEngine::cleanupEMResources() {
    /**
     * Cleanup EM-specific Vulkan resources
     * Called by destroyResources()
     *
     * NOTE: Component managers handle cleanup of resources they track:
     * - Pipelines: Cleaned up by _pipelineFactory
     * - Buffers/Memory: Cleaned up by _bufferManager
     * - Descriptor layouts: Cleaned up by _descriptorManager (if created through it)
     *
     * We only clean up resources created directly through Vulkan API
     */
    if (!_nova || !_nova->_architect) {
        return;
    }

    VkDevice device = _nova->_architect->logical_device;
    if (device == VK_NULL_HANDLE) {
        return;
    }

    // Reset pipeline handles (pipelines are managed by _pipelineFactory)
    _em_potentials_pipeline = VK_NULL_HANDLE;
    _em_field_strengths_pipeline = VK_NULL_HANDLE;
    _em_reduce_energy_pipeline = VK_NULL_HANDLE;

    // Destroy pipeline layouts (created directly, not through a manager)
    if (_em_potentials_layout != VK_NULL_HANDLE) {
        vkDestroyPipelineLayout(device, _em_potentials_layout, nullptr);
        _em_potentials_layout = VK_NULL_HANDLE;
    }
    if (_em_field_strengths_layout != VK_NULL_HANDLE) {
        vkDestroyPipelineLayout(device, _em_field_strengths_layout, nullptr);
        _em_field_strengths_layout = VK_NULL_HANDLE;
    }
    if (_em_reduce_energy_layout != VK_NULL_HANDLE) {
        vkDestroyPipelineLayout(device, _em_reduce_energy_layout, nullptr);
        _em_reduce_energy_layout = VK_NULL_HANDLE;
    }

    // Reset descriptor layout handles (layouts are managed by _descriptorManager)
    _em_potentials_desc_layout = VK_NULL_HANDLE;
    _em_field_strengths_desc_layout = VK_NULL_HANDLE;
    _em_reduce_energy_desc_layout = VK_NULL_HANDLE;

    // Reset buffer handles (buffers are managed by _bufferManager)
    _phi_buffer = VK_NULL_HANDLE;
    _A_x_buffer = VK_NULL_HANDLE;
    _A_y_buffer = VK_NULL_HANDLE;
    _E_x_buffer = VK_NULL_HANDLE;
    _E_y_buffer = VK_NULL_HANDLE;
    _B_z_buffer = VK_NULL_HANDLE;
    _em_energy_buffer = VK_NULL_HANDLE;
    _em_params_buffer = VK_NULL_HANDLE;
    _theta_previous_buffer = VK_NULL_HANDLE;

    // Reset memory handles
    _phi_memory = VK_NULL_HANDLE;
    _A_x_memory = VK_NULL_HANDLE;
    _A_y_memory = VK_NULL_HANDLE;
    _E_x_memory = VK_NULL_HANDLE;
    _E_y_memory = VK_NULL_HANDLE;
    _B_z_memory = VK_NULL_HANDLE;
    _em_energy_memory = VK_NULL_HANDLE;
    _em_params_memory = VK_NULL_HANDLE;
    _theta_previous_memory = VK_NULL_HANDLE;
}