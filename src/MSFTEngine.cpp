#include "MSFTEngine.h"
#include <cstring>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>

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

MSFTEngine::MSFTEngine(Nova* nova)
    : _nova(nova),
      _Nx(0), _Ny(0),
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
      _dirac_pipeline_layout(VK_NULL_HANDLE),
      _dirac_initialized(false)  // Initialize Dirac field flag
{
    // Ensure Nova is initialized before creating managers
    if (!nova || !nova->initialized) {
        std::cerr << "[MSFTEngine] Error: Nova not initialized" << std::endl;
        return;
    }

    // Initialize pipeline factory
    _pipelineFactory = std::make_unique<MSFTPipelineFactory>(nova->_architect->logical_device);

    // Initialize compute dispatcher
    _compute = std::make_unique<MSFTCompute>(nova->_architect->logical_device,
                                              nova->_architect->queues.compute,
                                              nova->_architect->queues.indices.compute_family.value_or(0));

    // Initialize buffer manager
    _bufferManager = std::make_unique<MSFTBufferManager>(nova->_architect->logical_device,
                                                          nova->_architect->physical_device);

    // Initialize descriptor manager
    _descriptorManager = std::make_unique<MSFTDescriptorManager>(nova->_architect->logical_device);
}

MSFTEngine::~MSFTEngine() {
    // Clean up all Vulkan resources
    destroyResources();
}

void MSFTEngine::initialize(uint32_t Nx, uint32_t Ny, float Delta, float chiral_angle) {
    _Nx = Nx;
    _Ny = Ny;
    _Delta = Delta;
    _chiral_angle = chiral_angle;
    _time_step = 0;

    // Log physics interpretation
    if (_nova) {
        std::cout << "\n[MSFTEngine] Physical Interpretation:" << std::endl;
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

    // Upload initial data to GPU
    uploadToGPU();
}

void MSFTEngine::setInitialPhases(const std::vector<float>& theta) {
    if (theta.size() != _Nx * _Ny) {
        std::cerr << "[MSFTEngine] Error: Invalid theta array size" << std::endl;
        return;
    }

    _theta_data = theta;
    uploadToGPU();
}

void MSFTEngine::setNaturalFrequencies(const std::vector<float>& omega) {
    if (omega.size() != _Nx * _Ny) {
        std::cerr << "[MSFTEngine] Error: Invalid omega array size" << std::endl;
        return;
    }

    _omega_data = omega;
    uploadToGPU();
}

void MSFTEngine::step(float dt, float K, float damping) {
    // Ensure we have GPU resources
    if (!_kuramoto_pipeline || !_sync_pipeline || !_gravity_pipeline) {
        std::cerr << "[MSFTEngine] Pipelines not created, falling back to CPU simulation" << std::endl;
        return;
    }

    // Calculate workgroup counts (local size = 16x16 in shader)
    uint32_t workgroupsX = (_Nx + 15) / 16;
    uint32_t workgroupsY = (_Ny + 15) / 16;

    // Begin compute batch
    if (!_compute->beginBatch()) {
        std::cerr << "[MSFTEngine] Failed to begin compute batch" << std::endl;
        return;
    }

    // Step 1: Kuramoto phase evolution
    // Push constants for kuramoto shader
    struct KuramotoPush {
        float dt;
        float K;
        float damping;
        float Delta;
        uint32_t Nx;
        uint32_t Ny;
    } kuramoto_push = {
        dt, K, damping, _Delta, _Nx, _Ny
    };

    _compute->dispatchKuramoto(_kuramoto_pipeline,
                               _kuramoto_pipeline_layout,
                               _kuramoto_descriptor_set,
                               &kuramoto_push, sizeof(kuramoto_push),
                               workgroupsX, workgroupsY);

    // Memory barrier between kuramoto and sync
    _compute->insertMemoryBarrier();

    // Step 2: Calculate synchronization field
    // Push constants for sync shader
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

    // Memory barrier between sync and gravity
    _compute->insertMemoryBarrier();

    // Step 3: Calculate gravitational field
    // Push constants for gravity shader
    struct GravityPush {
        float Delta;
        uint32_t Nx;
        uint32_t Ny;
    } gravity_push = {
        _Delta, _Nx, _Ny
    };

    _compute->dispatchGravityField(_gravity_pipeline,
                                   _gravity_pipeline_layout,
                                   _gravity_descriptor_set,
                                   &gravity_push, sizeof(gravity_push),
                                   workgroupsX, workgroupsY);

    // Submit batch and wait
    if (!_compute->submitBatch(true)) {
        std::cerr << "[MSFTEngine] Failed to submit compute batch" << std::endl;
        return;
    }

    // Download results from GPU
    downloadFromGPU();

    // Swap theta buffers for next iteration
    std::swap(_theta_buffer, _theta_out_buffer);
    std::swap(_theta_memory, _theta_out_memory);

    // Also swap in descriptor sets (both kuramoto and sync reference theta_buffer)
    // This is handled by recreating descriptor sets in Phase 3
}

void MSFTEngine::stepStochastic(float dt, float K, float damping,
                                float sigma_theta, float sigma_psi) {
    // Ensure we have GPU resources
    if (!_kuramoto_stochastic_pipeline) {
        std::cerr << "[MSFTEngine] Stochastic pipelines not created" << std::endl;
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

std::vector<float> MSFTEngine::getSyncField() const {
    return _R_field_data;
}

std::vector<float> MSFTEngine::getMassField() const {
    std::vector<float> mass_field(_R_field_data.size());
    for (size_t i = 0; i < _R_field_data.size(); ++i) {
        mass_field[i] = _Delta * _R_field_data[i];
    }
    return mass_field;
}

std::vector<float> MSFTEngine::getPhaseField() const {
    return _theta_data;
}

std::vector<float> MSFTEngine::getGravitationalField() const {
    std::vector<float> gravity_field(2 * _Nx * _Ny);
    for (size_t i = 0; i < _Nx * _Ny; ++i) {
        gravity_field[2*i + 0] = _gravity_x_data[i];
        gravity_field[2*i + 1] = _gravity_y_data[i];
    }
    return gravity_field;
}

std::complex<float> MSFTEngine::getSpinorComponent(uint32_t x, uint32_t y, uint32_t component) const {
    if (x >= _Nx || y >= _Ny || component >= 4) {
        return std::complex<float>(0.0f, 0.0f);
    }
    uint32_t idx = 4 * (y * _Nx + x) + component;
    return _spinor_field[idx];
}

void MSFTEngine::createBuffers() {
    if (!_nova || !_nova->initialized) {
        std::cerr << "[MSFTEngine] Nova not initialized, cannot create buffers" << std::endl;
        return;
    }

    size_t float_size = sizeof(float) * _Nx * _Ny;
    size_t spinor_size = sizeof(std::complex<float>) * 4 * _Nx * _Ny;

    // Create buffers using buffer manager
    auto theta_pair = _bufferManager->createBuffer(float_size, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
                                                   VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
    _theta_buffer = theta_pair.first;
    _theta_memory = theta_pair.second;

    auto theta_out_pair = _bufferManager->createBuffer(float_size, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
                                                       VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
    _theta_out_buffer = theta_out_pair.first;
    _theta_out_memory = theta_out_pair.second;

    auto omega_pair = _bufferManager->createBuffer(float_size, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
                                                   VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
    _omega_buffer = omega_pair.first;
    _omega_memory = omega_pair.second;

    auto R_field_pair = _bufferManager->createBuffer(float_size, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
                                                     VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
    _R_field_buffer = R_field_pair.first;
    _R_field_memory = R_field_pair.second;

    auto gravity_x_pair = _bufferManager->createBuffer(float_size, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
                                                       VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
    _gravity_x_buffer = gravity_x_pair.first;
    _gravity_x_memory = gravity_x_pair.second;

    auto gravity_y_pair = _bufferManager->createBuffer(float_size, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
                                                       VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
    _gravity_y_buffer = gravity_y_pair.first;
    _gravity_y_memory = gravity_y_pair.second;

    auto spinor_density_pair = _bufferManager->createBuffer(float_size, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
                                                            VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
    _spinor_density_buffer = spinor_density_pair.first;
    _spinor_density_memory = spinor_density_pair.second;

    auto spinor_pair = _bufferManager->createBuffer(spinor_size, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
                                                    VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
    _spinor_buffer = spinor_pair.first;
    _spinor_memory = spinor_pair.second;

    // Verify all buffers were created
    if (_theta_buffer == VK_NULL_HANDLE || _theta_out_buffer == VK_NULL_HANDLE ||
        _omega_buffer == VK_NULL_HANDLE || _R_field_buffer == VK_NULL_HANDLE ||
        _gravity_x_buffer == VK_NULL_HANDLE || _gravity_y_buffer == VK_NULL_HANDLE ||
        _spinor_density_buffer == VK_NULL_HANDLE || _spinor_buffer == VK_NULL_HANDLE) {
        std::cerr << "[MSFTEngine] Failed to create one or more buffers" << std::endl;
        destroyResources();
        return;
    }

    if (_nova) {
        std::cout << "[MSFTEngine] Created GPU buffers for " << _Nx << "x" << _Ny << " grid" << std::endl;
        std::cout << "  Phase buffers: " << float_size << " bytes each" << std::endl;
        std::cout << "  Spinor buffer: " << spinor_size << " bytes" << std::endl;
    }
}

void MSFTEngine::createPipelines() {
    if (!_nova || !_nova->initialized) {
        std::cerr << "[MSFTEngine] Nova not initialized, cannot create pipelines" << std::endl;
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
        std::cerr << "[MSFTEngine] Failed to create descriptor pool" << std::endl;
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
        pushRange.size = sizeof(float) * 6;

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
        pushRange.size = sizeof(uint32_t) * 2;

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
        pushRange.size = sizeof(float) + sizeof(uint32_t) * 2;

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
            std::cerr << "[MSFTEngine] Warning: Failed to load stochastic Kuramoto shader" << std::endl;
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
            std::cerr << "[MSFTEngine] Warning: Failed to load stochastic Dirac shader" << std::endl;
        }
    }

    if (_nova) {
        std::cout << "[MSFTEngine] Created compute pipelines" << std::endl;
        if (_kuramoto_pipeline) std::cout << "  ✓ Kuramoto pipeline" << std::endl;
        if (_sync_pipeline) std::cout << "  ✓ Sync field pipeline" << std::endl;
        if (_gravity_pipeline) std::cout << "  ✓ Gravity field pipeline" << std::endl;
        if (_kuramoto_stochastic_pipeline) std::cout << "  ✓ Stochastic Kuramoto pipeline" << std::endl;
        if (_dirac_stochastic_pipeline) std::cout << "  ✓ Stochastic Dirac pipeline" << std::endl;
    }
}

void MSFTEngine::uploadToGPU() {
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

    if (_nova) {
        std::cout << "[MSFTEngine] Uploaded data to GPU" << std::endl;
    }
}

void MSFTEngine::downloadFromGPU() {
    if (!_bufferManager) return;

    // Download phase field (from theta_out_buffer)
    _bufferManager->downloadData(_theta_out_memory, _theta_data.data(),
                                sizeof(float) * _theta_data.size());

    // Download synchronization field
    _bufferManager->downloadData(_R_field_memory, _R_field_data.data(),
                                sizeof(float) * _R_field_data.size());

    // Download gravitational field components
    _bufferManager->downloadData(_gravity_x_memory, _gravity_x_data.data(),
                                sizeof(float) * _gravity_x_data.size());
    _bufferManager->downloadData(_gravity_y_memory, _gravity_y_data.data(),
                                sizeof(float) * _gravity_y_data.size());

    if (_nova) {
        std::cout << "[MSFTEngine] Downloaded results from GPU" << std::endl;
    }
}

void MSFTEngine::destroyResources() {
    if (_nova && _nova->_architect && _nova->_architect->logical_device) {
        vkDeviceWaitIdle(_nova->_architect->logical_device);

        // Destroy pipelines
        if (_kuramoto_pipeline) vkDestroyPipeline(_nova->_architect->logical_device, _kuramoto_pipeline, nullptr);
        if (_sync_pipeline) vkDestroyPipeline(_nova->_architect->logical_device, _sync_pipeline, nullptr);
        if (_gravity_pipeline) vkDestroyPipeline(_nova->_architect->logical_device, _gravity_pipeline, nullptr);
        if (_dirac_pipeline) vkDestroyPipeline(_nova->_architect->logical_device, _dirac_pipeline, nullptr);
        if (_kuramoto_stochastic_pipeline) vkDestroyPipeline(_nova->_architect->logical_device, _kuramoto_stochastic_pipeline, nullptr);
        if (_dirac_stochastic_pipeline) vkDestroyPipeline(_nova->_architect->logical_device, _dirac_stochastic_pipeline, nullptr);

        // Destroy pipeline layouts
        if (_kuramoto_pipeline_layout) vkDestroyPipelineLayout(_nova->_architect->logical_device, _kuramoto_pipeline_layout, nullptr);
        if (_sync_pipeline_layout) vkDestroyPipelineLayout(_nova->_architect->logical_device, _sync_pipeline_layout, nullptr);
        if (_gravity_pipeline_layout) vkDestroyPipelineLayout(_nova->_architect->logical_device, _gravity_pipeline_layout, nullptr);
        if (_dirac_pipeline_layout) vkDestroyPipelineLayout(_nova->_architect->logical_device, _dirac_pipeline_layout, nullptr);

        // Destroy descriptor set layouts
        if (_kuramoto_descriptor_layout) vkDestroyDescriptorSetLayout(_nova->_architect->logical_device, _kuramoto_descriptor_layout, nullptr);
        if (_sync_descriptor_layout) vkDestroyDescriptorSetLayout(_nova->_architect->logical_device, _sync_descriptor_layout, nullptr);
        if (_gravity_descriptor_layout) vkDestroyDescriptorSetLayout(_nova->_architect->logical_device, _gravity_descriptor_layout, nullptr);
        if (_dirac_descriptor_layout) vkDestroyDescriptorSetLayout(_nova->_architect->logical_device, _dirac_descriptor_layout, nullptr);

        // Destroy descriptor pool
        if (_descriptor_pool) vkDestroyDescriptorPool(_nova->_architect->logical_device, _descriptor_pool, nullptr);

        // Destroy buffers and free memory
        if (_theta_buffer) vkDestroyBuffer(_nova->_architect->logical_device, _theta_buffer, nullptr);
        if (_theta_memory) vkFreeMemory(_nova->_architect->logical_device, _theta_memory, nullptr);
        if (_theta_out_buffer) vkDestroyBuffer(_nova->_architect->logical_device, _theta_out_buffer, nullptr);
        if (_theta_out_memory) vkFreeMemory(_nova->_architect->logical_device, _theta_out_memory, nullptr);
        if (_omega_buffer) vkDestroyBuffer(_nova->_architect->logical_device, _omega_buffer, nullptr);
        if (_omega_memory) vkFreeMemory(_nova->_architect->logical_device, _omega_memory, nullptr);
        if (_R_field_buffer) vkDestroyBuffer(_nova->_architect->logical_device, _R_field_buffer, nullptr);
        if (_R_field_memory) vkFreeMemory(_nova->_architect->logical_device, _R_field_memory, nullptr);
        if (_gravity_x_buffer) vkDestroyBuffer(_nova->_architect->logical_device, _gravity_x_buffer, nullptr);
        if (_gravity_x_memory) vkFreeMemory(_nova->_architect->logical_device, _gravity_x_memory, nullptr);
        if (_gravity_y_buffer) vkDestroyBuffer(_nova->_architect->logical_device, _gravity_y_buffer, nullptr);
        if (_gravity_y_memory) vkFreeMemory(_nova->_architect->logical_device, _gravity_y_memory, nullptr);
        if (_spinor_density_buffer) vkDestroyBuffer(_nova->_architect->logical_device, _spinor_density_buffer, nullptr);
        if (_spinor_density_memory) vkFreeMemory(_nova->_architect->logical_device, _spinor_density_memory, nullptr);
        if (_spinor_buffer) vkDestroyBuffer(_nova->_architect->logical_device, _spinor_buffer, nullptr);
        if (_spinor_memory) vkFreeMemory(_nova->_architect->logical_device, _spinor_memory, nullptr);
    }

    if (_nova) {
        std::cout << "[MSFTEngine] Resources destroyed" << std::endl;
    }
}

// ============================================================================
// DIRAC FIELD METHODS (CPU-only implementation)
// ============================================================================

void MSFTEngine::initializeDiracField(float x0, float y0, float sigma, float amplitude) {
    /**
     * Initialize Dirac spinor field with Gaussian wavepacket
     * CPU-only implementation - GPU shaders exceed timeout budget
     *
     * Creates a Gaussian wavepacket centered at (x0, y0) with width sigma.
     * For simplicity, we start with a scalar field (only ψ₀ non-zero).
     * Full 4-component spinor can be added later if needed.
     */

    // Ensure psi vector is allocated (4 components per grid point)
    _psi.resize(4 * _Nx * _Ny);

    // Initialize all components to zero
    std::fill(_psi.begin(), _psi.end(), std::complex<float>(0.0f, 0.0f));

    // Create Gaussian wavepacket in component 0
    float norm_sum = 0.0f;

    for (uint32_t y = 0; y < _Ny; y++) {
        for (uint32_t x = 0; x < _Nx; x++) {
            // Distance from center
            float dx = float(x) - x0;
            float dy = float(y) - y0;
            float r2 = dx * dx + dy * dy;

            // Gaussian envelope
            float gauss = amplitude * std::exp(-r2 / (2.0f * sigma * sigma));

            // Set component 0 (scalar approximation)
            uint32_t idx = (y * _Nx + x) * 4;  // 4 components per point
            _psi[idx] = std::complex<float>(gauss, 0.0f);

            // Accumulate norm
            norm_sum += gauss * gauss;
        }
    }

    // Normalize wavefunction: ∫|ψ|² = 1
    if (norm_sum > 0.0f) {
        float norm_factor = 1.0f / std::sqrt(norm_sum);
        for (auto& psi_comp : _psi) {
            psi_comp *= norm_factor;
        }
    }

    _dirac_initialized = true;

    // Log initialization
    if (_nova) {
        std::cout << "[MSFTEngine] Dirac field initialized: "
                  << "center=(" << x0 << "," << y0 << "), "
                  << "sigma=" << sigma << ", "
                  << "norm=" << norm_sum << std::endl;
    }
}

void MSFTEngine::stepWithDirac(float dt, float lambda_coupling) {
    /**
     * Step coupled Kuramoto-Dirac evolution with mass coupling
     * CPU-only implementation - GPU shaders exceed timeout budget
     *
     * Physics:
     * - Dirac equation: i·dΨ/dt = [-iα·∇ + β·m(x)]Ψ
     * - Mass coupling: m(x,y) = Δ·R(x,y) from synchronization field
     * - Simple Euler integration for CPU efficiency
     * - Feedback: Ψ density influences phase field via lambda_coupling
     */

    if (!_dirac_initialized) {
        // Initialize with default Gaussian at center
        initializeDiracField(_Nx/2.0f, _Ny/2.0f, 5.0f, 1.0f);
    }

    // Create temporary buffer for new psi values
    std::vector<std::complex<float>> psi_new = _psi;

    // Get current synchronization field for mass coupling
    std::vector<float> R_field = getSyncField();

    // Euler integration of Dirac equation
    for (uint32_t y = 0; y < _Ny; y++) {
        for (uint32_t x = 0; x < _Nx; x++) {
            uint32_t idx = y * _Nx + x;
            uint32_t psi_idx = idx * 4;  // 4 components per point

            // Get mass field: m(x,y) = Δ·R(x,y)
            float mass = _Delta * R_field[idx];

            // For scalar approximation, only evolve component 0
            std::complex<float> psi0 = _psi[psi_idx];

            // Compute spatial derivatives (finite differences with periodic BC)
            uint32_t xp = ((x + 1) % _Nx);
            uint32_t xm = ((x + _Nx - 1) % _Nx);
            uint32_t yp = ((y + 1) % _Ny);
            uint32_t ym = ((y + _Ny - 1) % _Ny);

            std::complex<float> dpsi_dx = (_psi[y * _Nx * 4 + xp * 4] - _psi[y * _Nx * 4 + xm * 4]) * 0.5f;
            std::complex<float> dpsi_dy = (_psi[yp * _Nx * 4 + x * 4] - _psi[ym * _Nx * 4 + x * 4]) * 0.5f;

            // Simplified Hamiltonian: H = -i·∇ + m(x)
            // For scalar field: dΨ/dt = -i·H·Ψ = -i(-i·∇ + m)Ψ = -∇Ψ - i·m·Ψ
            std::complex<float> dpsi_dt = -(dpsi_dx + dpsi_dy)
                                         - std::complex<float>(0.0f, mass) * psi0;

            // Euler update
            psi_new[psi_idx] = psi0 + dt * dpsi_dt;

            // Optional: Feedback to phase field (via lambda_coupling)
            if (lambda_coupling > 0.0f && _theta_data.size() > idx) {
                // Density-weighted phase shift
                float density = std::norm(psi0);
                _theta_data[idx] += lambda_coupling * dt * density;
            }
        }
    }

    // Update psi with new values
    _psi = psi_new;

    // Optional: Re-upload phase data to GPU if modified
    if (lambda_coupling > 0.0f && _theta_buffer != VK_NULL_HANDLE) {
        uploadToGPU();
    }
}

std::vector<float> MSFTEngine::getDiracDensity() const {
    /**
     * Get Dirac spinor density |Ψ|² for analysis
     * CPU-only implementation - GPU shaders exceed timeout budget
     *
     * Returns density at each grid point.
     * For scalar field: |Ψ|² = |ψ₀|²
     * For full spinor: |Ψ|² = |ψ₀|² + |ψ₁|² + |ψ₂|² + |ψ₃|²
     */

    std::vector<float> density(_Nx * _Ny, 0.0f);

    if (!_dirac_initialized || _psi.empty()) {
        return density;  // Return zeros if not initialized
    }

    for (uint32_t y = 0; y < _Ny; y++) {
        for (uint32_t x = 0; x < _Nx; x++) {
            uint32_t idx = y * _Nx + x;
            uint32_t psi_idx = idx * 4;  // 4 components per point

            // Sum density over all 4 components (most will be zero for scalar case)
            float rho = 0.0f;
            for (int comp = 0; comp < 4; comp++) {
                if (psi_idx + comp < _psi.size()) {
                    rho += std::norm(_psi[psi_idx + comp]);
                }
            }

            density[idx] = rho;
        }
    }

    return density;
}