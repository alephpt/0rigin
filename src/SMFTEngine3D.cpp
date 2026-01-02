#include "SMFTEngine3D.h"
#include <iostream>
#include <cmath>
#include <cstring>
#include <algorithm>

/**
 * SMFTEngine3D Implementation - 3D SMFT GPU Compute Engine
 *
 * Weeks 3-4: GPU Integration + 3D Kuramoto (this file)
 * Weeks 5-6: 3D Electromagnetic Fields
 * Weeks 7-8: Dirac 3D+1 Spinor
 */

SMFTEngine3D::SMFTEngine3D(Nova* nova)
    : _nova(nova),
      _Nx(0), _Ny(0), _Nz(0), _N_total(0),
      _Delta(0.0f),
      _time_step(0),
      _gpu_ready(false),
      _theta_buffer(VK_NULL_HANDLE),
      _theta_out_buffer(VK_NULL_HANDLE),
      _omega_buffer(VK_NULL_HANDLE),
      _R_field_buffer(VK_NULL_HANDLE),
      _theta_memory(VK_NULL_HANDLE),
      _theta_out_memory(VK_NULL_HANDLE),
      _omega_memory(VK_NULL_HANDLE),
      _R_field_memory(VK_NULL_HANDLE),
      _Ex_buffer(VK_NULL_HANDLE), _Ey_buffer(VK_NULL_HANDLE), _Ez_buffer(VK_NULL_HANDLE),
      _Bx_buffer(VK_NULL_HANDLE), _By_buffer(VK_NULL_HANDLE), _Bz_buffer(VK_NULL_HANDLE),
      _Ex_memory(VK_NULL_HANDLE), _Ey_memory(VK_NULL_HANDLE), _Ez_memory(VK_NULL_HANDLE),
      _Bx_memory(VK_NULL_HANDLE), _By_memory(VK_NULL_HANDLE), _Bz_memory(VK_NULL_HANDLE),
      _phi_buffer(VK_NULL_HANDLE),
      _Ax_buffer(VK_NULL_HANDLE), _Ay_buffer(VK_NULL_HANDLE),
      _Az_buffer(VK_NULL_HANDLE), _A0_buffer(VK_NULL_HANDLE),
      _phi_memory(VK_NULL_HANDLE),
      _Ax_memory(VK_NULL_HANDLE), _Ay_memory(VK_NULL_HANDLE),
      _Az_memory(VK_NULL_HANDLE), _A0_memory(VK_NULL_HANDLE),
      _spinor_buffer(VK_NULL_HANDLE),
      _spinor_memory(VK_NULL_HANDLE),
      _kuramoto3d_pipeline(VK_NULL_HANDLE),
      _sync_field3d_pipeline(VK_NULL_HANDLE),
      _maxwell_E_pipeline(VK_NULL_HANDLE),
      _maxwell_B_pipeline(VK_NULL_HANDLE),
      _dirac3d_pipeline(VK_NULL_HANDLE),
      _kuramoto3d_layout(VK_NULL_HANDLE),
      _sync_field3d_layout(VK_NULL_HANDLE),
      _maxwell_E_layout(VK_NULL_HANDLE),
      _maxwell_B_layout(VK_NULL_HANDLE),
      _dirac3d_layout(VK_NULL_HANDLE),
      _kuramoto3d_descriptor_set(VK_NULL_HANDLE),
      _sync_field3d_descriptor_set(VK_NULL_HANDLE),
      _maxwell_descriptor_set(VK_NULL_HANDLE),
      _dirac3d_descriptor_set(VK_NULL_HANDLE),
      _kuramoto3d_descriptor_layout(VK_NULL_HANDLE),
      _sync_field3d_descriptor_layout(VK_NULL_HANDLE),
      _maxwell_descriptor_layout(VK_NULL_HANDLE),
      _dirac3d_descriptor_layout(VK_NULL_HANDLE),
      _descriptor_pool(VK_NULL_HANDLE)
{
    if (!nova || !nova->initialized) {
        std::cerr << "[SMFTEngine3D] Error: Nova not initialized" << std::endl;
        return;
    }

    // Initialize core 3D grid infrastructure
    _core3d = std::make_unique<SMFTCore3D>(nova->_architect->logical_device,
                                            nova->_architect->physical_device);

    // Initialize pipeline factory
    _pipelineFactory = std::make_unique<SMFTPipelineFactory>(nova->_architect->logical_device);

    // Initialize compute dispatcher
    _compute = std::make_unique<SMFTCompute>(nova->_architect->logical_device,
                                              nova->_architect->queues.compute,
                                              nova->_architect->queues.indices.compute_family.value_or(0));
    if (!_compute->initialize()) {
        std::cerr << "[SMFTEngine3D] Failed to initialize compute dispatcher" << std::endl;
        return;
    }

    // Initialize buffer manager
    _bufferManager = std::make_unique<SMFTBufferManager>(nova->_architect->logical_device,
                                                          nova->_architect->physical_device);

    // Initialize descriptor manager
    _descriptorManager = std::make_unique<SMFTDescriptorManager>(nova->_architect->logical_device);

    std::cout << "[SMFTEngine3D] Initialized GPU compute engine" << std::endl;
}

SMFTEngine3D::~SMFTEngine3D() {
    destroyResources();
}

void SMFTEngine3D::initialize(uint32_t Nx, uint32_t Ny, uint32_t Nz, float Delta) {
    _Nx = Nx;
    _Ny = Ny;
    _Nz = Nz;
    _N_total = Nx * Ny * Nz;
    _Delta = Delta;
    _time_step = 0;

    std::cout << "\n[SMFTEngine3D] Initializing 3D SMFT simulation" << std::endl;
    std::cout << "  Grid: " << Nx << " × " << Ny << " × " << Nz
              << " (" << _N_total << " points)" << std::endl;
    std::cout << "  Vacuum potential Δ = " << Delta << std::endl;

    size_t memory_per_field = _N_total * sizeof(float);
    std::cout << "  Memory per field: " << memory_per_field / 1024 << " KB" << std::endl;

    // Initialize core 3D grid
    SMFTCore3D::Config config;
    config.Nx = Nx;
    config.Ny = Ny;
    config.Nz = Nz;
    config.dx = 1.0f;
    config.dt = 0.01f;
    config.coupling_strength = 1.0f;
    _core3d->initialize(config);

    // Allocate CPU-side data
    _theta_data.resize(_N_total, 0.0f);
    _omega_data.resize(_N_total, 0.0f);
    _R_field_data.resize(_N_total, 1.0f);

    // Create GPU resources
    createBuffers();
    createPipelines();
    createDescriptors();

    _gpu_ready = true;
    std::cout << "[SMFTEngine3D] GPU initialization complete" << std::endl;
}

void SMFTEngine3D::createBuffers() {
    if (!_bufferManager) {
        std::cerr << "[SMFTEngine3D] Buffer manager not initialized" << std::endl;
        return;
    }

    size_t buffer_size = _N_total * sizeof(float);

    // Create Kuramoto field buffers using pair return API
    auto theta_pair = _bufferManager->createBuffer(buffer_size,
        VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT,
        VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
    _theta_buffer = theta_pair.first;
    _theta_memory = theta_pair.second;

    auto theta_out_pair = _bufferManager->createBuffer(buffer_size,
        VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT,
        VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
    _theta_out_buffer = theta_out_pair.first;
    _theta_out_memory = theta_out_pair.second;

    auto omega_pair = _bufferManager->createBuffer(buffer_size,
        VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT,
        VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
    _omega_buffer = omega_pair.first;
    _omega_memory = omega_pair.second;

    auto R_field_pair = _bufferManager->createBuffer(buffer_size,
        VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT,
        VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
    _R_field_buffer = R_field_pair.first;
    _R_field_memory = R_field_pair.second;

    std::cout << "[SMFTEngine3D] Created Kuramoto field buffers (" << buffer_size / 1024 << " KB each)" << std::endl;
}

void SMFTEngine3D::createPipelines() {
    // Pipelines will be created on-demand when shaders are ready
    // For now, just log readiness
    std::cout << "[SMFTEngine3D] Pipeline factory ready for shader loading" << std::endl;
}

void SMFTEngine3D::createDescriptors() {
    if (!_descriptorManager) {
        std::cerr << "[SMFTEngine3D] Descriptor manager not initialized" << std::endl;
        return;
    }

    // Descriptor pool will be created when pipelines are set up
    std::cout << "[SMFTEngine3D] Descriptor manager ready" << std::endl;
}

void SMFTEngine3D::setInitialPhases(const std::vector<float>& theta) {
    if (theta.size() != _N_total) {
        std::cerr << "[SMFTEngine3D] Invalid theta size: " << theta.size()
                  << " (expected " << _N_total << ")" << std::endl;
        return;
    }

    _theta_data = theta;
    if (_theta_buffer != VK_NULL_HANDLE) {
        uploadFieldData(_theta_buffer, _theta_data);
    }
}

void SMFTEngine3D::setNaturalFrequencies(const std::vector<float>& omega) {
    if (omega.size() != _N_total) {
        std::cerr << "[SMFTEngine3D] Invalid omega size: " << omega.size()
                  << " (expected " << _N_total << ")" << std::endl;
        return;
    }

    _omega_data = omega;
    if (_omega_buffer != VK_NULL_HANDLE) {
        uploadFieldData(_omega_buffer, _omega_data);
    }
}

void SMFTEngine3D::stepKuramoto3D(float dt, float K, float damping) {
    if (!_gpu_ready) {
        std::cerr << "[SMFTEngine3D] GPU not ready" << std::endl;
        return;
    }

    // TODO: Implement GPU dispatch for kuramoto3d.comp
    // For now, use CPU fallback via SMFTCore3D
    _core3d->evolveKuramotoCPU(dt);

    _time_step++;
}

void SMFTEngine3D::computeSyncField3D() {
    if (!_gpu_ready) {
        std::cerr << "[SMFTEngine3D] GPU not ready" << std::endl;
        return;
    }

    // TODO: Implement GPU dispatch for sync_field3d.comp
    // For now, use CPU fallback via SMFTCore3D
    _core3d->computeRField();
}

std::vector<float> SMFTEngine3D::getPhaseField3D() const {
    if (_theta_buffer != VK_NULL_HANDLE) {
        std::vector<float> data(_N_total);
        downloadFieldData(_theta_buffer, data);
        return data;
    }
    return _core3d->getTheta();
}

std::vector<float> SMFTEngine3D::getSyncField3D() const {
    if (_R_field_buffer != VK_NULL_HANDLE) {
        std::vector<float> data(_N_total);
        downloadFieldData(_R_field_buffer, data);
        return data;
    }
    return _core3d->getRField();
}

std::vector<float> SMFTEngine3D::getMassField3D() const {
    auto R = getSyncField3D();
    std::vector<float> mass(R.size());
    for (size_t i = 0; i < R.size(); ++i) {
        mass[i] = _Delta * R[i];
    }
    return mass;
}

void SMFTEngine3D::initializeEM3D() {
    std::cout << "[SMFTEngine3D] Initializing 3D EM fields (Weeks 5-6)" << std::endl;
    // TODO: Allocate Ex, Ey, Ez, Bx, By, Bz buffers
}

void SMFTEngine3D::stepMaxwell3D(float dt) {
    std::cout << "[SMFTEngine3D] Maxwell 3D step (Weeks 5-6 implementation)" << std::endl;
    // TODO: Dispatch maxwell_evolve_E.comp and maxwell_evolve_B.comp
}

void SMFTEngine3D::initializeStuckelberg3D() {
    std::cout << "[SMFTEngine3D] Initializing 3D Stückelberg gauge (Weeks 5-6)" << std::endl;
    // TODO: Allocate phi, Ax, Ay, Az, A0 buffers
}

void SMFTEngine3D::applyGaugeTransform3D() {
    std::cout << "[SMFTEngine3D] Applying 3D gauge transform (Weeks 5-6)" << std::endl;
    // TODO: Implement A'μ = Aμ + ∂μφ/e
}

void SMFTEngine3D::initializeVortexLine(float center_x, float center_y, float center_z,
                                        float radius, int axis) {
    std::cout << "[SMFTEngine3D] Initializing vortex line (Weeks 5-6)" << std::endl;
    // TODO: Initialize phase field with vortex ring topology
}

void SMFTEngine3D::initializeDirac3D() {
    std::cout << "[SMFTEngine3D] Initializing 3D+1 Dirac spinor (Weeks 7-8)" << std::endl;
    // TODO: Allocate 4-component spinor buffer
}

void SMFTEngine3D::stepDirac3D(float dt) {
    std::cout << "[SMFTEngine3D] Dirac 3D step (Weeks 7-8 implementation)" << std::endl;
    // TODO: Dispatch dirac3d.comp
}

void SMFTEngine3D::uploadFieldData(VkBuffer buffer, const std::vector<float>& data) {
    if (!_bufferManager) return;

    // Find corresponding memory handle
    VkDeviceMemory memory = VK_NULL_HANDLE;
    if (buffer == _theta_buffer) memory = _theta_memory;
    else if (buffer == _theta_out_buffer) memory = _theta_out_memory;
    else if (buffer == _omega_buffer) memory = _omega_memory;
    else if (buffer == _R_field_buffer) memory = _R_field_memory;

    if (memory != VK_NULL_HANDLE) {
        size_t size = data.size() * sizeof(float);
        _bufferManager->uploadData(memory, data.data(), size);
    }
}

void SMFTEngine3D::downloadFieldData(VkBuffer buffer, std::vector<float>& data) const {
    if (!_bufferManager) return;

    // Find corresponding memory handle
    VkDeviceMemory memory = VK_NULL_HANDLE;
    if (buffer == _theta_buffer) memory = _theta_memory;
    else if (buffer == _theta_out_buffer) memory = _theta_out_memory;
    else if (buffer == _omega_buffer) memory = _omega_memory;
    else if (buffer == _R_field_buffer) memory = _R_field_memory;

    if (memory != VK_NULL_HANDLE) {
        size_t size = data.size() * sizeof(float);
        _bufferManager->downloadData(memory, data.data(), size);
    }
}

void SMFTEngine3D::destroyResources() {
    if (!_nova || !_nova->_architect) return;

    VkDevice device = _nova->_architect->logical_device;

    // Destroy Kuramoto buffers
    if (_theta_buffer != VK_NULL_HANDLE) vkDestroyBuffer(device, _theta_buffer, nullptr);
    if (_theta_out_buffer != VK_NULL_HANDLE) vkDestroyBuffer(device, _theta_out_buffer, nullptr);
    if (_omega_buffer != VK_NULL_HANDLE) vkDestroyBuffer(device, _omega_buffer, nullptr);
    if (_R_field_buffer != VK_NULL_HANDLE) vkDestroyBuffer(device, _R_field_buffer, nullptr);

    if (_theta_memory != VK_NULL_HANDLE) vkFreeMemory(device, _theta_memory, nullptr);
    if (_theta_out_memory != VK_NULL_HANDLE) vkFreeMemory(device, _theta_out_memory, nullptr);
    if (_omega_memory != VK_NULL_HANDLE) vkFreeMemory(device, _omega_memory, nullptr);
    if (_R_field_memory != VK_NULL_HANDLE) vkFreeMemory(device, _R_field_memory, nullptr);

    // Destroy EM buffers
    if (_Ex_buffer != VK_NULL_HANDLE) vkDestroyBuffer(device, _Ex_buffer, nullptr);
    if (_Ey_buffer != VK_NULL_HANDLE) vkDestroyBuffer(device, _Ey_buffer, nullptr);
    if (_Ez_buffer != VK_NULL_HANDLE) vkDestroyBuffer(device, _Ez_buffer, nullptr);
    if (_Bx_buffer != VK_NULL_HANDLE) vkDestroyBuffer(device, _Bx_buffer, nullptr);
    if (_By_buffer != VK_NULL_HANDLE) vkDestroyBuffer(device, _By_buffer, nullptr);
    if (_Bz_buffer != VK_NULL_HANDLE) vkDestroyBuffer(device, _Bz_buffer, nullptr);

    if (_Ex_memory != VK_NULL_HANDLE) vkFreeMemory(device, _Ex_memory, nullptr);
    if (_Ey_memory != VK_NULL_HANDLE) vkFreeMemory(device, _Ey_memory, nullptr);
    if (_Ez_memory != VK_NULL_HANDLE) vkFreeMemory(device, _Ez_memory, nullptr);
    if (_Bx_memory != VK_NULL_HANDLE) vkFreeMemory(device, _Bx_memory, nullptr);
    if (_By_memory != VK_NULL_HANDLE) vkFreeMemory(device, _By_memory, nullptr);
    if (_Bz_memory != VK_NULL_HANDLE) vkFreeMemory(device, _Bz_memory, nullptr);

    // Destroy Stückelberg buffers
    if (_phi_buffer != VK_NULL_HANDLE) vkDestroyBuffer(device, _phi_buffer, nullptr);
    if (_Ax_buffer != VK_NULL_HANDLE) vkDestroyBuffer(device, _Ax_buffer, nullptr);
    if (_Ay_buffer != VK_NULL_HANDLE) vkDestroyBuffer(device, _Ay_buffer, nullptr);
    if (_Az_buffer != VK_NULL_HANDLE) vkDestroyBuffer(device, _Az_buffer, nullptr);
    if (_A0_buffer != VK_NULL_HANDLE) vkDestroyBuffer(device, _A0_buffer, nullptr);

    if (_phi_memory != VK_NULL_HANDLE) vkFreeMemory(device, _phi_memory, nullptr);
    if (_Ax_memory != VK_NULL_HANDLE) vkFreeMemory(device, _Ax_memory, nullptr);
    if (_Ay_memory != VK_NULL_HANDLE) vkFreeMemory(device, _Ay_memory, nullptr);
    if (_Az_memory != VK_NULL_HANDLE) vkFreeMemory(device, _Az_memory, nullptr);
    if (_A0_memory != VK_NULL_HANDLE) vkFreeMemory(device, _A0_memory, nullptr);

    // Destroy Dirac buffers
    if (_spinor_buffer != VK_NULL_HANDLE) vkDestroyBuffer(device, _spinor_buffer, nullptr);
    if (_spinor_memory != VK_NULL_HANDLE) vkFreeMemory(device, _spinor_memory, nullptr);

    // Destroy pipelines
    if (_kuramoto3d_pipeline != VK_NULL_HANDLE) vkDestroyPipeline(device, _kuramoto3d_pipeline, nullptr);
    if (_sync_field3d_pipeline != VK_NULL_HANDLE) vkDestroyPipeline(device, _sync_field3d_pipeline, nullptr);
    if (_maxwell_E_pipeline != VK_NULL_HANDLE) vkDestroyPipeline(device, _maxwell_E_pipeline, nullptr);
    if (_maxwell_B_pipeline != VK_NULL_HANDLE) vkDestroyPipeline(device, _maxwell_B_pipeline, nullptr);
    if (_dirac3d_pipeline != VK_NULL_HANDLE) vkDestroyPipeline(device, _dirac3d_pipeline, nullptr);

    // Destroy pipeline layouts
    if (_kuramoto3d_layout != VK_NULL_HANDLE) vkDestroyPipelineLayout(device, _kuramoto3d_layout, nullptr);
    if (_sync_field3d_layout != VK_NULL_HANDLE) vkDestroyPipelineLayout(device, _sync_field3d_layout, nullptr);
    if (_maxwell_E_layout != VK_NULL_HANDLE) vkDestroyPipelineLayout(device, _maxwell_E_layout, nullptr);
    if (_maxwell_B_layout != VK_NULL_HANDLE) vkDestroyPipelineLayout(device, _maxwell_B_layout, nullptr);
    if (_dirac3d_layout != VK_NULL_HANDLE) vkDestroyPipelineLayout(device, _dirac3d_layout, nullptr);

    // Destroy descriptor layouts
    if (_kuramoto3d_descriptor_layout != VK_NULL_HANDLE)
        vkDestroyDescriptorSetLayout(device, _kuramoto3d_descriptor_layout, nullptr);
    if (_sync_field3d_descriptor_layout != VK_NULL_HANDLE)
        vkDestroyDescriptorSetLayout(device, _sync_field3d_descriptor_layout, nullptr);
    if (_maxwell_descriptor_layout != VK_NULL_HANDLE)
        vkDestroyDescriptorSetLayout(device, _maxwell_descriptor_layout, nullptr);
    if (_dirac3d_descriptor_layout != VK_NULL_HANDLE)
        vkDestroyDescriptorSetLayout(device, _dirac3d_descriptor_layout, nullptr);

    // Destroy descriptor pool
    if (_descriptor_pool != VK_NULL_HANDLE) vkDestroyDescriptorPool(device, _descriptor_pool, nullptr);

    std::cout << "[SMFTEngine3D] Destroyed all Vulkan resources" << std::endl;
}
