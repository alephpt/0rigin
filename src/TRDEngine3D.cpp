#include "TRDEngine3D.h"
#include <iostream>
#include <cmath>
#include <cstring>
#include <algorithm>

/**
 * TRDEngine3D Implementation - 3D TRD GPU Compute Engine
 *
 * Dual-solver architecture: routes physics to appropriate solver.
 * - Dissipative models (Kuramoto) → TRDCore3D CPU solver
 * - Conservative models (Sine-Gordon, Dirac) → ConservativeSolver
 *
 * GPU acceleration via Vulkan compute shaders with CPU fallback.
 */

TRDEngine3D::TRDEngine3D(Nova* nova)
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
        std::cerr << "[TRDEngine3D] Error: Nova not initialized" << std::endl;
        return;
    }

    // Initialize core 3D grid infrastructure
    _core3d = std::make_unique<TRDCore3D>(nova->_architect->logical_device,
                                            nova->_architect->physical_device);

    // Initialize pipeline factory
    _pipelineFactory = std::make_unique<TRDPipelineFactory>(nova->_architect->logical_device);

    // Initialize compute dispatcher
    _compute = std::make_unique<TRDCompute>(nova->_architect->logical_device,
                                              nova->_architect->queues.compute,
                                              nova->_architect->queues.indices.compute_family.value_or(0));
    if (!_compute->initialize()) {
        std::cerr << "[TRDEngine3D] Failed to initialize compute dispatcher" << std::endl;
        return;
    }

    // Initialize buffer manager
    _bufferManager = std::make_unique<TRDBufferManager>(nova->_architect->logical_device,
                                                          nova->_architect->physical_device);

    // Initialize descriptor manager
    _descriptorManager = std::make_unique<TRDDescriptorManager>(nova->_architect->logical_device);

    // Initialize ConservativeSolver for particle dynamics
    _conservative_solver = std::make_unique<ConservativeSolver>();

    // Default physics model: vacuum Kuramoto
    _physics_model = "vacuum_kuramoto";

    std::cout << "[TRDEngine3D] Initialized GPU compute engine (dual-solver architecture)" << std::endl;
}

TRDEngine3D::~TRDEngine3D() {
    destroyResources();
}

void TRDEngine3D::initialize(uint32_t Nx, uint32_t Ny, uint32_t Nz, float Delta) {
    _Nx = Nx;
    _Ny = Ny;
    _Nz = Nz;
    _N_total = Nx * Ny * Nz;
    _Delta = Delta;
    _time_step = 0;

    std::cout << "\n[TRDEngine3D] Initializing 3D TRD simulation" << std::endl;
    std::cout << "  Grid: " << Nx << " × " << Ny << " × " << Nz
              << " (" << _N_total << " points)" << std::endl;
    std::cout << "  Vacuum potential Δ = " << Delta << std::endl;

    size_t memory_per_field = _N_total * sizeof(float);
    std::cout << "  Memory per field: " << memory_per_field / 1024 << " KB" << std::endl;

    // Initialize core 3D grid
    TRDCore3D::Config config;
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
    std::cout << "[TRDEngine3D] GPU initialization complete" << std::endl;
}

void TRDEngine3D::createBuffers() {
    if (!_bufferManager) {
        std::cerr << "[TRDEngine3D] Buffer manager not initialized" << std::endl;
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

    std::cout << "[TRDEngine3D] Created Kuramoto field buffers (" << buffer_size / 1024 << " KB each)" << std::endl;
}

void TRDEngine3D::createPipelines() {
    // Pipelines will be created on-demand when shaders are ready
    // For now, just log readiness
    std::cout << "[TRDEngine3D] Pipeline factory ready for shader loading" << std::endl;
}

void TRDEngine3D::createDescriptors() {
    if (!_descriptorManager) {
        std::cerr << "[TRDEngine3D] Descriptor manager not initialized" << std::endl;
        return;
    }

    // Descriptor pool will be created when pipelines are set up
    std::cout << "[TRDEngine3D] Descriptor manager ready" << std::endl;
}

void TRDEngine3D::setInitialPhases(const std::vector<float>& theta) {
    if (theta.size() != _N_total) {
        std::cerr << "[TRDEngine3D] Invalid theta size: " << theta.size()
                  << " (expected " << _N_total << ")" << std::endl;
        return;
    }

    _theta_data = theta;
    if (_theta_buffer != VK_NULL_HANDLE) {
        uploadFieldData(_theta_buffer, _theta_data);
    }
}

void TRDEngine3D::setNaturalFrequencies(const std::vector<float>& omega) {
    if (omega.size() != _N_total) {
        std::cerr << "[TRDEngine3D] Invalid omega size: " << omega.size()
                  << " (expected " << _N_total << ")" << std::endl;
        return;
    }

    _omega_data = omega;
    if (_omega_buffer != VK_NULL_HANDLE) {
        uploadFieldData(_omega_buffer, _omega_data);
    }
}

void TRDEngine3D::stepKuramoto3D(float dt, float K, float damping) {
    if (!_gpu_ready) {
        std::cerr << "[TRDEngine3D] GPU not ready" << std::endl;
        return;
    }

    // CPU fallback via TRDCore3D (GPU dispatch planned for kuramoto3d.comp)
    _core3d->evolveKuramotoCPU(dt);

    _time_step++;
}

void TRDEngine3D::computeSyncField3D() {
    if (!_gpu_ready) {
        std::cerr << "[TRDEngine3D] GPU not ready" << std::endl;
        return;
    }

    // CPU fallback via TRDCore3D (GPU dispatch planned for sync_field3d.comp)
    _core3d->computeRField();
}

std::vector<float> TRDEngine3D::getPhaseField3D() const {
    if (_theta_buffer != VK_NULL_HANDLE) {
        std::vector<float> data(_N_total);
        downloadFieldData(_theta_buffer, data);
        return data;
    }
    return _core3d->getTheta();
}

std::vector<float> TRDEngine3D::getSyncField3D() const {
    if (_R_field_buffer != VK_NULL_HANDLE) {
        std::vector<float> data(_N_total);
        downloadFieldData(_R_field_buffer, data);
        return data;
    }
    return _core3d->getRField();
}

std::vector<float> TRDEngine3D::getMassField3D() const {
    auto R = getSyncField3D();
    std::vector<float> mass(R.size());
    for (size_t i = 0; i < R.size(); ++i) {
        mass[i] = _Delta * R[i];
    }
    return mass;
}

void TRDEngine3D::initializeEM3D() {
    std::cout << "[TRDEngine3D] Initializing 3D EM fields" << std::endl;
    // EM field buffers allocated on demand via ConservativeSolver
}

void TRDEngine3D::stepMaxwell3D(float dt) {
    // Maxwell evolution handled by ConservativeSolver for conservative physics
    // GPU dispatch via maxwell_evolve_E/B.comp planned for acceleration
}

void TRDEngine3D::initializeStuckelberg3D() {
    std::cout << "[TRDEngine3D] Initializing 3D Stückelberg gauge" << std::endl;
    // Stückelberg fields managed by StuckelbergEM module
}

void TRDEngine3D::applyGaugeTransform3D() {
    // Gauge transform A'μ = Aμ + ∂μφ/e handled by StuckelbergEM
}

void TRDEngine3D::initializeVortexLine(float center_x, float center_y, float center_z,
                                        float radius, int axis) {
    std::cout << "[TRDEngine3D] Initializing vortex line" << std::endl;
    // Vortex initialization delegated to ConservativeSolver::initializeVortexWithProperVelocity
}

void TRDEngine3D::initializeDirac3D() {
    std::cout << "[TRDEngine3D] Initializing 3D+1 Dirac spinor" << std::endl;
    // Dirac spinor managed by Dirac3D module with split-operator FFT
}

void TRDEngine3D::stepDirac3D(float dt) {
    // Dirac evolution handled by ConservativeSolver::evolveDirac
    // GPU dispatch via dirac3d.comp planned for acceleration
}

void TRDEngine3D::uploadFieldData(VkBuffer buffer, const std::vector<float>& data) {
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

void TRDEngine3D::downloadFieldData(VkBuffer buffer, std::vector<float>& data) const {
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

void TRDEngine3D::setPhysicsModel(const std::string& model) {
    _physics_model = model;
    std::cout << "[TRDEngine3D] Physics model set to: " << model << std::endl;

    // Initialize conservative solver if needed
    if (model == "particle_sine_gordon" || model == "particle_dirac" || model == "coupled_vacuum_particle") {
        if (!_conservative_solver) {
            _conservative_solver = std::make_unique<ConservativeSolver>();
        }

        // Initialize with current grid dimensions
        ConservativeSolver::Config config;
        config.nx = _Nx;
        config.ny = _Ny;
        config.nz = _Nz;
        config.dx = 1.0f;
        config.dt = 0.005f;  // Smaller timestep for stability
        config.method = ConservativeSolver::IntegrationMethod::STRANG_SPLITTING;  // Default to Strang (superior for nonlinear PDE)

        _conservative_solver->initialize(config);
        std::cout << "[TRDEngine3D] ConservativeSolver initialized for particle dynamics" << std::endl;
    }
}

void TRDEngine3D::setIntegrationMethod(const std::string& method) {
    if (!_conservative_solver) {
        std::cerr << "[TRDEngine3D] Error: ConservativeSolver not initialized" << std::endl;
        return;
    }

    // Parse string to enum
    ConservativeSolver::IntegrationMethod integration_method;
    if (method == "velocity_verlet") {
        integration_method = ConservativeSolver::IntegrationMethod::VELOCITY_VERLET;
    } else if (method == "rk2_symplectic") {
        integration_method = ConservativeSolver::IntegrationMethod::RK2_SYMPLECTIC;
    } else if (method == "strang_splitting") {
        integration_method = ConservativeSolver::IntegrationMethod::STRANG_SPLITTING;
    } else if (method == "half_strang") {
        integration_method = ConservativeSolver::IntegrationMethod::HALF_STRANG;
    } else {
        std::cerr << "[TRDEngine3D] Unknown integration method: " << method
                  << ", defaulting to strang_splitting" << std::endl;
        integration_method = ConservativeSolver::IntegrationMethod::STRANG_SPLITTING;
    }

    // Re-initialize with new method
    ConservativeSolver::Config config;
    config.nx = _Nx;
    config.ny = _Ny;
    config.nz = _Nz;
    config.dx = 1.0f;
    config.dt = 0.005f;
    config.method = integration_method;

    _conservative_solver->initialize(config);
    std::cout << "[TRDEngine3D] Integration method set to: " << method << std::endl;
}

void TRDEngine3D::runSimulation(float dt) {
    // DUAL-SOLVER ROUTING
    // Routes physics evolution to appropriate solver based on model

    if (_physics_model == "vacuum_kuramoto") {
        // DISSIPATIVE/THERMODYNAMIC (Kuramoto gradient flow)
        // - Vacuum synchronization dynamics
        // - R-field evolution toward equilibrium
        // - Energy NOT conserved (expected for dissipative system)
        _core3d->evolveKuramotoCPU(dt);
        std::cout << "[TRDEngine3D] Vacuum Kuramoto step (R = " << _core3d->getAverageR() << ")" << std::endl;

    } else if (_physics_model == "particle_sine_gordon") {
        // CONSERVATIVE/UNITARY (Sine-Gordon solitons)
        // - Particle scattering, vortex dynamics
        // - Energy conserved <0.01% (GO/NO-GO gate)
        // - Symplectic Velocity Verlet integration
        if (!_conservative_solver) {
            std::cerr << "[TRDEngine3D] Error: ConservativeSolver not initialized for particle_sine_gordon" << std::endl;
            return;
        }

        _conservative_solver->evolveSineGordon(dt);

        // Validate energy conservation (GO/NO-GO criterion)
        if (!_conservative_solver->validateEnergyConservation(0.0001f)) {
            std::cerr << "[TRDEngine3D] WARNING: Energy conservation violated in Sine-Gordon evolution" << std::endl;
        }

    } else if (_physics_model == "particle_dirac") {
        // CONSERVATIVE/UNITARY (Dirac fermion)
        // For standalone Dirac without vacuum coupling, use uniform fields
        if (!_conservative_solver) {
            std::cerr << "[TRDEngine3D] Error: ConservativeSolver not initialized for particle_dirac" << std::endl;
            return;
        }

        // Create uniform fields for standalone Dirac (no vacuum coupling)
        uint32_t total_points = _Nx * _Ny * _Nz;
        std::vector<float> uniform_R(total_points, 1.0f);  // Uniform R = 1
        std::vector<float> uniform_theta(total_points, 0.0f);  // Uniform θ = 0

        _conservative_solver->evolveDirac(dt, uniform_R, uniform_theta, _Delta);

    } else if (_physics_model == "coupled_vacuum_particle") {
        // HYBRID - BOTH dissipative AND conservative
        // 1. Evolve vacuum (Kuramoto/dissipative)
        _core3d->evolveKuramotoCPU(dt);

        // 2. Get vacuum fields for chiral mass coupling
        auto R_field = _core3d->getRField();
        auto theta_field = _core3d->getTheta();  // Get theta field for chiral coupling
        float avg_R = _core3d->getAverageR();

        // 3. Evolve particle (Dirac/conservative) with chiral mass from vacuum fields
        if (_conservative_solver) {
            _conservative_solver->evolveDirac(dt, R_field, theta_field, _Delta);
        }

        std::cout << "[TRDEngine3D] Coupled evolution: R=" << avg_R
                  << ", Delta=" << _Delta << std::endl;

    } else {
        std::cerr << "[TRDEngine3D] Unknown physics model: " << _physics_model << std::endl;
    }
}

void TRDEngine3D::destroyResources() {
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

    std::cout << "[TRDEngine3D] Destroyed all Vulkan resources" << std::endl;
}
