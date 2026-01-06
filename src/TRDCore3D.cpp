// src/TRDCore3D.cpp
#include "TRDCore3D.h"
#include <iostream>
#include <cmath>
#include <random>
#include <numeric>
#include <complex>

TRDCore3D::TRDCore3D(VkDevice device, VkPhysicalDevice physicalDevice)
    : _device(device), _physical_device(physicalDevice), _gpu_enabled(false) {
}

TRDCore3D::TRDCore3D()
    : _device(VK_NULL_HANDLE), _physical_device(VK_NULL_HANDLE), _gpu_enabled(false) {
}

void TRDCore3D::initialize(const Config& config) {
    _config = config;
    _Nx = config.Nx;
    _Ny = config.Ny;
    _Nz = config.Nz;
    _N_total = _Nx * _Ny * _Nz;

    // Allocate 3D field storage
    _theta_data.resize(_N_total, 0.0f);
    _omega_data.resize(_N_total, 0.0f);
    _R_field_data.resize(_N_total, 1.0f);  // Start fully synchronized

    // Log initialization
    std::cout << "[TRDCore3D] Initialized " << _Nx << "x" << _Ny << "x" << _Nz
              << " grid (" << _N_total << " points)" << std::endl;

    size_t memory_per_field = _N_total * sizeof(float);
    size_t total_memory = 3 * memory_per_field;
    std::cout << "[TRDCore3D] Memory usage: " << memory_per_field / 1024
              << " KB per field, " << total_memory / 1024 << " KB total" << std::endl;
}

TRDCore3D::Neighbors3D TRDCore3D::getNeighbors(uint32_t i, uint32_t j, uint32_t k) const {
    Neighbors3D neighbors;

    // X neighbors with periodic boundary
    neighbors.x_plus = index3D(wrapX(i + 1), j, k);
    neighbors.x_minus = index3D(wrapX(i - 1), j, k);

    // Y neighbors with periodic boundary
    neighbors.y_plus = index3D(i, wrapY(j + 1), k);
    neighbors.y_minus = index3D(i, wrapY(j - 1), k);

    // Z neighbors with periodic boundary
    neighbors.z_plus = index3D(i, j, wrapZ(k + 1));
    neighbors.z_minus = index3D(i, j, wrapZ(k - 1));

    return neighbors;
}

void TRDCore3D::initializeUniform(float phase) {
    std::fill(_theta_data.begin(), _theta_data.end(), phase);
    std::fill(_omega_data.begin(), _omega_data.end(), 0.0f);  // No intrinsic frequency
    std::fill(_R_field_data.begin(), _R_field_data.end(), 1.0f);  // Perfect sync for uniform

    std::cout << "[TRDCore3D] Initialized uniform phase field (θ = "
              << phase << " rad)" << std::endl;
}

void TRDCore3D::initializeRandom(uint32_t seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> phase_dist(-M_PI, M_PI);
    std::normal_distribution<float> omega_dist(0.0f, 0.1f);  // Small frequency spread

    for (uint32_t i = 0; i < _N_total; ++i) {
        _theta_data[i] = phase_dist(rng);
        _omega_data[i] = omega_dist(rng);
    }

    // Compute initial R field
    computeRField();

    std::cout << "[TRDCore3D] Initialized random phase field (seed = "
              << seed << ")" << std::endl;
}

float TRDCore3D::computeKuramotoCoupling(uint32_t idx) const {
    // Get 3D coordinates
    uint32_t i, j, k;
    coords3D(idx, i, j, k);

    // Get neighbors
    Neighbors3D neighbors = getNeighbors(i, j, k);

    // Current phase
    float theta_self = _theta_data[idx];

    // Kuramoto coupling: sum of sin(theta_neighbor - theta_self) over 6 neighbors
    float coupling = 0.0f;
    coupling += std::sin(_theta_data[neighbors.x_plus] - theta_self);
    coupling += std::sin(_theta_data[neighbors.x_minus] - theta_self);
    coupling += std::sin(_theta_data[neighbors.y_plus] - theta_self);
    coupling += std::sin(_theta_data[neighbors.y_minus] - theta_self);
    coupling += std::sin(_theta_data[neighbors.z_plus] - theta_self);
    coupling += std::sin(_theta_data[neighbors.z_minus] - theta_self);

    // Normalize by number of neighbors (6 in 3D)
    return coupling / 6.0f;
}

void TRDCore3D::evolveKuramotoCPU(float dt) {
    // Dispatch to appropriate integrator based on mode
    switch (_config.mode) {
        case IntegrationMode::EULER:
            evolveEulerCPU(dt);
            break;
        case IntegrationMode::SYMPLECTIC:
            evolveSymplecticCPU(dt);
            break;
        default:
            std::cerr << "[TRDCore3D] Unknown integration mode, defaulting to symplectic"
                      << std::endl;
            evolveSymplecticCPU(dt);
            break;
    }
}

void TRDCore3D::evolveSymplecticCPU(float dt) {
    static bool first_call = true;
    if (first_call) {
        std::cout << "[TRDCore3D] Using SYMPLECTIC integration (RK2 Midpoint Method)" << std::endl;
        #ifdef _OPENMP
        std::cout << "[TRDCore3D] OpenMP parallelization enabled" << std::endl;
        #endif
        first_call = false;
    }

    // RK2 Midpoint Method (symplectic for first-order systems)
    // For dθ/dt = f(θ), this method is:
    //   k1 = f(θ(t))
    //   k2 = f(θ(t) + k1·dt/2)
    //   θ(t+dt) = θ(t) + k2·dt

    // Step 1: Compute k1 = f(θ) at current state
    std::vector<float> k1(_N_total);
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (uint32_t idx = 0; idx < _N_total; ++idx) {
        float coupling = computeKuramotoCoupling(idx);
        k1[idx] = _omega_data[idx] + _config.coupling_strength * coupling;
    }

    // Step 2: Compute midpoint θ_mid = θ + k1·dt/2
    std::vector<float> theta_mid(_N_total);
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (uint32_t idx = 0; idx < _N_total; ++idx) {
        theta_mid[idx] = _theta_data[idx] + 0.5f * dt * k1[idx];

        // Wrap to [-π, π]
        while (theta_mid[idx] > M_PI) theta_mid[idx] -= 2 * M_PI;
        while (theta_mid[idx] < -M_PI) theta_mid[idx] += 2 * M_PI;
    }

    // Step 3: Temporarily update field to compute k2
    std::vector<float> theta_old = std::move(_theta_data);
    _theta_data = std::move(theta_mid);

    // Step 4: Compute k2 = f(θ_mid)
    std::vector<float> k2(_N_total);
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (uint32_t idx = 0; idx < _N_total; ++idx) {
        float coupling = computeKuramotoCoupling(idx);
        k2[idx] = _omega_data[idx] + _config.coupling_strength * coupling;
    }

    // Step 5: Final update θ(t+dt) = θ(t) + k2·dt
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (uint32_t idx = 0; idx < _N_total; ++idx) {
        _theta_data[idx] = theta_old[idx] + dt * k2[idx];

        // Wrap to [-π, π]
        while (_theta_data[idx] > M_PI) _theta_data[idx] -= 2 * M_PI;
        while (_theta_data[idx] < -M_PI) _theta_data[idx] += 2 * M_PI;
    }

    // Update synchronization field
    computeRField();
}

void TRDCore3D::evolveEulerCPU(float dt) {
    static bool first_call = true;
    if (first_call) {
        std::cout << "[TRDCore3D] Using EULER integration (dissipative)" << std::endl;
        first_call = false;
    }

    // Legacy Euler integration (dissipative)
    // Temporary storage for new phases
    std::vector<float> new_theta(_N_total);

    // Update each oscillator
    for (uint32_t idx = 0; idx < _N_total; ++idx) {
        // Kuramoto dynamics: dθ/dt = ω + K * sum(sin(θ_j - θ_i))
        float coupling = computeKuramotoCoupling(idx);
        float dtheta_dt = _omega_data[idx] + _config.coupling_strength * coupling;

        // Euler integration
        new_theta[idx] = _theta_data[idx] + dt * dtheta_dt;

        // Keep phase in [-π, π] range
        while (new_theta[idx] > M_PI) new_theta[idx] -= 2 * M_PI;
        while (new_theta[idx] < -M_PI) new_theta[idx] += 2 * M_PI;
    }

    // Update phase field
    _theta_data = std::move(new_theta);

    // Update synchronization field
    computeRField();
}

float TRDCore3D::computeEnergy() const {
    // Kuramoto model energy: E = -K * sum_{<i,j>} cos(θ_j - θ_i)
    // For each pair of neighbors, we have -K*cos(θ_j - θ_i)
    float E = 0.0f;

    for (uint32_t idx = 0; idx < _N_total; ++idx) {
        // Get 3D coordinates
        uint32_t i, j, k;
        coords3D(idx, i, j, k);

        // Get neighbors
        Neighbors3D neighbors = getNeighbors(i, j, k);

        float theta_self = _theta_data[idx];

        // Add coupling energy for each neighbor pair
        // We count each pair once by only summing forward neighbors (x+, y+, z+)
        float theta_xp = _theta_data[neighbors.x_plus];
        float theta_yp = _theta_data[neighbors.y_plus];
        float theta_zp = _theta_data[neighbors.z_plus];

        E -= _config.coupling_strength * std::cos(theta_xp - theta_self);
        E -= _config.coupling_strength * std::cos(theta_yp - theta_self);
        E -= _config.coupling_strength * std::cos(theta_zp - theta_self);
    }

    return E;
}

void TRDCore3D::computeRField() {
    // For now, compute global R and store it uniformly
    // In full implementation, compute local R in neighborhoods

    std::complex<float> order_parameter(0.0f, 0.0f);

    for (uint32_t idx = 0; idx < _N_total; ++idx) {
        float theta = _theta_data[idx];
        order_parameter += std::complex<float>(std::cos(theta), std::sin(theta));
    }

    order_parameter /= static_cast<float>(_N_total);
    float global_R = std::abs(order_parameter);

    // Store global R at all points (simplified for now)
    std::fill(_R_field_data.begin(), _R_field_data.end(), global_R);
}

float TRDCore3D::getAverageR() const {
    // Since we're storing global R uniformly, just return first value
    // In full implementation, would average local R values
    return _R_field_data.empty() ? 0.0f : _R_field_data[0];
}

void TRDCore3D::enableGPU() {
    if (_device == VK_NULL_HANDLE) {
        std::cerr << "[TRDCore3D] Error: Cannot enable GPU without valid Vulkan device"
                  << std::endl;
        return;
    }

    _gpu_enabled = true;
    std::cout << "[TRDCore3D] GPU compute enabled" << std::endl;

    // GPU buffer allocation and shader compilation will be added in Week 2
}