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