// src/TRDCore.cpp
#include "TRDCore.h"
#include <iostream>
#include <cmath>

TRDCore::TRDCore(VkDevice device, VkPhysicalDevice physicalDevice)
    : device_(device), physical_device_(physicalDevice) {
}

void TRDCore::initialize(const Config& config) {
    config_ = config;

    // Allocate TRD field storage
    size_t grid_size = config_.nx * config_.ny;
    theta_field_.resize(grid_size, 0.0f);
    R_field_.resize(grid_size, 1.0f);  // Start synchronized

    // Initialize with simple vortex configuration
    int cx = config_.nx / 2;
    int cy = config_.ny / 2;

    for (uint32_t j = 0; j < config_.ny; ++j) {
        for (uint32_t i = 0; i < config_.nx; ++i) {
            int idx = j * config_.nx + i;

            // Vortex phase: θ = atan2(y - cy, x - cx)
            float dx = static_cast<float>(i) - cx;
            float dy = static_cast<float>(j) - cy;
            theta_field_[idx] = std::atan2(dy, dx);

            // Gaussian synchronization profile
            float r_sq = dx*dx + dy*dy;
            float sigma = 10.0f;
            R_field_[idx] = std::exp(-r_sq / (2.0f * sigma * sigma));
        }
    }

    std::cout << "[TRDCore] Initialized " << config_.nx << "x" << config_.ny
              << " grid with vortex configuration" << std::endl;
}

void TRDCore::enableEM(float photon_mass_coupling) {
    em_enabled_ = true;

    // Create Stückelberg field evolution
    stuckelberg_em_ = std::make_unique<physics::StuckelbergEM>(
        config_.nx, config_.ny, config_.dx, photon_mass_coupling
    );

    std::cout << "[TRDCore] Stückelberg EM enabled (m_γ = "
              << photon_mass_coupling << ")" << std::endl;
    std::cout << "[TRDCore] Gauge-restored mechanism: A'_μ = A_μ + ∂_μφ/e (φ = θ)"
              << std::endl;
}

void TRDCore::evolveFields(float dt) {
    // Placeholder: Simple phase diffusion (would be full Kuramoto evolution)
    (void)dt;

    // For now, just maintain the vortex structure
    // In full implementation, this would be:
    //   1. Compute Kuramoto coupling
    //   2. Update phases: dθ/dt = ω + K·Σsin(θⱼ - θᵢ)
    //   3. Compute synchronization: R = |⟨e^(iθ)⟩|
}

void TRDCore::evolveEM(float dt) {
    if (!em_enabled_) {
        return;
    }

    // Stückelberg evolution
    stuckelberg_em_->computePotentials(
        theta_field_.data(),
        R_field_.data(),
        config_.nx,
        config_.ny,
        config_.dx,
        dt
    );

    stuckelberg_em_->computeFieldStrengths();
}

physics::FieldTensor TRDCore::getEMFieldAt(int i, int j) const {
    if (!em_enabled_) {
        return physics::FieldTensor{};
    }
    return stuckelberg_em_->getFieldAt(i, j);
}
