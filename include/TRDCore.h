// include/TRDCore.h
#pragma once

#include "physics/StuckelbergEM.h"
#include <vector>
#include <memory>
#include <vulkan/vulkan.h>

/**
 * TRDCore - Minimal demonstration core for TRD+EM coupled evolution
 *
 * This is a lightweight wrapper demonstrating the integration pattern
 * for coupling Stückelberg EM fields to TRD dynamics. In the full TRDEngine,
 * this functionality would be integrated directly.
 *
 * Physics:
 *   - TRD: θ(x,y) phase field, R(x,y) synchronization order parameter
 *   - Stückelberg EM: Gauge-restored massive photon A'_μ = A_μ + ∂_μφ/e
 *   - Coupling: Direct φ = θ coupling (gauge invariant!)
 */
class TRDCore {
public:
    struct Config {
        uint32_t nx = 64;
        uint32_t ny = 64;
        float dx = 1.0f;
        float dt = 0.01f;
    };

    TRDCore(VkDevice device, VkPhysicalDevice physicalDevice);
    ~TRDCore() = default;

    /**
     * Initialize simulation with grid configuration
     */
    void initialize(const Config& config);

    /**
     * Enable electromagnetic field evolution with Stückelberg mechanism
     * @param photon_mass_coupling: photon mass m_γ
     */
    void enableEM(float photon_mass_coupling);

    /**
     * Evolve TRD fields one timestep (placeholder - would be full Kuramoto)
     */
    void evolveFields(float dt);

    /**
     * Evolve EM fields one timestep using Stückelberg dynamics
     */
    void evolveEM(float dt);

    /**
     * Get electromagnetic field tensor at grid point (i,j)
     */
    physics::FieldTensor getEMFieldAt(int i, int j) const;

    /**
     * Get current phase field (for coupling to EM)
     */
    const std::vector<float>& getPhaseField() const { return theta_field_; }

    /**
     * Get current synchronization field (for photon mass)
     */
    const std::vector<float>& getSyncField() const { return R_field_; }

private:
    // Vulkan resources
    VkDevice device_;
    VkPhysicalDevice physical_device_;

    // Configuration
    Config config_;

    // TRD field storage
    std::vector<float> theta_field_;  // Phase field θ(x,y)
    std::vector<float> R_field_;      // Synchronization field R(x,y)

    // EM support
    bool em_enabled_ = false;
    std::unique_ptr<physics::StuckelbergEM> stuckelberg_em_;
};
