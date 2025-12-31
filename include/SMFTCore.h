// include/SMFTCore.h
#pragma once

#include "physics/ProcaEM.h"
#include "physics/ProcaVulkanPipeline.h"
#include <vector>
#include <memory>
#include <vulkan/vulkan.h>

/**
 * SMFTCore - Minimal demonstration core for SMFT+EM coupled evolution
 *
 * This is a lightweight wrapper demonstrating the integration pattern
 * for coupling Proca EM fields to SMFT dynamics. In the full MSFTEngine,
 * this functionality would be integrated directly.
 *
 * Physics:
 *   - SMFT: θ(x,y) phase field, R(x,y) synchronization order parameter
 *   - Proca EM: Massive photon with m_γ(x,y) = g·(1-R(x,y))
 *   - Coupling: SMFT Noether current j_μ sources EM field
 */
class SMFTCore {
public:
    struct Config {
        uint32_t nx = 64;
        uint32_t ny = 64;
        float dx = 1.0f;
        float dt = 0.01f;
    };

    SMFTCore(VkDevice device, VkPhysicalDevice physicalDevice);
    ~SMFTCore() = default;

    /**
     * Initialize simulation with grid configuration
     */
    void initialize(const Config& config);

    /**
     * Enable electromagnetic field evolution with Proca mechanism
     * @param photon_mass_coupling: g in m_γ = g(1-R)
     */
    void enableEM(float photon_mass_coupling);

    /**
     * Evolve SMFT fields one timestep (placeholder - would be full Kuramoto)
     */
    void evolveFields(float dt);

    /**
     * Evolve EM fields one timestep using Proca dynamics
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

    // SMFT field storage
    std::vector<float> theta_field_;  // Phase field θ(x,y)
    std::vector<float> R_field_;      // Synchronization field R(x,y)

    // EM support
    bool em_enabled_ = false;
    std::unique_ptr<physics::ProcaEM> proca_em_;
    std::unique_ptr<physics::ProcaVulkanPipeline> proca_pipeline_;
};
