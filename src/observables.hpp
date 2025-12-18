#pragma once

#include <vulkan/vulkan.h>
#include <vector>
#include <string>
#include <cmath>
#include "SMFT_buffers.hpp"

namespace SMFT {

/**
 * @file observables.hpp
 * @brief Comprehensive observables computation for SMFT+Dirac system
 *
 * OBSERVABLE CATEGORIES:
 * ======================
 * 1. Synchronization Metrics (Kuramoto order parameters)
 * 2. Spinor Properties (density, current, spin, chirality)
 * 3. Field Theory Observables (mass, energy, dispersion)
 * 4. Conservation Laws (probability, energy, norm)
 *
 * GPU ACCELERATION:
 * - All observables computed via Vulkan compute shaders
 * - Reduction kernels for global quantities
 * - Parallel computation for local field values
 */

/**
 * @brief Synchronization observables
 */
struct SyncObservables {
    float R_global;           // Global order parameter: |⟨e^{iθ}⟩|
    float phase_coherence;    // Phase coherence: ⟨cos(θ_i - θ_j)⟩
    float local_sync_mean;    // Mean local synchronization
    float local_sync_std;     // Std dev of local synchronization
    std::vector<float> R_field;          // Local R(x,y) field
    std::vector<float> sync_gradient;    // |∇R(x,y)|
};

/**
 * @brief Spinor observables
 */
struct SpinorObservables {
    float total_density;      // ∫ ρ d²x where ρ = |ψ|²
    float norm_squared;       // ||ψ||² = ∫ ψ†ψ d²x
    float chirality;          // ⟨ψ|γ⁵|ψ⟩

    std::vector<float> density;          // ρ(x,y) = |ψ|²
    std::vector<float> current_x;        // j^x = ψ†γ^x ψ
    std::vector<float> current_y;        // j^y = ψ†γ^y ψ
    std::vector<float> spin_density;     // S(x,y) = ψ†Σψ
    std::vector<float> chiral_density;   // χ(x,y) = ψ†γ⁵ψ
};

/**
 * @brief Field theory observables
 */
struct FieldTheoryObservables {
    float effective_mass_mean;   // ⟨m_eff⟩ where m_eff = Δ·R
    float effective_mass_max;    // max(m_eff)
    float hamiltonian;           // Total energy H
    float kinetic_energy;        // T = ∫ ψ†(-iγ·∇)ψ d²x
    float potential_energy;      // V = ∫ ψ†(m·𝟙)ψ d²x
    float localization;          // Localization measure: IPR

    std::vector<float> mass_field;       // m_eff(x,y) = Δ·R(x,y)
    std::vector<float> energy_density;   // H(x,y) local energy
    std::vector<float> mass_gradient;    // |∇m(x,y)|
};

/**
 * @brief Conservation law checks
 */
struct ConservationObservables {
    float probability_conservation;  // ∂_t ρ + ∇·j (should be ~0)
    float norm_deviation;            // |1 - ||ψ||²| (should be small)
    float energy_drift;              // ΔH/H over time
    float charge_total;              // Total charge Q = ∫ ρ d²x
};

/**
 * @brief Complete observable set
 */
struct AllObservables {
    SyncObservables sync;
    SpinorObservables spinor;
    FieldTheoryObservables field_theory;
    ConservationObservables conservation;

    float simulation_time;
    uint32_t timestep;
};

/**
 * @brief GPU-accelerated observables computation engine
 */
class ObservablesEngine {
public:
    /**
     * @brief Constructor
     */
    ObservablesEngine(
        VkDevice device,
        VkPhysicalDevice physical_device,
        VkCommandPool command_pool,
        VkQueue compute_queue
    );

    /**
     * @brief Destructor
     */
    ~ObservablesEngine();

    /**
     * @brief Initialize observables engine
     */
    void initialize(uint32_t grid_x, uint32_t grid_y, float dx = 1.0f);

    /**
     * @brief Compute all observables from current SMFT state
     *
     * @param buffers SMFT GPU buffers (theta, R, psi, density)
     * @param Delta Mass gap parameter
     * @param time Current simulation time
     * @param timestep Current timestep number
     * @return Complete observable set
     */
    AllObservables compute(
        const SMFTBuffers& buffers,
        float Delta,
        float time,
        uint32_t timestep
    );

    /**
     * @brief Compute synchronization observables only
     */
    SyncObservables computeSyncObservables(const SMFTBuffers& buffers);

    /**
     * @brief Compute spinor observables only
     */
    SpinorObservables computeSpinorObservables(const SMFTBuffers& buffers);

    /**
     * @brief Compute field theory observables only
     */
    FieldTheoryObservables computeFieldTheoryObservables(
        const SMFTBuffers& buffers,
        float Delta
    );

    /**
     * @brief Compute conservation law checks
     */
    ConservationObservables computeConservationObservables(
        const SMFTBuffers& buffers,
        const AllObservables& previous,
        float dt
    );

private:
    // Vulkan handles
    VkDevice device_;
    VkPhysicalDevice physical_device_;
    VkCommandPool command_pool_;
    VkQueue compute_queue_;

    // Grid parameters
    uint32_t grid_x_;
    uint32_t grid_y_;
    uint32_t N_total_;
    float dx_;

    // Compute pipelines
    VkPipeline sync_global_pipeline_;
    VkPipeline spinor_density_pipeline_;
    VkPipeline current_pipeline_;
    VkPipeline mass_field_pipeline_;
    VkPipeline energy_density_pipeline_;
    VkPipeline reduction_pipeline_;

    // Pipeline layout
    VkPipelineLayout pipeline_layout_;
    VkDescriptorSetLayout descriptor_set_layout_;
    VkDescriptorPool descriptor_pool_;
    VkDescriptorSet descriptor_set_;

    // Shader modules
    VkShaderModule sync_global_shader_;
    VkShaderModule spinor_density_shader_;
    VkShaderModule current_shader_;
    VkShaderModule mass_field_shader_;
    VkShaderModule energy_density_shader_;
    VkShaderModule reduction_shader_;

    // Reduction buffers (for parallel reduction)
    VkBuffer reduction_buffer_;
    VkDeviceMemory reduction_memory_;
    VkDeviceSize reduction_size_;

    // Command buffer
    VkCommandBuffer command_buffer_;

    // Helper methods
    void createShaderModules();
    void createDescriptorSetLayout();
    void createDescriptorPool();
    void createDescriptorSet();
    void createPipelines();
    void createCommandBuffer();
    void createReductionBuffer();

    void destroyPipelines();
    void destroyBuffers();

    VkShaderModule loadShaderModule(const char* filepath);

    // Reduction operations
    float parallelReduce(VkBuffer input, uint32_t count, const char* operation);
    float computeMean(const std::vector<float>& data);
    float computeStdDev(const std::vector<float>& data, float mean);
    float computeGradientMagnitude(const std::vector<float>& field, uint32_t x, uint32_t y);

    // Copy operations
    void copyFromBuffer(VkBuffer buffer, void* data, VkDeviceSize size) const;
    std::vector<float> readBuffer(VkBuffer buffer, VkDeviceSize size) const;
};

/**
 * @brief Observable time series tracker
 */
class ObservableTimeSeries {
public:
    struct TimePoint {
        float time;
        AllObservables obs;
    };

    /**
     * @brief Add observation at time point
     */
    void addObservation(float time, const AllObservables& obs);

    /**
     * @brief Get all observations
     */
    const std::vector<TimePoint>& getTimeSeries() const { return series_; }

    /**
     * @brief Get observable at specific index
     */
    const AllObservables& getObservable(size_t index) const {
        return series_[index].obs;
    }

    /**
     * @brief Get number of observations
     */
    size_t size() const { return series_.size(); }

    /**
     * @brief Clear all observations
     */
    void clear() { series_.clear(); }

    /**
     * @brief Export to CSV format
     */
    std::string toCSV() const;

    /**
     * @brief Export to JSON format
     */
    std::string toJSON() const;

private:
    std::vector<TimePoint> series_;
};

} // namespace SMFT
} // namespace Nova