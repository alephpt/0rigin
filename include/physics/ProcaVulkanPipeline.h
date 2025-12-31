// include/physics/ProcaVulkanPipeline.h
#pragma once
#include <vulkan/vulkan.h>
#include <vector>
#include <string>

namespace physics {

/**
 * ProcaVulkanPipeline - GPU compute pipeline for Proca massive photon evolution
 *
 * Implements GPU-accelerated Proca field evolution:
 *   (□ + m²)A_μ = j_μ
 *
 * Where:
 *   - A_μ = (φ, A_x, A_y, A_z) is the massive vector potential
 *   - m = m_γ(x,y) is spatially-varying photon mass
 *   - j_μ is the current source from SMFT Noether current
 */
class ProcaVulkanPipeline {
public:
    ProcaVulkanPipeline(VkDevice device, VkPhysicalDevice physicalDevice);
    ~ProcaVulkanPipeline();

    /**
     * Initialize GPU resources for grid of size nx × ny
     */
    void initialize(uint32_t nx, uint32_t ny);

    /**
     * Execute one timestep of Proca field evolution
     * Dispatches compute shader to evolve (□ + m²)A_μ = j_μ
     */
    void evolveProca(VkCommandBuffer cmd, float dt);

    /**
     * Compute field tensor F_μν from potentials A_μ
     * E = -∇φ - ∂A/∂t, B = ∇×A
     */
    void computeFieldTensor(VkCommandBuffer cmd);

    /**
     * Compute EM field energy: ∫(E² + B²)/2 d³x
     */
    void computeEnergy(VkCommandBuffer cmd);

    /**
     * Enforce Lorenz gauge condition: ∂_μ A^μ = 0
     */
    void enforceLorenzGauge(VkCommandBuffer cmd);

    // Buffer accessors
    VkBuffer getPhiBuffer() const { return phi_buffer_; }
    VkBuffer getAxBuffer() const { return Ax_buffer_; }
    VkBuffer getAyBuffer() const { return Ay_buffer_; }
    VkBuffer getBzBuffer() const { return Bz_buffer_; }

private:
    VkDevice device_;
    VkPhysicalDevice physical_device_;

    uint32_t nx_, ny_;

    // Compute pipelines
    VkPipeline evolve_pipeline_;
    VkPipeline field_tensor_pipeline_;
    VkPipeline energy_pipeline_;
    VkPipeline gauge_pipeline_;

    VkPipelineLayout pipeline_layout_;
    VkDescriptorSetLayout descriptor_layout_;
    VkDescriptorPool descriptor_pool_;
    VkDescriptorSet descriptor_set_;

    // Field buffers: A_μ = (φ, A_x, A_y, A_z)
    VkBuffer phi_buffer_, Ax_buffer_, Ay_buffer_;

    // Field strength buffers: F_μν → (E_x, E_y, B_z)
    VkBuffer Ex_buffer_, Ey_buffer_, Bz_buffer_;

    VkDeviceMemory buffer_memory_;

    void createShaderModule(const std::string& path, VkShaderModule& module);
    void createBuffer(VkDeviceSize size, VkBuffer& buffer, VkDeviceMemory& memory);
};

} // namespace physics
