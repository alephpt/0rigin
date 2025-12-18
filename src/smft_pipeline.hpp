#pragma once

#include <vulkan/vulkan.h>
#include <vector>
#include <memory>
#include <cmath>
#include "memory_barriers.hpp"
#include "SMFT_buffers.hpp"
#include "SMFT_params.hpp"

namespace Nova {
namespace SMFT {

class SMFTPipeline {
public:
    SMFTPipeline(
        VkDevice device,
        VkPhysicalDevice physical_device,
        VkCommandPool command_pool,
        VkQueue compute_queue
    );
    ~SMFTPipeline();

    void initialize(
        uint32_t grid_x,
        uint32_t grid_y,
        float delta = 1.0f,
        float chiral_angle = 0.0f
    );

    void setInitialPhases(const std::vector<float>& theta_data);
    void setNaturalFrequencies(const std::vector<float>& omega_data);
    void setInitialSpinorField(const std::vector<float>& psi_data);
    void step(float dt, float coupling = 1.0f, float damping = 0.1f);
    void evolve(int num_steps, float dt, float coupling = 1.0f, float damping = 0.1f);

    std::vector<float> getPhases() const;
    std::vector<float> getSyncField() const;
    std::vector<float> getSpinorField() const;
    std::vector<float> getSpinorDensity() const;
    
    const SMFTParams& getParams() const { return params_; }
    void setParams(const SMFTParams& params);
    float getTime() const { return current_time_; }
    void reset();

private:
    VkDevice device_;
    VkPhysicalDevice physical_device_;
    VkCommandPool command_pool_;
    VkQueue compute_queue_;

    VkPipeline kuramoto_pipeline_;
    VkPipeline sync_pipeline_;
    VkPipeline dirac_pipeline_;
    VkPipeline feedback_pipeline_;
    VkPipelineLayout pipeline_layout_;
    VkDescriptorSetLayout descriptor_set_layout_;
    VkDescriptorPool descriptor_pool_;
    VkDescriptorSet descriptor_set_;

    VkShaderModule kuramoto_shader_;
    VkShaderModule sync_shader_;
    VkShaderModule dirac_shader_;
    VkShaderModule feedback_shader_;

    SMFTBuffers buffers_;
    SMFTParams params_;
    float current_time_ = 0.0f;
    VkCommandBuffer command_buffer_;

    void createBuffers();
    void createShaderModules();
    void createDescriptorSetLayout();
    void createDescriptorPool();
    void createDescriptorSet();
    void createPipelines();
    void createCommandBuffer();
    void updateDescriptorSet();
    void recordCommandBuffer(float dt, float coupling, float damping);
    VkShaderModule loadShaderModule(const char* filepath);
    void copyToBuffer(VkBuffer buffer, const void* data, VkDeviceSize size);
    void copyFromBuffer(VkBuffer buffer, void* data, VkDeviceSize size) const;
    void destroyBuffers();
    void destroyPipelines();
};

} // namespace SMFT
} // namespace Nova