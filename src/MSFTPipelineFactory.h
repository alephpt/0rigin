#pragma once

#include <vulkan/vulkan.h>
#include <string>
#include <vector>
#include <stdexcept>

/**
 * MSFTPipelineFactory - Factory class for creating Vulkan compute pipelines
 *
 * Extracts pipeline creation logic from MSFTEngine to follow Single Responsibility Principle.
 * Each method is kept under 50 lines to meet code quality standards.
 *
 * This factory handles:
 * - Shader module loading from SPIR-V files
 * - Compute pipeline creation with appropriate layouts
 * - Pipeline lifecycle management (create/destroy)
 */
class MSFTPipelineFactory {
public:
    /**
     * Constructor
     * @param device Vulkan logical device for pipeline creation
     */
    explicit MSFTPipelineFactory(VkDevice device);

    /**
     * Destructor - ensures all created pipelines are destroyed
     */
    ~MSFTPipelineFactory();

    // Pipeline creation methods - each creates a specific compute pipeline

    /**
     * Create Kuramoto dynamics pipeline
     * @param shaderPath Path to kuramoto_step.comp.spv shader
     * @param pipelineLayout Pre-created pipeline layout for Kuramoto
     * @return Created pipeline or VK_NULL_HANDLE on failure
     */
    VkPipeline createKuramotoPipeline(const std::string& shaderPath,
                                      VkPipelineLayout pipelineLayout);

    /**
     * Create synchronization field calculation pipeline
     * @param shaderPath Path to sync_field.comp.spv shader
     * @param pipelineLayout Pre-created pipeline layout for sync field
     * @return Created pipeline or VK_NULL_HANDLE on failure
     */
    VkPipeline createSyncFieldPipeline(const std::string& shaderPath,
                                       VkPipelineLayout pipelineLayout);

    /**
     * Create gravitational field calculation pipeline
     * @param shaderPath Path to gravity_field.comp.spv shader
     * @param pipelineLayout Pre-created pipeline layout for gravity
     * @return Created pipeline or VK_NULL_HANDLE on failure
     */
    VkPipeline createGravityFieldPipeline(const std::string& shaderPath,
                                         VkPipelineLayout pipelineLayout);

    /**
     * Create stochastic Kuramoto evolution pipeline
     * @param shaderPath Path to kuramoto_stochastic.comp.spv shader
     * @param pipelineLayout Pre-created pipeline layout
     * @return Created pipeline or VK_NULL_HANDLE on failure
     */
    VkPipeline createKuramotoStochasticPipeline(const std::string& shaderPath,
                                                VkPipelineLayout pipelineLayout);

    /**
     * Create Dirac evolution pipeline
     * @param shaderPath Path to dirac.comp.spv shader
     * @param pipelineLayout Pre-created pipeline layout for Dirac
     * @return Created pipeline or VK_NULL_HANDLE on failure
     */
    VkPipeline createDiracPipeline(const std::string& shaderPath,
                                   VkPipelineLayout pipelineLayout);

    /**
     * Create stochastic Dirac evolution pipeline
     * @param shaderPath Path to dirac_stochastic.comp.spv shader
     * @param pipelineLayout Pre-created pipeline layout
     * @return Created pipeline or VK_NULL_HANDLE on failure
     */
    VkPipeline createDiracStochasticPipeline(const std::string& shaderPath,
                                             VkPipelineLayout pipelineLayout);

    /**
     * Destroy a specific pipeline
     * @param pipeline Pipeline to destroy
     */
    void destroyPipeline(VkPipeline pipeline);

    /**
     * Destroy all pipelines tracked by this factory
     */
    void destroyAllPipelines();

private:
    VkDevice _device;
    std::vector<VkPipeline> _createdPipelines;  // Track all created pipelines for cleanup

    /**
     * Load SPIR-V shader from file
     * @param path Path to .spv file
     * @return Shader bytecode or empty vector on failure
     */
    std::vector<uint32_t> loadShaderFile(const std::string& path);

    /**
     * Create shader module from SPIR-V bytecode
     * @param code SPIR-V bytecode
     * @return Shader module or VK_NULL_HANDLE on failure
     */
    VkShaderModule createShaderModule(const std::vector<uint32_t>& code);

    /**
     * Generic compute pipeline creation
     * @param shaderModule Compiled shader module
     * @param pipelineLayout Pipeline layout with descriptors and push constants
     * @return Created pipeline or VK_NULL_HANDLE on failure
     */
    VkPipeline createComputePipeline(VkShaderModule shaderModule,
                                    VkPipelineLayout pipelineLayout);

    /**
     * Helper to create pipeline from shader path
     * Combines loading, module creation, and pipeline creation
     * @param shaderPath Path to .spv file
     * @param pipelineLayout Pipeline layout to use
     * @param pipelineName Name for error reporting
     * @return Created pipeline or VK_NULL_HANDLE on failure
     */
    VkPipeline createPipelineFromShader(const std::string& shaderPath,
                                        VkPipelineLayout pipelineLayout,
                                        const std::string& pipelineName);
};