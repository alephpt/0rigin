#pragma once
#include <vulkan/vulkan.h>
#include <vector>

/**
 * MSFTDescriptorManager - Vulkan Descriptor Management for MSFT Engine
 *
 * Encapsulates all Vulkan descriptor-related operations:
 * - Descriptor pool creation and management
 * - Descriptor set layout creation
 * - Descriptor set allocation
 * - Descriptor set updates with buffer bindings
 *
 * This class handles the verbose Vulkan boilerplate for managing
 * descriptor sets across multiple compute pipelines (kuramoto, sync, gravity, dirac).
 */
class MSFTDescriptorManager {
public:
    /**
     * Constructor
     * @param device Vulkan logical device
     */
    explicit MSFTDescriptorManager(VkDevice device);

    /**
     * Destructor - cleans up all managed descriptor resources
     */
    ~MSFTDescriptorManager();

    // Disable copy/move semantics
    MSFTDescriptorManager(const MSFTDescriptorManager&) = delete;
    MSFTDescriptorManager& operator=(const MSFTDescriptorManager&) = delete;
    MSFTDescriptorManager(MSFTDescriptorManager&&) = delete;
    MSFTDescriptorManager& operator=(MSFTDescriptorManager&&) = delete;

    /**
     * Create descriptor pool for allocating descriptor sets
     * @param maxSets Maximum number of descriptor sets to allocate
     * @param poolSizes Vector of descriptor type/count pairs
     * @return Created descriptor pool handle
     */
    VkDescriptorPool createDescriptorPool(
        uint32_t maxSets,
        const std::vector<VkDescriptorPoolSize>& poolSizes);

    /**
     * Create descriptor set layout for compute shader bindings
     * @param bindings Vector of descriptor set layout bindings
     * @return Created descriptor set layout handle
     */
    VkDescriptorSetLayout createDescriptorSetLayout(
        const std::vector<VkDescriptorSetLayoutBinding>& bindings);

    /**
     * Allocate descriptor set from pool
     * @param pool Descriptor pool to allocate from
     * @param layout Descriptor set layout to use
     * @return Allocated descriptor set handle
     */
    VkDescriptorSet allocateDescriptorSet(
        VkDescriptorPool pool,
        VkDescriptorSetLayout layout);

    /**
     * Update descriptor set with buffer bindings
     * @param descriptorSet Descriptor set to update
     * @param buffers Vector of buffer handles to bind
     */
    void updateDescriptorSet(
        VkDescriptorSet descriptorSet,
        const std::vector<VkBuffer>& buffers);

    /**
     * Create binding structure for descriptor set layout
     * @param binding Binding index
     * @param type Descriptor type (e.g., STORAGE_BUFFER)
     * @return Configured descriptor set layout binding
     */
    VkDescriptorSetLayoutBinding createBinding(
        uint32_t binding,
        VkDescriptorType type);

    /**
     * Destroy descriptor pool
     * @param pool Pool to destroy
     */
    void destroyDescriptorPool(VkDescriptorPool pool);

    /**
     * Destroy descriptor set layout
     * @param layout Layout to destroy
     */
    void destroyDescriptorSetLayout(VkDescriptorSetLayout layout);

    /**
     * Cleanup all managed resources
     */
    void destroyAll();

private:
    VkDevice _device;

    // Track all created resources for cleanup
    std::vector<VkDescriptorPool> _pools;
    std::vector<VkDescriptorSetLayout> _layouts;
};
