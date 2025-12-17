#pragma once
#include <vulkan/vulkan.h>
#include <vector>
#include <utility>
#include <cstring>

/**
 * MSFTBufferManager - Manages Vulkan buffer creation and memory allocation
 *
 * Extracts buffer management logic from MSFTEngine to reduce God Object complexity.
 * Handles creation, memory allocation, data transfer, and cleanup of Vulkan buffers.
 */
class MSFTBufferManager {
public:
    /**
     * Constructor
     * @param device Vulkan logical device handle
     * @param physicalDevice Vulkan physical device handle
     */
    MSFTBufferManager(VkDevice device, VkPhysicalDevice physicalDevice);

    /**
     * Destructor - automatically cleans up all managed buffers
     */
    ~MSFTBufferManager();

    /**
     * Create a buffer with specified size and usage
     * @param size Buffer size in bytes
     * @param usage Vulkan buffer usage flags
     * @param properties Memory property flags
     * @return Pair of (buffer, memory) handles
     */
    std::pair<VkBuffer, VkDeviceMemory> createBuffer(
        VkDeviceSize size,
        VkBufferUsageFlags usage,
        VkMemoryPropertyFlags properties
    );

    /**
     * Create storage buffer for compute shader access
     * @param size Buffer size in bytes
     * @return Pair of (buffer, memory) handles
     */
    std::pair<VkBuffer, VkDeviceMemory> createStorageBuffer(VkDeviceSize size);

    /**
     * Upload data to a buffer
     * @param memory Device memory handle
     * @param data Pointer to source data
     * @param size Size in bytes to copy
     */
    void uploadData(VkDeviceMemory memory, const void* data, VkDeviceSize size);

    /**
     * Download data from a buffer
     * @param memory Device memory handle
     * @param data Pointer to destination
     * @param size Size in bytes to copy
     */
    void downloadData(VkDeviceMemory memory, void* data, VkDeviceSize size);

    /**
     * Map memory for direct access
     * @param memory Device memory handle
     * @param size Size to map
     * @return Pointer to mapped memory
     */
    void* mapMemory(VkDeviceMemory memory, VkDeviceSize size);

    /**
     * Unmap previously mapped memory
     * @param memory Device memory handle
     */
    void unmapMemory(VkDeviceMemory memory);

    /**
     * Initialize buffer with zeros
     * @param memory Device memory handle
     * @param size Size in bytes
     */
    void clearBuffer(VkDeviceMemory memory, VkDeviceSize size);

    /**
     * Destroy a specific buffer and free its memory
     * @param buffer Buffer handle
     * @param memory Memory handle
     */
    void destroyBuffer(VkBuffer buffer, VkDeviceMemory memory);

    /**
     * Destroy all managed buffers
     */
    void destroyAllBuffers();

    /**
     * Get device handle
     * @return VkDevice handle
     */
    VkDevice getDevice() const { return _device; }

private:
    VkDevice _device;
    VkPhysicalDevice _physicalDevice;

    // Track allocated resources for cleanup
    std::vector<std::pair<VkBuffer, VkDeviceMemory>> _managedBuffers;

    /**
     * Find suitable memory type index
     * @param typeFilter Memory type filter bits
     * @param properties Required memory properties
     * @return Memory type index
     */
    uint32_t findMemoryType(uint32_t typeFilter, VkMemoryPropertyFlags properties);

    /**
     * Allocate memory for a buffer
     * @param buffer Buffer to allocate memory for
     * @param properties Memory property flags
     * @return Allocated device memory handle
     */
    VkDeviceMemory allocateBufferMemory(VkBuffer buffer, VkMemoryPropertyFlags properties);
};