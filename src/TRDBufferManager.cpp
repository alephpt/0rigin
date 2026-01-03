#include "TRDBufferManager.h"
#include <stdexcept>
#include <algorithm>

TRDBufferManager::TRDBufferManager(VkDevice device, VkPhysicalDevice physicalDevice)
    : _device(device), _physicalDevice(physicalDevice) {
}

TRDBufferManager::~TRDBufferManager() {
    destroyAllBuffers();
}

std::pair<VkBuffer, VkDeviceMemory> TRDBufferManager::createBuffer(
    VkDeviceSize size,
    VkBufferUsageFlags usage,
    VkMemoryPropertyFlags properties) {

    VkBuffer buffer = VK_NULL_HANDLE;
    VkDeviceMemory bufferMemory = VK_NULL_HANDLE;

    // Create buffer
    VkBufferCreateInfo bufferInfo{};
    bufferInfo.sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;
    bufferInfo.size = size;
    bufferInfo.usage = usage;
    bufferInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;

    VkResult result = vkCreateBuffer(_device, &bufferInfo, nullptr, &buffer);
    if (result != VK_SUCCESS) {
        throw std::runtime_error("Failed to create buffer");
    }

    // Allocate memory
    bufferMemory = allocateBufferMemory(buffer, properties);

    // Bind memory to buffer
    vkBindBufferMemory(_device, buffer, bufferMemory, 0);

    // Track for cleanup
    _managedBuffers.push_back({buffer, bufferMemory});

    return {buffer, bufferMemory};
}

std::pair<VkBuffer, VkDeviceMemory> TRDBufferManager::createStorageBuffer(VkDeviceSize size) {
    // Standard storage buffer for compute shader access
    VkBufferUsageFlags usage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT |
                               VK_BUFFER_USAGE_TRANSFER_SRC_BIT |
                               VK_BUFFER_USAGE_TRANSFER_DST_BIT;

    VkMemoryPropertyFlags properties = VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT |
                                       VK_MEMORY_PROPERTY_HOST_COHERENT_BIT;

    return createBuffer(size, usage, properties);
}

VkDeviceMemory TRDBufferManager::allocateBufferMemory(
    VkBuffer buffer,
    VkMemoryPropertyFlags properties) {

    // Get memory requirements
    VkMemoryRequirements memRequirements;
    vkGetBufferMemoryRequirements(_device, buffer, &memRequirements);

    // Allocate memory
    VkMemoryAllocateInfo allocInfo{};
    allocInfo.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
    allocInfo.allocationSize = memRequirements.size;
    allocInfo.memoryTypeIndex = findMemoryType(memRequirements.memoryTypeBits, properties);

    VkDeviceMemory bufferMemory;
    VkResult result = vkAllocateMemory(_device, &allocInfo, nullptr, &bufferMemory);
    if (result != VK_SUCCESS) {
        vkDestroyBuffer(_device, buffer, nullptr);
        throw std::runtime_error("Failed to allocate buffer memory");
    }

    return bufferMemory;
}

uint32_t TRDBufferManager::findMemoryType(
    uint32_t typeFilter,
    VkMemoryPropertyFlags properties) {

    VkPhysicalDeviceMemoryProperties memProperties;
    vkGetPhysicalDeviceMemoryProperties(_physicalDevice, &memProperties);

    for (uint32_t i = 0; i < memProperties.memoryTypeCount; i++) {
        bool typeMatches = (typeFilter & (1 << i));
        bool propertiesMatch = (memProperties.memoryTypes[i].propertyFlags & properties) == properties;

        if (typeMatches && propertiesMatch) {
            return i;
        }
    }

    throw std::runtime_error("Failed to find suitable memory type");
}

void TRDBufferManager::uploadData(
    VkDeviceMemory memory,
    const void* data,
    VkDeviceSize size) {

    void* mappedMemory = mapMemory(memory, size);
    memcpy(mappedMemory, data, size);
    unmapMemory(memory);
}

void TRDBufferManager::downloadData(
    VkDeviceMemory memory,
    void* data,
    VkDeviceSize size) {

    void* mappedMemory = mapMemory(memory, size);
    memcpy(data, mappedMemory, size);
    unmapMemory(memory);
}

void* TRDBufferManager::mapMemory(VkDeviceMemory memory, VkDeviceSize size) {
    void* mappedMemory = nullptr;
    VkResult result = vkMapMemory(_device, memory, 0, size, 0, &mappedMemory);
    if (result != VK_SUCCESS) {
        throw std::runtime_error("Failed to map memory");
    }
    return mappedMemory;
}

void TRDBufferManager::unmapMemory(VkDeviceMemory memory) {
    vkUnmapMemory(_device, memory);
}

void TRDBufferManager::clearBuffer(VkDeviceMemory memory, VkDeviceSize size) {
    void* mappedMemory = mapMemory(memory, size);
    memset(mappedMemory, 0, size);
    unmapMemory(memory);
}

std::vector<std::pair<VkBuffer, VkDeviceMemory>> TRDBufferManager::createAccumulatorBuffers(VkDeviceSize size) {
    std::vector<std::pair<VkBuffer, VkDeviceMemory>> accumulators;

    // Create theta_sum accumulator buffer
    auto theta_sum = createStorageBuffer(size);
    clearBuffer(theta_sum.second, size);  // Initialize to zero
    accumulators.push_back(theta_sum);

    // Create R_sum accumulator buffer
    auto R_sum = createStorageBuffer(size);
    clearBuffer(R_sum.second, size);  // Initialize to zero
    accumulators.push_back(R_sum);

    return accumulators;
}

void TRDBufferManager::destroyBuffer(VkBuffer buffer, VkDeviceMemory memory) {
    // Remove from managed buffers list
    auto it = std::remove_if(_managedBuffers.begin(), _managedBuffers.end(),
        [buffer, memory](const std::pair<VkBuffer, VkDeviceMemory>& pair) {
            return pair.first == buffer && pair.second == memory;
        });
    _managedBuffers.erase(it, _managedBuffers.end());

    // Destroy Vulkan resources
    if (buffer != VK_NULL_HANDLE) {
        vkDestroyBuffer(_device, buffer, nullptr);
    }
    if (memory != VK_NULL_HANDLE) {
        vkFreeMemory(_device, memory, nullptr);
    }
}

void TRDBufferManager::destroyAllBuffers() {
    for (const auto& [buffer, memory] : _managedBuffers) {
        if (buffer != VK_NULL_HANDLE) {
            vkDestroyBuffer(_device, buffer, nullptr);
        }
        if (memory != VK_NULL_HANDLE) {
            vkFreeMemory(_device, memory, nullptr);
        }
    }
    _managedBuffers.clear();
}