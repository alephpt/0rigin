#include "SMFTDescriptorManager.h"
#include <stdexcept>
#include <algorithm>

SMFTDescriptorManager::SMFTDescriptorManager(VkDevice device)
    : _device(device) {
}

SMFTDescriptorManager::~SMFTDescriptorManager() {
    destroyAll();
}

VkDescriptorPool SMFTDescriptorManager::createDescriptorPool(
    uint32_t maxSets,
    const std::vector<VkDescriptorPoolSize>& poolSizes) {

    VkDescriptorPoolCreateInfo poolInfo{};
    poolInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
    poolInfo.poolSizeCount = static_cast<uint32_t>(poolSizes.size());
    poolInfo.pPoolSizes = poolSizes.data();
    poolInfo.maxSets = maxSets;

    VkDescriptorPool pool = VK_NULL_HANDLE;
    if (vkCreateDescriptorPool(_device, &poolInfo, nullptr, &pool) != VK_SUCCESS) {
        throw std::runtime_error("Failed to create descriptor pool");
    }

    _pools.push_back(pool);
    return pool;
}

VkDescriptorSetLayout SMFTDescriptorManager::createDescriptorSetLayout(
    const std::vector<VkDescriptorSetLayoutBinding>& bindings) {

    VkDescriptorSetLayoutCreateInfo layoutInfo{};
    layoutInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
    layoutInfo.bindingCount = static_cast<uint32_t>(bindings.size());
    layoutInfo.pBindings = bindings.data();

    VkDescriptorSetLayout layout = VK_NULL_HANDLE;
    if (vkCreateDescriptorSetLayout(_device, &layoutInfo, nullptr, &layout) != VK_SUCCESS) {
        throw std::runtime_error("Failed to create descriptor set layout");
    }

    _layouts.push_back(layout);
    return layout;
}

VkDescriptorSet SMFTDescriptorManager::allocateDescriptorSet(
    VkDescriptorPool pool,
    VkDescriptorSetLayout layout) {

    VkDescriptorSetAllocateInfo allocInfo{};
    allocInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO;
    allocInfo.descriptorPool = pool;
    allocInfo.descriptorSetCount = 1;
    allocInfo.pSetLayouts = &layout;

    VkDescriptorSet descriptorSet = VK_NULL_HANDLE;
    if (vkAllocateDescriptorSets(_device, &allocInfo, &descriptorSet) != VK_SUCCESS) {
        throw std::runtime_error("Failed to allocate descriptor set");
    }

    return descriptorSet;
}

void SMFTDescriptorManager::updateDescriptorSet(
    VkDescriptorSet descriptorSet,
    const std::vector<VkBuffer>& buffers) {

    // Create buffer info structures
    std::vector<VkDescriptorBufferInfo> bufferInfos;
    bufferInfos.reserve(buffers.size());
    for (const auto& buffer : buffers) {
        VkDescriptorBufferInfo info{};
        info.buffer = buffer;
        info.offset = 0;
        info.range = VK_WHOLE_SIZE;
        bufferInfos.push_back(info);
    }

    // Create write descriptor structures
    std::vector<VkWriteDescriptorSet> writes;
    writes.reserve(buffers.size());
    for (size_t i = 0; i < buffers.size(); ++i) {
        VkWriteDescriptorSet write{};
        write.sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
        write.dstSet = descriptorSet;
        write.dstBinding = static_cast<uint32_t>(i);
        write.dstArrayElement = 0;
        write.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        write.descriptorCount = 1;
        write.pBufferInfo = &bufferInfos[i];
        writes.push_back(write);
    }

    vkUpdateDescriptorSets(_device, static_cast<uint32_t>(writes.size()),
                          writes.data(), 0, nullptr);
}

VkDescriptorSetLayoutBinding SMFTDescriptorManager::createBinding(
    uint32_t binding,
    VkDescriptorType type) {

    VkDescriptorSetLayoutBinding bindingInfo{};
    bindingInfo.binding = binding;
    bindingInfo.descriptorType = type;
    bindingInfo.descriptorCount = 1;
    bindingInfo.stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
    bindingInfo.pImmutableSamplers = nullptr;

    return bindingInfo;
}

void SMFTDescriptorManager::destroyDescriptorPool(VkDescriptorPool pool) {
    if (pool != VK_NULL_HANDLE) {
        vkDestroyDescriptorPool(_device, pool, nullptr);

        // Remove from tracking vector
        auto it = std::find(_pools.begin(), _pools.end(), pool);
        if (it != _pools.end()) {
            _pools.erase(it);
        }
    }
}

void SMFTDescriptorManager::destroyDescriptorSetLayout(VkDescriptorSetLayout layout) {
    if (layout != VK_NULL_HANDLE) {
        vkDestroyDescriptorSetLayout(_device, layout, nullptr);

        // Remove from tracking vector
        auto it = std::find(_layouts.begin(), _layouts.end(), layout);
        if (it != _layouts.end()) {
            _layouts.erase(it);
        }
    }
}

void SMFTDescriptorManager::destroyAll() {
    // Destroy all tracked pools
    for (auto pool : _pools) {
        if (pool != VK_NULL_HANDLE) {
            vkDestroyDescriptorPool(_device, pool, nullptr);
        }
    }
    _pools.clear();

    // Destroy all tracked layouts
    for (auto layout : _layouts) {
        if (layout != VK_NULL_HANDLE) {
            vkDestroyDescriptorSetLayout(_device, layout, nullptr);
        }
    }
    _layouts.clear();
}
