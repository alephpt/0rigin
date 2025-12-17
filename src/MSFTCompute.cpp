#include "MSFTCompute.h"
#include <stdexcept>
#include <cstring>

MSFTCompute::MSFTCompute(VkDevice device, VkQueue queue, uint32_t queueFamilyIndex)
    : _device(device)
    , _queue(queue)
    , _queueFamilyIndex(queueFamilyIndex)
    , _commandPool(VK_NULL_HANDLE)
    , _commandBuffer(VK_NULL_HANDLE)
    , _fence(VK_NULL_HANDLE)
    , _isRecording(false) {
}

MSFTCompute::~MSFTCompute() {
    cleanup();
}

bool MSFTCompute::initialize() {
    // Create command pool
    VkCommandPoolCreateInfo poolInfo{};
    poolInfo.sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
    poolInfo.queueFamilyIndex = _queueFamilyIndex;
    poolInfo.flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;

    if (vkCreateCommandPool(_device, &poolInfo, nullptr, &_commandPool) != VK_SUCCESS) {
        return false;
    }

    // Allocate command buffer
    VkCommandBufferAllocateInfo allocInfo{};
    allocInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;
    allocInfo.commandPool = _commandPool;
    allocInfo.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
    allocInfo.commandBufferCount = 1;

    if (vkAllocateCommandBuffers(_device, &allocInfo, &_commandBuffer) != VK_SUCCESS) {
        vkDestroyCommandPool(_device, _commandPool, nullptr);
        _commandPool = VK_NULL_HANDLE;
        return false;
    }

    // Create fence for synchronization
    VkFenceCreateInfo fenceInfo{};
    fenceInfo.sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO;

    if (vkCreateFence(_device, &fenceInfo, nullptr, &_fence) != VK_SUCCESS) {
        cleanup();
        return false;
    }

    return true;
}

bool MSFTCompute::beginBatch() {
    if (_isRecording || _commandBuffer == VK_NULL_HANDLE) {
        return false;
    }

    VkCommandBufferBeginInfo beginInfo{};
    beginInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    beginInfo.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;

    if (vkBeginCommandBuffer(_commandBuffer, &beginInfo) != VK_SUCCESS) {
        return false;
    }

    _isRecording = true;
    return true;
}

void MSFTCompute::dispatchCompute(VkPipeline pipeline,
                                 VkPipelineLayout layout,
                                 VkDescriptorSet descriptorSet,
                                 const void* pushConstants,
                                 size_t pushSize,
                                 uint32_t groupCountX,
                                 uint32_t groupCountY) {
    if (!_isRecording) return;

    // Bind pipeline
    vkCmdBindPipeline(_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline);

    // Bind descriptor set
    if (descriptorSet != VK_NULL_HANDLE) {
        vkCmdBindDescriptorSets(_commandBuffer,
                              VK_PIPELINE_BIND_POINT_COMPUTE,
                              layout, 0, 1, &descriptorSet, 0, nullptr);
    }

    // Push constants
    if (pushConstants && pushSize > 0) {
        vkCmdPushConstants(_commandBuffer, layout,
                          VK_SHADER_STAGE_COMPUTE_BIT,
                          0, pushSize, pushConstants);
    }

    // Dispatch compute
    vkCmdDispatch(_commandBuffer, groupCountX, groupCountY, 1);
}

void MSFTCompute::dispatchKuramoto(VkPipeline pipeline,
                                  VkPipelineLayout layout,
                                  VkDescriptorSet descriptorSet,
                                  const void* pushConstants,
                                  size_t pushSize,
                                  uint32_t groupCountX,
                                  uint32_t groupCountY) {
    dispatchCompute(pipeline, layout, descriptorSet,
                   pushConstants, pushSize, groupCountX, groupCountY);
}

void MSFTCompute::dispatchSyncField(VkPipeline pipeline,
                                   VkPipelineLayout layout,
                                   VkDescriptorSet descriptorSet,
                                   const void* pushConstants,
                                   size_t pushSize,
                                   uint32_t groupCountX,
                                   uint32_t groupCountY) {
    dispatchCompute(pipeline, layout, descriptorSet,
                   pushConstants, pushSize, groupCountX, groupCountY);
}

void MSFTCompute::dispatchGravityField(VkPipeline pipeline,
                                      VkPipelineLayout layout,
                                      VkDescriptorSet descriptorSet,
                                      const void* pushConstants,
                                      size_t pushSize,
                                      uint32_t groupCountX,
                                      uint32_t groupCountY) {
    dispatchCompute(pipeline, layout, descriptorSet,
                   pushConstants, pushSize, groupCountX, groupCountY);
}

void MSFTCompute::dispatchDiracEvolution(VkPipeline pipeline,
                                        VkPipelineLayout layout,
                                        VkDescriptorSet descriptorSet,
                                        const void* pushConstants,
                                        size_t pushSize,
                                        uint32_t groupCountX,
                                        uint32_t groupCountY) {
    dispatchCompute(pipeline, layout, descriptorSet,
                   pushConstants, pushSize, groupCountX, groupCountY);
}

void MSFTCompute::copyBuffer(VkBuffer src, VkBuffer dst, VkDeviceSize size) {
    if (!_isRecording) return;

    VkBufferCopy copyRegion{};
    copyRegion.size = size;
    vkCmdCopyBuffer(_commandBuffer, src, dst, 1, &copyRegion);
}

void MSFTCompute::insertMemoryBarrier() {
    if (!_isRecording) return;

    VkMemoryBarrier barrier{};
    barrier.sType = VK_STRUCTURE_TYPE_MEMORY_BARRIER;
    barrier.srcAccessMask = VK_ACCESS_SHADER_WRITE_BIT;
    barrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT;

    vkCmdPipelineBarrier(_commandBuffer,
                        VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                        VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                        0, 1, &barrier, 0, nullptr, 0, nullptr);
}

void MSFTCompute::insertBufferBarrier(VkBuffer buffer,
                                     VkAccessFlags srcAccess,
                                     VkAccessFlags dstAccess) {
    if (!_isRecording) return;

    VkBufferMemoryBarrier barrier{};
    barrier.sType = VK_STRUCTURE_TYPE_BUFFER_MEMORY_BARRIER;
    barrier.srcAccessMask = srcAccess;
    barrier.dstAccessMask = dstAccess;
    barrier.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    barrier.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    barrier.buffer = buffer;
    barrier.offset = 0;
    barrier.size = VK_WHOLE_SIZE;

    vkCmdPipelineBarrier(_commandBuffer,
                        VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                        VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                        0, 0, nullptr, 1, &barrier, 0, nullptr);
}

bool MSFTCompute::submitBatch(bool waitForCompletion) {
    if (!_isRecording) return false;

    // End command buffer recording
    if (vkEndCommandBuffer(_commandBuffer) != VK_SUCCESS) {
        _isRecording = false;
        return false;
    }

    _isRecording = false;

    // Reset fence
    vkResetFences(_device, 1, &_fence);

    // Submit command buffer
    VkSubmitInfo submitInfo{};
    submitInfo.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;
    submitInfo.commandBufferCount = 1;
    submitInfo.pCommandBuffers = &_commandBuffer;

    VkResult result = vkQueueSubmit(_queue, 1, &submitInfo,
                                   waitForCompletion ? _fence : VK_NULL_HANDLE);

    if (result != VK_SUCCESS) {
        return false;
    }

    // Wait for completion if requested
    if (waitForCompletion) {
        vkWaitForFences(_device, 1, &_fence, VK_TRUE, UINT64_MAX);
    }

    return true;
}

void MSFTCompute::calculateWorkgroups(uint32_t width, uint32_t height,
                                     uint32_t& workgroupsX,
                                     uint32_t& workgroupsY,
                                     uint32_t localSizeX,
                                     uint32_t localSizeY) {
    workgroupsX = (width + localSizeX - 1) / localSizeX;
    workgroupsY = (height + localSizeY - 1) / localSizeY;
}

void MSFTCompute::cleanup() {
    if (_device == VK_NULL_HANDLE) return;

    if (_fence != VK_NULL_HANDLE) {
        vkDestroyFence(_device, _fence, nullptr);
        _fence = VK_NULL_HANDLE;
    }

    if (_commandBuffer != VK_NULL_HANDLE && _commandPool != VK_NULL_HANDLE) {
        vkFreeCommandBuffers(_device, _commandPool, 1, &_commandBuffer);
        _commandBuffer = VK_NULL_HANDLE;
    }

    if (_commandPool != VK_NULL_HANDLE) {
        vkDestroyCommandPool(_device, _commandPool, nullptr);
        _commandPool = VK_NULL_HANDLE;
    }
}