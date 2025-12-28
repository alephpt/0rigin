#pragma once

#include <vulkan/vulkan.h>
#include <memory>
#include <vector>

/**
 * SMFTCompute - Handles GPU compute dispatch operations for SMFT simulation
 *
 * Responsibilities:
 * - Command buffer creation and recording
 * - Pipeline binding and descriptor set binding
 * - Push constant management
 * - Compute shader dispatch
 * - Synchronization between compute operations
 * - Queue submission and fence management
 */
class SMFTCompute {
public:
    SMFTCompute(VkDevice device, VkQueue queue, uint32_t queueFamilyIndex);
    ~SMFTCompute();

    // Initialize command pool
    bool initialize();

    // Begin a new compute batch
    bool beginBatch();

    // Dispatch operations for various compute shaders
    void dispatchKuramoto(VkPipeline pipeline,
                         VkPipelineLayout layout,
                         VkDescriptorSet descriptorSet,
                         const void* pushConstants,
                         size_t pushSize,
                         uint32_t groupCountX,
                         uint32_t groupCountY);

    void dispatchSyncField(VkPipeline pipeline,
                          VkPipelineLayout layout,
                          VkDescriptorSet descriptorSet,
                          const void* pushConstants,
                          size_t pushSize,
                          uint32_t groupCountX,
                          uint32_t groupCountY);

    void dispatchGravityField(VkPipeline pipeline,
                             VkPipelineLayout layout,
                             VkDescriptorSet descriptorSet,
                             const void* pushConstants,
                             size_t pushSize,
                             uint32_t groupCountX,
                             uint32_t groupCountY);

    void dispatchDiracEvolution(VkPipeline pipeline,
                               VkPipelineLayout layout,
                               VkDescriptorSet descriptorSet,
                               const void* pushConstants,
                               size_t pushSize,
                               uint32_t groupCountX,
                               uint32_t groupCountY);

    void dispatchAccumulation(VkPipeline pipeline,
                             VkPipelineLayout layout,
                             VkDescriptorSet descriptorSet,
                             const void* pushConstants,
                             size_t pushSize,
                             uint32_t groupCountX,
                             uint32_t groupCountY);

    // EM Field compute dispatches (Phase 5 - Sprint 3)
    void dispatchEMPotentials(VkPipeline pipeline,
                             VkPipelineLayout layout,
                             VkDescriptorSet descriptorSet,
                             const void* pushConstants,
                             size_t pushSize,
                             uint32_t groupCountX,
                             uint32_t groupCountY);

    void dispatchEMFieldStrengths(VkPipeline pipeline,
                                 VkPipelineLayout layout,
                                 VkDescriptorSet descriptorSet,
                                 const void* pushConstants,
                                 size_t pushSize,
                                 uint32_t groupCountX,
                                 uint32_t groupCountY);

    void dispatchEMReduceEnergy(VkPipeline pipeline,
                               VkPipelineLayout layout,
                               VkDescriptorSet descriptorSet,
                               const void* pushConstants,
                               size_t pushSize,
                               uint32_t groupCountX,
                               uint32_t groupCountY);

    // Buffer operations
    void copyBuffer(VkBuffer src, VkBuffer dst, VkDeviceSize size);
    void fillBuffer(VkBuffer buffer, uint32_t value, VkDeviceSize size);

    // Synchronization
    void insertMemoryBarrier();
    void insertBufferBarrier(VkBuffer buffer,
                           VkAccessFlags srcAccess,
                           VkAccessFlags dstAccess);

    // End batch and submit to queue
    bool submitBatch(bool waitForCompletion = true);

    // Utility to calculate workgroup counts
    static void calculateWorkgroups(uint32_t width, uint32_t height,
                                   uint32_t& workgroupsX,
                                   uint32_t& workgroupsY,
                                   uint32_t localSizeX = 16,
                                   uint32_t localSizeY = 16);

private:
    VkDevice _device;
    VkQueue _queue;
    uint32_t _queueFamilyIndex;

    VkCommandPool _commandPool;
    VkCommandBuffer _commandBuffer;
    VkFence _fence;

    bool _isRecording;

    // Helper to dispatch any compute shader
    void dispatchCompute(VkPipeline pipeline,
                        VkPipelineLayout layout,
                        VkDescriptorSet descriptorSet,
                        const void* pushConstants,
                        size_t pushSize,
                        uint32_t groupCountX,
                        uint32_t groupCountY);

    void cleanup();
};