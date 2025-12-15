#pragma once
#include <vulkan/vulkan.h>
#include "MSFT_buffers.hpp"
#include "MSFT_params.hpp"

namespace MSFT {

class MSFTMemoryBarriers {
public:
    static void recordMSFTPipeline(
        VkCommandBuffer commandBuffer,
        VkPipeline kuramoto_pipeline,
        VkPipeline sync_pipeline,
        VkPipeline dirac_pipeline,
        VkPipeline feedback_pipeline,
        VkPipelineLayout pipeline_layout,
        VkDescriptorSet descriptor_set,
        const MSFTBuffers& buffers,
        const MSFTParams& params,
        uint32_t Nx,
        uint32_t Ny
    ) {
        VkMemoryBarrier barrier = {};
        barrier.sType = VK_STRUCTURE_TYPE_MEMORY_BARRIER;
        barrier.srcAccessMask = VK_ACCESS_SHADER_WRITE_BIT;
        barrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT;

        uint32_t groupCountX = (Nx + 15) / 16;
        uint32_t groupCountY = (Ny + 15) / 16;

        vkCmdPushConstants(
            commandBuffer,
            pipeline_layout,
            VK_SHADER_STAGE_COMPUTE_BIT,
            0,
            sizeof(MSFTParams),
            &params
        );

        // 1. Kuramoto step
        vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, kuramoto_pipeline);
        vkCmdBindDescriptorSets(commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline_layout, 0, 1, &descriptor_set, 0, nullptr);
        vkCmdDispatch(commandBuffer, groupCountX, groupCountY, 1);

        vkCmdPipelineBarrier(commandBuffer,
            VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
            0, 1, &barrier, 0, nullptr, 0, nullptr);

        // 2. Sync field
        vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, sync_pipeline);
        vkCmdDispatch(commandBuffer, groupCountX, groupCountY, 1);

        vkCmdPipelineBarrier(commandBuffer,
            VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
            0, 1, &barrier, 0, nullptr, 0, nullptr);

        // 3. Dirac RK4
        vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, dirac_pipeline);
        vkCmdDispatch(commandBuffer, groupCountX, groupCountY, 1);

        vkCmdPipelineBarrier(commandBuffer,
            VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
            0, 1, &barrier, 0, nullptr, 0, nullptr);
        
        // 4. Spinor feedback
        vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, feedback_pipeline);
        vkCmdDispatch(commandBuffer, groupCountX, groupCountY, 1);

        vkCmdPipelineBarrier(commandBuffer,
            VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
            0, 1, &barrier, 0, nullptr, 0, nullptr);
    }
};

} // namespace MSFT