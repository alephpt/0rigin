#include "SMFT_pipeline.hpp"
#include <fstream>
#include <iostream>
#include <cstring>
#include <stdexcept>
#include <algorithm>

namespace Nova {
namespace SMFT {

// ============================================================================
// SMFTPipeline Implementation
// ============================================================================

SMFTPipeline::SMFTPipeline(
    VkDevice device,
    VkPhysicalDevice physical_device,
    VkCommandPool command_pool,
    VkQueue compute_queue
) : device_(device),
    physical_device_(physical_device),
    command_pool_(command_pool),
    compute_queue_(compute_queue),
    kuramoto_pipeline_(VK_NULL_HANDLE),
    sync_pipeline_(VK_NULL_HANDLE),
    dirac_pipeline_(VK_NULL_HANDLE),
    feedback_pipeline_(VK_NULL_HANDLE),
    pipeline_layout_(VK_NULL_HANDLE),
    descriptor_set_layout_(VK_NULL_HANDLE),
    descriptor_pool_(VK_NULL_HANDLE),
    descriptor_set_(VK_NULL_HANDLE),
    kuramoto_shader_(VK_NULL_HANDLE),
    sync_shader_(VK_NULL_HANDLE),
    dirac_shader_(VK_NULL_HANDLE),
    feedback_shader_(VK_NULL_HANDLE),
    command_buffer_(VK_NULL_HANDLE),
    current_time_(0.0f) {
}

SMFTPipeline::~SMFTPipeline() {
    destroyPipelines();
    destroyBuffers();

    if (command_buffer_ != VK_NULL_HANDLE) {
        vkFreeCommandBuffers(device_, command_pool_, 1, &command_buffer_);
    }

    if (descriptor_pool_ != VK_NULL_HANDLE) {
        vkDestroyDescriptorPool(device_, descriptor_pool_, nullptr);
    }

    if (descriptor_set_layout_ != VK_NULL_HANDLE) {
        vkDestroyDescriptorSetLayout(device_, descriptor_set_layout_, nullptr);
    }
}

void SMFTPipeline::initialize(
    uint32_t grid_x,
    uint32_t grid_y,
    float delta,
    float chiral_angle
) {
    params_.Nx = grid_x;
    params_.Ny = grid_y;
    params_.N_total = grid_x * grid_y;
    params_.Delta = delta;
    params_.chiral_angle = chiral_angle;

    // Auto-configure operator splitting (Dirac sub-stepping)
    params_.auto_configure_substeps();

    std::cout << "[SMFT] Operator splitting: Dirac sub-steps=" << params_.dirac_substeps
              << " (dt_dirac=" << (params_.dt / params_.dirac_substeps) << ")" << std::endl;

    // Validate timestep
    if (!params_.validate_dt()) {
        std::cout << "[SMFT] Warning: dt=" << params_.dt
                  << " may be unsafe for Δ=" << delta
                  << ". Recommended: dt<=" << params_.compute_safe_dt() << std::endl;
    }

    createBuffers();
    createShaderModules();
    createDescriptorSetLayout();
    createDescriptorPool();
    createDescriptorSet();
    createPipelines();
    createCommandBuffer();

    std::cout << "[SMFT] Pipeline initialized: " << grid_x << "x" << grid_y
              << " grid, Δ=" << delta << ", dt=" << params_.dt << std::endl;
}

void SMFTPipeline::createBuffers() {
    VkDeviceSize theta_size = params_.N_total * sizeof(float);
    VkDeviceSize R_size = params_.N_total * sizeof(float);
    VkDeviceSize psi_size = params_.N_total * 4 * 2 * sizeof(float);  // 4 components × 2 (complex) × float
    VkDeviceSize density_size = params_.N_total * sizeof(float);

    // Helper lambda to create buffer
    auto createBuffer = [this](VkDeviceSize size, VkBufferUsageFlags usage,
                                VkBuffer& buffer, VkDeviceMemory& memory) {
        VkBufferCreateInfo bufferInfo{};
        bufferInfo.sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;
        bufferInfo.size = size;
        bufferInfo.usage = usage;
        bufferInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;

        if (vkCreateBuffer(device_, &bufferInfo, nullptr, &buffer) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create buffer");
        }

        VkMemoryRequirements memRequirements;
        vkGetBufferMemoryRequirements(device_, buffer, &memRequirements);

        VkMemoryAllocateInfo allocInfo{};
        allocInfo.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
        allocInfo.allocationSize = memRequirements.size;

        // Find memory type (host visible for CPU access)
        VkPhysicalDeviceMemoryProperties memProperties;
        vkGetPhysicalDeviceMemoryProperties(physical_device_, &memProperties);

        uint32_t memoryTypeIndex = UINT32_MAX;
        for (uint32_t i = 0; i < memProperties.memoryTypeCount; i++) {
            if ((memRequirements.memoryTypeBits & (1 << i)) &&
                (memProperties.memoryTypes[i].propertyFlags &
                 (VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT))) {
                memoryTypeIndex = i;
                break;
            }
        }

        if (memoryTypeIndex == UINT32_MAX) {
            throw std::runtime_error("Failed to find suitable memory type");
        }

        allocInfo.memoryTypeIndex = memoryTypeIndex;

        if (vkAllocateMemory(device_, &allocInfo, nullptr, &memory) != VK_SUCCESS) {
            throw std::runtime_error("Failed to allocate buffer memory");
        }

        vkBindBufferMemory(device_, buffer, memory, 0);
    };

    // Create all buffers
    VkBufferUsageFlags storageUsage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT;

    createBuffer(theta_size, storageUsage, buffers_.theta, buffers_.theta_memory);
    createBuffer(theta_size, storageUsage, buffers_.theta_out, buffers_.theta_out_memory);
    createBuffer(theta_size, storageUsage, buffers_.omega, buffers_.omega_memory);
    createBuffer(R_size, storageUsage, buffers_.R_field, buffers_.R_memory);
    createBuffer(psi_size, storageUsage, buffers_.psi, buffers_.psi_memory);
    createBuffer(density_size, storageUsage, buffers_.spinor_density, buffers_.density_memory);

    buffers_.theta_size = theta_size;
    buffers_.omega_size = theta_size;
    buffers_.R_size = R_size;
    buffers_.psi_size = psi_size;
    buffers_.density_size = density_size;

    std::cout << "[SMFT] Created buffers: theta=" << theta_size << "B, R=" << R_size
              << "B, psi=" << psi_size << "B, ρ=" << density_size << "B" << std::endl;
}

void SMFTPipeline::createShaderModules() {
    kuramoto_shader_ = loadShaderModule("shaders/SMFT/kuramoto_step.comp.spv");
    sync_shader_ = loadShaderModule("shaders/SMFT/sync_field.comp.spv");
    dirac_shader_ = loadShaderModule("shaders/SMFT/dirac_rk4.comp.spv");
    feedback_shader_ = loadShaderModule("shaders/SMFT/spinor_feedback.comp.spv");

    std::cout << "[SMFT] Loaded shader modules" << std::endl;
}

VkShaderModule SMFTPipeline::loadShaderModule(const char* filepath) {
    std::ifstream file(filepath, std::ios::ate | std::ios::binary);

    if (!file.is_open()) {
        throw std::runtime_error(std::string("Failed to open shader file: ") + filepath);
    }

    size_t fileSize = (size_t)file.tellg();
    std::vector<char> buffer(fileSize);

    file.seekg(0);
    file.read(buffer.data(), fileSize);
    file.close();

    VkShaderModuleCreateInfo createInfo{};
    createInfo.sType = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO;
    createInfo.codeSize = buffer.size();
    createInfo.pCode = reinterpret_cast<const uint32_t*>(buffer.data());

    VkShaderModule shaderModule;
    if (vkCreateShaderModule(device_, &createInfo, nullptr, &shaderModule) != VK_SUCCESS) {
        throw std::runtime_error("Failed to create shader module");
    }

    return shaderModule;
}

void SMFTPipeline::createDescriptorSetLayout() {
    // Descriptor set layout for SMFT pipeline
    // All shaders share same layout for consistency
    std::vector<VkDescriptorSetLayoutBinding> bindings = {
        // Binding 0: theta / psi (read/write depending on shader)
        {0, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
        // Binding 1: theta_out / R_field / spinor_density (write)
        {1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
        // Binding 2: omega / R_field (read)
        {2, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
        // Binding 3: spinor_density (read for kuramoto feedback)
        {3, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
    };

    VkDescriptorSetLayoutCreateInfo layoutInfo{};
    layoutInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
    layoutInfo.bindingCount = static_cast<uint32_t>(bindings.size());
    layoutInfo.pBindings = bindings.data();

    if (vkCreateDescriptorSetLayout(device_, &layoutInfo, nullptr, &descriptor_set_layout_) != VK_SUCCESS) {
        throw std::runtime_error("Failed to create descriptor set layout");
    }
}

void SMFTPipeline::createDescriptorPool() {
    std::vector<VkDescriptorPoolSize> poolSizes = {
        {VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 16}  // 4 bindings × 4 shaders
    };

    VkDescriptorPoolCreateInfo poolInfo{};
    poolInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
    poolInfo.poolSizeCount = static_cast<uint32_t>(poolSizes.size());
    poolInfo.pPoolSizes = poolSizes.data();
    poolInfo.maxSets = 4;  // One set per shader

    if (vkCreateDescriptorPool(device_, &poolInfo, nullptr, &descriptor_pool_) != VK_SUCCESS) {
        throw std::runtime_error("Failed to create descriptor pool");
    }
}

void SMFTPipeline::createDescriptorSet() {
    VkDescriptorSetAllocateInfo allocInfo{};
    allocInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO;
    allocInfo.descriptorPool = descriptor_pool_;
    allocInfo.descriptorSetCount = 1;
    allocInfo.pSetLayouts = &descriptor_set_layout_;

    if (vkAllocateDescriptorSets(device_, &allocInfo, &descriptor_set_) != VK_SUCCESS) {
        throw std::runtime_error("Failed to allocate descriptor set");
    }

    updateDescriptorSet();
}

void SMFTPipeline::updateDescriptorSet() {
    // Update descriptor set with buffer bindings
    std::vector<VkWriteDescriptorSet> descriptorWrites(4);

    VkDescriptorBufferInfo thetaInfo{};
    thetaInfo.buffer = buffers_.theta;
    thetaInfo.offset = 0;
    thetaInfo.range = buffers_.theta_size;

    descriptorWrites[0].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
    descriptorWrites[0].dstSet = descriptor_set_;
    descriptorWrites[0].dstBinding = 0;
    descriptorWrites[0].dstArrayElement = 0;
    descriptorWrites[0].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
    descriptorWrites[0].descriptorCount = 1;
    descriptorWrites[0].pBufferInfo = &thetaInfo;

    VkDescriptorBufferInfo thetaOutInfo{};
    thetaOutInfo.buffer = buffers_.theta_out;
    thetaOutInfo.offset = 0;
    thetaOutInfo.range = buffers_.theta_size;

    descriptorWrites[1].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
    descriptorWrites[1].dstSet = descriptor_set_;
    descriptorWrites[1].dstBinding = 1;
    descriptorWrites[1].dstArrayElement = 0;
    descriptorWrites[1].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
    descriptorWrites[1].descriptorCount = 1;
    descriptorWrites[1].pBufferInfo = &thetaOutInfo;

    VkDescriptorBufferInfo omegaInfo{};
    omegaInfo.buffer = buffers_.omega;
    omegaInfo.offset = 0;
    omegaInfo.range = buffers_.omega_size;

    descriptorWrites[2].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
    descriptorWrites[2].dstSet = descriptor_set_;
    descriptorWrites[2].dstBinding = 2;
    descriptorWrites[2].dstArrayElement = 0;
    descriptorWrites[2].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
    descriptorWrites[2].descriptorCount = 1;
    descriptorWrites[2].pBufferInfo = &omegaInfo;

    VkDescriptorBufferInfo densityInfo{};
    densityInfo.buffer = buffers_.spinor_density;
    densityInfo.offset = 0;
    densityInfo.range = buffers_.density_size;

    descriptorWrites[3].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
    descriptorWrites[3].dstSet = descriptor_set_;
    descriptorWrites[3].dstBinding = 3;
    descriptorWrites[3].dstArrayElement = 0;
    descriptorWrites[3].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
    descriptorWrites[3].descriptorCount = 1;
    descriptorWrites[3].pBufferInfo = &densityInfo;

    vkUpdateDescriptorSets(device_, static_cast<uint32_t>(descriptorWrites.size()),
                           descriptorWrites.data(), 0, nullptr);
}

void SMFTPipeline::createPipelines() {
    // Create pipeline layout with push constants
    VkPushConstantRange pushConstantRange{};
    pushConstantRange.stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
    pushConstantRange.offset = 0;
    pushConstantRange.size = sizeof(SMFTParams);

    VkPipelineLayoutCreateInfo pipelineLayoutInfo{};
    pipelineLayoutInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
    pipelineLayoutInfo.setLayoutCount = 1;
    pipelineLayoutInfo.pSetLayouts = &descriptor_set_layout_;
    pipelineLayoutInfo.pushConstantRangeCount = 1;
    pipelineLayoutInfo.pPushConstantRanges = &pushConstantRange;

    if (vkCreatePipelineLayout(device_, &pipelineLayoutInfo, nullptr, &pipeline_layout_) != VK_SUCCESS) {
        throw std::runtime_error("Failed to create pipeline layout");
    }

    // Helper to create compute pipeline
    auto createComputePipeline = [this](VkShaderModule shader, VkPipeline& pipeline) {
        VkPipelineShaderStageCreateInfo shaderStageInfo{};
        shaderStageInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
        shaderStageInfo.stage = VK_SHADER_STAGE_COMPUTE_BIT;
        shaderStageInfo.module = shader;
        shaderStageInfo.pName = "main";

        VkComputePipelineCreateInfo pipelineInfo{};
        pipelineInfo.sType = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO;
        pipelineInfo.stage = shaderStageInfo;
        pipelineInfo.layout = pipeline_layout_;

        if (vkCreateComputePipelines(device_, VK_NULL_HANDLE, 1, &pipelineInfo, nullptr, &pipeline) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create compute pipeline");
        }
    };

    createComputePipeline(kuramoto_shader_, kuramoto_pipeline_);
    createComputePipeline(sync_shader_, sync_pipeline_);
    createComputePipeline(dirac_shader_, dirac_pipeline_);
    createComputePipeline(feedback_shader_, feedback_pipeline_);

    std::cout << "[SMFT] Created compute pipelines" << std::endl;
}

void SMFTPipeline::createCommandBuffer() {
    VkCommandBufferAllocateInfo allocInfo{};
    allocInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;
    allocInfo.commandPool = command_pool_;
    allocInfo.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
    allocInfo.commandBufferCount = 1;

    if (vkAllocateCommandBuffers(device_, &allocInfo, &command_buffer_) != VK_SUCCESS) {
        throw std::runtime_error("Failed to allocate command buffer");
    }
}

void SMFTPipeline::setInitialPhases(const std::vector<float>& theta_data) {
    if (theta_data.size() != params_.N_total) {
        throw std::runtime_error("Phase data size mismatch");
    }
    copyToBuffer(buffers_.theta, theta_data.data(), buffers_.theta_size);
}

void SMFTPipeline::setNaturalFrequencies(const std::vector<float>& omega_data) {
    if (omega_data.size() != params_.N_total) {
        throw std::runtime_error("Omega data size mismatch");
    }
    copyToBuffer(buffers_.omega, omega_data.data(), buffers_.omega_size);
}

void SMFTPipeline::setInitialSpinorField(const std::vector<float>& psi_data) {
    if (psi_data.size() != params_.N_total * 8) {  // 4 components × 2 (complex)
        throw std::runtime_error("Spinor data size mismatch");
    }
    copyToBuffer(buffers_.psi, psi_data.data(), buffers_.psi_size);
}

void SMFTPipeline::step(float dt, float coupling, float damping) {
    params_.dt = dt;
    params_.K = coupling;
    params_.damping = damping;

    recordCommandBuffer(dt, coupling, damping);

    // Submit command buffer
    VkSubmitInfo submitInfo{};
    submitInfo.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;
    submitInfo.commandBufferCount = 1;
    submitInfo.pCommandBuffers = &command_buffer_;

    vkQueueSubmit(compute_queue_, 1, &submitInfo, VK_NULL_HANDLE);
    vkQueueWaitIdle(compute_queue_);

    current_time_ += dt;
}

void SMFTPipeline::recordCommandBuffer(float dt, float coupling, float damping) {
    vkResetCommandBuffer(command_buffer_, 0);

    VkCommandBufferBeginInfo beginInfo{};
    beginInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    beginInfo.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;

    vkBeginCommandBuffer(command_buffer_, &beginInfo);

    // Use SMFTMemoryBarriers to record full pipeline with barriers
    SMFTMemoryBarriers::recordSMFTPipeline(
        command_buffer_,
        kuramoto_pipeline_,
        sync_pipeline_,
        dirac_pipeline_,
        feedback_pipeline_,
        pipeline_layout_,
        descriptor_set_,
        buffers_,
        params_,
        params_.Nx,
        params_.Ny
    );

    vkEndCommandBuffer(command_buffer_);
}

void SMFTPipeline::evolve(int num_steps, float dt, float coupling, float damping) {
    for (int i = 0; i < num_steps; i++) {
        step(dt, coupling, damping);
    }
}

std::vector<float> SMFTPipeline::getPhases() const {
    std::vector<float> theta_data(params_.N_total);
    copyFromBuffer(buffers_.theta, theta_data.data(), buffers_.theta_size);
    return theta_data;
}

std::vector<float> SMFTPipeline::getSyncField() const {
    std::vector<float> R_data(params_.N_total);
    copyFromBuffer(buffers_.R_field, R_data.data(), buffers_.R_size);
    return R_data;
}

std::vector<float> SMFTPipeline::getSpinorField() const {
    std::vector<float> psi_data(params_.N_total * 8);
    copyFromBuffer(buffers_.psi, psi_data.data(), buffers_.psi_size);
    return psi_data;
}

std::vector<float> SMFTPipeline::getSpinorDensity() const {
    std::vector<float> density_data(params_.N_total);
    copyFromBuffer(buffers_.spinor_density, density_data.data(), buffers_.density_size);
    return density_data;
}

void SMFTPipeline::setParams(const SMFTParams& params) {
    params_ = params;
}

void SMFTPipeline::reset() {
    current_time_ = 0.0f;
}

void SMFTPipeline::copyToBuffer(VkBuffer buffer, const void* data, VkDeviceSize size) {
    void* mapped;
    VkDeviceMemory memory;

    // Find corresponding memory for buffer
    if (buffer == buffers_.theta) memory = buffers_.theta_memory;
    else if (buffer == buffers_.theta_out) memory = buffers_.theta_out_memory;
    else if (buffer == buffers_.omega) memory = buffers_.omega_memory;
    else if (buffer == buffers_.R_field) memory = buffers_.R_memory;
    else if (buffer == buffers_.psi) memory = buffers_.psi_memory;
    else if (buffer == buffers_.spinor_density) memory = buffers_.density_memory;
    else throw std::runtime_error("Unknown buffer");

    vkMapMemory(device_, memory, 0, size, 0, &mapped);
    memcpy(mapped, data, size);
    vkUnmapMemory(device_, memory);
}

void SMFTPipeline::copyFromBuffer(VkBuffer buffer, void* data, VkDeviceSize size) const {
    void* mapped;
    VkDeviceMemory memory;

    if (buffer == buffers_.theta) memory = buffers_.theta_memory;
    else if (buffer == buffers_.theta_out) memory = buffers_.theta_out_memory;
    else if (buffer == buffers_.omega) memory = buffers_.omega_memory;
    else if (buffer == buffers_.R_field) memory = buffers_.R_memory;
    else if (buffer == buffers_.psi) memory = buffers_.psi_memory;
    else if (buffer == buffers_.spinor_density) memory = buffers_.density_memory;
    else throw std::runtime_error("Unknown buffer");

    vkMapMemory(device_, memory, 0, size, 0, &mapped);
    memcpy(data, mapped, size);
    vkUnmapMemory(device_, memory);
}

void SMFTPipeline::destroyBuffers() {
    auto destroyBuffer = [this](VkBuffer& buffer, VkDeviceMemory& memory) {
        if (buffer != VK_NULL_HANDLE) {
            vkDestroyBuffer(device_, buffer, nullptr);
            buffer = VK_NULL_HANDLE;
        }
        if (memory != VK_NULL_HANDLE) {
            vkFreeMemory(device_, memory, nullptr);
            memory = VK_NULL_HANDLE;
        }
    };

    destroyBuffer(buffers_.theta, buffers_.theta_memory);
    destroyBuffer(buffers_.theta_out, buffers_.theta_out_memory);
    destroyBuffer(buffers_.omega, buffers_.omega_memory);
    destroyBuffer(buffers_.R_field, buffers_.R_memory);
    destroyBuffer(buffers_.psi, buffers_.psi_memory);
    destroyBuffer(buffers_.spinor_density, buffers_.density_memory);
}

void SMFTPipeline::destroyPipelines() {
    if (kuramoto_pipeline_ != VK_NULL_HANDLE) {
        vkDestroyPipeline(device_, kuramoto_pipeline_, nullptr);
        kuramoto_pipeline_ = VK_NULL_HANDLE;
    }
    if (sync_pipeline_ != VK_NULL_HANDLE) {
        vkDestroyPipeline(device_, sync_pipeline_, nullptr);
        sync_pipeline_ = VK_NULL_HANDLE;
    }
    if (dirac_pipeline_ != VK_NULL_HANDLE) {
        vkDestroyPipeline(device_, dirac_pipeline_, nullptr);
        dirac_pipeline_ = VK_NULL_HANDLE;
    }
    if (feedback_pipeline_ != VK_NULL_HANDLE) {
        vkDestroyPipeline(device_, feedback_pipeline_, nullptr);
        feedback_pipeline_ = VK_NULL_HANDLE;
    }
    if (pipeline_layout_ != VK_NULL_HANDLE) {
        vkDestroyPipelineLayout(device_, pipeline_layout_, nullptr);
        pipeline_layout_ = VK_NULL_HANDLE;
    }
    if (kuramoto_shader_ != VK_NULL_HANDLE) {
        vkDestroyShaderModule(device_, kuramoto_shader_, nullptr);
        kuramoto_shader_ = VK_NULL_HANDLE;
    }
    if (sync_shader_ != VK_NULL_HANDLE) {
        vkDestroyShaderModule(device_, sync_shader_, nullptr);
        sync_shader_ = VK_NULL_HANDLE;
    }
    if (dirac_shader_ != VK_NULL_HANDLE) {
        vkDestroyShaderModule(device_, dirac_shader_, nullptr);
        dirac_shader_ = VK_NULL_HANDLE;
    }
    if (feedback_shader_ != VK_NULL_HANDLE) {
        vkDestroyShaderModule(device_, feedback_shader_, nullptr);
        feedback_shader_ = VK_NULL_HANDLE;
    }
}

} // namespace SMFT
} // namespace Nova