#include "observables.hpp"
#include <stdexcept>
#include <iostream>

namespace SMFT {

ObservablesEngine::ObservablesEngine(
    VkDevice device,
    VkPhysicalDevice physical_device,
    VkCommandPool command_pool,
    VkQueue compute_queue
) : device_(device),
    physical_device_(physical_device),
    command_pool_(command_pool),
    compute_queue_(compute_queue) {
}

ObservablesEngine::~ObservablesEngine() {
    destroyPipelines();
    destroyBuffers();
}

void ObservablesEngine::initialize(uint32_t grid_x, uint32_t grid_y, float dx) {
    grid_x_ = grid_x;
    grid_y_ = grid_y;
    N_total_ = grid_x * grid_y;
    dx_ = dx;

    createShaderModules();
    createDescriptorSetLayout();
    createDescriptorPool();
    createDescriptorSet();
    createPipelines();
    createCommandBuffer();
    createReductionBuffer();
}

AllObservables ObservablesEngine::compute(
    const SMFTBuffers& buffers,
    float Delta,
    float time,
    uint32_t timestep
) {
    AllObservables obs;
    obs.simulation_time = time;
    obs.timestep = timestep;

    // ...
    
    return obs;
}
void ObservablesEngine::createShaderModules(){}
void ObservablesEngine::createDescriptorSetLayout(){}
void ObservablesEngine::createDescriptorPool(){}
void ObservablesEngine::createDescriptorSet(){}
void ObservablesEngine::createPipelines(){}
void ObservablesEngine::createCommandBuffer(){}
void ObservablesEngine::createReductionBuffer(){}
void ObservablesEngine::destroyPipelines(){}
void ObservablesEngine::destroyBuffers(){}
VkShaderModule ObservablesEngine::loadShaderModule(const char* filepath) { return VK_NULL_HANDLE; }

} // namespace SMFT
