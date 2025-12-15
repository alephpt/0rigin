#pragma once
#include <vulkan/vulkan.h>

namespace Nova {
namespace MSFT {

struct MSFTBuffers {
    // Buffers
    VkBuffer theta;
    VkBuffer theta_out;
    VkBuffer omega;
    VkBuffer R_field;
    VkBuffer psi;
    VkBuffer spinor_density;

    // Memory
    VkDeviceMemory theta_memory;
    VkDeviceMemory theta_out_memory;
    VkDeviceMemory omega_memory;
    VkDeviceMemory R_memory;
    VkDeviceMemory psi_memory;
    VkDeviceMemory density_memory;

    // Sizes
    VkDeviceSize theta_size;
    VkDeviceSize omega_size;
    VkDeviceSize R_size;
    VkDeviceSize psi_size;
    VkDeviceSize density_size;
};

} // namespace MSFT
} // namespace Nova