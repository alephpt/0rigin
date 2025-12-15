/**
 * @file queue_diagnostics.hpp
 * @brief GPU queue diagnostics for MSFT pipeline
 */

#pragma once

#include <vulkan/vulkan.h>
#include <iostream>

namespace Nova {
namespace Core {

/**
 * @brief Diagnostics utility for Vulkan queue capabilities
 */
class QueueDiagnostics {
public:
    QueueDiagnostics(VkPhysicalDevice physicalDevice)
        : physicalDevice_(physicalDevice) {}

    /**
     * @brief Print GPU compute capabilities
     */
    void printCapabilities() {
        VkPhysicalDeviceProperties props;
        vkGetPhysicalDeviceProperties(physicalDevice_, &props);

        VkPhysicalDeviceFeatures features;
        vkGetPhysicalDeviceFeatures(physicalDevice_, &features);

        std::cout << "[Queue] Compute capabilities:" << std::endl;
        std::cout << "  Max workgroup size: "
                  << props.limits.maxComputeWorkGroupSize[0] << " x "
                  << props.limits.maxComputeWorkGroupSize[1] << " x "
                  << props.limits.maxComputeWorkGroupSize[2] << std::endl;
        std::cout << "  Max workgroup invocations: "
                  << props.limits.maxComputeWorkGroupInvocations << std::endl;
        std::cout << "  FP64 support: "
                  << (features.shaderFloat64 ? "Yes" : "No") << std::endl;
    }

private:
    VkPhysicalDevice physicalDevice_;
};

} // namespace Core
} // namespace Nova