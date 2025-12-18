/**
 * Headless SMFT Physics Test
 * Tests SMFT compute pipeline without graphics/window
 * Validates physics simulation independently
 */

#include <vulkan/vulkan.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <random>
#include <filesystem>
#include "../src/SMFTEngine.h"

// Helper function to create minimal Vulkan instance
VkInstance createInstance() {
    VkApplicationInfo appInfo{};
    appInfo.sType = VK_STRUCTURE_TYPE_APPLICATION_INFO;
    appInfo.pApplicationName = "SMFT Headless Test";
    appInfo.applicationVersion = VK_MAKE_VERSION(1, 0, 0);
    appInfo.pEngineName = "No Engine";
    appInfo.engineVersion = VK_MAKE_VERSION(1, 0, 0);
    appInfo.apiVersion = VK_API_VERSION_1_3;

    VkInstanceCreateInfo createInfo{};
    createInfo.sType = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO;
    createInfo.pApplicationInfo = &appInfo;

    // Enable validation layers in debug
    const std::vector<const char*> validationLayers = {
        "VK_LAYER_KHRONOS_validation"
    };

    #ifdef DEBUG
    createInfo.enabledLayerCount = static_cast<uint32_t>(validationLayers.size());
    createInfo.ppEnabledLayerNames = validationLayers.data();
    #else
    createInfo.enabledLayerCount = 0;
    #endif

    VkInstance instance;
    VkResult result = vkCreateInstance(&createInfo, nullptr, &instance);
    if (result != VK_SUCCESS) {
        std::cerr << "Failed to create Vulkan instance: " << result << std::endl;
        return VK_NULL_HANDLE;
    }

    return instance;
}

// Helper function to select physical device
VkPhysicalDevice selectPhysicalDevice(VkInstance instance) {
    uint32_t deviceCount = 0;
    vkEnumeratePhysicalDevices(instance, &deviceCount, nullptr);

    if (deviceCount == 0) {
        std::cerr << "No Vulkan devices found!" << std::endl;
        return VK_NULL_HANDLE;
    }

    std::vector<VkPhysicalDevice> devices(deviceCount);
    vkEnumeratePhysicalDevices(instance, &deviceCount, devices.data());

    // Select first device with compute capability
    for (const auto& device : devices) {
        VkPhysicalDeviceProperties deviceProperties;
        vkGetPhysicalDeviceProperties(device, &deviceProperties);

        std::cout << "GPU: " << deviceProperties.deviceName << std::endl;
        std::cout << "  Type: ";
        switch(deviceProperties.deviceType) {
            case VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU:
                std::cout << "Discrete GPU" << std::endl;
                break;
            case VK_PHYSICAL_DEVICE_TYPE_INTEGRATED_GPU:
                std::cout << "Integrated GPU" << std::endl;
                break;
            default:
                std::cout << "Other" << std::endl;
        }

        // Check for compute queue
        uint32_t queueFamilyCount = 0;
        vkGetPhysicalDeviceQueueFamilyProperties(device, &queueFamilyCount, nullptr);

        std::vector<VkQueueFamilyProperties> queueFamilies(queueFamilyCount);
        vkGetPhysicalDeviceQueueFamilyProperties(device, &queueFamilyCount, queueFamilies.data());

        for (uint32_t i = 0; i < queueFamilyCount; i++) {
            if (queueFamilies[i].queueFlags & VK_QUEUE_COMPUTE_BIT) {
                std::cout << "  Compute queue family found at index " << i << std::endl;
                return device;
            }
        }
    }

    std::cerr << "No suitable GPU found!" << std::endl;
    return VK_NULL_HANDLE;
}

// Helper to find compute queue family index
uint32_t findComputeQueueFamily(VkPhysicalDevice device) {
    uint32_t queueFamilyCount = 0;
    vkGetPhysicalDeviceQueueFamilyProperties(device, &queueFamilyCount, nullptr);

    std::vector<VkQueueFamilyProperties> queueFamilies(queueFamilyCount);
    vkGetPhysicalDeviceQueueFamilyProperties(device, &queueFamilyCount, queueFamilies.data());

    for (uint32_t i = 0; i < queueFamilyCount; i++) {
        if (queueFamilies[i].queueFlags & VK_QUEUE_COMPUTE_BIT) {
            return i;
        }
    }

    return UINT32_MAX;
}

// Helper function to create logical device
VkDevice createLogicalDevice(VkPhysicalDevice physicalDevice, uint32_t queueFamilyIndex, VkQueue& computeQueue) {
    float queuePriority = 1.0f;

    VkDeviceQueueCreateInfo queueCreateInfo{};
    queueCreateInfo.sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO;
    queueCreateInfo.queueFamilyIndex = queueFamilyIndex;
    queueCreateInfo.queueCount = 1;
    queueCreateInfo.pQueuePriorities = &queuePriority;

    VkPhysicalDeviceFeatures deviceFeatures{};

    VkDeviceCreateInfo createInfo{};
    createInfo.sType = VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO;
    createInfo.queueCreateInfoCount = 1;
    createInfo.pQueueCreateInfos = &queueCreateInfo;
    createInfo.pEnabledFeatures = &deviceFeatures;
    createInfo.enabledExtensionCount = 0;

    VkDevice device;
    VkResult result = vkCreateDevice(physicalDevice, &createInfo, nullptr, &device);
    if (result != VK_SUCCESS) {
        std::cerr << "Failed to create logical device: " << result << std::endl;
        return VK_NULL_HANDLE;
    }

    vkGetDeviceQueue(device, queueFamilyIndex, 0, &computeQueue);

    return device;
}

// Minimal Nova structure for headless operation
struct NovaHeadless {
    struct Architect {
        VkDevice logical_device;
        VkPhysicalDevice physical_device;
        struct Queues {
            VkQueue compute;
            VkQueue graphics; // Same as compute for headless
            struct Indices {
                std::optional<uint32_t> compute_family;
                std::optional<uint32_t> graphics_family;
            } indices;
        } queues;
    };

    Architect* _architect;

    NovaHeadless(VkDevice device, VkPhysicalDevice physicalDevice, VkQueue queue, uint32_t queueFamily) {
        _architect = new Architect();
        _architect->logical_device = device;
        _architect->physical_device = physicalDevice;
        _architect->queues.compute = queue;
        _architect->queues.graphics = queue; // Use same queue
        _architect->queues.indices.compute_family = queueFamily;
        _architect->queues.indices.graphics_family = queueFamily;
    }

    ~NovaHeadless() {
        delete _architect;
    }
};

// Write field data to file
void writeFieldToFile(const std::string& filename, const std::vector<float>& data, uint32_t Nx, uint32_t Ny) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    file << std::fixed << std::setprecision(6);
    for (uint32_t y = 0; y < Ny; y++) {
        for (uint32_t x = 0; x < Nx; x++) {
            file << data[y * Nx + x];
            if (x < Nx - 1) file << " ";
        }
        file << "\n";
    }

    file.close();
    std::cout << "Wrote " << filename << std::endl;
}

int main() {
    std::cout << "=== SMFT Headless Physics Test ===" << std::endl;
    std::cout << "Testing compute-only physics simulation without graphics" << std::endl << std::endl;

    // 1. Initialize Vulkan
    VkInstance instance = createInstance();
    if (instance == VK_NULL_HANDLE) {
        return -1;
    }

    VkPhysicalDevice physicalDevice = selectPhysicalDevice(instance);
    if (physicalDevice == VK_NULL_HANDLE) {
        vkDestroyInstance(instance, nullptr);
        return -1;
    }

    uint32_t queueFamilyIndex = findComputeQueueFamily(physicalDevice);
    if (queueFamilyIndex == UINT32_MAX) {
        std::cerr << "No compute queue family found!" << std::endl;
        vkDestroyInstance(instance, nullptr);
        return -1;
    }

    VkQueue computeQueue;
    VkDevice device = createLogicalDevice(physicalDevice, queueFamilyIndex, computeQueue);
    if (device == VK_NULL_HANDLE) {
        vkDestroyInstance(instance, nullptr);
        return -1;
    }

    std::cout << "\nVulkan initialized successfully" << std::endl;
    std::cout << "================================" << std::endl << std::endl;

    // 2. Create headless Nova wrapper
    NovaHeadless* nova = new NovaHeadless(device, physicalDevice, computeQueue, queueFamilyIndex);

    // 3. Initialize SMFTEngine with test parameters
    SMFTEngine* smftEngine = new SMFTEngine(reinterpret_cast<Nova*>(nova));

    const uint32_t Nx = 128;
    const uint32_t Ny = 128;
    const float Delta = 2.5f;     // Vacuum potential from SMFT.cpp
    const float chiral_angle = 0.0f;

    std::cout << "Initializing SMFT Engine:" << std::endl;
    std::cout << "  Grid size: " << Nx << " x " << Ny << std::endl;
    std::cout << "  Delta (vacuum potential): " << Delta << std::endl;
    std::cout << "  Chiral angle: " << chiral_angle << std::endl << std::endl;

    smftEngine->initialize(Nx, Ny, Delta, chiral_angle);

    // 4. Set initial conditions (from SMFT.cpp)
    std::vector<float> theta_init(Nx * Ny);
    std::vector<float> omega_init(Nx * Ny);

    std::mt19937 rng(42); // Deterministic seed
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);

    for (uint32_t y = 0; y < Ny; y++) {
        for (uint32_t x = 0; x < Nx; x++) {
            uint32_t idx = y * Nx + x;
            theta_init[idx] = 2.0f * M_PI * dist(rng);
            omega_init[idx] = 0.5f * (dist(rng) - 0.5f);

            // Add spatial wave modulation
            float wave_x = std::cos(2.0f * M_PI * x / Nx * 2.0f);
            float wave_y = std::cos(2.0f * M_PI * y / Ny * 2.0f);
            theta_init[idx] += 0.2f * wave_x * wave_y;
        }
    }

    smftEngine->setInitialPhases(theta_init);
    smftEngine->setNaturalFrequencies(omega_init);

    // 5. Run simulation
    const int numSteps = 1000;
    const int outputInterval = 100;
    const float dt = 0.01f;  // From SMFT.cpp line 89
    const float K = 1.0f;    // From SMFT.cpp line 90 (user adjustable)
    const float damping = 0.1f;

    std::cout << "Running simulation:" << std::endl;
    std::cout << "  Steps: " << numSteps << std::endl;
    std::cout << "  dt: " << dt << std::endl;
    std::cout << "  Coupling strength K: " << K << std::endl;
    std::cout << "  Damping: " << damping << std::endl << std::endl;

    auto startTime = std::chrono::high_resolution_clock::now();

    for (int step = 0; step < numSteps; step++) {
        auto stepStart = std::chrono::high_resolution_clock::now();

        // Run physics step with timeout protection
        try {
            smftEngine->step(dt, K, damping);
        } catch (const std::exception& e) {
            std::cerr << "Error at step " << step << ": " << e.what() << std::endl;
            break;
        }

        // Output every 100 steps
        if (step % outputInterval == 0 || step == numSteps - 1) {
            std::vector<float> R_field = smftEngine->getSyncField();

            // Calculate statistics
            float R_sum = 0.0f;
            float R_min = 1.0f;
            float R_max = 0.0f;

            for (float r : R_field) {
                R_sum += r;
                R_min = std::min(R_min, r);
                R_max = std::max(R_max, r);
            }

            float R_avg = R_sum / R_field.size();

            auto stepEnd = std::chrono::high_resolution_clock::now();
            auto stepDuration = std::chrono::duration_cast<std::chrono::milliseconds>(stepEnd - stepStart);

            std::cout << "Step " << std::setw(4) << step
                      << ": R_avg = " << std::fixed << std::setprecision(4) << R_avg
                      << ", R_min = " << R_min
                      << ", R_max = " << R_max
                      << " (step time: " << stepDuration.count() << " ms)" << std::endl;
        }
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    auto totalDuration = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime);

    std::cout << "\nSimulation completed in " << totalDuration.count() << " seconds" << std::endl;

    // 6. Write final field data to files
    std::cout << "\nWriting output files..." << std::endl;

    // Create output directory
    std::filesystem::create_directories("output");

    // Get final field data
    std::vector<float> R_field = smftEngine->getSyncField();
    std::vector<float> theta_field = smftEngine->getPhaseField();
    std::vector<float> gravity_field = smftEngine->getGravitationalField();

    // Separate gravity field into x and y components
    std::vector<float> gravity_x(Nx * Ny);
    std::vector<float> gravity_y(Nx * Ny);
    for (uint32_t i = 0; i < Nx * Ny; i++) {
        gravity_x[i] = gravity_field[2*i];
        gravity_y[i] = gravity_field[2*i + 1];
    }

    // Write field data
    writeFieldToFile("output/R_field.dat", R_field, Nx, Ny);
    writeFieldToFile("output/theta.dat", theta_field, Nx, Ny);
    writeFieldToFile("output/gravity_x.dat", gravity_x, Nx, Ny);
    writeFieldToFile("output/gravity_y.dat", gravity_y, Nx, Ny);

    // 7. Cleanup
    delete smftEngine;
    delete nova;

    vkDestroyDevice(device, nullptr);
    vkDestroyInstance(instance, nullptr);

    std::cout << "\nTest completed successfully!" << std::endl;

    return 0;
}