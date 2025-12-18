/**
 * Minimal MSFT Compute Test - Direct Vulkan without Nova
 * Tests physics compute without any graphics dependencies
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
#include <cassert>

struct ComputeContext {
    VkInstance instance;
    VkPhysicalDevice physicalDevice;
    VkDevice device;
    VkQueue computeQueue;
    uint32_t queueFamilyIndex;
};

// Initialize Vulkan for compute only
bool initVulkan(ComputeContext& ctx) {
    // Create instance
    VkApplicationInfo appInfo{};
    appInfo.sType = VK_STRUCTURE_TYPE_APPLICATION_INFO;
    appInfo.pApplicationName = "MSFT Compute Test";
    appInfo.applicationVersion = VK_MAKE_VERSION(1, 0, 0);
    appInfo.apiVersion = VK_API_VERSION_1_2;

    VkInstanceCreateInfo instanceInfo{};
    instanceInfo.sType = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO;
    instanceInfo.pApplicationInfo = &appInfo;

    if (vkCreateInstance(&instanceInfo, nullptr, &ctx.instance) != VK_SUCCESS) {
        std::cerr << "Failed to create Vulkan instance" << std::endl;
        return false;
    }

    // Select physical device
    uint32_t deviceCount = 0;
    vkEnumeratePhysicalDevices(ctx.instance, &deviceCount, nullptr);
    if (deviceCount == 0) {
        std::cerr << "No Vulkan devices found" << std::endl;
        return false;
    }

    std::vector<VkPhysicalDevice> devices(deviceCount);
    vkEnumeratePhysicalDevices(ctx.instance, &deviceCount, devices.data());

    // Select first device with compute
    for (auto device : devices) {
        VkPhysicalDeviceProperties props;
        vkGetPhysicalDeviceProperties(device, &props);

        std::cout << "GPU: " << props.deviceName << std::endl;

        // Find compute queue
        uint32_t queueFamilyCount = 0;
        vkGetPhysicalDeviceQueueFamilyProperties(device, &queueFamilyCount, nullptr);

        std::vector<VkQueueFamilyProperties> queueFamilies(queueFamilyCount);
        vkGetPhysicalDeviceQueueFamilyProperties(device, &queueFamilyCount, queueFamilies.data());

        for (uint32_t i = 0; i < queueFamilyCount; i++) {
            if (queueFamilies[i].queueFlags & VK_QUEUE_COMPUTE_BIT) {
                ctx.physicalDevice = device;
                ctx.queueFamilyIndex = i;
                break;
            }
        }

        if (ctx.physicalDevice != VK_NULL_HANDLE) break;
    }

    if (ctx.physicalDevice == VK_NULL_HANDLE) {
        std::cerr << "No suitable GPU found" << std::endl;
        return false;
    }

    // Create logical device
    float queuePriority = 1.0f;

    VkDeviceQueueCreateInfo queueInfo{};
    queueInfo.sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO;
    queueInfo.queueFamilyIndex = ctx.queueFamilyIndex;
    queueInfo.queueCount = 1;
    queueInfo.pQueuePriorities = &queuePriority;

    VkDeviceCreateInfo deviceInfo{};
    deviceInfo.sType = VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO;
    deviceInfo.queueCreateInfoCount = 1;
    deviceInfo.pQueueCreateInfos = &queueInfo;

    if (vkCreateDevice(ctx.physicalDevice, &deviceInfo, nullptr, &ctx.device) != VK_SUCCESS) {
        std::cerr << "Failed to create logical device" << std::endl;
        return false;
    }

    vkGetDeviceQueue(ctx.device, ctx.queueFamilyIndex, 0, &ctx.computeQueue);

    std::cout << "Vulkan initialized successfully" << std::endl;
    return true;
}

// Simple CPU implementation of Kuramoto dynamics
void cpuKuramotoStep(std::vector<float>& theta, const std::vector<float>& omega,
                     uint32_t Nx, uint32_t Ny, float K, float dt) {
    std::vector<float> theta_new(Nx * Ny);

    for (uint32_t y = 0; y < Ny; y++) {
        for (uint32_t x = 0; x < Nx; x++) {
            uint32_t idx = y * Nx + x;

            // Calculate coupling with neighbors
            float coupling = 0.0f;

            // 4-neighbor coupling (von Neumann neighborhood)
            uint32_t neighbors[4][2] = {
                {(x + 1) % Nx, y},           // Right
                {(x + Nx - 1) % Nx, y},      // Left
                {x, (y + 1) % Ny},           // Down
                {x, (y + Ny - 1) % Ny}       // Up
            };

            for (int n = 0; n < 4; n++) {
                uint32_t nx = neighbors[n][0];
                uint32_t ny = neighbors[n][1];
                uint32_t nidx = ny * Nx + nx;
                coupling += std::sin(theta[nidx] - theta[idx]);
            }

            // Kuramoto equation: dθ/dt = ω + K * coupling
            theta_new[idx] = theta[idx] + dt * (omega[idx] + K * coupling / 4.0f);
        }
    }

    theta = theta_new;
}

// Calculate synchronization field R
std::vector<float> calculateSyncField(const std::vector<float>& theta,
                                      uint32_t Nx, uint32_t Ny) {
    std::vector<float> R_field(Nx * Ny);

    for (uint32_t y = 0; y < Ny; y++) {
        for (uint32_t x = 0; x < Nx; x++) {
            uint32_t idx = y * Nx + x;

            // Local averaging over neighborhood
            float sum_cos = 0.0f;
            float sum_sin = 0.0f;
            int count = 0;

            // Include self and 4 neighbors
            int radius = 1;
            for (int dy = -radius; dy <= radius; dy++) {
                for (int dx = -radius; dx <= radius; dx++) {
                    // Periodic boundaries
                    int nx = (x + dx + Nx) % Nx;
                    int ny = (y + dy + Ny) % Ny;
                    uint32_t nidx = ny * Nx + nx;

                    sum_cos += std::cos(theta[nidx]);
                    sum_sin += std::sin(theta[nidx]);
                    count++;
                }
            }

            // R = |⟨e^(iθ)⟩|
            float avg_cos = sum_cos / count;
            float avg_sin = sum_sin / count;
            R_field[idx] = std::sqrt(avg_cos * avg_cos + avg_sin * avg_sin);
        }
    }

    return R_field;
}

// Calculate gravitational field from R
void calculateGravityField(const std::vector<float>& R_field,
                          std::vector<float>& gravity_x,
                          std::vector<float>& gravity_y,
                          uint32_t Nx, uint32_t Ny, float Delta) {
    for (uint32_t y = 0; y < Ny; y++) {
        for (uint32_t x = 0; x < Nx; x++) {
            uint32_t idx = y * Nx + x;

            // Central differences for gradient
            uint32_t x_plus = (x + 1) % Nx;
            uint32_t x_minus = (x + Nx - 1) % Nx;
            uint32_t y_plus = (y + 1) % Ny;
            uint32_t y_minus = (y + Ny - 1) % Ny;

            float dR_dx = (R_field[y * Nx + x_plus] - R_field[y * Nx + x_minus]) / 2.0f;
            float dR_dy = (R_field[y_plus * Nx + x] - R_field[y_minus * Nx + x]) / 2.0f;

            // g = -Δ·∇R
            gravity_x[idx] = -Delta * dR_dx;
            gravity_y[idx] = -Delta * dR_dy;
        }
    }
}

// Write field data to file
void writeFieldToFile(const std::string& filename, const std::vector<float>& data,
                     uint32_t Nx, uint32_t Ny) {
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
    std::cout << "=== MSFT Compute-Only Test ===" << std::endl;
    std::cout << "CPU-based physics simulation for validation" << std::endl << std::endl;

    // Initialize Vulkan (just for device query)
    ComputeContext ctx{};
    if (!initVulkan(ctx)) {
        return -1;
    }

    // Parameters
    const uint32_t Nx = 128;
    const uint32_t Ny = 128;
    const float Delta = 2.5f;
    const float dt = 0.01f;
    const float K = 1.0f;
    const float damping = 0.1f;

    std::cout << "\nSimulation parameters:" << std::endl;
    std::cout << "  Grid: " << Nx << " x " << Ny << std::endl;
    std::cout << "  Delta: " << Delta << std::endl;
    std::cout << "  dt: " << dt << std::endl;
    std::cout << "  K: " << K << std::endl << std::endl;

    // Initialize fields
    std::vector<float> theta(Nx * Ny);
    std::vector<float> omega(Nx * Ny);

    std::mt19937 rng(42);
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);

    for (uint32_t y = 0; y < Ny; y++) {
        for (uint32_t x = 0; x < Nx; x++) {
            uint32_t idx = y * Nx + x;
            theta[idx] = 2.0f * M_PI * dist(rng);
            omega[idx] = 0.5f * (dist(rng) - 0.5f);

            // Add spatial modulation
            float wave_x = std::cos(2.0f * M_PI * x / Nx * 2.0f);
            float wave_y = std::cos(2.0f * M_PI * y / Ny * 2.0f);
            theta[idx] += 0.2f * wave_x * wave_y;
        }
    }

    // Run simulation
    const int numSteps = 10000;  // Extended to 10k steps (t_max = 100)
    const int outputInterval = 500;

    // Timeseries storage for Phase 0 analysis
    std::vector<float> timeseries_R_avg(numSteps);
    std::vector<float> timeseries_R_min(numSteps);
    std::vector<float> timeseries_R_max(numSteps);

    auto startTime = std::chrono::high_resolution_clock::now();

    for (int step = 0; step < numSteps; step++) {
        // Update phases
        cpuKuramotoStep(theta, omega, Nx, Ny, K, dt);

        // Calculate sync field and stats at EVERY step (for Phase 0 analysis)
        std::vector<float> R_field = calculateSyncField(theta, Nx, Ny);

        float R_sum = 0.0f, R_min = 1.0f, R_max = 0.0f;
        for (float r : R_field) {
            R_sum += r;
            R_min = std::min(R_min, r);
            R_max = std::max(R_max, r);
        }
        float R_avg = R_sum / R_field.size();

        // Store timeseries
        timeseries_R_avg[step] = R_avg;
        timeseries_R_min[step] = R_min;
        timeseries_R_max[step] = R_max;

        // Print stats (only every 100 steps for readability)
        if (step % outputInterval == 0 || step == numSteps - 1) {
            std::cout << "Step " << std::setw(4) << step
                      << ": R_avg = " << std::fixed << std::setprecision(4) << R_avg
                      << ", R_min = " << R_min
                      << ", R_max = " << R_max << std::endl;
        }

        // Save final results
        if (step == numSteps - 1) {
            system("mkdir -p output");

            // Calculate gravity field
            std::vector<float> gravity_x(Nx * Ny);
            std::vector<float> gravity_y(Nx * Ny);
            calculateGravityField(R_field, gravity_x, gravity_y, Nx, Ny, Delta);

            // Write spatial fields (final state)
            writeFieldToFile("output/R_field.dat", R_field, Nx, Ny);
            writeFieldToFile("output/theta.dat", theta, Nx, Ny);
            writeFieldToFile("output/gravity_x.dat", gravity_x, Nx, Ny);
            writeFieldToFile("output/gravity_y.dat", gravity_y, Nx, Ny);

            // Write timeseries data for Phase 0 analysis
            std::ofstream ts_avg("output/timeseries_R_avg.dat");
            std::ofstream ts_min("output/timeseries_R_min.dat");
            std::ofstream ts_max("output/timeseries_R_max.dat");

            ts_avg << std::fixed << std::setprecision(6);
            ts_min << std::fixed << std::setprecision(6);
            ts_max << std::fixed << std::setprecision(6);

            for (int i = 0; i < numSteps; i++) {
                ts_avg << i << " " << timeseries_R_avg[i] << "\n";
                ts_min << i << " " << timeseries_R_min[i] << "\n";
                ts_max << i << " " << timeseries_R_max[i] << "\n";
            }

            ts_avg.close();
            ts_min.close();
            ts_max.close();

            std::cout << "Wrote output/timeseries_R_avg.dat" << std::endl;
            std::cout << "Wrote output/timeseries_R_min.dat" << std::endl;
            std::cout << "Wrote output/timeseries_R_max.dat" << std::endl;
        }
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime);

    std::cout << "\nSimulation completed in " << duration.count() << " seconds" << std::endl;

    // Cleanup
    vkDestroyDevice(ctx.device, nullptr);
    vkDestroyInstance(ctx.instance, nullptr);

    std::cout << "Test completed successfully!" << std::endl;

    return 0;
}