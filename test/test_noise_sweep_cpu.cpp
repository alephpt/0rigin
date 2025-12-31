/**
 * Rigorous CPU-Based Noise Sweep Experiment
 *
 * Answers methodological questions from immediate.md:
 * Q1: Simulation duration - 10,000 steps (t_max = 100 time units)
 * Q2: Initial conditions - Pre-synchronized (R₀ ≈ 0.99) to test stability
 * Q3: Damping - γ = 0.1 active throughout
 * Q4: Noise implementation - Proper Euler-Maruyama with σ·√(dt)·N(0,1)
 * Q5: Temporal evolution - Full R(t) timeseries saved for all σ
 * Q6: Grid resolution - 128×128 (can be changed for convergence test)
 */

#include <vulkan/vulkan.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <random>
#include <string>
#include <sys/stat.h>

// Simple Vulkan context (only for device query, no compute)
struct ComputeContext {
    VkInstance instance;
    VkPhysicalDevice physicalDevice;
    VkDevice device;
    VkQueue computeQueue;
    uint32_t queueFamilyIndex;
};

bool initVulkan(ComputeContext& ctx) {
    VkApplicationInfo appInfo{};
    appInfo.sType = VK_STRUCTURE_TYPE_APPLICATION_INFO;
    appInfo.pApplicationName = "MSFT CPU Noise Sweep";
    appInfo.applicationVersion = VK_MAKE_VERSION(1, 0, 0);
    appInfo.apiVersion = VK_API_VERSION_1_2;

    VkInstanceCreateInfo instanceInfo{};
    instanceInfo.sType = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO;
    instanceInfo.pApplicationInfo = &appInfo;

    if (vkCreateInstance(&instanceInfo, nullptr, &ctx.instance) != VK_SUCCESS) {
        std::cerr << "Failed to create Vulkan instance" << std::endl;
        return false;
    }

    uint32_t deviceCount = 0;
    vkEnumeratePhysicalDevices(ctx.instance, &deviceCount, nullptr);
    if (deviceCount == 0) return false;

    std::vector<VkPhysicalDevice> devices(deviceCount);
    vkEnumeratePhysicalDevices(ctx.instance, &deviceCount, devices.data());

    for (auto device : devices) {
        VkPhysicalDeviceProperties props;
        vkGetPhysicalDeviceProperties(device, &props);
        std::cout << "GPU: " << props.deviceName << std::endl;

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

    if (ctx.physicalDevice == VK_NULL_HANDLE) return false;

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
        return false;
    }

    vkGetDeviceQueue(ctx.device, ctx.queueFamilyIndex, 0, &ctx.computeQueue);
    return true;
}

// Kuramoto step with proper Euler-Maruyama noise
void kuramotoStepWithNoise(std::vector<float>& theta, const std::vector<float>& omega,
                           uint32_t Nx, uint32_t Ny, float K, float damping,
                           float dt, float sigma, std::mt19937& rng) {
    std::vector<float> theta_new(Nx * Ny);
    std::normal_distribution<float> noise(0.0f, 1.0f);

    for (uint32_t y = 0; y < Ny; y++) {
        for (uint32_t x = 0; x < Nx; x++) {
            uint32_t idx = y * Nx + x;

            // Coupling term (4-neighbor von Neumann)
            float coupling = 0.0f;
            uint32_t neighbors[4][2] = {
                {(x + 1) % Nx, y},
                {(x + Nx - 1) % Nx, y},
                {x, (y + 1) % Ny},
                {x, (y + Ny - 1) % Ny}
            };

            for (int n = 0; n < 4; n++) {
                uint32_t nx = neighbors[n][0];
                uint32_t ny = neighbors[n][1];
                uint32_t nidx = ny * Nx + nx;
                coupling += std::sin(theta[nidx] - theta[idx]);
            }

            // Damping term: -γ·sin(θ)
            float damping_force = -damping * std::sin(theta[idx]);

            // Deterministic drift
            float drift = omega[idx] + (K / 4.0f) * coupling + damping_force;

            // Stochastic term: σ·√(dt)·N(0,1)
            // CRITICAL: Proper Euler-Maruyama scaling
            float noise_term = sigma * std::sqrt(dt) * noise(rng);

            // Update: θ(t+dt) = θ(t) + drift·dt + noise
            theta_new[idx] = theta[idx] + drift * dt + noise_term;
        }
    }

    theta = theta_new;
}

// Compute global Kuramoto order parameter R_global = |⟨e^(iθ)⟩|
float computeGlobalR(const std::vector<float>& theta) {
    double sum_real = 0.0;
    double sum_imag = 0.0;

    for (float t : theta) {
        sum_real += std::cos(t);
        sum_imag += std::sin(t);
    }

    sum_real /= theta.size();
    sum_imag /= theta.size();

    return std::sqrt(sum_real * sum_real + sum_imag * sum_imag);
}

// Compute local R field
std::vector<float> computeLocalRField(const std::vector<float>& theta,
                                       uint32_t Nx, uint32_t Ny, int radius = 1) {
    std::vector<float> R_field(Nx * Ny);

    for (uint32_t y = 0; y < Ny; y++) {
        for (uint32_t x = 0; x < Nx; x++) {
            uint32_t idx = y * Nx + x;

            float sum_cos = 0.0f;
            float sum_sin = 0.0f;
            int count = 0;

            for (int dy = -radius; dy <= radius; dy++) {
                for (int dx = -radius; dx <= radius; dx++) {
                    int nx = (x + dx + Nx) % Nx;
                    int ny = (y + dy + Ny) % Ny;
                    uint32_t nidx = ny * Nx + nx;

                    sum_cos += std::cos(theta[nidx]);
                    sum_sin += std::sin(theta[nidx]);
                    count++;
                }
            }

            float avg_cos = sum_cos / count;
            float avg_sin = sum_sin / count;
            R_field[idx] = std::sqrt(avg_cos * avg_cos + avg_sin * avg_sin);
        }
    }

    return R_field;
}

// Compute localization L = ∫ R⁴ dA
float computeLocalization(const std::vector<float>& R_field) {
    float L = 0.0f;
    for (float R : R_field) {
        L += R * R * R * R;
    }
    return L;
}

int main() {
    std::cout << "=== MSFT Rigorous CPU Noise Sweep ===" << std::endl;
    std::cout << "Methodology aligned with immediate.md requirements\n" << std::endl;

    ComputeContext ctx{};
    if (!initVulkan(ctx)) {
        std::cerr << "Vulkan init failed (non-critical for CPU test)" << std::endl;
    }

    // === PARAMETERS (Answering immediate.md questions) ===

    // Q6: Grid resolution
    const uint32_t Nx = 128;
    const uint32_t Ny = 128;
    std::cout << "Grid: " << Nx << " × " << Ny << std::endl;

    // Physics parameters
    const float Delta = 2.5f;
    const float dt = 0.01f;
    const float K = 1.0f;

    // Q3: Damping parameter (ACTIVE)
    const float damping = 0.1f;
    std::cout << "Damping γ = " << damping << " (ACTIVE)" << std::endl;

    // Q1: Simulation duration
    const int warmup_steps = 5000;      // Warmup to reach R₀ ≈ 0.99
    const int measurement_steps = 10000; // Long measurement for steady state
    std::cout << "Warmup steps: " << warmup_steps << " (t = " << warmup_steps * dt << ")" << std::endl;
    std::cout << "Measurement steps: " << measurement_steps << " (t = " << measurement_steps * dt << ")" << std::endl;

    // Q4: Noise levels for sweep (CORRECTED - covering actual transition)
    std::vector<float> sigma_values = {
        0.0f,       // Deterministic baseline
        1e-5f,      // Falsification threshold
        1e-4f, 1e-3f, 1e-2f,  // Sub-threshold
        0.05f, 0.1f, 0.15f, 0.2f, 0.25f, 0.3f, 0.35f, 0.4f,  // Transition region
        0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f,  // Critical region
        1.5f, 2.0f, 3.0f  // Above threshold
    };

    std::cout << "\nNoise sweep σ values: ";
    for (float s : sigma_values) {
        std::cout << s << " ";
    }
    std::cout << "\n" << std::endl;

    // Create output directory
    system("mkdir -p /home/persist/neotec/0rigin/output/noise_sweep");

    // Random number generator
    std::mt19937 rng(42);  // Fixed seed for reproducibility

    // === MAIN EXPERIMENT LOOP ===

    std::ofstream summary("/home/persist/neotec/0rigin/output/noise_sweep/summary.dat");
    summary << "# Rigorous CPU Noise Sweep Results\n";
    summary << "# Methodology: immediate.md compliant\n";
    summary << "# Grid: " << Nx << "×" << Ny << ", dt = " << dt << ", K = " << K << ", γ = " << damping << "\n";
    summary << "# Warmup: " << warmup_steps << " steps, Measurement: " << measurement_steps << " steps\n";
    summary << "# Noise: Proper Euler-Maruyama σ·√(dt)·N(0,1)\n";
    summary << "# IC: Pre-synchronized (R₀ ≈ 0.99)\n";
    summary << "#\n";
    summary << "# sigma  R_initial  R_final  R_mean  R_std  L_final  L_mean\n";
    summary << std::fixed << std::setprecision(6);

    for (float sigma : sigma_values) {
        std::cout << "\n=== σ = " << std::scientific << sigma << " ===" << std::endl;

        // === Q2: Initial conditions - PRE-SYNCHRONIZED ===
        std::vector<float> theta(Nx * Ny);
        std::vector<float> omega(Nx * Ny, 0.0f);  // Zero natural frequencies for pure sync test

        // Initialize with small perturbation around θ = 0 (synchronized)
        std::uniform_real_distribution<float> init_dist(-0.1f, 0.1f);
        for (uint32_t i = 0; i < Nx * Ny; i++) {
            theta[i] = init_dist(rng);
        }

        float R_initial = computeGlobalR(theta);
        std::cout << "R_initial = " << std::fixed << std::setprecision(4) << R_initial << std::endl;

        // WARMUP PHASE (no noise, reach full synchronization)
        std::cout << "Warmup phase (σ = 0)..." << std::flush;
        for (int step = 0; step < warmup_steps; step++) {
            kuramotoStepWithNoise(theta, omega, Nx, Ny, K, damping, dt, 0.0f, rng);
        }
        float R_warmup = computeGlobalR(theta);
        std::cout << " R_warmup = " << R_warmup << std::endl;

        if (R_warmup < 0.95f) {
            std::cerr << "WARNING: Warmup failed to synchronize (R = " << R_warmup << ")" << std::endl;
        }

        // MEASUREMENT PHASE (with noise σ)
        std::cout << "Measurement phase (σ = " << std::scientific << sigma << ")..." << std::endl;

        // Timeseries storage
        std::vector<float> R_global_series(measurement_steps);
        std::vector<float> L_series(measurement_steps);

        for (int step = 0; step < measurement_steps; step++) {
            // Evolution with noise
            kuramotoStepWithNoise(theta, omega, Nx, Ny, K, damping, dt, sigma, rng);

            // Measure observables
            R_global_series[step] = computeGlobalR(theta);

            std::vector<float> R_field = computeLocalRField(theta, Nx, Ny);
            L_series[step] = computeLocalization(R_field);

            // Progress indicator
            if (step % 2000 == 0 || step == measurement_steps - 1) {
                std::cout << "  Step " << std::setw(5) << step
                          << ": R = " << std::fixed << std::setprecision(4) << R_global_series[step]
                          << ", L = " << std::setw(8) << (int)L_series[step] << std::endl;
            }
        }

        // Compute statistics (last 50% of data for steady state)
        int steady_start = measurement_steps / 2;
        double R_sum = 0.0, R_sum_sq = 0.0;
        double L_sum = 0.0;

        for (int i = steady_start; i < measurement_steps; i++) {
            R_sum += R_global_series[i];
            R_sum_sq += R_global_series[i] * R_global_series[i];
            L_sum += L_series[i];
        }

        int N_samples = measurement_steps - steady_start;
        float R_mean = R_sum / N_samples;
        float R_std = std::sqrt(R_sum_sq / N_samples - R_mean * R_mean);
        float L_mean = L_sum / N_samples;

        float R_final = R_global_series.back();
        float L_final = L_series.back();

        std::cout << "FINAL: R_final = " << std::fixed << std::setprecision(4) << R_final
                  << ", R_mean = " << R_mean << " ± " << R_std
                  << ", L_final = " << (int)L_final << std::endl;

        // Write summary
        summary << std::scientific << sigma << "  "
                << std::fixed << R_warmup << "  " << R_final << "  "
                << R_mean << "  " << R_std << "  "
                << (int)L_final << "  " << (int)L_mean << "\n";
        summary.flush();

        // === Q5: Save full timeseries at critical σ ===
        std::string sigma_str = std::to_string((int)(sigma * 1e7));  // Convert to integer for filename
        std::string ts_filename = "/home/persist/neotec/0rigin/output/noise_sweep/timeseries_sigma_" + sigma_str + ".dat";

        std::ofstream ts_file(ts_filename);
        ts_file << "# Timeseries for σ = " << std::scientific << sigma << "\n";
        ts_file << "# step  t  R_global  L\n";
        ts_file << std::fixed << std::setprecision(6);

        for (int i = 0; i < measurement_steps; i++) {
            ts_file << i << "  " << (i * dt) << "  "
                    << R_global_series[i] << "  " << L_series[i] << "\n";
        }
        ts_file.close();
        std::cout << "Wrote " << ts_filename << std::endl;
    }

    summary.close();
    std::cout << "\n=== EXPERIMENT COMPLETE ===" << std::endl;
    std::cout << "Summary: /home/persist/neotec/0rigin/output/noise_sweep/summary.dat" << std::endl;

    // Cleanup
    if (ctx.device != VK_NULL_HANDLE) {
        vkDestroyDevice(ctx.device, nullptr);
        vkDestroyInstance(ctx.instance, nullptr);
    }

    return 0;
}
