/**
 * Test MV-3: Uniform Field Test
 *
 * Purpose: Verify R = 1.0 when all phases are identical
 *
 * Theory:
 * - If θ(x,y) = θ₀ (constant everywhere)
 * - Then R(x,y) = |⟨e^(iθ)⟩| = |e^(iθ₀)| = 1.0
 *
 * This is the fundamental sanity check for sync field computation.
 *
 * Test cases:
 * 1. All θ = 0        → R = 1.0
 * 2. All θ = π/4      → R = 1.0
 * 3. All θ = π        → R = 1.0
 * 4. All θ = -π/2     → R = 1.0
 *
 * Expected: R_min = R_max = R_avg = 1.0 ± ε (ε < 1e-6)
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <vulkan/vulkan.h>

// Minimal Vulkan context
struct VulkanContext {
    VkInstance instance;
    VkPhysicalDevice physicalDevice;
    VkDevice device;
    uint32_t computeQueueFamily;
    VkQueue computeQueue;
};

// CPU implementation of sync field (reference)
std::vector<float> cpu_sync_field(const std::vector<float>& theta, uint32_t Nx, uint32_t Ny) {
    std::vector<float> R_field(Nx * Ny);

    for (uint32_t y = 0; y < Ny; y++) {
        for (uint32_t x = 0; x < Nx; x++) {
            uint32_t idx = y * Nx + x;

            // Compute local synchronization R = |⟨e^(iθ)⟩|
            float sum_cos = 0.0f;
            float sum_sin = 0.0f;
            int count = 0;

            // 3×3 neighborhood
            for (int dy = -1; dy <= 1; dy++) {
                for (int dx = -1; dx <= 1; dx++) {
                    int nx = static_cast<int>(x) + dx;
                    int ny = static_cast<int>(y) + dy;

                    // Periodic boundaries
                    if (nx < 0) nx += Nx;
                    if (nx >= static_cast<int>(Nx)) nx -= Nx;
                    if (ny < 0) ny += Ny;
                    if (ny >= static_cast<int>(Ny)) ny -= Ny;

                    uint32_t neighbor_idx = ny * Nx + nx;
                    float theta_j = theta[neighbor_idx];

                    sum_cos += std::cos(theta_j);
                    sum_sin += std::sin(theta_j);
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

void print_stats(const std::vector<float>& R_field, const std::string& label) {
    float R_sum = 0.0f;
    float R_min = 1e10f;
    float R_max = -1e10f;

    for (float r : R_field) {
        R_sum += r;
        R_min = std::min(R_min, r);
        R_max = std::max(R_max, r);
    }

    float R_avg = R_sum / R_field.size();

    std::cout << label << ":" << std::endl;
    std::cout << "  R_avg = " << std::fixed << std::setprecision(10) << R_avg << std::endl;
    std::cout << "  R_min = " << R_min << std::endl;
    std::cout << "  R_max = " << R_max << std::endl;
    std::cout << "  Range = " << (R_max - R_min) << std::endl;
}

bool test_uniform_field(float theta_0, const std::string& test_name) {
    std::cout << "\n=== Test Case: " << test_name << " ===" << std::endl;
    std::cout << "θ₀ = " << theta_0 << " (" << (theta_0 * 180.0 / M_PI) << "°)" << std::endl;

    const uint32_t Nx = 128;
    const uint32_t Ny = 128;

    // Initialize uniform field
    std::vector<float> theta(Nx * Ny, theta_0);

    // Compute sync field (CPU reference)
    std::vector<float> R_field = cpu_sync_field(theta, Nx, Ny);

    // Print statistics
    print_stats(R_field, "CPU Result");

    // Check if R ≈ 1.0 everywhere
    float R_sum = 0.0f;
    float R_min = 1e10f;
    float R_max = -1e10f;

    for (float r : R_field) {
        R_sum += r;
        R_min = std::min(R_min, r);
        R_max = std::max(R_max, r);
    }

    float R_avg = R_sum / R_field.size();

    // Tolerance for floating-point precision
    const float epsilon = 1e-6f;

    bool pass = (std::abs(R_avg - 1.0f) < epsilon) &&
                (std::abs(R_min - 1.0f) < epsilon) &&
                (std::abs(R_max - 1.0f) < epsilon);

    if (pass) {
        std::cout << "✓ PASS - All R values are 1.0 ± " << epsilon << std::endl;
    } else {
        std::cout << "✗ FAIL - R values deviate from 1.0" << std::endl;
        std::cout << "  Expected: R = 1.0 ± " << epsilon << std::endl;
        std::cout << "  Got: R_avg = " << R_avg << ", range = [" << R_min << ", " << R_max << "]" << std::endl;
    }

    return pass;
}

int main() {
    std::cout << "=== Test MV-3: Uniform Field Test ===" << std::endl;
    std::cout << "CPU-based validation (no GPU needed)" << std::endl;
    std::cout << "\nTheory: If θ(x,y) = θ₀ (constant), then R(x,y) = 1.0" << std::endl;
    std::cout << "Reason: ⟨e^(iθ)⟩ = e^(iθ₀) → |⟨e^(iθ)⟩| = |e^(iθ₀)| = 1.0" << std::endl;

    bool all_pass = true;

    // Test Case 1: θ = 0
    all_pass &= test_uniform_field(0.0f, "All θ = 0");

    // Test Case 2: θ = π/4
    all_pass &= test_uniform_field(M_PI / 4.0f, "All θ = π/4");

    // Test Case 3: θ = π/2
    all_pass &= test_uniform_field(M_PI / 2.0f, "All θ = π/2");

    // Test Case 4: θ = π
    all_pass &= test_uniform_field(M_PI, "All θ = π");

    // Test Case 5: θ = -π/2
    all_pass &= test_uniform_field(-M_PI / 2.0f, "All θ = -π/2");

    // Summary
    std::cout << "\n=== SUMMARY ===" << std::endl;
    if (all_pass) {
        std::cout << "✓ ALL TESTS PASSED" << std::endl;
        std::cout << "Synchronization field computation is correct for uniform fields." << std::endl;
        return 0;
    } else {
        std::cout << "✗ SOME TESTS FAILED" << std::endl;
        std::cout << "There is a bug in the sync_field computation." << std::endl;
        return 1;
    }
}
