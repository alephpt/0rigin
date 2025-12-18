/**
 * @file test_prng_verify.cpp
 * @brief Verification Test 1: PRNG Quality (from Directive.md)
 *
 * Tests PCG + Box-Muller implementation:
 * - Mean should be ≈ 0
 * - Variance should be ≈ 1
 *
 * Success Criteria:
 * - |mean| < 0.1
 * - 0.9 < variance < 1.1
 *
 * If FAILS: PRNG is broken → do NOT proceed with noise sweep
 */

#include "../lib/Nova/Nova.h"
#include "../src/MSFTEngine.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

int main() {
    std::cout << "=== PRNG Verification Test ===" << std::endl;
    std::cout << "Testing PCG + Box-Muller implementation\n" << std::endl;

    // Initialize Nova (compute-only mode)
    NovaConfig config = {
        .name = "PRNG Verification",
        .screen = {800, 600},
        .debug_level = "error",  // Reduce debug spam
        .dimensions = "2D",
        .camera_type = "orthographic",
        .compute = true,
    };
    Nova nova(config);
    nova.initialized = true;

    // Small grid for quick test
    const int Nx = 64;
    const int Ny = 64;
    const int N_samples = Nx * Ny;

    std::cout << "Grid: " << Nx << " x " << Ny << " (" << N_samples << " samples)" << std::endl;

    // Initialize MSFT engine
    MSFTEngine engine(&nova);
    engine.initialize(Nx, Ny, 2.5f, 0.0f);

    // Set random initial phases
    std::vector<float> theta_init(N_samples);
    srand(42);
    for (int i = 0; i < N_samples; i++) {
        theta_init[i] = (float(rand()) / RAND_MAX) * 2.0f * M_PI - M_PI;
    }
    engine.setInitialPhases(theta_init);

    // Run ONE stochastic step with sigma=1.0 to generate random numbers
    std::cout << "\nRunning single stochastic step (sigma=1.0, dt=1.0)..." << std::endl;
    engine.stepStochastic(1.0f, 1.0f, 1.0f);

    // Get the phase field (which should have noise added)
    std::vector<float> theta_after = engine.getPhaseField();

    // Compute differences (noise added)
    std::vector<float> noise_samples(N_samples);
    for (int i = 0; i < N_samples; i++) {
        // Unwrap angle difference
        float diff = theta_after[i] - theta_init[i];
        while (diff > M_PI) diff -= 2.0f * M_PI;
        while (diff < -M_PI) diff += 2.0f * M_PI;
        noise_samples[i] = diff;
    }

    // Compute statistics
    float sum = 0.0f;
    float sum_sq = 0.0f;
    for (float x : noise_samples) {
        sum += x;
        sum_sq += x * x;
    }

    float mean = sum / N_samples;
    float variance = (sum_sq / N_samples) - (mean * mean);
    float std_dev = sqrt(variance);

    std::cout << "\n=== Results ===" << std::endl;
    std::cout << "Mean: " << mean << std::endl;
    std::cout << "Variance: " << variance << std::endl;
    std::cout << "Std Dev: " << std_dev << std::endl;

    // Success criteria
    bool mean_ok = fabs(mean) < 0.1f;
    bool var_ok = (variance > 0.9f) && (variance < 1.1f);

    std::cout << "\n=== Verification ===" << std::endl;
    std::cout << "Mean check (|mean| < 0.1): " << (mean_ok ? "✓ PASS" : "✗ FAIL") << std::endl;
    std::cout << "Variance check (0.9 < var < 1.1): " << (var_ok ? "✓ PASS" : "✗ FAIL") << std::endl;

    if (mean_ok && var_ok) {
        std::cout << "\n✓✓✓ PRNG VERIFICATION PASSED ✓✓✓" << std::endl;
        std::cout << "Proceed to noise sweep experiment." << std::endl;
        return 0;
    } else {
        std::cout << "\n✗✗✗ PRNG VERIFICATION FAILED ✗✗✗" << std::endl;
        std::cout << "DO NOT proceed with noise sweep." << std::endl;
        std::cout << "Debug PRNG implementation first." << std::endl;
        return 1;
    }
}
