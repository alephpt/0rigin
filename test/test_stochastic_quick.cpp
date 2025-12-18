/**
 * @file test_stochastic_quick.cpp
 * @brief Quick test to verify stochastic pipeline loads and runs
 */

#include "../lib/Nova/Nova.h"
#include "../src/SMFTEngine.h"
#include <iostream>
#include <vector>
#include <cmath>

int main() {
    std::cout << "=== Quick Stochastic Pipeline Test ===" << std::endl;

    // Small grid for quick test
    const int Nx = 32;
    const int Ny = 32;

    // Initialize phase field
    std::vector<float> theta_init(Nx * Ny);
    srand(42);
    for (int i = 0; i < Nx * Ny; i++) {
        theta_init[i] = (float(rand()) / RAND_MAX) * 2.0f * M_PI - M_PI;
    }

    // Initialize Nova (headless mode)
    NovaConfig config = {
        .name = "Stochastic Quick Test",
        .screen = {800, 600},
        .debug_level = "error",
        .dimensions = "2D",
        .camera_type = "orthographic",
        .compute = true,
    };
    Nova nova(config);
    nova.initialized = true;

    // Initialize SMFT engine
    SMFTEngine engine(&nova);
    engine.initialize(Nx, Ny, 2.5f, 0.0f);
    engine.setInitialPhases(theta_init);

    std::cout << "Running ONE stochastic step (sigma=1e-3)..." << std::endl;

    // Run one stochastic step
    engine.stepStochastic(0.01f, 1.0f, 1e-3f);

    // Get results
    std::vector<float> theta_after = engine.getPhaseField();

    // Compute mean phase change
    float mean_change = 0.0f;
    for (int i = 0; i < Nx * Ny; i++) {
        float diff = theta_after[i] - theta_init[i];
        while (diff > M_PI) diff -= 2.0f * M_PI;
        while (diff < -M_PI) diff += 2.0f * M_PI;
        mean_change += std::abs(diff);
    }
    mean_change /= (Nx * Ny);

    std::cout << "\n=== Results ===" << std::endl;
    std::cout << "Mean phase change: " << mean_change << std::endl;

    if (mean_change > 0.0f && mean_change < 1.0f) {
        std::cout << "✓ Stochastic pipeline appears to be working!" << std::endl;
        std::cout << "  (phases changed by reasonable amount)" << std::endl;
        return 0;
    } else {
        std::cerr << "✗ Suspicious result - check pipeline" << std::endl;
        return 1;
    }
}
