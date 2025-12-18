/**
 * @file test_stochastic_R_check.cpp
 * @brief Verify that stochastic mode computes R field correctly
 */

#include "../lib/Nova/Nova.h"
#include "../src/SMFTEngine.h"
#include <iostream>
#include <vector>
#include <cmath>

int main() {
    std::cout << "=== Stochastic R Field Verification ===" << std::endl;

    const int Nx = 64;
    const int Ny = 64;

    NovaConfig config = {
        .name = "Stochastic R Check",
        .screen = {800, 600},
        .debug_level = "error",
        .dimensions = "2D",
        .camera_type = "orthographic",
        .compute = true,
    };
    Nova nova(config);
    nova.initialized = true;

    SMFTEngine engine(&nova);
    engine.initialize(Nx, Ny, 2.5f, 0.0f);

    // Set random initial phases
    std::vector<float> theta_init(Nx * Ny);
    srand(42);
    for (int i = 0; i < Nx * Ny; i++) {
        theta_init[i] = (float(rand()) / RAND_MAX) * 2.0f * M_PI - M_PI;
    }
    engine.setInitialPhases(theta_init);

    // Run 10 stochastic steps
    std::cout << "Running 10 stochastic steps (sigma=1e-4)..." << std::endl;
    for (int step = 0; step < 10; step++) {
        engine.stepStochastic(0.01f, 1.0f, 1e-4f);
    }

    // Get R field
    std::vector<float> R_field = engine.getSyncField();

    // Check if R has non-zero values
    float sum_R = 0.0f;
    float max_R = 0.0f;
    int non_zero_count = 0;

    for (float r : R_field) {
        sum_R += r;
        if (r > 0.0f) non_zero_count++;
        if (r > max_R) max_R = r;
    }

    float mean_R = sum_R / R_field.size();

    std::cout << "\n=== R Field Statistics ===" << std::endl;
    std::cout << "Mean R: " << mean_R << std::endl;
    std::cout << "Max R: " << max_R << std::endl;
    std::cout << "Non-zero values: " << non_zero_count << " / " << R_field.size() << std::endl;

    if (mean_R > 0.0f && non_zero_count > 0) {
        std::cout << "\n✓ SUCCESS: R field is being computed!" << std::endl;
        std::cout << "  Stochastic pipeline is working correctly." << std::endl;
        return 0;
    } else {
        std::cerr << "\n✗ FAILURE: R field is all zeros!" << std::endl;
        std::cerr << "  Stochastic pipeline NOT computing R correctly." << std::endl;
        return 1;
    }
}
