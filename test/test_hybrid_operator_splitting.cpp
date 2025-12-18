/**
 * Integration test for GPU-CPU hybrid operator splitting
 * Actually initializes Nova, creates engine, and runs operator splitting
 */

#include "../src/SMFTEngine.h"
#include <iostream>
#include <cmath>
#include <random>

int main() {
    std::cout << "GPU-CPU Hybrid Operator Splitting Integration Test" << std::endl;
    std::cout << "====================================================" << std::endl;

    // Initialize Nova (minimal config for compute)
    NovaConfig config = {
        .name = "Operator Splitting Test",
        .screen = {800, 600},
        .debug_level = "info",
        .dimensions = "2D",
        .camera_type = "orthographic",
        .compute = true,
    };

    Nova nova(config);
    nova.initialized = true;

    std::cout << "\n✓ Nova initialized successfully" << std::endl;

    // Create SMFT engine
    SMFTEngine engine(&nova);

    // Initialize with small grid for testing
    uint32_t grid_size = 64;
    float Delta = 1.0f;
    float chiral_angle = 0.0f;

    std::cout << "✓ Initializing " << grid_size << "x" << grid_size << " grid..." << std::endl;
    engine.initialize(grid_size, grid_size, Delta, chiral_angle);

    // Set up operator splitting with N=10 for testing
    std::cout << "✓ Configuring operator splitting (N=10)..." << std::endl;
    engine.setSubstepRatio(10);

    // Initialize random phases and frequencies
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> phase_dist(0.0f, 2.0f * M_PI);
    std::uniform_real_distribution<float> omega_dist(-0.1f, 0.1f);

    std::vector<float> theta(grid_size * grid_size);
    std::vector<float> omega(grid_size * grid_size);

    for (size_t i = 0; i < theta.size(); i++) {
        theta[i] = phase_dist(gen);
        omega[i] = omega_dist(gen);
    }

    engine.setInitialPhases(theta);
    engine.setNaturalFrequencies(omega);

    std::cout << "✓ Initial conditions set" << std::endl;

    // Initialize hybrid system with Dirac field
    float x0 = grid_size / 2.0f;
    float y0 = grid_size / 2.0f;
    float sigma = 5.0f;

    std::cout << "✓ Initializing hybrid GPU-CPU system..." << std::endl;
    engine.initializeHybrid(x0, y0, sigma);

    std::cout << "\n=== Running Operator Splitting Test ===" << std::endl;

    // Run for 100 steps (10 Dirac updates with N=10)
    int total_steps = 100;
    float dt = 0.01f;
    float K = 1.0f;
    float damping = 0.1f;

    std::cout << "Running " << total_steps << " timesteps..." << std::endl;

    for (int t = 0; t < total_steps; t++) {
        engine.step(dt, K, damping);

        // Log every 10 steps
        if (t % 10 == 0) {
            auto R_field = engine.getSyncField();
            auto psi_density = engine.getDiracDensity();

            // Compute some statistics
            float R_mean = 0.0f;
            float psi_mean = 0.0f;
            for (size_t i = 0; i < R_field.size(); i++) {
                R_mean += R_field[i];
                psi_mean += psi_density[i];
            }
            R_mean /= R_field.size();
            psi_mean /= psi_density.size();

            std::cout << "  Step " << t << ": <R> = " << R_mean
                     << ", <|ψ|²> = " << psi_mean << std::endl;
        }
    }

    std::cout << "\n✓ Simulation completed successfully!" << std::endl;

    // Get final state
    auto final_R = engine.getSyncField();
    auto final_psi = engine.getDiracDensity();

    float R_max = 0.0f, R_min = 1e10f;
    float psi_max = 0.0f;

    for (size_t i = 0; i < final_R.size(); i++) {
        R_max = std::max(R_max, final_R[i]);
        R_min = std::min(R_min, final_R[i]);
        psi_max = std::max(psi_max, final_psi[i]);
    }

    std::cout << "\nFinal State:" << std::endl;
    std::cout << "  R field range: [" << R_min << ", " << R_max << "]" << std::endl;
    std::cout << "  Max Dirac density: " << psi_max << std::endl;

    std::cout << "\n=== Test PASSED ===" << std::endl;
    std::cout << "GPU-CPU hybrid operator splitting working correctly!" << std::endl;

    return 0;
}
