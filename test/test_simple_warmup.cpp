/**
 * Simple warmup test - minimal diagnostic
 * Tests ONLY: Initialize → Warmup → Check R_global
 */

#include "../lib/Nova/Nova.h"
#include "../src/MSFTEngine.h"
#include "../src/MSFTCommon.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>

int main() {
    std::cout << "=== Simple Warmup Test ===" << std::endl;

    const int Nx = 64;
    const int Ny = 64;
    const int N_warmup = 1000;
    const float dt = 0.01f;
    const float K = 27.21f;

    NovaConfig config = {
        .name = "Simple Warmup",
        .screen = {800, 600},
        .debug_level = "error",  // Minimal output
        .dimensions = "2D",
        .camera_type = "orthographic",
        .compute = true,
    };
    Nova nova(config);
    nova.initialized = true;

    MSFTEngine engine(&nova);
    engine.initialize(Nx, Ny, 2.5f, 0.0f);

    // Set random initial phases
    std::vector<float> theta_init(Nx * Ny);
    srand(42);
    for (int i = 0; i < Nx * Ny; i++) {
        theta_init[i] = (float(rand()) / RAND_MAX) * 2.0f * M_PI - M_PI;
    }
    engine.setInitialPhases(theta_init);

    // Check IC upload
    std::vector<float> theta_check = engine.getPhaseField();
    float R_initial = MSFT::compute_global_R(theta_check);
    std::cout << "R_initial (after IC upload): " << R_initial << std::endl;
    std::cout << "First 10 phases: ";
    for (int i = 0; i < 10; i++) {
        std::cout << theta_check[i] << " ";
    }
    std::cout << std::endl;

    // Warmup
    std::cout << "\nRunning warmup (" << N_warmup << " steps)..." << std::endl;
    for (int step = 0; step < N_warmup; step++) {
        engine.step(dt, K, 0.0f);  // Deterministic, damping=0

        if (step % 200 == 0) {
            std::vector<float> theta = engine.getPhaseField();
            float R = MSFT::compute_global_R(theta);
            std::cout << "  Step " << step << ": R_global = " << R;
            std::cout << ", first_phase = " << theta[0];
            std::cout << ", mean_phase = " << (std::accumulate(theta.begin(), theta.end(), 0.0) / theta.size());
            std::cout << std::endl;
        }
    }

    // Final check
    std::vector<float> theta_final = engine.getPhaseField();
    float R_final = MSFT::compute_global_R(theta_final);
    std::cout << "\n=== RESULT ===" << std::endl;
    std::cout << "R_final (after warmup): " << R_final << std::endl;
    std::cout << "First 10 phases: ";
    for (int i = 0; i < 10; i++) {
        std::cout << theta_final[i] << " ";
    }
    std::cout << std::endl;

    if (R_final > 0.8) {
        std::cout << "\n✓ SUCCESS: System synchronized (R > 0.8)" << std::endl;
        return 0;
    } else if (std::isnan(R_final)) {
        std::cout << "\n✗ FAILURE: R is NaN (GPU compute broken)" << std::endl;
        return 2;
    } else {
        std::cout << "\n✗ FAILURE: System did not synchronize (R < 0.8)" << std::endl;
        return 1;
    }
}
