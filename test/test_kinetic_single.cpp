/**
 * test_kinetic_single.cpp
 *
 * Test single kinetic step (not split) to isolate the issue
 */

#include "Dirac3D.h"
#include <iostream>
#include <iomanip>
#include <cmath>

// Custom test that applies just kinetic evolution
void testSingleKineticStep() {
    const uint32_t N = 16;
    const float dt = 0.01f;

    Dirac3D dirac(N, N, N);
    dirac.initializeGaussian(0.0f, 0.0f, 0.0f, 2.0f);

    float norm_initial = dirac.getNorm();
    std::cout << "Initial norm: " << std::scientific << norm_initial << std::endl;

    // Apply ONE kinetic step (full dt, not half-step)
    dirac.applyKineticHalfStep(dt);  // Despite name, we pass full dt

    float norm_after = dirac.getNorm();
    float drift = (norm_after - norm_initial) / norm_initial;

    std::cout << "After one kinetic step (dt=" << dt << "):" << std::endl;
    std::cout << "  Norm: " << norm_after << std::endl;
    std::cout << "  Drift: " << drift << " (" << drift*100 << "%)" << std::endl;

    if (std::abs(drift) < 1e-8) {
        std::cout << "✓ Single kinetic step is perfectly unitary" << std::endl;
    } else {
        std::cout << "✗ Single kinetic step has unitarity error!" << std::endl;
    }
}

int main() {
    std::cout << "=== Single Kinetic Step Test ===" << std::endl;
    testSingleKineticStep();

    std::cout << "\n=== Multiple Single Steps Test ===" << std::endl;
    const uint32_t N = 16;
    const float dt = 0.01f;
    const int num_steps = 100;

    Dirac3D dirac(N, N, N);
    dirac.initializeGaussian(0.0f, 0.0f, 0.0f, 2.0f);

    float norm_initial = dirac.getNorm();

    std::cout << std::setw(6) << "Step"
              << std::setw(15) << "Norm"
              << std::setw(15) << "Drift"
              << std::setw(20) << "Drift/step"
              << std::endl;
    std::cout << std::string(56, '-') << std::endl;

    for (int step = 0; step <= num_steps; ++step) {
        if (step % 10 == 0) {
            float norm = dirac.getNorm();
            float drift = (norm - norm_initial) / norm_initial;
            float drift_per_step = (step > 0) ? drift / step : 0;

            std::cout << std::setw(6) << step
                      << std::scientific << std::setprecision(6)
                      << std::setw(15) << norm
                      << std::setw(15) << drift
                      << std::setw(20) << drift_per_step
                      << std::endl;
        }

        if (step < num_steps) {
            dirac.applyKineticHalfStep(dt);
        }
    }

    float norm_final = dirac.getNorm();
    float total_drift = (norm_final - norm_initial) / norm_initial;
    std::cout << "\nTotal drift after " << num_steps << " steps: "
              << total_drift << " (" << total_drift*100 << "%)" << std::endl;
    std::cout << "Drift per step: " << total_drift/num_steps << std::endl;

    return 0;
}