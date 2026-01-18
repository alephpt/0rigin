/**
 * test_dirac_kinetic_only.cpp
 *
 * Test kinetic evolution only (no mass) to isolate the problem
 */

#include "Dirac3D.h"
#include <iostream>
#include <iomanip>
#include <cmath>

int main() {
    std::cout << "=== Dirac Kinetic-Only Evolution Test ===" << std::endl;
    std::cout << "Testing free particle evolution (no mass terms)\n" << std::endl;

    const uint32_t N = 16;
    const float dt = 0.01f;
    const int num_steps = 1000;

    Dirac3D dirac(N, N, N);
    dirac.initializeGaussian(0.0f, 0.0f, 0.0f, 2.0f);

    float norm_initial = dirac.getNorm();

    std::cout << "Initial norm: " << std::scientific << norm_initial << std::endl;
    std::cout << "Time step: " << dt << std::endl;
    std::cout << "Total steps: " << num_steps << "\n" << std::endl;

    // Zero mass field for pure kinetic evolution
    std::vector<float> zero_mass(N*N*N, 0.0f);

    std::cout << std::setw(6) << "Step"
              << std::setw(15) << "Norm"
              << std::setw(15) << "ΔN/N"
              << std::endl;
    std::cout << std::string(36, '-') << std::endl;

    for (int step = 0; step <= num_steps; ++step) {
        if (step % 100 == 0) {
            float norm = dirac.getNorm();
            float drift = (norm - norm_initial) / norm_initial;
            std::cout << std::setw(6) << step
                      << std::scientific << std::setprecision(6)
                      << std::setw(15) << norm
                      << std::setw(15) << drift
                      << std::endl;
        }

        if (step < num_steps) {
            dirac.step(zero_mass, dt);
        }
    }

    float norm_final = dirac.getNorm();
    float total_drift = std::abs(norm_final - norm_initial) / norm_initial;

    std::cout << "\n=== Results ===" << std::endl;
    std::cout << "Initial norm: " << norm_initial << std::endl;
    std::cout << "Final norm: " << norm_final << std::endl;
    std::cout << "Total drift: " << total_drift << " (" << total_drift*100 << "%)" << std::endl;

    if (total_drift < 1e-6) {
        std::cout << "✓ Kinetic evolution is perfectly unitary (drift < 1e-6)" << std::endl;
    } else if (total_drift < 1e-4) {
        std::cout << "✓ Kinetic evolution is unitary (drift < 0.01%)" << std::endl;
    } else {
        std::cout << "✗ Kinetic evolution has unitarity issues!" << std::endl;
    }

    return 0;
}