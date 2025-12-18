/**
 * Find actual critical noise threshold
 * Sweep from σ = 10⁻⁴ to σ = 10
 */

#include "../src/MSFTCommon.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>

using namespace MSFT;

int main() {
    std::cout << "=== Finding Critical Noise Threshold ===" << std::endl << std::endl;

    const uint32_t Nx = 64;
    const uint32_t Ny = 64;
    const float dt = 0.01f;
    const float K = 1.0f;
    const float damping = 0.1f;

    std::mt19937 rng(42);

    // Broader sweep
    std::vector<float> sigma_values = {
        1e-4f, 3e-4f,
        1e-3f, 3e-3f,
        1e-2f, 3e-2f,
        1e-1f, 3e-1f,
        1.0f, 3.0f
    };

    std::cout << "Grid: " << Nx << "×" << Ny << std::endl;
    std::cout << "K = " << K << ", γ = " << damping << ", dt = " << dt << std::endl;
    std::cout << "Simulation: 2000 steps per σ value" << std::endl << std::endl;

    std::cout << std::setw(12) << "sigma" << std::setw(12) << "R_final" << std::setw(15) << "interpretation" << std::endl;
    std::cout << std::string(40, '-') << std::endl;

    for (float sigma : sigma_values) {
        std::vector<float> theta(Nx * Ny, 0.0f);  // Start synchronized
        std::vector<float> omega(Nx * Ny, 0.0f);

        // Run simulation
        for (int step = 0; step < 2000; step++) {
            stepKuramotoWithNoise(theta, omega, dt, K, damping, sigma, Nx, Ny, rng);
        }

        float R_final = computeGlobalR(theta);

        std::cout << std::scientific << std::setprecision(2) << std::setw(12) << sigma
                  << std::fixed << std::setprecision(4) << std::setw(12) << R_final;

        if (R_final > 0.95f) {
            std::cout << "    synchronized";
        } else if (R_final > 0.5f) {
            std::cout << "    partial sync";
        } else if (R_final > 0.2f) {
            std::cout << "    weak sync";
        } else {
            std::cout << "    thermal gas";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
    std::cout << "Critical threshold σ_c is where R drops below 0.95" << std::endl;

    return 0;
}
