/**
 * test_dirac_basic_stability.cpp
 *
 * Basic stability test for Dirac3D integrator
 * Tests with uniform constant mass to isolate numerical issues
 */

#include "Dirac3D.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

int main() {
    std::cout << "=== Dirac3D Basic Stability Test ===" << std::endl;

    const uint32_t Nx = 16, Ny = 16, Nz = 16;
    const float dt = 0.01f;
    const float Delta = 2.5f;
    const int num_steps = 100;

    std::cout << "Grid: " << Nx << "×" << Ny << "×" << Nz << std::endl;
    std::cout << "Time step: dt = " << dt << std::endl;
    std::cout << "Evolution steps: " << num_steps << std::endl;

    // Initialize Dirac solver
    Dirac3D dirac(Nx, Ny, Nz);

    // Initialize Gaussian wavepacket using built-in method
    dirac.initializeGaussian(Nx/2.0f, Ny/2.0f, Nz/2.0f, 3.0f);

    // Test 1: Uniform R=1, θ=0 (pure scalar mass)
    std::cout << "\n--- Test 1: Uniform scalar mass (R=1, θ=0) ---" << std::endl;

    std::vector<float> R_field(Nx * Ny * Nz, 1.0f);
    std::vector<float> theta_field(Nx * Ny * Nz, 0.0f);

    float norm0 = dirac.getNorm();
    std::cout << "Initial norm: " << std::fixed << std::setprecision(6) << norm0 << std::endl;

    // Evolve
    for (int step = 1; step <= num_steps; ++step) {
        dirac.stepWithChiralMass(R_field, theta_field, Delta, dt);

        if (step % 20 == 0) {
            float norm = dirac.getNorm();
            float dNorm = std::abs(norm - norm0);
            std::cout << "Step " << std::setw(3) << step
                     << ": norm = " << norm
                     << ", Δnorm = " << std::scientific << dNorm << std::endl;

            if (dNorm > 0.1) {
                std::cout << "❌ UNSTABLE: Norm deviation > 0.1" << std::endl;
                return 1;
            }
        }
    }

    float norm_final = dirac.getNorm();
    float norm_drift = std::abs(norm_final - norm0);

    std::cout << "\nFinal norm: " << std::fixed << norm_final << std::endl;
    std::cout << "Norm drift: " << std::scientific << norm_drift << std::endl;

    if (norm_drift < 1e-3) {
        std::cout << "✓ PASS: Basic stability achieved" << std::endl;
        return 0;
    } else {
        std::cout << "❌ FAIL: Norm drift exceeds threshold" << std::endl;
        return 1;
    }
}