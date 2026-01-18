/**
 * test_dirac_debug.cpp
 *
 * Minimal test to diagnose Dirac norm conservation issue
 */

#include "include/Dirac3D.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

int main() {
    std::cout << "=== Dirac Norm Conservation Debug Test ===" << std::endl;

    // Very small grid and timestep for debugging
    const uint32_t N = 8;
    const float dt = 0.001f;  // Very small timestep
    const float Delta = 0.1f;  // Small mass scale
    const int num_steps = 100;

    std::cout << "Parameters:" << std::endl;
    std::cout << "  Grid: " << N << "^3" << std::endl;
    std::cout << "  dt: " << dt << std::endl;
    std::cout << "  Delta: " << Delta << std::endl;
    std::cout << "  Steps: " << num_steps << std::endl << std::endl;

    // Create solver
    Dirac3D dirac(N, N, N);

    // Initialize with Gaussian
    dirac.initializeGaussian(0.0f, 0.0f, 0.0f, 1.0f);

    // Create uniform fields for simplicity
    std::vector<float> R_field(N*N*N, 1.0f);  // Uniform R=1
    std::vector<float> theta_field(N*N*N, 0.5f);  // Uniform theta=0.5

    float norm_initial = dirac.getNorm();
    std::cout << "Initial norm: " << std::scientific << norm_initial << std::endl << std::endl;

    // Test different evolution methods
    std::cout << "Testing stepWithChiralMass (Strang + VV hybrid):" << std::endl;

    for (int step = 0; step <= num_steps; ++step) {
        if (step % 10 == 0) {
            float norm = dirac.getNorm();
            float drift = (norm - norm_initial) / norm_initial;
            std::cout << "  Step " << std::setw(3) << step
                      << ": norm=" << std::scientific << norm
                      << ", drift=" << std::setprecision(2) << drift * 100 << "%"
                      << std::endl;
        }

        if (step < num_steps) {
            dirac.stepWithChiralMass(R_field, theta_field, Delta, dt);
        }
    }

    float norm_final = dirac.getNorm();
    float total_drift = std::abs(norm_final - norm_initial) / norm_initial;

    std::cout << std::endl;
    std::cout << "Final norm: " << norm_final << std::endl;
    std::cout << "Total drift: " << total_drift * 100 << "%" << std::endl;

    // Test with scalar mass only for comparison
    std::cout << std::endl;
    std::cout << "Testing step() with scalar mass only:" << std::endl;

    Dirac3D dirac2(N, N, N);
    dirac2.initializeGaussian(0.0f, 0.0f, 0.0f, 1.0f);

    // Scalar mass field
    std::vector<float> mass_field(N*N*N);
    for (size_t i = 0; i < N*N*N; ++i) {
        mass_field[i] = Delta * R_field[i] * std::cos(theta_field[i]);
    }

    float norm2_initial = dirac2.getNorm();

    for (int step = 0; step < 100; ++step) {
        dirac2.step(mass_field, dt);
    }

    float norm2_final = dirac2.getNorm();
    float drift2 = std::abs(norm2_final - norm2_initial) / norm2_initial;

    std::cout << "  Initial: " << norm2_initial << std::endl;
    std::cout << "  Final: " << norm2_final << std::endl;
    std::cout << "  Drift: " << drift2 * 100 << "%" << std::endl;

    std::cout << std::endl;
    std::cout << "=== Analysis ===" << std::endl;

    if (total_drift > 0.01) {
        std::cout << "✗ CRITICAL: Norm not conserved (drift > 1%)" << std::endl;
        std::cout << "  This indicates non-unitary evolution" << std::endl;

        // Compute expected oscillation frequency
        float omega_max = Delta * (1.0f + std::cos(0.5f));  // Max frequency
        float T_min = 2 * M_PI / omega_max;  // Min period
        float dt_crit = T_min / 20;  // Need ~20 points per period

        std::cout << "  Max frequency: " << omega_max << " rad/time" << std::endl;
        std::cout << "  Min period: " << T_min << " time units" << std::endl;
        std::cout << "  Critical dt: " << dt_crit << std::endl;

        if (dt > dt_crit) {
            std::cout << "  Current dt=" << dt << " > dt_crit=" << dt_crit << std::endl;
            std::cout << "  → Timestep too large for stability!" << std::endl;
        }
    } else {
        std::cout << "✓ Norm conserved to within 1%" << std::endl;
    }

    return (total_drift > 0.01) ? 1 : 0;
}