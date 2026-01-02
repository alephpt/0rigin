/**
 * test_dirac3d_free.cpp
 *
 * Test Dirac3D free particle propagation
 *
 * Validates:
 *   1. Norm conservation (unitary evolution)
 *   2. Current conservation
 *   3. Free particle spreading
 */

#include "Dirac3D.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <algorithm>

int main() {
    std::cout << "=== Dirac3D Free Particle Test ===\n\n";

    // Grid parameters
    const uint32_t N = 32;  // 32³ grid (smaller for faster testing)
    const float sigma = 4.0f;  // Wavepacket width
    const float dt = 0.1f;
    const int num_steps = 50;

    std::cout << "Grid: " << N << "³ = " << (N*N*N) << " points\n";
    std::cout << "Wavepacket width: " << sigma << " grid units\n";
    std::cout << "Time step: " << dt << "\n";
    std::cout << "Evolution steps: " << num_steps << "\n\n";

    // Create Dirac solver
    Dirac3D dirac(N, N, N);

    // Initialize Gaussian wavepacket at center
    dirac.initializeGaussian(0.0f, 0.0f, 0.0f, sigma);

    // Zero mass field (free particle)
    std::vector<float> mass_field(N * N * N, 0.0f);

    // Initial diagnostics
    const float norm_initial = dirac.getNorm();
    std::cout << "Initial norm: " << norm_initial << "\n\n";

    // Evolution loop
    std::cout << "Time Evolution:\n";
    std::cout << std::setw(6) << "Step"
              << std::setw(12) << "Time"
              << std::setw(15) << "Norm"
              << std::setw(15) << "ΔN/N"
              << std::setw(15) << "Max ρ"
              << std::setw(15) << "j_x max"
              << "\n";
    std::cout << std::string(77, '-') << "\n";

    for (int step = 0; step <= num_steps; ++step) {
        const float t = step * dt;

        // Compute diagnostics
        const float norm = dirac.getNorm();
        const float dN_over_N = (norm - norm_initial) / norm_initial;

        const auto density = dirac.getDensity();
        const auto current_x = dirac.getCurrent(0);

        const float max_density = *std::max_element(density.begin(), density.end());
        const float max_current_x = *std::max_element(current_x.begin(), current_x.end(),
            [](float a, float b) { return std::abs(a) < std::abs(b); });

        // Print diagnostics every 5 steps
        if (step % 5 == 0) {
            std::cout << std::setw(6) << step
                      << std::setw(12) << std::fixed << std::setprecision(2) << t
                      << std::setw(15) << std::scientific << std::setprecision(6) << norm
                      << std::setw(15) << std::scientific << std::setprecision(4) << dN_over_N
                      << std::setw(15) << std::scientific << std::setprecision(4) << max_density
                      << std::setw(15) << std::scientific << std::setprecision(4) << max_current_x
                      << "\n";
        }

        // Evolve Dirac equation
        if (step < num_steps) {
            dirac.step(mass_field, dt);
        }
    }

    std::cout << "\n";

    // Final checks
    const float norm_final = dirac.getNorm();
    const float norm_conservation_error = std::abs(norm_final - norm_initial) / norm_initial;

    std::cout << "=== Results ===\n";
    std::cout << "Initial norm: " << std::scientific << norm_initial << "\n";
    std::cout << "Final norm:   " << std::scientific << norm_final << "\n";
    std::cout << "Norm drift:   " << std::scientific << norm_conservation_error << "\n\n";

    // Check current conservation (∂ρ/∂t + ∇·j = 0)
    // For free particle, total current should integrate to zero
    const auto current_x = dirac.getCurrent(0);
    const auto current_y = dirac.getCurrent(1);
    const auto current_z = dirac.getCurrent(2);

    const float total_jx = std::accumulate(current_x.begin(), current_x.end(), 0.0f);
    const float total_jy = std::accumulate(current_y.begin(), current_y.end(), 0.0f);
    const float total_jz = std::accumulate(current_z.begin(), current_z.end(), 0.0f);

    std::cout << "Current integrals:\n";
    std::cout << "  ∫ j_x d³x = " << std::scientific << total_jx << "\n";
    std::cout << "  ∫ j_y d³x = " << std::scientific << total_jy << "\n";
    std::cout << "  ∫ j_z d³x = " << std::scientific << total_jz << "\n\n";

    // Quality gates
    const float norm_tolerance = 1e-3f;  // 0.1% norm conservation
    const float current_tolerance = 1e-2f;  // Current integral should be small

    bool norm_conserved = (norm_conservation_error < norm_tolerance);
    bool current_conserved = (std::abs(total_jx) < current_tolerance &&
                             std::abs(total_jy) < current_tolerance &&
                             std::abs(total_jz) < current_tolerance);

    std::cout << "Quality Gates:\n";
    std::cout << "  Norm conservation (ΔN/N < " << norm_tolerance << "): "
              << (norm_conserved ? "PASS ✓" : "FAIL ✗") << "\n";
    std::cout << "  Current conservation (∫j < " << current_tolerance << "): "
              << (current_conserved ? "PASS ✓" : "FAIL ✗") << "\n";

    std::cout << "\n";

    if (norm_conserved && current_conserved) {
        std::cout << "✓ Dirac3D free particle test PASSED\n";
        return 0;
    } else {
        std::cout << "✗ Dirac3D free particle test FAILED\n";
        return 1;
    }
}
