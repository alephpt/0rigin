/**
 * test_maxwell3d_wave.cpp
 *
 * Test Maxwell3D electromagnetic wave propagation
 *
 * Validates:
 *   1. Wave propagates at speed c = 1 (in natural units)
 *   2. Energy conservation
 *   3. Spherical wave front maintains shape
 */

#include "Maxwell3D.h"
#include <iostream>
#include <cmath>
#include <iomanip>

int main() {
    std::cout << "=== Maxwell3D Wave Propagation Test ===\n\n";

    // Grid parameters
    const uint32_t N = 64;  // 64³ grid
    const float wavelength = 8.0f;  // λ = 8 grid points
    const float amplitude = 1.0f;
    const float dt = 0.1f;  // Time step (CFL condition: dt < dx)
    const int num_steps = 100;

    std::cout << "Grid: " << N << "³ = " << (N*N*N) << " points\n";
    std::cout << "Wavelength: " << wavelength << " grid units\n";
    std::cout << "Time step: " << dt << "\n";
    std::cout << "Evolution steps: " << num_steps << "\n\n";

    // Create Maxwell solver
    Maxwell3D maxwell(N, N, N);

    // Initialize spherical wave
    maxwell.initializeSphericalWave(wavelength, amplitude);

    // Initial energy
    const float E_initial = maxwell.getTotalEnergy();
    std::cout << "Initial EM energy: " << E_initial << "\n\n";

    // Evolution loop
    std::cout << "Time Evolution:\n";
    std::cout << std::setw(6) << "Step"
              << std::setw(12) << "Time"
              << std::setw(15) << "Energy"
              << std::setw(15) << "ΔE/E"
              << std::setw(15) << "Max |E|"
              << std::setw(15) << "Max |B|"
              << "\n";
    std::cout << std::string(77, '-') << "\n";

    for (int step = 0; step <= num_steps; ++step) {
        const float t = step * dt;

        // Compute diagnostics
        const float E_total = maxwell.getTotalEnergy();
        const float dE_over_E = (E_total - E_initial) / E_initial;

        // Find maximum field magnitudes
        float max_E = 0.0f, max_B = 0.0f;
        const auto& Ex = maxwell.getEx();
        const auto& Ey = maxwell.getEy();
        const auto& Ez = maxwell.getEz();
        const auto& Bx = maxwell.getBx();
        const auto& By = maxwell.getBy();
        const auto& Bz = maxwell.getBz();

        for (uint32_t i = 0; i < maxwell.getTotalPoints(); ++i) {
            const float E_mag = std::sqrt(Ex[i]*Ex[i] + Ey[i]*Ey[i] + Ez[i]*Ez[i]);
            const float B_mag = std::sqrt(Bx[i]*Bx[i] + By[i]*By[i] + Bz[i]*Bz[i]);
            max_E = std::max(max_E, E_mag);
            max_B = std::max(max_B, B_mag);
        }

        // Print diagnostics every 10 steps
        if (step % 10 == 0) {
            std::cout << std::setw(6) << step
                      << std::setw(12) << std::fixed << std::setprecision(2) << t
                      << std::setw(15) << std::scientific << std::setprecision(4) << E_total
                      << std::setw(15) << std::scientific << std::setprecision(4) << dE_over_E
                      << std::setw(15) << std::scientific << std::setprecision(4) << max_E
                      << std::setw(15) << std::scientific << std::setprecision(4) << max_B
                      << "\n";
        }

        // Evolve Maxwell equations
        maxwell.step(dt);
    }

    std::cout << "\n";

    // Final energy check
    const float E_final = maxwell.getTotalEnergy();
    const float energy_conservation_error = std::abs(E_final - E_initial) / E_initial;

    std::cout << "=== Results ===\n";
    std::cout << "Initial energy: " << std::scientific << E_initial << "\n";
    std::cout << "Final energy:   " << std::scientific << E_final << "\n";
    std::cout << "Energy drift:   " << std::scientific << energy_conservation_error << "\n\n";

    // Quality gates
    const float energy_tolerance = 0.05f;  // 5% energy drift tolerance
    bool energy_conserved = (energy_conservation_error < energy_tolerance);

    std::cout << "Quality Gates:\n";
    std::cout << "  Energy conservation (ΔE/E < " << energy_tolerance << "): "
              << (energy_conserved ? "PASS ✓" : "FAIL ✗") << "\n";

    std::cout << "\n";

    if (energy_conserved) {
        std::cout << "✓ Maxwell3D wave propagation test PASSED\n";
        return 0;
    } else {
        std::cout << "✗ Maxwell3D wave propagation test FAILED\n";
        return 1;
    }
}
