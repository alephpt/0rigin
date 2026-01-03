/**
 * test_3d_full_smft.cpp
 *
 * Week 9-10: Full 3D TRD Integration Test
 *
 * Combines:
 *   - TRDCore3D (Kuramoto synchronization)
 *   - Maxwell3D (electromagnetic fields)
 *   - Dirac3D (fermion dynamics)
 *
 * Validates full TRD evolution cycle in 3D
 */

#include "TRDCore3D.h"
#include "Maxwell3D.h"
#include "Dirac3D.h"
#include <iostream>
#include <iomanip>
#include <cmath>

int main() {
    std::cout << "=== Full 3D TRD Integration Test ===\n\n";

    // Grid parameters (32³ for testing)
    const uint32_t N = 32;
    const float dt = 0.01f;
    const int num_steps = 100;

    std::cout << "Grid: " << N << "³ = " << (N*N*N) << " points\n";
    std::cout << "Time step: " << dt << "\n";
    std::cout << "Evolution steps: " << num_steps << "\n\n";

    // ========================================================================
    // Initialize Kuramoto (TRDCore3D)
    // ========================================================================

    std::cout << "Initializing Kuramoto phase field...\n";
    TRDCore3D smft_core;

    TRDCore3D::Config config;
    config.Nx = N;
    config.Ny = N;
    config.Nz = N;
    config.dx = 1.0f;
    config.dt = dt;
    config.coupling_strength = 1.0f;

    smft_core.initialize(config);
    smft_core.initializeRandom(42);

    const float R_initial = smft_core.getAverageR();
    std::cout << "  Initial synchronization R: " << R_initial << "\n";

    // ========================================================================
    // Initialize Maxwell (electromagnetic fields)
    // ========================================================================

    std::cout << "Initializing electromagnetic fields...\n";
    Maxwell3D maxwell(N, N, N);

    // Initialize with small-amplitude wave
    maxwell.initializeSphericalWave(8.0f, 0.1f);

    const float EM_energy_initial = maxwell.getTotalEnergy();
    std::cout << "  Initial EM energy: " << EM_energy_initial << "\n";

    // ========================================================================
    // Initialize Dirac (fermion field)
    // ========================================================================

    std::cout << "Initializing Dirac spinor field...\n";
    Dirac3D dirac(N, N, N);

    // Gaussian wavepacket
    dirac.initializeGaussian(0.0f, 0.0f, 0.0f, 4.0f);

    const float dirac_norm_initial = dirac.getNorm();
    std::cout << "  Initial Dirac norm: " << dirac_norm_initial << "\n\n";

    // ========================================================================
    // Coupled evolution loop
    // ========================================================================

    std::cout << "Starting coupled 3D TRD evolution...\n";
    std::cout << std::setw(6) << "Step"
              << std::setw(10) << "Time"
              << std::setw(12) << "R_avg"
              << std::setw(15) << "EM_energy"
              << std::setw(15) << "Dirac_norm"
              << "\n";
    std::cout << std::string(58, '-') << "\n";

    for (int step = 0; step <= num_steps; ++step) {
        const float t = step * dt;

        // Compute diagnostics
        if (step % 10 == 0) {
            smft_core.computeRField();
            const float R_avg = smft_core.getAverageR();
            const float EM_energy = maxwell.getTotalEnergy();
            const float dirac_norm = dirac.getNorm();

            std::cout << std::setw(6) << step
                      << std::setw(10) << std::fixed << std::setprecision(3) << t
                      << std::setw(12) << std::fixed << std::setprecision(6) << R_avg
                      << std::setw(15) << std::scientific << std::setprecision(4) << EM_energy
                      << std::setw(15) << std::scientific << std::setprecision(6) << dirac_norm
                      << "\n";
        }

        // Evolve all components
        if (step < num_steps) {
            // 1. Evolve Kuramoto phase field
            smft_core.evolveKuramotoCPU(dt);

            // 2. Evolve Maxwell fields
            maxwell.step(dt);

            // 3. Evolve Dirac spinor
            // Use Kuramoto R-field as effective mass
            const auto& R_field = smft_core.getRField();
            std::vector<float> mass_field(N * N * N);
            for (uint32_t i = 0; i < N * N * N; ++i) {
                mass_field[i] = 1.0f - R_field[i]; // Mass ~ (1 - R)
            }
            dirac.step(mass_field, dt);
        }
    }

    std::cout << "\n";

    // ========================================================================
    // Final validation
    // ========================================================================

    std::cout << "=== Final State ===\n";

    smft_core.computeRField();
    const float R_final = smft_core.getAverageR();
    const float EM_energy_final = maxwell.getTotalEnergy();
    const float dirac_norm_final = dirac.getNorm();

    std::cout << "Kuramoto synchronization:\n";
    std::cout << "  Initial R: " << R_initial << "\n";
    std::cout << "  Final R:   " << R_final << "\n";
    std::cout << "  ΔR:        " << (R_final - R_initial) << "\n\n";

    std::cout << "EM energy:\n";
    std::cout << "  Initial: " << std::scientific << EM_energy_initial << "\n";
    std::cout << "  Final:   " << std::scientific << EM_energy_final << "\n";
    std::cout << "  Drift:   " << std::scientific
              << std::abs(EM_energy_final - EM_energy_initial) / EM_energy_initial << "\n\n";

    std::cout << "Dirac norm:\n";
    std::cout << "  Initial: " << std::scientific << dirac_norm_initial << "\n";
    std::cout << "  Final:   " << std::scientific << dirac_norm_final << "\n";
    std::cout << "  Drift:   " << std::scientific
              << std::abs(dirac_norm_final - dirac_norm_initial) / dirac_norm_initial << "\n\n";

    // Quality gates
    const float R_tolerance = 0.2f;  // R can increase significantly
    const float EM_tolerance = 0.1f;  // 10% EM energy drift
    const float dirac_tolerance = 1e-2f;  // 1% Dirac norm drift

    bool R_valid = (R_final >= 0.0f && R_final <= 1.0f);
    bool EM_conserved = (std::abs(EM_energy_final - EM_energy_initial) / EM_energy_initial < EM_tolerance);
    bool dirac_conserved = (std::abs(dirac_norm_final - dirac_norm_initial) / dirac_norm_initial < dirac_tolerance);

    std::cout << "Quality Gates:\n";
    std::cout << "  R in valid range [0,1]: " << (R_valid ? "PASS ✓" : "FAIL ✗") << "\n";
    std::cout << "  EM energy conserved (< " << EM_tolerance << "): "
              << (EM_conserved ? "PASS ✓" : "FAIL ✗") << "\n";
    std::cout << "  Dirac norm conserved (< " << dirac_tolerance << "): "
              << (dirac_conserved ? "PASS ✓" : "FAIL ✗") << "\n";

    std::cout << "\n";

    if (R_valid && EM_conserved && dirac_conserved) {
        std::cout << "✓ Full 3D TRD integration test PASSED\n";
        std::cout << "\n🎉 3D MIGRATION COMPLETE (Weeks 5-10)\n";
        return 0;
    } else {
        std::cout << "✗ Full 3D TRD integration test FAILED\n";
        return 1;
    }
}
