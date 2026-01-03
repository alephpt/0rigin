// test/test_smft_em_integration.cpp
// Test integration of Stückelberg EM into TRDEngine

#include "TRDEngine.h"
#include "DiracEvolution.h"
#include "Nova.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

/**
 * Test: Verify Stückelberg EM integration in TRDEngine
 *
 * Expected behavior:
 * 1. TRDEngine initializes with StuckelbergEM
 * 2. Phase vortex generates magnetic field B_z ≠ 0
 * 3. Dirac evolution includes EM minimal coupling
 * 4. EM observables accessible via getEM_Bz()
 */
int main(int argc, char* argv[]) {
    std::cout << "=== TEST: TRDEngine + Stückelberg EM Integration ===" << std::endl;
    std::cout << std::endl;

    // Initialize Nova
    NovaConfig config{};
    config.name = "TRD EM Test";
    config.screen = {800, 600};
    config.debug_level = "none";
    config.dimensions = "2D";
    config.camera_type = "orbit";
    config.compute = true;

    Nova nova(config);
    if (!nova.initialized) {
        std::cerr << "Failed to initialize Nova" << std::endl;
        return 1;
    }

    // Create TRDEngine
    TRDEngine engine(&nova);

    // Initialize with 64x64 grid
    const int Nx = 64, Ny = 64;
    const float Delta = 1.0f; // Mass gap
    const float chiral_angle = 0.0f;

    engine.initialize(Nx, Ny, Delta, chiral_angle);

    // Create phase vortex (same as test_stuckelberg_vortex_bfield)
    std::vector<float> theta(Nx * Ny);
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            float x = i - Nx/2.0f;
            float y = j - Ny/2.0f;
            theta[j * Nx + i] = std::atan2(y, x);
        }
    }
    engine.setInitialPhases(theta);

    // Set uniform natural frequencies
    std::vector<float> omega(Nx * Ny, 0.0f);
    engine.setNaturalFrequencies(omega);

    // Initialize Dirac field
    engine.initializeDiracField(Nx/2.0f, Ny/2.0f, 5.0f, 1.0f);

    std::cout << "Initial setup:" << std::endl;
    std::cout << "  Grid: " << Nx << " x " << Ny << std::endl;
    std::cout << "  Phase vortex centered at (" << Nx/2 << "," << Ny/2 << ")" << std::endl;
    std::cout << "  Dirac wavepacket at center" << std::endl;
    std::cout << std::endl;

    // Evolve with EM coupling
    const float dt = 0.01f;
    const float K = 1.0f; // Kuramoto coupling
    const float damping = 0.1f;
    const float lambda = 0.1f; // Dirac feedback
    const int substeps = 10; // Operator splitting ratio

    std::cout << "Evolution parameters:" << std::endl;
    std::cout << "  dt = " << dt << std::endl;
    std::cout << "  K = " << K << std::endl;
    std::cout << "  damping = " << damping << std::endl;
    std::cout << "  lambda = " << lambda << std::endl;
    std::cout << "  substeps = " << substeps << std::endl;
    std::cout << std::endl;

    std::cout << "Evolving with Stückelberg EM..." << std::endl;
    std::cout << "Step    max|B_z|    EM_energy   Dirac_norm" << std::endl;
    std::cout << "--------------------------------------------" << std::endl;

    for (int step = 0; step <= 10; step++) {
        if (step > 0) {
            engine.stepWithDirac(dt, lambda, substeps, K, damping);
        }

        // Get EM observables
        const auto& B_field = engine.getEM_Bz();
        float max_B = 0.0f;
        if (!B_field.empty()) {
            max_B = *std::max_element(B_field.begin(), B_field.end());
        }
        float em_energy = engine.getEM_Energy();

        // Get Dirac norm
        float dirac_norm = 0.0f;
        if (engine.getDiracEvolution()) {
            dirac_norm = engine.getDiracEvolution()->getNorm();
        }

        printf("%4d    %9.6f   %9.6f   %9.6f\n",
               step, max_B, em_energy, dirac_norm);
    }

    std::cout << std::endl;

    // Final verification
    const auto& B_field_final = engine.getEM_Bz();
    float max_B_final = *std::max_element(B_field_final.begin(), B_field_final.end());

    std::cout << "=== VERIFICATION ===" << std::endl;

    // Test 1: EM field generated
    bool em_test = (max_B_final > 0.01f);
    std::cout << "EM field generation: "
              << (em_test ? "✓ PASS" : "✗ FAIL")
              << " (B_max = " << max_B_final << ")" << std::endl;

    // Test 2: Dirac norm preserved (unitary evolution)
    float final_norm = engine.getDiracEvolution()->getNorm();
    bool norm_test = (std::abs(final_norm - 1.0f) < 1e-6f);
    std::cout << "Dirac norm preservation: "
              << (norm_test ? "✓ PASS" : "✗ FAIL")
              << " (norm = " << final_norm << ")" << std::endl;

    // Test 3: EM energy computed
    float final_energy = engine.getEM_Energy();
    bool energy_test = (final_energy > 0.0f);
    std::cout << "EM energy computation: "
              << (energy_test ? "✓ PASS" : "✗ FAIL")
              << " (energy = " << final_energy << ")" << std::endl;

    std::cout << std::endl;

    // Overall result
    bool all_pass = em_test && norm_test && energy_test;
    if (all_pass) {
        std::cout << "✅ SUCCESS: Stückelberg EM fully integrated into TRDEngine!" << std::endl;
        std::cout << "   - Phase vortex generates B field via direct φ=θ coupling" << std::endl;
        std::cout << "   - Dirac evolution includes minimal EM coupling" << std::endl;
        std::cout << "   - EM observables accessible and valid" << std::endl;
    } else {
        std::cout << "❌ FAILURE: Integration incomplete" << std::endl;
    }

    return all_pass ? 0 : 1;
}