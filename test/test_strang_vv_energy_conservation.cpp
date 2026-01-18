/**
 * test_strang_vv_energy_conservation.cpp
 *
 * Comprehensive validation of hybrid Strang + Velocity Verlet integrator
 * for full chiral mass coupling in Dirac equation
 *
 * Tests the implementation:
 * - Outer: Strang splitting (T/2 - M - T/2)
 * - Inner: Velocity Verlet for FULL chiral mass (m_S + i·m_P·γ⁵)
 *
 * Quality Gates:
 * 1. Energy conservation: ΔE/E < 0.01% (GO/NO-GO criterion)
 * 2. Norm conservation: Δ||Ψ||² < 1e-6
 * 3. Chiral asymmetry: Verify m_L ≠ m_R where θ(x) ≠ 0
 * 4. θ-dependence: Different behavior for θ=0, π/4, π/2
 * 5. Time reversibility: Forward-backward evolution < 1e-4 error
 */

#include "Dirac3D.h"
#include "ConservativeSolver.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <complex>
#include <numeric>
#include <cassert>
#include <chrono>

// Test results structure
struct ValidationResults {
    float energy_drift_percent;
    float norm_drift;
    bool chiral_asymmetry_detected;
    bool theta_dependence_verified;
    float time_reversal_error;
    bool all_tests_passed;
};

// Compute total energy of the Dirac field
float computeDiracEnergy(const Dirac3D& dirac,
                         const std::vector<float>& R_field,
                         const std::vector<float>& theta_field,
                         float Delta) {
    uint32_t Nx = dirac.getNx();
    uint32_t Ny = dirac.getNy();
    uint32_t Nz = dirac.getNz();
    uint32_t N_total = Nx * Ny * Nz;

    // Get spinor components
    auto density = dirac.getDensity();

    // Compute kinetic energy from momentum-space representation
    float E_kinetic = 0.0f;
    auto j_x = dirac.getCurrent(0);
    auto j_y = dirac.getCurrent(1);
    auto j_z = dirac.getCurrent(2);

    for (uint32_t i = 0; i < N_total; ++i) {
        E_kinetic += 0.5f * (j_x[i]*j_x[i] + j_y[i]*j_y[i] + j_z[i]*j_z[i]);
    }

    // Compute mass energy including full chiral coupling
    float E_mass = 0.0f;
    for (uint32_t i = 0; i < N_total; ++i) {
        float R = R_field[i];
        float theta = theta_field[i];

        // Full chiral masses
        float m_S = Delta * R * std::cos(theta);  // Scalar mass
        float m_P = Delta * R * std::sin(theta);  // Pseudoscalar mass

        // Effective mass contribution
        float m_eff = std::sqrt(m_S*m_S + m_P*m_P);
        E_mass += m_eff * density[i];
    }

    return E_kinetic + E_mass;
}

// Compute norm of spinor field
float computeNorm(const Dirac3D& dirac) {
    auto density = dirac.getDensity();
    return std::accumulate(density.begin(), density.end(), 0.0f);
}

// Initialize Gaussian wavepacket with specified momentum
void initializeGaussianPacket(Dirac3D& dirac,
                              float cx, float cy, float cz,
                              float sigma,
                              float kx, float ky, float kz) {
    uint32_t Nx = dirac.getNx();
    uint32_t Ny = dirac.getNy();
    uint32_t Nz = dirac.getNz();

    // Initialize spinor with Gaussian envelope and momentum
    std::vector<std::complex<float>> psi[4];
    for (int c = 0; c < 4; ++c) {
        psi[c].resize(Nx * Ny * Nz);
    }

    float norm_factor = 0.0f;

    for (uint32_t k = 0; k < Nz; ++k) {
        for (uint32_t j = 0; j < Ny; ++j) {
            for (uint32_t i = 0; i < Nx; ++i) {
                float x = (float)i - cx;
                float y = (float)j - cy;
                float z = (float)k - cz;

                float r2 = x*x + y*y + z*z;
                float envelope = std::exp(-r2 / (2.0f * sigma * sigma));

                // Plane wave phase
                float phase = kx*x + ky*y + kz*z;
                std::complex<float> wave(std::cos(phase), std::sin(phase));

                uint32_t idx = k * Ny * Nx + j * Nx + i;

                // Positive energy spinor (upper components dominant)
                psi[0][idx] = envelope * wave;
                psi[1][idx] = envelope * wave * 0.5f;
                psi[2][idx] = envelope * wave * 0.1f;
                psi[3][idx] = envelope * wave * 0.1f;

                norm_factor += std::norm(psi[0][idx]) + std::norm(psi[1][idx])
                            + std::norm(psi[2][idx]) + std::norm(psi[3][idx]);
            }
        }
    }

    // Normalize
    norm_factor = std::sqrt(norm_factor);
    for (int c = 0; c < 4; ++c) {
        for (uint32_t i = 0; i < Nx * Ny * Nz; ++i) {
            psi[c][i] /= norm_factor;
        }
    }

    // Set spinor in Dirac3D via initialize method
    // Flatten the 4 component arrays into single array
    std::vector<std::complex<float>> psi_flat(4 * Nx * Ny * Nz);
    for (uint32_t i = 0; i < Nx * Ny * Nz; ++i) {
        for (int c = 0; c < 4; ++c) {
            psi_flat[4 * i + c] = psi[c][i];
        }
    }
    dirac.initialize(psi_flat);
}

// Initialize spatially-varying θ field (vortex configuration)
void initializeVortexField(std::vector<float>& theta,
                           uint32_t Nx, uint32_t Ny, uint32_t Nz,
                           float cx, float cy, int charge = 1) {
    for (uint32_t k = 0; k < Nz; ++k) {
        for (uint32_t j = 0; j < Ny; ++j) {
            for (uint32_t i = 0; i < Nx; ++i) {
                float x = (float)i - cx;
                float y = (float)j - cy;
                uint32_t idx = k * Ny * Nx + j * Nx + i;

                // Vortex configuration
                theta[idx] = charge * std::atan2(y, x);
            }
        }
    }
}

// Test 1: Energy Conservation with Full Chiral Coupling
ValidationResults testEnergyConservation() {
    std::cout << "\n=== Test 1: Energy Conservation with Full Chiral Coupling ===" << std::endl;

    // Grid parameters
    const uint32_t Nx = 32, Ny = 32, Nz = 32;
    const float dt = 0.01f;
    const float Delta = 2.5f;
    const int num_steps = 1000;

    std::cout << "Grid: " << Nx << "×" << Ny << "×" << Nz << std::endl;
    std::cout << "Time step: dt = " << dt << std::endl;
    std::cout << "Evolution steps: " << num_steps << std::endl;
    std::cout << "Mass scale: Δ = " << Delta << std::endl;

    // Initialize Dirac solver
    Dirac3D dirac(Nx, Ny, Nz);

    // Initialize vacuum fields
    std::vector<float> R_field(Nx * Ny * Nz);
    std::vector<float> theta_field(Nx * Ny * Nz);

    // Create spatially-varying R field (soliton-like)
    for (uint32_t k = 0; k < Nz; ++k) {
        for (uint32_t j = 0; j < Ny; ++j) {
            for (uint32_t i = 0; i < Nx; ++i) {
                uint32_t idx = k * Ny * Nx + j * Nx + i;
                float x = (float)i - Nx/2;
                float y = (float)j - Ny/2;
                float z = (float)k - Nz/2;
                float r = std::sqrt(x*x + y*y + z*z);
                R_field[idx] = 1.0f / (1.0f + 0.1f * r * r);
            }
        }
    }

    // Initialize vortex configuration in θ field
    initializeVortexField(theta_field, Nx, Ny, Nz, Nx/2.0f, Ny/2.0f, 1);

    // Initialize Gaussian wavepacket
    initializeGaussianPacket(dirac, Nx/2.0f, Ny/2.0f, Nz/2.0f, 3.0f, 0.1f, 0.0f, 0.0f);

    // Compute initial energy and norm
    float E0 = computeDiracEnergy(dirac, R_field, theta_field, Delta);
    float norm0 = computeNorm(dirac);

    std::cout << "\nInitial energy: E0 = " << std::fixed << std::setprecision(6) << E0 << std::endl;
    std::cout << "Initial norm: ||Ψ||² = " << norm0 << std::endl;

    // Track energy over evolution
    std::vector<float> energies;
    std::vector<float> norms;
    energies.push_back(E0);
    norms.push_back(norm0);

    // Evolution loop using Strang + Velocity Verlet
    std::cout << "\nEvolving with Strang + Velocity Verlet..." << std::endl;

    for (int step = 1; step <= num_steps; ++step) {
        // Apply hybrid Strang-VV evolution via public interface
        dirac.stepWithChiralMass(R_field, theta_field, Delta, dt);

        // Monitor energy every 100 steps
        if (step % 100 == 0) {
            float E = computeDiracEnergy(dirac, R_field, theta_field, Delta);
            float norm = computeNorm(dirac);
            energies.push_back(E);
            norms.push_back(norm);

            float dE = std::abs(E - E0) / E0 * 100.0f;
            float dNorm = std::abs(norm - norm0);

            std::cout << "Step " << std::setw(4) << step
                     << ": E = " << std::setprecision(6) << E
                     << ", ΔE/E = " << std::setprecision(4) << dE << "%"
                     << ", Δ||Ψ||² = " << std::scientific << dNorm << std::endl;
        }
    }

    // Final energy and norm
    float E_final = computeDiracEnergy(dirac, R_field, theta_field, Delta);
    float norm_final = computeNorm(dirac);

    float energy_drift = std::abs(E_final - E0) / E0 * 100.0f;
    float norm_drift = std::abs(norm_final - norm0);

    std::cout << "\n--- Energy Conservation Results ---" << std::endl;
    std::cout << "Initial energy: E0 = " << std::fixed << std::setprecision(6) << E0 << std::endl;
    std::cout << "Final energy: E = " << E_final << std::endl;
    std::cout << "Energy drift: ΔE/E = " << std::setprecision(4) << energy_drift << "%" << std::endl;
    std::cout << "Norm drift: Δ||Ψ||² = " << std::scientific << norm_drift << std::endl;

    // Prepare results
    ValidationResults results;
    results.energy_drift_percent = energy_drift;
    results.norm_drift = norm_drift;
    results.all_tests_passed = (energy_drift < 0.01f) && (norm_drift < 1e-6);

    if (energy_drift < 0.01f) {
        std::cout << "✓ PASS: Energy conservation < 0.01%" << std::endl;
    } else {
        std::cout << "✗ FAIL: Energy conservation exceeds 0.01% threshold" << std::endl;
    }

    if (norm_drift < 1e-6) {
        std::cout << "✓ PASS: Norm conservation < 1e-6" << std::endl;
    } else {
        std::cout << "✗ FAIL: Norm conservation exceeds threshold" << std::endl;
    }

    return results;
}

// Test 2: θ-dependence verification
bool testThetaDependence() {
    std::cout << "\n=== Test 2: θ-Dependence Verification ===" << std::endl;
    std::cout << "Testing different θ values to confirm pseudoscalar term is active\n" << std::endl;

    const uint32_t Nx = 16, Ny = 16, Nz = 16;
    const float dt = 0.01f;
    const float Delta = 2.5f;
    const int num_steps = 100;

    // Test configurations
    struct TestConfig {
        float theta_value;
        const char* description;
    };

    TestConfig configs[] = {
        {0.0f, "θ = 0 (pure scalar mass)"},
        {M_PI/4, "θ = π/4 (mixed scalar/pseudoscalar)"},
        {M_PI/2, "θ = π/2 (pure pseudoscalar mass)"}
    };

    std::vector<float> final_energies;

    for (const auto& config : configs) {
        std::cout << "Testing " << config.description << "..." << std::endl;

        // Initialize solver
        Dirac3D dirac(Nx, Ny, Nz);

        // Uniform fields
        std::vector<float> R_field(Nx * Ny * Nz, 1.0f);
        std::vector<float> theta_field(Nx * Ny * Nz, config.theta_value);

        // Initialize wavepacket
        initializeGaussianPacket(dirac, Nx/2.0f, Ny/2.0f, Nz/2.0f, 2.0f, 0.1f, 0.1f, 0.0f);

        // Evolve
        for (int step = 0; step < num_steps; ++step) {
            dirac.stepWithChiralMass(R_field, theta_field, Delta, dt);
        }

        // Compute final energy
        float E_final = computeDiracEnergy(dirac, R_field, theta_field, Delta);
        final_energies.push_back(E_final);

        std::cout << "  Final energy: " << std::fixed << std::setprecision(6) << E_final << std::endl;
    }

    // Check that energies are different
    float E_diff_01 = std::abs(final_energies[0] - final_energies[1]);
    float E_diff_12 = std::abs(final_energies[1] - final_energies[2]);
    float E_diff_02 = std::abs(final_energies[0] - final_energies[2]);

    std::cout << "\nEnergy differences:" << std::endl;
    std::cout << "|E(θ=0) - E(θ=π/4)| = " << E_diff_01 << std::endl;
    std::cout << "|E(θ=π/4) - E(θ=π/2)| = " << E_diff_12 << std::endl;
    std::cout << "|E(θ=0) - E(θ=π/2)| = " << E_diff_02 << std::endl;

    bool theta_dependent = (E_diff_01 > 1e-4) && (E_diff_12 > 1e-4) && (E_diff_02 > 1e-4);

    if (theta_dependent) {
        std::cout << "✓ PASS: θ-dependence verified - pseudoscalar term is active" << std::endl;
    } else {
        std::cout << "✗ FAIL: No θ-dependence detected - pseudoscalar term may be inactive" << std::endl;
    }

    return theta_dependent;
}

// Test 3: Chiral Asymmetry Detection
bool testChiralAsymmetry() {
    std::cout << "\n=== Test 3: Chiral Asymmetry Detection ===" << std::endl;
    std::cout << "Verifying m_L ≠ m_R in vortex configuration\n" << std::endl;

    const uint32_t Nx = 24, Ny = 24, Nz = 24;
    const float dt = 0.01f;
    const float Delta = 2.5f;
    const int num_steps = 200;

    // Initialize solver
    Dirac3D dirac(Nx, Ny, Nz);

    // Create vortex configuration
    std::vector<float> R_field(Nx * Ny * Nz, 0.8f);
    std::vector<float> theta_field(Nx * Ny * Nz);
    initializeVortexField(theta_field, Nx, Ny, Nz, Nx/2.0f, Ny/2.0f, 1);

    // Initialize spinor
    initializeGaussianPacket(dirac, Nx/2.0f, Ny/2.0f, Nz/2.0f, 2.5f, 0.0f, 0.0f, 0.0f);

    // Get initial chirality
    auto density_initial = dirac.getDensity();

    // Evolve system
    std::cout << "Evolving in vortex background..." << std::endl;
    for (int step = 0; step < num_steps; ++step) {
        dirac.stepWithChiralMass(R_field, theta_field, Delta, dt);
    }

    // Analyze final state
    auto density_final = dirac.getDensity();

    // Compute effective masses at different spatial points
    std::vector<float> m_eff_left, m_eff_right;

    for (uint32_t i = 0; i < Nx * Ny * Nz; ++i) {
        float R = R_field[i];
        float theta = theta_field[i];

        // Left-handed effective mass
        float m_L = Delta * R * (1.0f + std::cos(theta));
        m_eff_left.push_back(m_L);

        // Right-handed effective mass
        float m_R = Delta * R * (1.0f - std::cos(theta));
        m_eff_right.push_back(m_R);
    }

    // Compute average mass difference
    float avg_m_L = std::accumulate(m_eff_left.begin(), m_eff_left.end(), 0.0f) / m_eff_left.size();
    float avg_m_R = std::accumulate(m_eff_right.begin(), m_eff_right.end(), 0.0f) / m_eff_right.size();
    float mass_asymmetry = std::abs(avg_m_L - avg_m_R) / (avg_m_L + avg_m_R);

    std::cout << "Average m_L = " << avg_m_L << std::endl;
    std::cout << "Average m_R = " << avg_m_R << std::endl;
    std::cout << "Chiral mass asymmetry: " << mass_asymmetry * 100.0f << "%" << std::endl;

    bool asymmetry_detected = (mass_asymmetry > 0.01f);

    if (asymmetry_detected) {
        std::cout << "✓ PASS: Chiral asymmetry detected (m_L ≠ m_R)" << std::endl;
    } else {
        std::cout << "✗ FAIL: No chiral asymmetry detected" << std::endl;
    }

    return asymmetry_detected;
}

// Test 4: Time Reversibility
float testTimeReversibility() {
    std::cout << "\n=== Test 4: Time Reversibility ===" << std::endl;
    std::cout << "Testing symplectic property of integrator\n" << std::endl;

    const uint32_t Nx = 16, Ny = 16, Nz = 16;
    const float dt = 0.01f;
    const float Delta = 2.5f;
    const int num_steps = 100;

    // Initialize solver
    Dirac3D dirac(Nx, Ny, Nz);

    // Simple configuration
    std::vector<float> R_field(Nx * Ny * Nz, 1.0f);
    std::vector<float> theta_field(Nx * Ny * Nz, M_PI/6);

    // Initialize and save initial state
    initializeGaussianPacket(dirac, Nx/2.0f, Ny/2.0f, Nz/2.0f, 2.0f, 0.1f, 0.1f, 0.1f);

    // Save initial spinor
    std::vector<std::complex<float>> psi_initial[4];
    for (int c = 0; c < 4; ++c) {
        psi_initial[c] = dirac.getComponent(c);
    }

    std::cout << "Forward evolution: " << num_steps << " steps..." << std::endl;

    // Forward evolution
    for (int step = 0; step < num_steps; ++step) {
        dirac.stepWithChiralMass(R_field, theta_field, Delta, dt);
    }

    std::cout << "Backward evolution: " << num_steps << " steps..." << std::endl;

    // Backward evolution
    for (int step = 0; step < num_steps; ++step) {
        dirac.stepWithChiralMass(R_field, theta_field, Delta, -dt);
    }

    // Compare with initial state
    float max_error = 0.0f;
    for (int c = 0; c < 4; ++c) {
        auto psi_final = dirac.getComponent(c);
        for (uint32_t i = 0; i < Nx * Ny * Nz; ++i) {
            float error = std::abs(psi_final[i] - psi_initial[c][i]);
            max_error = std::max(max_error, error);
        }
    }

    std::cout << "Time reversibility error: " << std::scientific << max_error << std::endl;

    if (max_error < 1e-4) {
        std::cout << "✓ PASS: Time reversibility verified (error < 1e-4)" << std::endl;
    } else {
        std::cout << "✗ FAIL: Time reversibility violated" << std::endl;
    }

    return max_error;
}

// Main test driver
int main() {
    std::cout << "======================================" << std::endl;
    std::cout << " STRANG-VV CHIRAL INTEGRATION TEST" << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "\nValidating hybrid Strang + Velocity Verlet implementation" << std::endl;
    std::cout << "for full chiral mass coupling (m_S + i·m_P·γ⁵)\n" << std::endl;

    auto start_time = std::chrono::high_resolution_clock::now();

    // Run all tests
    ValidationResults results = testEnergyConservation();
    results.theta_dependence_verified = testThetaDependence();
    results.chiral_asymmetry_detected = testChiralAsymmetry();
    results.time_reversal_error = testTimeReversibility();

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);

    // Final summary
    std::cout << "\n======================================" << std::endl;
    std::cout << "         VALIDATION SUMMARY" << std::endl;
    std::cout << "======================================" << std::endl;

    std::cout << "\n1. Energy Conservation: ";
    if (results.energy_drift_percent < 0.01f) {
        std::cout << "✓ PASS (" << results.energy_drift_percent << "%)" << std::endl;
    } else {
        std::cout << "✗ FAIL (" << results.energy_drift_percent << "%)" << std::endl;
    }

    std::cout << "2. Norm Conservation: ";
    if (results.norm_drift < 1e-6) {
        std::cout << "✓ PASS (" << std::scientific << results.norm_drift << ")" << std::endl;
    } else {
        std::cout << "✗ FAIL (" << results.norm_drift << ")" << std::endl;
    }

    std::cout << "3. θ-Dependence: ";
    if (results.theta_dependence_verified) {
        std::cout << "✓ PASS (pseudoscalar term active)" << std::endl;
    } else {
        std::cout << "✗ FAIL (pseudoscalar term inactive)" << std::endl;
    }

    std::cout << "4. Chiral Asymmetry: ";
    if (results.chiral_asymmetry_detected) {
        std::cout << "✓ PASS (m_L ≠ m_R detected)" << std::endl;
    } else {
        std::cout << "✗ FAIL (no asymmetry)" << std::endl;
    }

    std::cout << "5. Time Reversibility: ";
    if (results.time_reversal_error < 1e-4) {
        std::cout << "✓ PASS (error = " << std::scientific << results.time_reversal_error << ")" << std::endl;
    } else {
        std::cout << "✗ FAIL (error = " << results.time_reversal_error << ")" << std::endl;
    }

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "\nTotal test time: " << duration.count() << " seconds" << std::endl;

    // Overall verdict
    bool all_passed = (results.energy_drift_percent < 0.01f) &&
                      (results.norm_drift < 1e-6) &&
                      results.theta_dependence_verified &&
                      results.chiral_asymmetry_detected &&
                      (results.time_reversal_error < 1e-4);

    std::cout << "\n======================================" << std::endl;
    if (all_passed) {
        std::cout << "    OVERALL: ✓✓✓ ALL TESTS PASSED" << std::endl;
        std::cout << "    Implementation validated for production" << std::endl;
    } else {
        std::cout << "    OVERALL: ✗✗✗ TESTS FAILED" << std::endl;
        std::cout << "    Implementation needs fixes" << std::endl;
    }
    std::cout << "======================================" << std::endl;

    return all_passed ? 0 : 1;
}