/**
 * test_chiral_asymmetry.cpp
 *
 * Focused test on chiral asymmetry detection in the Strang-VV implementation
 * Verifies that the full chiral coupling creates distinct left/right-handed masses
 *
 * Physics:
 * - Chiral mass operator: M = Δ·R·e^{iθγ⁵} = m_S + i·m_P·γ⁵
 * - Left-handed mass: m_L = Δ·R·(1 + cos(θ))
 * - Right-handed mass: m_R = Δ·R·(1 - cos(θ))
 * - When θ ≠ 0, π: m_L ≠ m_R (chiral asymmetry)
 *
 * Test scenarios:
 * 1. Uniform θ field - verify global asymmetry
 * 2. Vortex configuration - spatially-varying asymmetry
 * 3. Domain wall - sharp transition in chirality
 */

#include "Dirac3D.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>
#include <numeric>
#include <fstream>

// Compute chirality expectation value
float computeChirality(const Dirac3D& dirac) {
    uint32_t N_total = dirac.getNx() * dirac.getNy() * dirac.getNz();

    // γ⁵ eigenvalues: +1 for upper components, -1 for lower
    float chirality = 0.0f;

    for (int c = 0; c < 4; ++c) {
        auto psi_c = dirac.getComponent(c);
        float gamma5_eigenvalue = (c < 2) ? 1.0f : -1.0f;

        for (uint32_t i = 0; i < N_total; ++i) {
            chirality += gamma5_eigenvalue * std::norm(psi_c[i]);
        }
    }

    return chirality;
}

// Compute left/right-handed densities
std::pair<float, float> computeChiralDensities(const Dirac3D& dirac) {
    uint32_t N_total = dirac.getNx() * dirac.getNy() * dirac.getNz();

    float density_left = 0.0f;
    float density_right = 0.0f;

    // Upper components (c=0,1) are left-handed in chiral basis
    // Lower components (c=2,3) are right-handed
    for (uint32_t i = 0; i < N_total; ++i) {
        auto psi0 = dirac.getComponent(0)[i];
        auto psi1 = dirac.getComponent(1)[i];
        auto psi2 = dirac.getComponent(2)[i];
        auto psi3 = dirac.getComponent(3)[i];

        density_left += std::norm(psi0) + std::norm(psi1);
        density_right += std::norm(psi2) + std::norm(psi3);
    }

    return {density_left, density_right};
}

// Initialize domain wall configuration
void initializeDomainWall(std::vector<float>& theta,
                          uint32_t Nx, uint32_t Ny, uint32_t Nz,
                          float wall_position, float wall_width) {
    for (uint32_t k = 0; k < Nz; ++k) {
        for (uint32_t j = 0; j < Ny; ++j) {
            for (uint32_t i = 0; i < Nx; ++i) {
                uint32_t idx = k * Ny * Nx + j * Nx + i;
                float x = (float)i;

                // Smooth domain wall: θ transitions from 0 to π
                theta[idx] = M_PI * 0.5f * (1.0f + std::tanh((x - wall_position) / wall_width));
            }
        }
    }
}

// Test 1: Uniform θ field asymmetry
void testUniformFieldAsymmetry() {
    std::cout << "\n=== Test 1: Uniform θ Field Chiral Asymmetry ===" << std::endl;

    const uint32_t Nx = 24, Ny = 24, Nz = 24;
    const float dt = 0.01f;
    const float Delta = 2.5f;
    const int num_steps = 200;

    // Test different uniform θ values
    std::vector<float> theta_values = {0.0f, M_PI/6, M_PI/4, M_PI/3, M_PI/2, 2*M_PI/3, M_PI};

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "\nθ value | m_L/Δ | m_R/Δ | Asymmetry | Chirality" << std::endl;
    std::cout << "--------|-------|-------|-----------|----------" << std::endl;

    for (float theta_val : theta_values) {
        // Initialize solver
        Dirac3D dirac(Nx, Ny, Nz);

        // Uniform fields
        std::vector<float> R_field(Nx * Ny * Nz, 1.0f);
        std::vector<float> theta_field(Nx * Ny * Nz, theta_val);

        // Theoretical masses
        float m_L = Delta * (1.0f + std::cos(theta_val));
        float m_R = Delta * (1.0f - std::cos(theta_val));

        // Initialize spinor (symmetric initially)
        std::vector<std::complex<float>> psi[4];
        for (int c = 0; c < 4; ++c) {
            psi[c].resize(Nx * Ny * Nz);
        }

        // Create localized wavepacket
        for (uint32_t k = 0; k < Nz; ++k) {
            for (uint32_t j = 0; j < Ny; ++j) {
                for (uint32_t i = 0; i < Nx; ++i) {
                    float x = (float)i - Nx/2;
                    float y = (float)j - Ny/2;
                    float z = (float)k - Nz/2;
                    float r2 = x*x + y*y + z*z;
                    float envelope = std::exp(-r2 / 50.0f);

                    uint32_t idx = k * Ny * Nx + j * Nx + i;

                    // Equal mix of all components initially
                    psi[0][idx] = envelope * 0.5f;
                    psi[1][idx] = envelope * 0.5f;
                    psi[2][idx] = envelope * 0.5f;
                    psi[3][idx] = envelope * 0.5f;
                }
            }
        }

        // Initialize via flat array
        std::vector<std::complex<float>> psi_flat(4 * Nx * Ny * Nz);
        for (uint32_t i = 0; i < Nx * Ny * Nz; ++i) {
            for (int c = 0; c < 4; ++c) {
                psi_flat[4 * i + c] = psi[c][i];
            }
        }
        dirac.initialize(psi_flat);

        // Evolve system
        for (int step = 0; step < num_steps; ++step) {
            dirac.stepWithChiralMass(R_field, theta_field, Delta, dt);
        }

        // Measure chirality
        float chirality = computeChirality(dirac);
        auto [rho_L, rho_R] = computeChiralDensities(dirac);
        float asymmetry = (rho_L - rho_R) / (rho_L + rho_R);

        std::cout << std::setw(7) << theta_val/M_PI << "π |"
                 << std::setw(6) << m_L/Delta << " |"
                 << std::setw(6) << m_R/Delta << " |"
                 << std::setw(10) << asymmetry << " |"
                 << std::setw(9) << chirality << std::endl;
    }

    std::cout << "\nNote: Asymmetry = (ρ_L - ρ_R)/(ρ_L + ρ_R)" << std::endl;
    std::cout << "      Positive = left-handed preference" << std::endl;
    std::cout << "      Negative = right-handed preference" << std::endl;
}

// Test 2: Vortex configuration asymmetry
void testVortexAsymmetry() {
    std::cout << "\n=== Test 2: Vortex Configuration Chiral Asymmetry ===" << std::endl;

    const uint32_t Nx = 32, Ny = 32, Nz = 32;
    const float dt = 0.01f;
    const float Delta = 2.5f;
    const int num_steps = 300;

    // Initialize solver
    Dirac3D dirac(Nx, Ny, Nz);

    // R field with radial profile
    std::vector<float> R_field(Nx * Ny * Nz);
    std::vector<float> theta_field(Nx * Ny * Nz);

    // Create vortex with charge n=1
    for (uint32_t k = 0; k < Nz; ++k) {
        for (uint32_t j = 0; j < Ny; ++j) {
            for (uint32_t i = 0; i < Nx; ++i) {
                float x = (float)i - Nx/2;
                float y = (float)j - Ny/2;
                float z = (float)k - Nz/2;
                float r_xy = std::sqrt(x*x + y*y);

                uint32_t idx = k * Ny * Nx + j * Nx + i;

                // Vortex core profile
                R_field[idx] = std::tanh(r_xy / 5.0f);

                // Azimuthal angle
                theta_field[idx] = std::atan2(y, x);
            }
        }
    }

    // Initialize spinor at vortex core
    std::vector<std::complex<float>> psi[4];
    for (int c = 0; c < 4; ++c) {
        psi[c].resize(Nx * Ny * Nz);
    }

    for (uint32_t k = 0; k < Nz; ++k) {
        for (uint32_t j = 0; j < Ny; ++j) {
            for (uint32_t i = 0; i < Nx; ++i) {
                float x = (float)i - Nx/2;
                float y = (float)j - Ny/2;
                float z = (float)k - Nz/2;
                float r2 = x*x + y*y + z*z;
                float envelope = std::exp(-r2 / 25.0f);

                uint32_t idx = k * Ny * Nx + j * Nx + i;

                psi[0][idx] = envelope * 0.7f;
                psi[1][idx] = envelope * 0.3f;
                psi[2][idx] = envelope * 0.3f;
                psi[3][idx] = envelope * 0.2f;
            }
        }
    }

    // Initialize via flat array
    std::vector<std::complex<float>> psi_flat(4 * Nx * Ny * Nz);
    for (uint32_t i = 0; i < Nx * Ny * Nz; ++i) {
        for (int c = 0; c < 4; ++c) {
            psi_flat[4 * i + c] = psi[c][i];
        }
    }
    dirac.initialize(psi_flat);

    // Initial measurements
    auto [rho_L_init, rho_R_init] = computeChiralDensities(dirac);
    float asymmetry_init = (rho_L_init - rho_R_init) / (rho_L_init + rho_R_init);

    std::cout << "\nInitial state:" << std::endl;
    std::cout << "  Left-handed density: " << rho_L_init << std::endl;
    std::cout << "  Right-handed density: " << rho_R_init << std::endl;
    std::cout << "  Asymmetry: " << asymmetry_init << std::endl;

    // Evolve in vortex background
    std::cout << "\nEvolving in vortex background..." << std::endl;

    std::vector<float> asymmetry_history;
    asymmetry_history.push_back(asymmetry_init);

    for (int step = 1; step <= num_steps; ++step) {
        dirac.stepWithChiralMass(R_field, theta_field, Delta, dt);

        if (step % 50 == 0) {
            auto [rho_L, rho_R] = computeChiralDensities(dirac);
            float asymmetry = (rho_L - rho_R) / (rho_L + rho_R);
            asymmetry_history.push_back(asymmetry);

            std::cout << "  Step " << std::setw(3) << step
                     << ": Asymmetry = " << std::fixed << std::setprecision(4)
                     << asymmetry << std::endl;
        }
    }

    // Final measurements
    auto [rho_L_final, rho_R_final] = computeChiralDensities(dirac);
    float asymmetry_final = (rho_L_final - rho_R_final) / (rho_L_final + rho_R_final);

    std::cout << "\nFinal state:" << std::endl;
    std::cout << "  Left-handed density: " << rho_L_final << std::endl;
    std::cout << "  Right-handed density: " << rho_R_final << std::endl;
    std::cout << "  Asymmetry: " << asymmetry_final << std::endl;

    // Analyze spatial distribution
    std::cout << "\nSpatial chirality distribution:" << std::endl;

    // Sample along radial direction
    uint32_t j = Ny/2, k = Nz/2;
    std::cout << "  r  |  θ  | m_L/Δ | m_R/Δ | Local ρ" << std::endl;
    std::cout << "-----|-----|-------|-------|--------" << std::endl;

    for (uint32_t i = Nx/2; i < Nx-2; i += 2) {
        float x = (float)i - Nx/2;
        float r = std::abs(x);
        uint32_t idx = k * Ny * Nx + j * Nx + i;

        float theta = theta_field[idx];
        float R = R_field[idx];
        float m_L = Delta * R * (1.0f + std::cos(theta));
        float m_R = Delta * R * (1.0f - std::cos(theta));

        // Local density
        float local_rho = 0.0f;
        for (int c = 0; c < 4; ++c) {
            auto psi_c = dirac.getSpinorComponent(c);
            local_rho += std::norm(psi_c[idx]);
        }

        std::cout << std::setw(4) << r << " |"
                 << std::setw(4) << theta/M_PI << "π |"
                 << std::setw(6) << m_L/Delta << " |"
                 << std::setw(6) << m_R/Delta << " |"
                 << std::setw(7) << local_rho << std::endl;
    }
}

// Test 3: Domain wall asymmetry
void testDomainWallAsymmetry() {
    std::cout << "\n=== Test 3: Domain Wall Chiral Asymmetry ===" << std::endl;

    const uint32_t Nx = 64, Ny = 16, Nz = 16;
    const float dt = 0.01f;
    const float Delta = 2.5f;
    const int num_steps = 400;

    // Initialize solver
    Dirac3D dirac(Nx, Ny, Nz);

    // Create domain wall
    std::vector<float> R_field(Nx * Ny * Nz, 1.0f);
    std::vector<float> theta_field(Nx * Ny * Nz);

    float wall_position = Nx/2.0f;
    float wall_width = 3.0f;
    initializeDomainWall(theta_field, Nx, Ny, Nz, wall_position, wall_width);

    std::cout << "\nDomain wall configuration:" << std::endl;
    std::cout << "  Position: x = " << wall_position << std::endl;
    std::cout << "  Width: " << wall_width << std::endl;
    std::cout << "  θ transitions from 0 to π" << std::endl;

    // Initialize spinor near wall
    std::vector<std::complex<float>> psi[4];
    for (int c = 0; c < 4; ++c) {
        psi[c].resize(Nx * Ny * Nz);
    }

    for (uint32_t k = 0; k < Nz; ++k) {
        for (uint32_t j = 0; j < Ny; ++j) {
            for (uint32_t i = 0; i < Nx; ++i) {
                float x = (float)i - wall_position;
                float y = (float)j - Ny/2;
                float z = (float)k - Nz/2;
                float r2 = x*x + y*y + z*z;
                float envelope = std::exp(-r2 / 50.0f);

                uint32_t idx = k * Ny * Nx + j * Nx + i;

                // Initial mix
                psi[0][idx] = envelope * 0.5f;
                psi[1][idx] = envelope * 0.5f;
                psi[2][idx] = envelope * 0.5f;
                psi[3][idx] = envelope * 0.5f;
            }
        }
    }

    // Initialize via flat array
    std::vector<std::complex<float>> psi_flat(4 * Nx * Ny * Nz);
    for (uint32_t i = 0; i < Nx * Ny * Nz; ++i) {
        for (int c = 0; c < 4; ++c) {
            psi_flat[4 * i + c] = psi[c][i];
        }
    }
    dirac.initialize(psi_flat);

    std::cout << "\nEvolving near domain wall..." << std::endl;

    // Evolution
    for (int step = 0; step < num_steps; ++step) {
        dirac.stepWithChiralMass(R_field, theta_field, Delta, dt);

        if (step % 100 == 0 && step > 0) {
            std::cout << "  Step " << step << " completed" << std::endl;
        }
    }

    // Analyze chirality across the wall
    std::cout << "\nChirality profile across domain wall:" << std::endl;
    std::cout << "   x   |   θ/π  | m_L/Δ | m_R/Δ | Local chirality" << std::endl;
    std::cout << "-------|--------|-------|-------|----------------" << std::endl;

    uint32_t j = Ny/2, k = Nz/2;

    for (uint32_t i = Nx/2 - 10; i <= Nx/2 + 10; i += 2) {
        uint32_t idx = k * Ny * Nx + j * Nx + i;

        float theta = theta_field[idx];
        float m_L = Delta * (1.0f + std::cos(theta));
        float m_R = Delta * (1.0f - std::cos(theta));

        // Local chirality
        float local_chirality = 0.0f;
        for (int c = 0; c < 4; ++c) {
            auto psi_c = dirac.getSpinorComponent(c);
            float gamma5_eigenvalue = (c < 2) ? 1.0f : -1.0f;
            local_chirality += gamma5_eigenvalue * std::norm(psi_c[idx]);
        }

        std::cout << std::setw(6) << i << " |"
                 << std::setw(7) << theta/M_PI << " |"
                 << std::setw(6) << m_L/Delta << " |"
                 << std::setw(6) << m_R/Delta << " |"
                 << std::setw(15) << local_chirality << std::endl;
    }

    // Compute chirality in each domain
    float chirality_left = 0.0f;   // θ ≈ 0 domain
    float chirality_right = 0.0f;  // θ ≈ π domain

    for (uint32_t i = 0; i < Nx; ++i) {
        for (uint32_t j = 0; j < Ny; ++j) {
            for (uint32_t k = 0; k < Nz; ++k) {
                uint32_t idx = k * Ny * Nx + j * Nx + i;

                float local_chirality = 0.0f;
                for (int c = 0; c < 4; ++c) {
                    auto psi_c = dirac.getSpinorComponent(c);
                    float gamma5_eigenvalue = (c < 2) ? 1.0f : -1.0f;
                    local_chirality += gamma5_eigenvalue * std::norm(psi_c[idx]);
                }

                if (i < Nx/2 - 5) {
                    chirality_left += local_chirality;
                } else if (i > Nx/2 + 5) {
                    chirality_right += local_chirality;
                }
            }
        }
    }

    std::cout << "\nDomain chiralities:" << std::endl;
    std::cout << "  Left domain (θ ≈ 0): " << chirality_left << std::endl;
    std::cout << "  Right domain (θ ≈ π): " << chirality_right << std::endl;
    std::cout << "  Difference: " << (chirality_left - chirality_right) << std::endl;
}

// Main driver
int main() {
    std::cout << "==========================================" << std::endl;
    std::cout << "    CHIRAL ASYMMETRY VALIDATION SUITE" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << "\nValidating full chiral coupling:" << std::endl;
    std::cout << "M = Δ·R·e^{iθγ⁵} = m_S + i·m_P·γ⁵" << std::endl;
    std::cout << "where m_S = Δ·R·cos(θ), m_P = Δ·R·sin(θ)\n" << std::endl;

    // Run all asymmetry tests
    testUniformFieldAsymmetry();
    testVortexAsymmetry();
    testDomainWallAsymmetry();

    std::cout << "\n==========================================" << std::endl;
    std::cout << "        ASYMMETRY TESTS COMPLETE" << std::endl;
    std::cout << "==========================================" << std::endl;

    return 0;
}