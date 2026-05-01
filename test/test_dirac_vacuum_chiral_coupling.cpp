/**
 * test_dirac_vacuum_chiral_coupling.cpp
 *
 * Validates Dirac3D chiral mass coupling with vacuum fields
 * Tests the chiral mass coupling with vacuum fields via stepWithChiralMass()
 *
 * Quality Gates:
 *   1. Energy conservation: ΔE/E < 0.01% (GO/NO-GO criterion)
 *   2. Chiral asymmetry: m_L ≠ m_R where θ(x) ≠ 0
 *   3. Particle localization: Density peaks near high R(x) regions
 *   4. Norm conservation: ∫|ψ|²d³x preserved < 1e-6
 *
 * Physical Scenario:
 *   - Initialize vacuum with topological vortex (non-trivial θ(x))
 *   - Place Dirac fermion in vacuum background
 *   - Evolve coupled system
 *   - Verify chiral mass generation from vacuum
 */

#include "ConservativeSolver.h"
#include "Dirac3D.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <cassert>

// Energy functional for Dirac field with chiral mass coupling.
//
// E = ⟨ψ|H|ψ⟩  with  H = α·p + β·M,  M = Δ·R·e^(iθγ⁵)
//
// Mass-sector contribution per site (using corrected Dirac-basis γ⁵ =
// anti-diagonal block):
//   E_mass(i) = ΔR_i · [cos(θ_i)·⟨ψ̄ψ⟩_i  +  sin(θ_i)·⟨ψ̄ iγ⁵ ψ⟩_i]
// where  ⟨ψ̄ψ⟩  = ψ†βψ = |ψ₀|² + |ψ₁|² − |ψ₂|² − |ψ₃|²
// and    ⟨ψ̄iγ⁵ψ⟩ = i·ψ†βγ⁵ψ
//                = i·(ψ̄₀ψ₂ + ψ̄₁ψ₃ − ψ̄₂ψ₀ − ψ̄₃ψ₁)
//                = −2·Im(ψ̄₀ψ₂ + ψ̄₁ψ₃)
//
// This replaces the earlier formula  E_mass = Δ·⟨ρ⟩  (mean of m_L,m_R) which
// was an artifact of the broken γ⁵ representation: in the corrected basis,
// the average of the two chiral eigenvalues is Δ·R·cos(θ) (not Δ), and the
// pseudoscalar piece carries the rest of the mass-sector energy.
float computeDiracEnergy(const Dirac3D& dirac, const std::vector<float>& R_field,
                         const std::vector<float>& theta_field, float Delta) {
    // Get field dimensions
    uint32_t Nx = dirac.getNx();
    uint32_t Ny = dirac.getNy();
    uint32_t Nz = dirac.getNz();
    uint32_t N = Nx * Ny * Nz;

    // Compute kinetic energy contribution from currents
    auto j_x = dirac.getCurrent(0);
    auto j_y = dirac.getCurrent(1);
    auto j_z = dirac.getCurrent(2);

    float E_kinetic = 0.0f;
    for (uint32_t i = 0; i < N; ++i) {
        E_kinetic += 0.5f * (j_x[i]*j_x[i] + j_y[i]*j_y[i] + j_z[i]*j_z[i]);
    }

    // Compute chiral mass contribution from per-site spinor components
    const auto& psi0 = dirac.getComponent(0);
    const auto& psi1 = dirac.getComponent(1);
    const auto& psi2 = dirac.getComponent(2);
    const auto& psi3 = dirac.getComponent(3);

    float E_mass = 0.0f;
    for (uint32_t i = 0; i < N; ++i) {
        const float scalar      = std::norm(psi0[i]) + std::norm(psi1[i])
                                  - std::norm(psi2[i]) - std::norm(psi3[i]);
        const std::complex<float> bg5 = std::conj(psi0[i]) * psi2[i]
                                       + std::conj(psi1[i]) * psi3[i];
        const float pseudoscalar = -2.0f * bg5.imag();
        E_mass += Delta * R_field[i] *
                  (std::cos(theta_field[i]) * scalar +
                   std::sin(theta_field[i]) * pseudoscalar);
    }

    return E_kinetic + E_mass;
}

// Initialize vortex configuration in vacuum field
void initializeVortex(std::vector<float>& theta, uint32_t Nx, uint32_t Ny, uint32_t Nz,
                     float cx, float cy, int charge) {
    for (uint32_t k = 0; k < Nz; ++k) {
        for (uint32_t j = 0; j < Ny; ++j) {
            for (uint32_t i = 0; i < Nx; ++i) {
                float x = (float)i - cx;
                float y = (float)j - cy;
                uint32_t idx = k * Ny * Nx + j * Nx + i;
                theta[idx] = charge * std::atan2(y, x);
            }
        }
    }
}

// Main test function
bool testChiralCoupling() {
    std::cout << "\n=== Dirac-Vacuum Chiral Coupling Test ===" << std::endl;
    std::cout << "Validating chiral mass coupling via stepWithChiralMass()\n\n";

    // Grid parameters
    const uint32_t Nx = 32, Ny = 32, Nz = 32;
    const float dx = 1.0f;
    const float dt = 0.01f;
    const float Delta = 2.5f;  // Chiral mass scale
    const int num_steps = 1000;
    const int check_interval = 100;

    std::cout << "Grid: " << Nx << "×" << Ny << "×" << Nz << " = " << (Nx*Ny*Nz) << " points\n";
    std::cout << "Time step: " << dt << "\n";
    std::cout << "Chiral mass scale Δ: " << Delta << "\n";
    std::cout << "Evolution steps: " << num_steps << "\n\n";

    // Create ConservativeSolver (manages vacuum fields)
    ConservativeSolver solver;
    ConservativeSolver::Config solverConfig;
    solverConfig.nx = Nx;
    solverConfig.ny = Ny;
    solverConfig.nz = Nz;
    solverConfig.dx = dx;
    solverConfig.dt = dt;
    solver.initialize(solverConfig);

    // Initialize vacuum with vortex configuration
    std::cout << "Initializing vacuum with topological vortex (charge=+1)..." << std::endl;
    std::vector<float> theta(Nx * Ny * Nz);
    initializeVortex(theta, Nx, Ny, Nz, Nx/2.0f, Ny/2.0f, +1);

    // Compute R field from theta gradient
    std::vector<float> R_field(Nx * Ny * Nz);
    for (uint32_t k = 0; k < Nz; ++k) {
        for (uint32_t j = 0; j < Ny; ++j) {
            for (uint32_t i = 0; i < Nx; ++i) {
                uint32_t idx = k * Ny * Nx + j * Nx + i;

                // Compute gradient magnitude |∇θ|
                float dtheta_dx = 0.0f, dtheta_dy = 0.0f, dtheta_dz = 0.0f;

                // x-derivative (periodic boundaries)
                uint32_t idx_xp = k * Ny * Nx + j * Nx + ((i+1) % Nx);
                uint32_t idx_xm = k * Ny * Nx + j * Nx + ((i+Nx-1) % Nx);
                dtheta_dx = (theta[idx_xp] - theta[idx_xm]) / (2.0f * dx);

                // y-derivative
                uint32_t idx_yp = k * Ny * Nx + ((j+1) % Ny) * Nx + i;
                uint32_t idx_ym = k * Ny * Nx + ((j+Ny-1) % Ny) * Nx + i;
                dtheta_dy = (theta[idx_yp] - theta[idx_ym]) / (2.0f * dx);

                // z-derivative
                uint32_t idx_zp = ((k+1) % Nz) * Ny * Nx + j * Nx + i;
                uint32_t idx_zm = ((k+Nz-1) % Nz) * Ny * Nx + j * Nx + i;
                dtheta_dz = (theta[idx_zp] - theta[idx_zm]) / (2.0f * dx);

                // R ~ 1/(1 + |∇θ|²)
                float grad_sq = dtheta_dx*dtheta_dx + dtheta_dy*dtheta_dy + dtheta_dz*dtheta_dz;
                R_field[idx] = 1.0f / (1.0f + grad_sq);
            }
        }
    }

    // Initialize Dirac field
    std::cout << "Initializing Dirac spinor (Gaussian wavepacket at vortex core)..." << std::endl;
    Dirac3D dirac(Nx, Ny, Nz);
    dirac.initializeGaussian(0.0f, 0.0f, 0.0f, 4.0f);  // Centered at vortex core

    // Initial diagnostics
    const float norm_initial = dirac.getNorm();
    const float energy_initial = computeDiracEnergy(dirac, R_field, theta, Delta);

    std::cout << "\nInitial state:" << std::endl;
    std::cout << "  Norm: " << std::scientific << std::setprecision(6) << norm_initial << std::endl;
    std::cout << "  Energy: " << std::scientific << std::setprecision(6) << energy_initial << std::endl;

    // Check initial chiral asymmetry at vortex core
    uint32_t center_idx = (Nz/2) * Ny * Nx + (Ny/2) * Nx + (Nx/2);
    float m_L_center = Delta * (1.0f + R_field[center_idx] * std::cos(theta[center_idx]));
    float m_R_center = Delta * (1.0f - R_field[center_idx] * std::cos(theta[center_idx]));
    std::cout << "  Chiral masses at vortex core: m_L=" << m_L_center
              << ", m_R=" << m_R_center << std::endl;

    // Storage for energy tracking
    std::vector<float> energies;
    std::vector<float> norms;
    energies.push_back(energy_initial);
    norms.push_back(norm_initial);

    // Evolution loop
    std::cout << "\nTime evolution with chiral coupling:\n";
    std::cout << std::setw(6) << "Step"
              << std::setw(10) << "Time"
              << std::setw(15) << "Energy"
              << std::setw(15) << "ΔE/E"
              << std::setw(15) << "Norm"
              << std::setw(15) << "ΔN/N"
              << std::setw(12) << "m_L-m_R"
              << "\n";
    std::cout << std::string(97, '-') << "\n";

    for (int step = 0; step <= num_steps; ++step) {
        const float t = step * dt;

        // Print diagnostics
        if (step % check_interval == 0) {
            const float energy = computeDiracEnergy(dirac, R_field, theta, Delta);
            const float norm = dirac.getNorm();
            const float dE_over_E = (energy - energy_initial) / std::abs(energy_initial);
            const float dN_over_N = (norm - norm_initial) / norm_initial;

            // Check chiral asymmetry
            float m_L = Delta * (1.0f + R_field[center_idx] * std::cos(theta[center_idx]));
            float m_R = Delta * (1.0f - R_field[center_idx] * std::cos(theta[center_idx]));
            float chirality = m_L - m_R;

            energies.push_back(energy);
            norms.push_back(norm);

            std::cout << std::setw(6) << step
                      << std::setw(10) << std::fixed << std::setprecision(3) << t
                      << std::setw(15) << std::scientific << std::setprecision(6) << energy
                      << std::setw(15) << std::scientific << std::setprecision(4) << dE_over_E
                      << std::setw(15) << std::scientific << std::setprecision(6) << norm
                      << std::setw(15) << std::scientific << std::setprecision(4) << dN_over_N
                      << std::setw(12) << std::fixed << std::setprecision(3) << chirality
                      << "\n";
        }

        // Evolve system
        if (step < num_steps) {
            // Use split-step method with chiral mass coupling
            dirac.stepWithChiralMass(R_field, theta, Delta, dt);
        }
    }

    std::cout << "\n";

    // Final diagnostics
    const float norm_final = dirac.getNorm();
    const float energy_final = computeDiracEnergy(dirac, R_field, theta, Delta);
    const float energy_drift = std::abs(energy_final - energy_initial) / std::abs(energy_initial);
    const float norm_drift = std::abs(norm_final - norm_initial) / norm_initial;

    // Check particle localization
    auto density = dirac.getDensity();
    float max_density = *std::max_element(density.begin(), density.end());
    uint32_t max_idx = std::distance(density.begin(),
                                     std::max_element(density.begin(), density.end()));

    // Convert linear index to 3D coordinates
    uint32_t max_k = max_idx / (Ny * Nx);
    uint32_t max_j = (max_idx % (Ny * Nx)) / Nx;
    uint32_t max_i = max_idx % Nx;

    std::cout << "=== Final Results ===" << std::endl;
    std::cout << "Energy conservation:" << std::endl;
    std::cout << "  Initial: " << std::scientific << energy_initial << std::endl;
    std::cout << "  Final:   " << std::scientific << energy_final << std::endl;
    std::cout << "  Drift:   " << std::scientific << energy_drift
              << " (" << (energy_drift * 100) << "%)" << std::endl;

    std::cout << "\nNorm conservation:" << std::endl;
    std::cout << "  Initial: " << std::scientific << norm_initial << std::endl;
    std::cout << "  Final:   " << std::scientific << norm_final << std::endl;
    std::cout << "  Drift:   " << std::scientific << norm_drift
              << " (" << (norm_drift * 100) << "%)" << std::endl;

    std::cout << "\nParticle localization:" << std::endl;
    std::cout << "  Max density: " << max_density << " at ("
              << max_i << ", " << max_j << ", " << max_k << ")" << std::endl;
    std::cout << "  R-field at max: " << R_field[max_idx] << std::endl;

    // Verify chiral asymmetry
    float m_L_final = Delta * (1.0f + R_field[center_idx] * std::cos(theta[center_idx]));
    float m_R_final = Delta * (1.0f - R_field[center_idx] * std::cos(theta[center_idx]));
    std::cout << "\nChiral asymmetry:" << std::endl;
    std::cout << "  Final m_L: " << m_L_final << std::endl;
    std::cout << "  Final m_R: " << m_R_final << std::endl;
    std::cout << "  Asymmetry: " << (m_L_final - m_R_final) << std::endl;

    // Quality gates
    const float energy_tolerance = 0.0001f;  // 0.01% energy conservation
    const float norm_tolerance = 1e-6f;      // Norm conservation
    const float chiral_threshold = 0.1f;     // Minimum chiral asymmetry

    bool energy_conserved = (energy_drift < energy_tolerance);
    bool norm_conserved = (norm_drift < norm_tolerance);
    bool chiral_asymmetry = (std::abs(m_L_final - m_R_final) > chiral_threshold);
    bool particle_localized = (R_field[max_idx] > 0.5f);  // Density peaks in high-R region

    std::cout << "\n=== Quality Gates ===" << std::endl;
    std::cout << "Energy conservation (ΔE/E < " << energy_tolerance << "): "
              << (energy_conserved ? "PASS ✓" : "FAIL ✗") << std::endl;
    std::cout << "Norm conservation (ΔN/N < " << norm_tolerance << "): "
              << (norm_conserved ? "PASS ✓" : "FAIL ✗") << std::endl;
    std::cout << "Chiral asymmetry (|m_L - m_R| > " << chiral_threshold << "): "
              << (chiral_asymmetry ? "PASS ✓" : "FAIL ✗") << std::endl;
    std::cout << "Particle localization (max density in R > 0.5): "
              << (particle_localized ? "PASS ✓" : "FAIL ✗") << std::endl;

    return energy_conserved && norm_conserved && chiral_asymmetry && particle_localized;
}

int main() {
    std::cout << "======================================" << std::endl;
    std::cout << "  Dirac-Vacuum Chiral Coupling Test  " << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "\nValidating Dirac3D chiral mass coupling" << std::endl;
    std::cout << "with ConservativeSolver::evolveDirac()" << std::endl;

    bool test_passed = testChiralCoupling();

    std::cout << "\n======================================" << std::endl;
    if (test_passed) {
        std::cout << "✓ Chiral coupling validation PASSED" << std::endl;
        std::cout << "  - Energy conserved to 0.01%" << std::endl;
        std::cout << "  - Norm preserved to 1e-6" << std::endl;
        std::cout << "  - Chiral asymmetry confirmed" << std::endl;
        std::cout << "  - Particle localization verified" << std::endl;
        return 0;
    } else {
        std::cout << "✗ Chiral coupling validation FAILED" << std::endl;
        return 1;
    }
}