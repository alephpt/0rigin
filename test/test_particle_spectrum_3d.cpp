/**
 * test_particle_spectrum_3d.cpp
 *
 * B1 CRITICAL TEST: Particle Mass Spectrum Derivation
 *
 * Goal: Validate TRD's core prediction that particle masses emerge from
 *       topological vortex excitations in the Kuramoto phase field
 *
 * Physics:
 *   - Kuramoto phase field θ(x,y,z) supports topological defects
 *   - Vortex types: point vortex (Q=1), double vortex (Q=2), knots (Q≥3)
 *   - Effective mass: m_eff = E_vortex / c²
 *   - Vortex energy: E = ∫(|∇θ|² + V(R)) d³x
 *   - Topological charge: Q = (1/2π)∮∇θ·dl
 *   - Expected: E_n ~ n·E_fundamental (n = winding number)
 *
 * Critical Quality Gate:
 *   m_2/m_1 ≈ 206.768 ± 50% (electron/muon mass ratio)
 *
 * Memory Reference:
 *   "Bekenstein-Hawking determines clustering via gravitational surface tension.
 *    Δ = √(ℏc/G) = Planck Mass"
 *
 * This test is MAKE-OR-BREAK for TRD theory validation.
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include <functional>

const float PI = 3.14159265358979323846f;
const float HBAR = 1.0f;  // Natural units
const float C_LIGHT = 1.0f;  // Natural units

/**
 * 3D Grid for Kuramoto phase field
 */
class KuramotoGrid3D {
private:
    int N;
    float dx;
    float K;  // Kuramoto coupling strength
    std::vector<float> theta;  // Phase field
    std::vector<float> R;      // Synchronization field

public:
    KuramotoGrid3D(int size, float spacing, float coupling)
        : N(size), dx(spacing), K(coupling),
          theta(size * size * size, 0.0f),
          R(size * size * size, 1.0f) {}

    int getSize() const { return N; }
    float getSpacing() const { return dx; }
    float getCoupling() const { return K; }

    float& theta_at(int ix, int iy, int iz) {
        return theta[ix + N * (iy + N * iz)];
    }

    const float& theta_at(int ix, int iy, int iz) const {
        return theta[ix + N * (iy + N * iz)];
    }

    float& R_at(int ix, int iy, int iz) {
        return R[ix + N * (iy + N * iz)];
    }

    const float& R_at(int ix, int iy, int iz) const {
        return R[ix + N * (iy + N * iz)];
    }

    // Periodic boundary conditions
    int wrap(int i) const {
        return (i + N) % N;
    }
};

/**
 * Initialize single vortex (Q=1) - Fundamental excitation
 * Phase: θ(x,y,z) = atan2(y-y₀, x-x₀) (constant along z)
 */
void initSingleVortex(KuramotoGrid3D& grid, float x0, float y0) {
    const int N = grid.getSize();
    const float dx = grid.getSpacing();

    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            for (int iz = 0; iz < N; ++iz) {
                float x = (ix - N/2) * dx;
                float y = (iy - N/2) * dx;

                float theta = std::atan2(y - y0, x - x0);
                grid.theta_at(ix, iy, iz) = theta;
            }
        }
    }

    std::cout << "  Initialized: Single vortex (Q=1) at (" << x0 << ", " << y0 << ")\n";
}

/**
 * Initialize double vortex (Q=2) - Second excitation level
 * Two vortices with same winding direction
 */
void initDoubleVortex(KuramotoGrid3D& grid, float separation) {
    const int N = grid.getSize();
    const float dx = grid.getSpacing();

    // Two vortices at ±separation/2 along x-axis
    float x1 = -separation / 2.0f;
    float x2 = +separation / 2.0f;
    float y0 = 0.0f;

    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            for (int iz = 0; iz < N; ++iz) {
                float x = (ix - N/2) * dx;
                float y = (iy - N/2) * dx;

                // Superposition of two vortex phases
                float theta1 = std::atan2(y - y0, x - x1);
                float theta2 = std::atan2(y - y0, x - x2);

                // Combined phase (additive for same winding)
                grid.theta_at(ix, iy, iz) = theta1 + theta2;
            }
        }
    }

    std::cout << "  Initialized: Double vortex (Q=2) with separation " << separation << "\n";
}

/**
 * Initialize triple vortex (Q=3) - Third excitation level
 * Three vortices in triangular configuration
 */
void initTripleVortex(KuramotoGrid3D& grid, float radius) {
    const int N = grid.getSize();
    const float dx = grid.getSpacing();

    // Three vortices at 120° intervals
    std::array<float, 3> x_pos = {
        radius * std::cos(0.0f),
        radius * std::cos(2.0f * PI / 3.0f),
        radius * std::cos(4.0f * PI / 3.0f)
    };
    std::array<float, 3> y_pos = {
        radius * std::sin(0.0f),
        radius * std::sin(2.0f * PI / 3.0f),
        radius * std::sin(4.0f * PI / 3.0f)
    };

    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            for (int iz = 0; iz < N; ++iz) {
                float x = (ix - N/2) * dx;
                float y = (iy - N/2) * dx;

                float theta_total = 0.0f;
                for (int v = 0; v < 3; ++v) {
                    float theta_v = std::atan2(y - y_pos[v], x - x_pos[v]);
                    theta_total += theta_v;
                }

                grid.theta_at(ix, iy, iz) = theta_total;
            }
        }
    }

    std::cout << "  Initialized: Triple vortex (Q=3) with radius " << radius << "\n";
}

/**
 * Compute total energy of vortex configuration
 * E = ∫[(∇θ)² + V(R)] d³x
 *
 * Where:
 *   - Gradient energy: (∇θ)²
 *   - Potential energy: V(R) = K·R²·(1 - cos(θ_i - θ_j))
 */
float computeVortexEnergy(const KuramotoGrid3D& grid) {
    const int N = grid.getSize();
    const float dx = grid.getSpacing();
    const float K = grid.getCoupling();

    float total_energy = 0.0f;

    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            for (int iz = 0; iz < N; ++iz) {
                // Central difference for gradient
                int ix_p = grid.wrap(ix + 1);
                int ix_m = grid.wrap(ix - 1);
                int iy_p = grid.wrap(iy + 1);
                int iy_m = grid.wrap(iy - 1);
                int iz_p = grid.wrap(iz + 1);
                int iz_m = grid.wrap(iz - 1);

                float theta_c = grid.theta_at(ix, iy, iz);

                // Gradient energy: (∇θ)²
                float dtheta_dx = (grid.theta_at(ix_p, iy, iz) - grid.theta_at(ix_m, iy, iz)) / (2.0f * dx);
                float dtheta_dy = (grid.theta_at(ix, iy_p, iz) - grid.theta_at(ix, iy_m, iz)) / (2.0f * dx);
                float dtheta_dz = (grid.theta_at(ix, iy, iz_p) - grid.theta_at(ix, iy, iz_m)) / (2.0f * dx);

                float grad_sq = dtheta_dx * dtheta_dx + dtheta_dy * dtheta_dy + dtheta_dz * dtheta_dz;

                // Potential energy from Kuramoto coupling
                float R = grid.R_at(ix, iy, iz);
                float potential = 0.0f;

                // Average over 6 neighbors
                std::array<float, 6> neighbors = {
                    grid.theta_at(ix_p, iy, iz),
                    grid.theta_at(ix_m, iy, iz),
                    grid.theta_at(ix, iy_p, iz),
                    grid.theta_at(ix, iy_m, iz),
                    grid.theta_at(ix, iy, iz_p),
                    grid.theta_at(ix, iy, iz_m)
                };

                for (float theta_n : neighbors) {
                    potential += K * R * R * (1.0f - std::cos(theta_c - theta_n));
                }
                potential /= 6.0f;

                // Total energy density
                float energy_density = 0.5f * grad_sq + potential;
                total_energy += energy_density * dx * dx * dx;
            }
        }
    }

    return total_energy;
}

/**
 * Compute topological charge Q
 * Q = (1/2π) ∮ ∇θ·dl
 *
 * For vortex: measure phase winding around contour
 */
float computeTopologicalCharge(const KuramotoGrid3D& grid) {
    const int N = grid.getSize();
    const int mid = N / 2;
    const int radius = N / 4;  // Contour radius

    float phase_sum = 0.0f;
    int num_points = 100;

    for (int i = 0; i < num_points; ++i) {
        float angle = 2.0f * PI * i / num_points;
        int ix = mid + static_cast<int>(radius * std::cos(angle));
        int iy = mid + static_cast<int>(radius * std::sin(angle));
        int iz = mid;

        ix = grid.wrap(ix);
        iy = grid.wrap(iy);

        if (i > 0) {
            float angle_prev = 2.0f * PI * (i - 1) / num_points;
            int ix_prev = mid + static_cast<int>(radius * std::cos(angle_prev));
            int iy_prev = mid + static_cast<int>(radius * std::sin(angle_prev));

            ix_prev = grid.wrap(ix_prev);
            iy_prev = grid.wrap(iy_prev);

            float theta_curr = grid.theta_at(ix, iy, iz);
            float theta_prev = grid.theta_at(ix_prev, iy_prev, iz);

            // Unwrap phase difference
            float dtheta = theta_curr - theta_prev;
            while (dtheta > PI) dtheta -= 2.0f * PI;
            while (dtheta < -PI) dtheta += 2.0f * PI;

            phase_sum += dtheta;
        }
    }

    // Close the loop
    float theta_final = grid.theta_at(mid + radius, mid, mid);
    float theta_initial = grid.theta_at(mid + radius, mid, mid);
    float dtheta = theta_final - theta_initial;
    while (dtheta > PI) dtheta -= 2.0f * PI;
    while (dtheta < -PI) dtheta += 2.0f * PI;
    phase_sum += dtheta;

    return phase_sum / (2.0f * PI);
}

/**
 * Test 1: Single vortex (Q=1) → Fundamental mass m₁
 */
bool testSingleVortexMass() {
    std::cout << "\n=== Test 1: Single Vortex (Q=1) - Fundamental Mass ===\n";

    const int N = 32;
    const float dx = 0.5f;
    const float K = 1.0f;

    KuramotoGrid3D grid(N, dx, K);
    initSingleVortex(grid, 0.0f, 0.0f);

    float energy = computeVortexEnergy(grid);
    float charge = computeTopologicalCharge(grid);
    float mass = energy / (C_LIGHT * C_LIGHT);

    std::cout << "  Topological charge Q: " << charge << "\n";
    std::cout << "  Vortex energy E₁: " << energy << "\n";
    std::cout << "  Effective mass m₁: " << mass << "\n";

    // Quality gate: Q should be close to 1
    bool charge_valid = std::abs(charge - 1.0f) < 0.2f;

    std::cout << "  Topological charge (≈1.0): " << (charge_valid ? "PASS ✓" : "FAIL ✗") << "\n";

    return charge_valid;
}

/**
 * Test 2: Double vortex (Q=2) → Second mass m₂
 */
bool testDoubleVortexMass() {
    std::cout << "\n=== Test 2: Double Vortex (Q=2) - Second Mass ===\n";

    const int N = 32;
    const float dx = 0.5f;
    const float K = 1.0f;

    KuramotoGrid3D grid(N, dx, K);
    initDoubleVortex(grid, 4.0f);  // Separation = 4.0

    float energy = computeVortexEnergy(grid);
    float charge = computeTopologicalCharge(grid);
    float mass = energy / (C_LIGHT * C_LIGHT);

    std::cout << "  Topological charge Q: " << charge << "\n";
    std::cout << "  Vortex energy E₂: " << energy << "\n";
    std::cout << "  Effective mass m₂: " << mass << "\n";

    // Quality gate: Q should be close to 2
    bool charge_valid = std::abs(charge - 2.0f) < 0.4f;

    std::cout << "  Topological charge (≈2.0): " << (charge_valid ? "PASS ✓" : "FAIL ✗") << "\n";

    return charge_valid;
}

/**
 * Test 3: Mass ratio m₂/m₁ → Compare to electron/muon (206.768)
 */
bool testMassRatio() {
    std::cout << "\n=== Test 3: Mass Ratio m₂/m₁ → Electron/Muon ===\n";

    const int N = 32;
    const float dx = 0.5f;
    const float K = 1.0f;

    // Measure m₁
    KuramotoGrid3D grid1(N, dx, K);
    initSingleVortex(grid1, 0.0f, 0.0f);
    float E1 = computeVortexEnergy(grid1);
    float m1 = E1 / (C_LIGHT * C_LIGHT);

    // Measure m₂
    KuramotoGrid3D grid2(N, dx, K);
    initDoubleVortex(grid2, 4.0f);
    float E2 = computeVortexEnergy(grid2);
    float m2 = E2 / (C_LIGHT * C_LIGHT);

    // Mass ratio
    float ratio_measured = m2 / m1;
    float ratio_expected = 206.768f;  // m_muon / m_electron
    float error = std::abs(ratio_measured - ratio_expected) / ratio_expected;

    std::cout << "\n  Fundamental mass m₁: " << m1 << "\n";
    std::cout << "  Second mass m₂: " << m2 << "\n";
    std::cout << "  Measured ratio m₂/m₁: " << ratio_measured << "\n";
    std::cout << "  Expected ratio (muon/electron): " << ratio_expected << "\n";
    std::cout << "  Relative error: " << (error * 100.0f) << "%\n\n";

    // CRITICAL QUALITY GATE: Within factor of 2 (50% error)
    bool ratio_valid = error < 0.5f;

    std::cout << "CRITICAL QUALITY GATE:\n";
    std::cout << "  Mass ratio within factor 2 (error < 50%): "
              << (ratio_valid ? "PASS ✓✓✓" : "FAIL ✗✗✗") << "\n";

    if (!ratio_valid) {
        std::cout << "\n  ⚠️  WARNING: Mass ratio does not match particle physics!\n";
        std::cout << "  This suggests:\n";
        std::cout << "    1. Vortex excitation model needs refinement\n";
        std::cout << "    2. Additional physics (radial modes, coupling) required\n";
        std::cout << "    3. Topological stability analysis needed\n";
    } else {
        std::cout << "\n  🎉 SUCCESS: TRD predicts particle mass spectrum!\n";
        std::cout << "  Topological excitations → Particle masses confirmed!\n";
    }

    return ratio_valid;
}

/**
 * Test 4: Triple vortex (Q=3) → Third mass m₃ (tau lepton?)
 */
bool testTripleVortexMass() {
    std::cout << "\n=== Test 4: Triple Vortex (Q=3) - Third Mass (Tau?) ===\n";

    const int N = 32;
    const float dx = 0.5f;
    const float K = 1.0f;

    KuramotoGrid3D grid(N, dx, K);
    initTripleVortex(grid, 3.0f);  // Radius = 3.0

    float energy = computeVortexEnergy(grid);
    float charge = computeTopologicalCharge(grid);
    float mass = energy / (C_LIGHT * C_LIGHT);

    std::cout << "  Topological charge Q: " << charge << "\n";
    std::cout << "  Vortex energy E₃: " << energy << "\n";
    std::cout << "  Effective mass m₃: " << mass << "\n";

    // Quality gate: Q should be close to 3
    bool charge_valid = std::abs(charge - 3.0f) < 0.6f;

    std::cout << "  Topological charge (≈3.0): " << (charge_valid ? "PASS ✓" : "FAIL ✗") << "\n";

    // Compare to tau lepton mass ratio
    // m_tau / m_electron ≈ 3477
    std::cout << "\n  Note: For tau lepton, expect m₃/m₁ ≈ 3477\n";
    std::cout << "  (Testing requires more sophisticated vortex configurations)\n";

    return charge_valid;
}

/**
 * Advanced Analysis: Energy scaling with topological charge
 */
void analyzeEnergyScaling() {
    std::cout << "\n=== Advanced Analysis: Energy Scaling E(Q) ===\n";

    const int N = 32;
    const float dx = 0.5f;
    const float K = 1.0f;

    std::vector<float> energies;
    std::vector<int> charges = {1, 2, 3};

    // Q=1
    KuramotoGrid3D grid1(N, dx, K);
    initSingleVortex(grid1, 0.0f, 0.0f);
    energies.push_back(computeVortexEnergy(grid1));

    // Q=2
    KuramotoGrid3D grid2(N, dx, K);
    initDoubleVortex(grid2, 4.0f);
    energies.push_back(computeVortexEnergy(grid2));

    // Q=3
    KuramotoGrid3D grid3(N, dx, K);
    initTripleVortex(grid3, 3.0f);
    energies.push_back(computeVortexEnergy(grid3));

    std::cout << "\n  Q    Energy       E/E₁\n";
    std::cout << "  ------------------------\n";
    for (size_t i = 0; i < charges.size(); ++i) {
        float ratio = energies[i] / energies[0];
        std::cout << "  " << charges[i] << "    "
                  << std::setw(10) << energies[i] << "  "
                  << std::setw(6) << std::setprecision(3) << ratio << "\n";
    }

    std::cout << "\n  Expected: E(Q) ~ Q·E₁ (linear scaling)\n";
    std::cout << "  Observed: Check if E₂/E₁ ≈ 2, E₃/E₁ ≈ 3\n";
}

/**
 * Main test runner
 */
int runParticleSpectrum3DTest() {
    std::cout << "========================================\n";
    std::cout << "  B1 PARTICLE SPECTRUM DERIVATION\n";
    std::cout << "========================================\n";
    std::cout << "Theory: Particles are topological vortex excitations\n";
    std::cout << "Model: m_eff = E_vortex / c²\n";
    std::cout << "Critical Gate: m₂/m₁ ≈ 206.768 (muon/electron)\n\n";

    bool all_pass = true;

    all_pass &= testSingleVortexMass();
    all_pass &= testDoubleVortexMass();
    all_pass &= testMassRatio();  // CRITICAL TEST
    all_pass &= testTripleVortexMass();

    analyzeEnergyScaling();

    std::cout << "\n========================================\n";
    std::cout << "FINAL VERDICT: " << (all_pass ? "PASS ✓✓✓" : "FAIL ✗✗✗") << "\n";
    std::cout << "========================================\n";

    if (!all_pass) {
        std::cout << "\n⚠️  CRITICAL: TRD particle spectrum derivation failed!\n";
        std::cout << "Theory requires refinement or additional physics.\n";
    } else {
        std::cout << "\n🎉 BREAKTHROUGH: TRD successfully predicts particle masses!\n";
        std::cout << "Topological vortex excitations → Particle spectrum confirmed!\n";
    }

    return all_pass ? 0 : 1;
}
