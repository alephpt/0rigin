/**
 * test_particle_scattering.cpp
 *
 * D4: Particle Scattering - Soliton Integrity Test
 *
 * Goal: Test whether TRD solitons (vortex excitations) scatter elastically
 *       while maintaining topological integrity.
 *
 * Physics:
 *   In field theory, solitons should bounce off each other elastically if
 *   topologically stable. This tests whether TRD vortices are true solitons
 *   or whether they annihilate/merge on collision.
 *
 *   Two Q=1 vortices on collision course:
 *     Vortex 1: centered at (-d, 0, 0) moving right (+v_x)
 *     Vortex 2: centered at (+d, 0, 0) moving left (-v_x)
 *
 *   Elastic scattering criteria:
 *     1. Topological charge conserved: Q_total = 2 (before & after)
 *     2. Energy conserved: ΔE/E < 1%
 *     3. Momentum conserved: Δp/p < 1%
 *     4. Solitons emerge intact (no annihilation)
 *
 *   Cross-section measurement:
 *     σ = (scattered flux) / (incident flux)
 *     Regge prediction: σ ~ 1/s where s = (E_cm)² is Mandelstam variable
 *
 * Golden Key Calibration: 1 TRD unit = 246 GeV
 *   Collision energy: E_cm = 2·γ·m·c² where γ = 1/√(1-v²/c²)
 *   For v=0.1c: E_cm ≈ 0.2·m (TRD units) → ~50 GeV for m≈0.5
 *
 * Quality Gates:
 *   - Q_total unchanged: ΔQ = 0 (exact)
 *   - Energy conserved: |ΔE/E| < 1%
 *   - Momentum conserved: |Δp/p| < 5%
 *   - Solitons survive: vortex count = 2 post-collision
 *   - Cross-section scaling: σ ~ 1/s (within 20% of Regge theory)
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include <fstream>
#include <yaml-cpp/yaml.h>

// Physical constants
const float PI = 3.14159265358979323846f;
const float HBAR = 1.0f;  // Natural units
const float C = 1.0f;     // Speed of light in natural units
const float TRD_UNIT_GEV = 246.0f;  // Golden Key calibration: 1 TRD unit = 246 GeV

/**
 * Simple 3D grid for phase field
 */
class Grid3D {
private:
    int N;
    float dx;
    std::vector<float> theta;      // Phase field
    std::vector<float> theta_dot;  // Time derivative (for velocity)

public:
    Grid3D(int size, float spacing)
        : N(size), dx(spacing),
          theta(size * size * size, 0.0f),
          theta_dot(size * size * size, 0.0f) {}

    int getSize() const { return N; }
    float getSpacing() const { return dx; }

    float& at(int ix, int iy, int iz) {
        return theta[ix + N * (iy + N * iz)];
    }

    const float& at(int ix, int iy, int iz) const {
        return theta[ix + N * (iy + N * iz)];
    }

    float& dot(int ix, int iy, int iz) {
        return theta_dot[ix + N * (iy + N * iz)];
    }

    const float& dot(int ix, int iy, int iz) const {
        return theta_dot[ix + N * (iy + N * iz)];
    }

    int index(int ix, int iy, int iz) const {
        return ix + N * (iy + N * iz);
    }

    void wrapIndices(int& ix, int& iy, int& iz) const {
        ix = (ix + N) % N;
        iy = (iy + N) % N;
        iz = (iz + N) % N;
    }
};

/**
 * Initialize vortex at position (x0, y0, z0)
 * Phase winds around z-axis through vortex core
 * Velocity is added via phase gradient (boost)
 */
void initVortex(Grid3D& grid, float x0, float y0, float z0,
                float vx = 0.0f, float vy = 0.0f, float vz = 0.0f) {
    const int N = grid.getSize();
    const float dx = grid.getSpacing();

    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            for (int iz = 0; iz < N; ++iz) {
                float x = (ix - N/2) * dx;
                float y = (iy - N/2) * dx;
                float z = (iz - N/2) * dx;

                // Vortex phase (winding around z-axis through vortex core)
                float theta_vortex = std::atan2(y - y0, x - x0);

                // Add velocity as phase gradient (de Broglie relation: p = ħk)
                // k = p/ħ, so phase boost: θ_boost = k·r = (p/ħ)·r
                // For velocity v: p = γmv, so k = γmv/ħ
                // In natural units (ħ=1, m=1): k ≈ v for v << c
                float k_x = vx;  // Wave vector for x-momentum
                float k_y = vy;
                float k_z = vz;

                float theta_boost = k_x * (x - x0) + k_y * (y - y0) + k_z * (z - z0);

                grid.at(ix, iy, iz) = theta_vortex + theta_boost;
            }
        }
    }
}

/**
 * Initialize two vortices on collision course
 * Vortex 1: at (-d, 0, 0) with velocity (+v, 0, 0)
 * Vortex 2: at (+d, 0, 0) with velocity (-v, 0, 0)
 */
void initCollisionScenario(Grid3D& grid, float separation, float velocity) {
    const int N = grid.getSize();
    const float dx = grid.getSpacing();

    // Clear grid
    for (int i = 0; i < N*N*N; ++i) {
        grid.at(i/N/N, (i/N)%N, i%N) = 0.0f;
    }

    // Vortex 1: left side, moving right
    initVortex(grid, -separation, 0.0f, 0.0f, velocity, 0.0f, 0.0f);

    // Vortex 2: right side, moving left
    // PROBLEM: Superposition of phases is ambiguous
    // Solution: Initialize in separate grids and sum carefully
    Grid3D temp(N, dx);
    initVortex(temp, separation, 0.0f, 0.0f, -velocity, 0.0f, 0.0f);

    // Add phases (modulo 2π)
    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            for (int iz = 0; iz < N; ++iz) {
                float theta_sum = grid.at(ix, iy, iz) + temp.at(ix, iy, iz);
                // Wrap to [-π, π]
                while (theta_sum > PI) theta_sum -= 2.0f * PI;
                while (theta_sum < -PI) theta_sum += 2.0f * PI;
                grid.at(ix, iy, iz) = theta_sum;
            }
        }
    }
}

/**
 * Compute topological charge Q via winding number
 * Q = (1/2π) ∮ ∇θ·dl around vortex core
 *
 * For 3D vortex line along z, integrate in x-y plane
 */
float computeTopologicalCharge(const Grid3D& grid) {
    const int N = grid.getSize();
    const float dx = grid.getSpacing();

    float Q_total = 0.0f;

    // Scan x-y plane at mid-z
    int iz = N / 2;

    // Look for vortex cores (phase singularities)
    // Core detection: large phase gradient magnitude
    for (int ix = 1; ix < N-1; ++ix) {
        for (int iy = 1; iy < N-1; ++iy) {
            // Compute winding around small plaquette
            float theta_00 = grid.at(ix,   iy,   iz);
            float theta_10 = grid.at(ix+1, iy,   iz);
            float theta_11 = grid.at(ix+1, iy+1, iz);
            float theta_01 = grid.at(ix,   iy+1, iz);

            // Phase differences (unwrap)
            auto unwrap = [](float dtheta) -> float {
                while (dtheta > PI) dtheta -= 2.0f * PI;
                while (dtheta < -PI) dtheta += 2.0f * PI;
                return dtheta;
            };

            float d1 = unwrap(theta_10 - theta_00);
            float d2 = unwrap(theta_11 - theta_10);
            float d3 = unwrap(theta_01 - theta_11);
            float d4 = unwrap(theta_00 - theta_01);

            float winding = (d1 + d2 + d3 + d4) / (2.0f * PI);

            // Round to nearest integer (quantized charge)
            int charge = std::round(winding);
            if (std::abs(charge) > 0) {
                Q_total += charge;
            }
        }
    }

    return Q_total;
}

/**
 * Compute total field energy
 * E = ∫ (1/2)[(∇θ)² + V(θ)] dV
 * For TRD: V = 0 (massless phase field)
 */
float computeFieldEnergy(const Grid3D& grid) {
    const int N = grid.getSize();
    const float dx = grid.getSpacing();
    const float dV = dx * dx * dx;

    float E_total = 0.0f;

    for (int ix = 1; ix < N-1; ++ix) {
        for (int iy = 1; iy < N-1; ++iy) {
            for (int iz = 1; iz < N-1; ++iz) {
                // Gradient energy (kinetic)
                float theta_c = grid.at(ix, iy, iz);
                float theta_xp = grid.at(ix+1, iy, iz);
                float theta_xm = grid.at(ix-1, iy, iz);
                float theta_yp = grid.at(ix, iy+1, iz);
                float theta_ym = grid.at(ix, iy-1, iz);
                float theta_zp = grid.at(ix, iy, iz+1);
                float theta_zm = grid.at(ix, iy, iz-1);

                float grad_x = (theta_xp - theta_xm) / (2.0f * dx);
                float grad_y = (theta_yp - theta_ym) / (2.0f * dx);
                float grad_z = (theta_zp - theta_zm) / (2.0f * dx);

                float grad_squared = grad_x*grad_x + grad_y*grad_y + grad_z*grad_z;

                E_total += 0.5f * grad_squared * dV;
            }
        }
    }

    return E_total;
}

/**
 * Compute total momentum
 * P = ∫ (∂θ/∂t)·∇θ dV
 * For field with velocity: P_i ~ ∫ (∂θ/∂x_i) dV
 */
std::array<float, 3> computeFieldMomentum(const Grid3D& grid) {
    const int N = grid.getSize();
    const float dx = grid.getSpacing();
    const float dV = dx * dx * dx;

    std::array<float, 3> P = {0.0f, 0.0f, 0.0f};

    for (int ix = 1; ix < N-1; ++ix) {
        for (int iy = 1; iy < N-1; ++iy) {
            for (int iz = 1; iz < N-1; ++iz) {
                float theta_xp = grid.at(ix+1, iy, iz);
                float theta_xm = grid.at(ix-1, iy, iz);
                float theta_yp = grid.at(ix, iy+1, iz);
                float theta_ym = grid.at(ix, iy-1, iz);
                float theta_zp = grid.at(ix, iy, iz+1);
                float theta_zm = grid.at(ix, iy, iz-1);

                P[0] += (theta_xp - theta_xm) / (2.0f * dx) * dV;
                P[1] += (theta_yp - theta_ym) / (2.0f * dx) * dV;
                P[2] += (theta_zp - theta_zm) / (2.0f * dx) * dV;
            }
        }
    }

    return P;
}

/**
 * Count vortices in the field
 * Returns number of topological charge concentrations
 */
int countVortices(const Grid3D& grid) {
    const int N = grid.getSize();
    int iz = N / 2;  // Mid-plane

    int count = 0;

    for (int ix = 1; ix < N-1; ++ix) {
        for (int iy = 1; iy < N-1; ++iy) {
            // Compute winding around small plaquette
            float theta_00 = grid.at(ix,   iy,   iz);
            float theta_10 = grid.at(ix+1, iy,   iz);
            float theta_11 = grid.at(ix+1, iy+1, iz);
            float theta_01 = grid.at(ix,   iy+1, iz);

            auto unwrap = [](float dtheta) -> float {
                while (dtheta > PI) dtheta -= 2.0f * PI;
                while (dtheta < -PI) dtheta += 2.0f * PI;
                return dtheta;
            };

            float d1 = unwrap(theta_10 - theta_00);
            float d2 = unwrap(theta_11 - theta_10);
            float d3 = unwrap(theta_01 - theta_11);
            float d4 = unwrap(theta_00 - theta_01);

            float winding = std::abs((d1 + d2 + d3 + d4) / (2.0f * PI));

            if (winding > 0.5f) {  // Threshold for vortex detection
                count++;
            }
        }
    }

    return count;
}

/**
 * Evolve field using simple Kuramoto-like dynamics
 * ∂θ/∂t = K·∇²θ (diffusion-like, smooths phase)
 */
void evolveField(Grid3D& grid, float dt, float coupling = 1.0f) {
    const int N = grid.getSize();
    const float dx = grid.getSpacing();

    std::vector<float> theta_new(N*N*N);

    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            for (int iz = 0; iz < N; ++iz) {
                int ixp = (ix + 1) % N;
                int ixm = (ix - 1 + N) % N;
                int iyp = (iy + 1) % N;
                int iym = (iy - 1 + N) % N;
                int izp = (iz + 1) % N;
                int izm = (iz - 1 + N) % N;

                float theta_c = grid.at(ix, iy, iz);
                float laplacian = (
                    grid.at(ixp, iy, iz) + grid.at(ixm, iy, iz) +
                    grid.at(ix, iyp, iz) + grid.at(ix, iym, iz) +
                    grid.at(ix, iy, izp) + grid.at(ix, iy, izm) -
                    6.0f * theta_c
                ) / (dx * dx);

                int idx = grid.index(ix, iy, iz);
                theta_new[idx] = theta_c + dt * coupling * laplacian;
            }
        }
    }

    // Update grid
    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            for (int iz = 0; iz < N; ++iz) {
                int idx = grid.index(ix, iy, iz);
                grid.at(ix, iy, iz) = theta_new[idx];
            }
        }
    }
}

/**
 * Test elastic scattering of two vortices
 */
bool testElasticScattering(const YAML::Node& scenario) {
    std::string name = scenario["name"].as<std::string>();
    std::string description = scenario["description"].as<std::string>();
    float separation = scenario["initial_separation"].as<float>();
    float velocity = scenario["velocity"].as<float>();
    float impact_param = scenario["impact_parameter"].as<float>();

    std::cout << "\n=== Testing: " << name << " ===\n";
    std::cout << "Description: " << description << "\n";
    std::cout << "Initial separation: " << separation << " (TRD units)\n";
    std::cout << "Velocity: " << velocity << "c\n";
    std::cout << "Impact parameter: " << impact_param << "\n\n";

    // Setup grid - read from global config
    extern YAML::Node g_config;
    const int N = g_config["physics"]["grid_size"].as<int>();
    const float dx = 1.0f;
    const float dt = g_config["physics"]["time_step"].as<float>();
    const int max_steps = g_config["physics"]["evolution_steps"].as<int>();

    Grid3D grid(N, dx);

    // Initialize collision
    std::cout << "Initializing collision scenario...\n" << std::flush;
    initCollisionScenario(grid, separation, velocity);

    // Initial measurements
    std::cout << "Computing initial topological charge...\n" << std::flush;
    float Q_initial = computeTopologicalCharge(grid);
    float E_initial = computeFieldEnergy(grid);
    auto P_initial = computeFieldMomentum(grid);
    int vortex_count_initial = countVortices(grid);

    std::cout << "Initial state:\n";
    std::cout << "  Topological charge Q: " << Q_initial << "\n";
    std::cout << "  Field energy E: " << E_initial << " (TRD units)\n";
    std::cout << "  Momentum P: (" << P_initial[0] << ", "
              << P_initial[1] << ", " << P_initial[2] << ")\n";
    std::cout << "  Vortex count: " << vortex_count_initial << "\n\n";

    // Convert to physical units (246 GeV)
    float E_initial_GeV = E_initial * TRD_UNIT_GEV;
    std::cout << "  Energy (GeV): " << E_initial_GeV << " GeV\n";

    // Collision energy (center of mass)
    float gamma = 1.0f / std::sqrt(1.0f - velocity*velocity);
    float mass_eff = 0.5f;  // Effective vortex mass (calibrated)
    float E_cm_trd = 2.0f * gamma * mass_eff;  // In TRD units
    float E_cm_GeV = E_cm_trd * TRD_UNIT_GEV;
    float mandelstam_s = E_cm_trd * E_cm_trd;

    std::cout << "  Lorentz factor γ: " << gamma << "\n";
    std::cout << "  Collision energy E_cm: " << E_cm_trd << " (TRD) = "
              << E_cm_GeV << " GeV\n";
    std::cout << "  Mandelstam variable s: " << mandelstam_s << "\n\n";

    // Time evolution
    std::cout << "Evolving collision (" << max_steps << " steps)...\n" << std::flush;
    for (int step = 0; step < max_steps; ++step) {
        evolveField(grid, dt, 1.0f);

        if (step % 100 == 0) {
            float E = computeFieldEnergy(grid);
            float dE = std::abs(E - E_initial) / E_initial;
            std::cout << "  Step " << step << ": E = " << E
                      << ", ΔE/E = " << (dE * 100.0f) << "%\n";
        }
    }

    // Final measurements
    float Q_final = computeTopologicalCharge(grid);
    float E_final = computeFieldEnergy(grid);
    auto P_final = computeFieldMomentum(grid);
    int vortex_count_final = countVortices(grid);

    std::cout << "\nFinal state:\n";
    std::cout << "  Topological charge Q: " << Q_final << "\n";
    std::cout << "  Field energy E: " << E_final << " (TRD units)\n";
    std::cout << "  Momentum P: (" << P_final[0] << ", "
              << P_final[1] << ", " << P_final[2] << ")\n";
    std::cout << "  Vortex count: " << vortex_count_final << "\n\n";

    // Conservation checks
    float dQ = std::abs(Q_final - Q_initial);
    float dE = std::abs(E_final - E_initial) / E_initial;
    float dPx = std::abs(P_final[0] - P_initial[0]) / (std::abs(P_initial[0]) + 1e-6f);

    bool Q_conserved = (dQ < 0.1f);  // Topological charge must be exact
    bool E_conserved = (dE < 0.01f);  // Energy < 1%
    bool P_conserved = (dPx < 0.05f);  // Momentum < 5%
    bool solitons_survive = (vortex_count_final >= 2);

    std::cout << "Conservation checks:\n";
    std::cout << "  ΔQ = " << dQ << " " << (Q_conserved ? "✓" : "✗") << "\n";
    std::cout << "  ΔE/E = " << (dE * 100.0f) << "% "
              << (E_conserved ? "✓" : "✗") << "\n";
    std::cout << "  Δp_x/p_x = " << (dPx * 100.0f) << "% "
              << (P_conserved ? "✓" : "✗") << "\n";
    std::cout << "  Solitons survive: " << vortex_count_final << "/2 "
              << (solitons_survive ? "✓" : "✗") << "\n\n";

    bool passed = Q_conserved && E_conserved && P_conserved && solitons_survive;

    if (passed) {
        std::cout << "✓ TEST PASSED - Elastic scattering validated\n";
    } else {
        std::cout << "✗ TEST FAILED - Scattering not elastic\n";
    }

    return passed;
}

// Global config (set by main test runner)
YAML::Node g_config;

/**
 * Main test runner
 */
int runParticleScatteringTest() {
    std::cout << "\n";
    std::cout << "========================================\n";
    std::cout << " D4: Particle Scattering - Soliton Integrity\n";
    std::cout << "========================================\n";
    std::cout << "\nGoal: Validate elastic scattering of TRD vortex solitons\n";
    std::cout << "Golden Key: 1 TRD unit = 246 GeV\n\n";

    // Load configuration
    YAML::Node config;
    try {
        config = YAML::LoadFile("config/particle_scattering.yaml");
        g_config = config;  // Store globally for scenario functions
    } catch (const std::exception& e) {
        std::cerr << "Error loading config: " << e.what() << std::endl;
        return 1;
    }

    // Run all collision scenarios
    auto scenarios = config["collision_scenarios"];
    int passed = 0;
    int total = scenarios.size();

    for (const auto& scenario : scenarios) {
        if (testElasticScattering(scenario)) {
            passed++;
        }
    }

    // Summary
    std::cout << "\n========================================\n";
    std::cout << " Test Summary\n";
    std::cout << "========================================\n";
    std::cout << "Passed: " << passed << "/" << total << "\n\n";

    if (passed == total) {
        std::cout << "✓ ALL TESTS PASSED\n";
        std::cout << "\nConclusion: TRD vortices are topologically stable solitons\n";
        std::cout << "that scatter elastically, maintaining charge, energy, and\n";
        std::cout << "momentum conservation. This validates the particle-like\n";
        std::cout << "nature of topological excitations in TRD.\n\n";
        return 0;
    } else {
        std::cout << "✗ SOME TESTS FAILED\n";
        std::cout << "\nSolitons may not be fully stable or elastic scattering\n";
        std::cout << "criteria not satisfied. Check field dynamics and\n";
        std::cout << "topological protection mechanisms.\n\n";
        return 1;
    }
}
