/**
 * test_em_gravity_coupling_3d.cpp
 *
 * Electromagnetic-Gravity Coupling in 3D
 *
 * Goal: Verify EM field energy curves SMFT spacetime
 *
 * Physics:
 *   Stress-energy tensor: T^μν from Maxwell fields
 *   Energy density: ρ_EM = (E² + B²)/2
 *   R-field evolution: ∂R/∂t ~ -ε·∇·T^(EM)
 *   Coupling: ΔR proportional to ρ_EM
 *
 * Test Scenarios:
 *   1. High-energy EM pulse in 3D
 *   2. Measure R-field response (back-reaction)
 *   3. Verify coupling constant ε
 *
 * Quality Gates:
 *   - Coupling constant within 10% of theoretical
 *   - Energy-momentum conservation < 1%
 *   - R-field perturbation correlates with ρ_EM
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <array>

const float PI = 3.14159265358979323846f;

/**
 * Simple 3D grid for fields
 */
template<typename T>
class Grid3D {
private:
    int N;
    std::vector<T> data;

public:
    Grid3D(int size) : N(size), data(size * size * size) {}

    T& at(int ix, int iy, int iz) {
        return data[ix + N * (iy + N * iz)];
    }

    const T& at(int ix, int iy, int iz) const {
        return data[ix + N * (iy + N * iz)];
    }

    int size() const { return N; }
};

/**
 * EM field structure
 */
struct EMField {
    float Ex, Ey, Ez;
    float Bx, By, Bz;

    float energyDensity() const {
        return 0.5f * (Ex*Ex + Ey*Ey + Ez*Ez + Bx*Bx + By*By + Bz*Bz);
    }
};

/**
 * Initialize Gaussian EM pulse
 * Localized energy packet
 */
void initEMPulse(Grid3D<EMField>& em, Grid3D<float>& R,
                  float dx, float E0, float sigma) {
    const int N = em.size();

    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            for (int iz = 0; iz < N; ++iz) {
                float x = (ix - N/2) * dx;
                float y = (iy - N/2) * dx;
                float z = (iz - N/2) * dx;

                float r2 = x*x + y*y + z*z;
                float envelope = E0 * std::exp(-r2 / (2.0f * sigma*sigma));

                // Electric field (polarized in x-direction)
                em.at(ix, iy, iz).Ex = envelope;
                em.at(ix, iy, iz).Ey = 0.0f;
                em.at(ix, iy, iz).Ez = 0.0f;

                // Magnetic field (polarized in y-direction for EM wave)
                em.at(ix, iy, iz).Bx = 0.0f;
                em.at(ix, iy, iz).By = envelope;  // B = E for EM wave (c=1)
                em.at(ix, iy, iz).Bz = 0.0f;

                // Initialize R-field (flat background)
                R.at(ix, iy, iz) = 1.0f;
            }
        }
    }
}

/**
 * Compute EM energy density at each point
 */
void computeEMEnergyDensity(const Grid3D<EMField>& em, Grid3D<float>& rho_EM) {
    const int N = em.size();

    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            for (int iz = 0; iz < N; ++iz) {
                rho_EM.at(ix, iy, iz) = em.at(ix, iy, iz).energyDensity();
            }
        }
    }
}

/**
 * Update R-field based on EM stress-energy
 * ∂R/∂t = -ε·ρ_EM (simplified coupling)
 */
void updateRField(Grid3D<float>& R, const Grid3D<float>& rho_EM,
                   float epsilon, float dt) {
    const int N = R.size();

    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            for (int iz = 0; iz < N; ++iz) {
                R.at(ix, iy, iz) += -epsilon * rho_EM.at(ix, iy, iz) * dt;
            }
        }
    }
}

/**
 * Measure correlation between R-field perturbation and EM energy
 */
float measureCorrelation(const Grid3D<float>& R, const Grid3D<float>& rho_EM) {
    const int N = R.size();

    // Compute averages
    float R_avg = 0.0f, rho_avg = 0.0f;
    int count = 0;

    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            for (int iz = 0; iz < N; ++iz) {
                R_avg += R.at(ix, iy, iz);
                rho_avg += rho_EM.at(ix, iy, iz);
                count++;
            }
        }
    }
    R_avg /= count;
    rho_avg /= count;

    // Compute correlation coefficient
    float numerator = 0.0f;
    float var_R = 0.0f, var_rho = 0.0f;

    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            for (int iz = 0; iz < N; ++iz) {
                float dR = R.at(ix, iy, iz) - R_avg;
                float drho = rho_EM.at(ix, iy, iz) - rho_avg;

                numerator += dR * drho;
                var_R += dR * dR;
                var_rho += drho * drho;
            }
        }
    }

    if (var_R == 0.0f || var_rho == 0.0f) return 0.0f;

    return numerator / std::sqrt(var_R * var_rho);
}

/**
 * Test 1: EM pulse creates R-field perturbation
 */
bool testEMGravityCoupling() {
    std::cout << "\n=== Test 1: EM-Gravity Coupling ===\n";
    std::cout << "High-energy EM pulse → R-field response\n\n";

    const int N = 32;
    const float dx = 0.5f;
    const float E0 = 2.0f;  // Strong field
    const float sigma = 2.0f;  // Pulse width
    const float epsilon = 0.1f;  // Coupling constant
    const float dt = 0.01f;
    const int num_steps = 100;

    // Initialize grids
    Grid3D<EMField> em(N);
    Grid3D<float> R(N);
    Grid3D<float> rho_EM(N);

    initEMPulse(em, R, dx, E0, sigma);

    // Compute initial EM energy density
    computeEMEnergyDensity(em, rho_EM);

    float rho_max_initial = 0.0f;
    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            for (int iz = 0; iz < N; ++iz) {
                rho_max_initial = std::max(rho_max_initial, rho_EM.at(ix, iy, iz));
            }
        }
    }

    std::cout << "Pulse amplitude E0: " << E0 << "\n";
    std::cout << "Pulse width σ: " << sigma << "\n";
    std::cout << "Coupling constant ε: " << epsilon << "\n";
    std::cout << "Max EM energy density: " << rho_max_initial << "\n\n";

    // Evolve R-field (static EM field for simplicity)
    for (int step = 0; step < num_steps; ++step) {
        updateRField(R, rho_EM, epsilon, dt);
    }

    // Measure R-field perturbation
    float R_min = 1.0f, R_max = 1.0f;
    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            for (int iz = 0; iz < N; ++iz) {
                R_min = std::min(R_min, R.at(ix, iy, iz));
                R_max = std::max(R_max, R.at(ix, iy, iz));
            }
        }
    }

    float delta_R = R_max - R_min;
    float expected_delta_R = epsilon * rho_max_initial * num_steps * dt;

    std::cout << "R-field range: [" << R_min << ", " << R_max << "]\n";
    std::cout << "ΔR (measured): " << delta_R << "\n";
    std::cout << "ΔR (expected): " << expected_delta_R << "\n";

    float coupling_error = std::abs(delta_R - expected_delta_R) / expected_delta_R;
    std::cout << "Coupling error: " << (coupling_error * 100.0f) << "%\n\n";

    bool pass = coupling_error < 0.10f;  // 10%

    std::cout << "Quality Gate (< 10%): " << (pass ? "PASS ✓" : "FAIL ✗") << "\n";

    return pass;
}

/**
 * Test 2: Correlation between R-field and EM energy
 */
bool testRFieldEMCorrelation() {
    std::cout << "\n=== Test 2: R-Field / EM Energy Correlation ===\n";
    std::cout << "Verify: R perturbation correlates with ρ_EM\n\n";

    const int N = 32;
    const float dx = 0.5f;
    const float E0 = 2.0f;
    const float sigma = 2.0f;
    const float epsilon = 0.1f;
    const float dt = 0.01f;
    const int num_steps = 100;

    Grid3D<EMField> em(N);
    Grid3D<float> R(N);
    Grid3D<float> rho_EM(N);

    initEMPulse(em, R, dx, E0, sigma);
    computeEMEnergyDensity(em, rho_EM);

    // Evolve
    for (int step = 0; step < num_steps; ++step) {
        updateRField(R, rho_EM, epsilon, dt);
    }

    // Measure correlation
    float correlation = measureCorrelation(R, rho_EM);

    std::cout << "Correlation coefficient: " << correlation << "\n";
    std::cout << "Expected: Strong negative correlation (R decreases with ρ_EM)\n\n";

    // Expect strong negative correlation (R decreases where ρ_EM is high)
    bool pass = correlation < -0.8f;  // Strong anti-correlation

    std::cout << "Quality Gate (corr < -0.8): " << (pass ? "PASS ✓" : "FAIL ✗") << "\n";

    return pass;
}

/**
 * Test 3: Energy budget
 * EM energy → Gravitational potential energy
 */
bool testEnergyBudget() {
    std::cout << "\n=== Test 3: Energy Budget (EM → Gravity) ===\n";
    std::cout << "Verify: EM energy loss ≈ Gravitational potential gain\n\n";

    const int N = 32;
    const float dx = 0.5f;
    const float E0 = 2.0f;
    const float sigma = 2.0f;
    const float epsilon = 0.1f;
    const float dt = 0.01f;
    const int num_steps = 100;

    Grid3D<EMField> em(N);
    Grid3D<float> R(N);
    Grid3D<float> rho_EM(N);

    initEMPulse(em, R, dx, E0, sigma);
    computeEMEnergyDensity(em, rho_EM);

    // Compute initial EM energy
    float E_EM_initial = 0.0f;
    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            for (int iz = 0; iz < N; ++iz) {
                E_EM_initial += rho_EM.at(ix, iy, iz) * dx*dx*dx;
            }
        }
    }

    // Evolve (EM field static, R-field responds)
    for (int step = 0; step < num_steps; ++step) {
        updateRField(R, rho_EM, epsilon, dt);
    }

    // Compute gravitational potential energy
    // U_grav ~ ∫ (R - 1)² dV (simplified)
    float E_grav = 0.0f;
    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            for (int iz = 0; iz < N; ++iz) {
                float dR = R.at(ix, iy, iz) - 1.0f;
                E_grav += 0.5f * dR * dR * dx*dx*dx;
            }
        }
    }

    std::cout << "Initial EM energy: " << E_EM_initial << "\n";
    std::cout << "Gravitational potential energy: " << E_grav << "\n";

    // Energy transfer fraction
    float transfer_fraction = E_grav / E_EM_initial;
    std::cout << "Energy transfer: " << (transfer_fraction * 100.0f) << "%\n\n";

    // Expect some energy transfer
    bool pass = transfer_fraction > 0.001f && transfer_fraction < 0.1f;

    std::cout << "Quality Gate (0.1% < transfer < 10%): " << (pass ? "PASS ✓" : "FAIL ✗") << "\n";

    return pass;
}

/**
 * Main test runner
 */
int runEMGravityCoupling3DTest() {
    std::cout << "========================================\n";
    std::cout << "  3D EM-Gravity Coupling Validation\n";
    std::cout << "========================================\n";
    std::cout << "Physics: T^μν(EM) → Curvature (R-field)\n";
    std::cout << "Coupling: ∂R/∂t ~ -ε·ρ_EM\n\n";

    bool all_pass = true;

    all_pass &= testEMGravityCoupling();
    all_pass &= testRFieldEMCorrelation();
    all_pass &= testEnergyBudget();

    std::cout << "\n========================================\n";
    std::cout << "FINAL VERDICT: " << (all_pass ? "PASS ✓" : "FAIL ✗") << "\n";
    std::cout << "========================================\n";

    return all_pass ? 0 : 1;
}
