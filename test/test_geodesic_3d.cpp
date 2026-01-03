/**
 * test_geodesic_3d.cpp
 *
 * Geodesic Equation Validation in 3D
 *
 * Goal: Verify particles follow geodesics in curved TRD spacetime
 *
 * Physics:
 *   Metric: g_μν = R²(x,y,z)·η_μν (conformal to Minkowski)
 *   Christoffel symbols: Γ^μ_αβ = ∂_α ln(R)·δ^μ_β + ∂_β ln(R)·δ^μ_α - η^μν·η_αβ·∂_ν ln(R)
 *   Geodesic equation: d²x^μ/dτ² = -Γ^μ_αβ (dx^α/dτ)(dx^β/dτ)
 *
 * Test Scenarios:
 *   1. Flat space (R=1) → Verify straight-line motion
 *   2. Gaussian R-peak → Verify curved trajectory (deflection)
 *   3. Two-body system → Geodesic deviation (tidal forces)
 *
 * Integration: RK4 for geodesic equation
 *
 * Quality Gates:
 *   - Trajectory deviation < 1% from analytical (where available)
 *   - Energy along geodesic conserved < 0.1%
 *   - Numerical stability over long timescales
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <array>
#include <functional>

const float PI = 3.14159265358979323846f;

/**
 * Geodesic particle structure
 */
struct GeodesicParticle {
    // Position (spatial)
    float x, y, z;

    // 4-velocity (u^μ = dx^μ/dτ)
    float u_t, u_x, u_y, u_z;

    // Proper time
    float tau;

    // Trajectory history
    std::vector<std::array<float, 8>> history;  // {τ, x, y, z, u_t, u_x, u_y, u_z}

    void recordState() {
        history.push_back({tau, x, y, z, u_t, u_x, u_y, u_z});
    }

    float getEnergy(float R) const {
        // E = -u_t (timelike Killing vector)
        return -u_t * R * R;  // Account for conformal factor
    }
};

/**
 * Compute R-field and its gradient at position (x, y, z)
 */
struct RFieldData {
    float R;
    float dR_dx, dR_dy, dR_dz;
};

/**
 * Flat R-field: R = 1 everywhere
 */
RFieldData flatRField(float x, float y, float z) {
    return {1.0f, 0.0f, 0.0f, 0.0f};
}

/**
 * Gaussian R-peak: R(x,y,z) = 1 + A·exp(-r²/2σ²)
 * Represents localized mass concentration
 */
RFieldData gaussianRPeak(float x, float y, float z, float A, float sigma) {
    float r2 = x*x + y*y + z*z;
    float exp_factor = std::exp(-r2 / (2.0f * sigma*sigma));

    float R = 1.0f + A * exp_factor;

    // Gradient: ∇R = -A·(r/σ²)·exp(-r²/2σ²)
    float grad_factor = -A * exp_factor / (sigma*sigma);
    float dR_dx = grad_factor * x;
    float dR_dy = grad_factor * y;
    float dR_dz = grad_factor * z;

    return {R, dR_dx, dR_dy, dR_dz};
}

/**
 * Compute Christoffel symbols for conformal metric g_μν = R²·η_μν
 *
 * Γ^μ_αβ = (∂_α ln R)δ^μ_β + (∂_β ln R)δ^μ_α - η^μν η_αβ (∂_ν ln R)
 */
void computeChristoffel(const RFieldData& data, float Gamma[4][4][4]) {
    // Clear all components
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 4; ++k)
                Gamma[i][j][k] = 0.0f;

    float R = data.R;
    float d_ln_R[4] = {0.0f, data.dR_dx/R, data.dR_dy/R, data.dR_dz/R};

    // Minkowski metric η_μν (signature -,+,+,+)
    float eta[4][4] = {
        {-1.0f, 0.0f, 0.0f, 0.0f},
        {0.0f, 1.0f, 0.0f, 0.0f},
        {0.0f, 0.0f, 1.0f, 0.0f},
        {0.0f, 0.0f, 0.0f, 1.0f}
    };

    // Compute Γ^μ_αβ
    for (int mu = 0; mu < 4; ++mu) {
        for (int alpha = 0; alpha < 4; ++alpha) {
            for (int beta = 0; beta < 4; ++beta) {
                // Term 1: (∂_α ln R)·δ^μ_β
                if (mu == beta) {
                    Gamma[mu][alpha][beta] += d_ln_R[alpha];
                }

                // Term 2: (∂_β ln R)·δ^μ_α
                if (mu == alpha) {
                    Gamma[mu][alpha][beta] += d_ln_R[beta];
                }

                // Term 3: -η^μν·η_αβ·(∂_ν ln R)
                for (int nu = 0; nu < 4; ++nu) {
                    Gamma[mu][alpha][beta] -= eta[mu][nu] * eta[alpha][beta] * d_ln_R[nu];
                }
            }
        }
    }
}

/**
 * RK4 integration step for geodesic equation
 *
 * dx^μ/dτ = u^μ
 * du^μ/dτ = -Γ^μ_αβ u^α u^β
 */
void geodesicRK4Step(GeodesicParticle& p,
                      std::function<RFieldData(float, float, float)> getRField,
                      float dtau) {
    // State vector: (x, y, z, u_x, u_y, u_z)
    // Note: u_t determined by normalization condition

    auto computeDerivatives = [&](float x, float y, float z,
                                   float ux, float uy, float uz) -> std::array<float, 6> {
        // Get R-field
        auto data = getRField(x, y, z);

        // Compute Christoffel symbols
        float Gamma[4][4][4];
        computeChristoffel(data, Gamma);

        // Reconstruct u_t from normalization: g_μν u^μ u^ν = -1 (timelike)
        // For conformal metric: R²(-u_t² + u_x² + u_y² + u_z²) = -1
        float u_spatial_sq = ux*ux + uy*uy + uz*uz;
        float ut = std::sqrt(1.0f / (data.R*data.R) + u_spatial_sq);

        // 4-velocity
        float u[4] = {ut, ux, uy, uz};

        // Compute accelerations: du^μ/dτ = -Γ^μ_αβ u^α u^β
        float du[4] = {0.0f, 0.0f, 0.0f, 0.0f};
        for (int mu = 0; mu < 4; ++mu) {
            for (int alpha = 0; alpha < 4; ++alpha) {
                for (int beta = 0; beta < 4; ++beta) {
                    du[mu] -= Gamma[mu][alpha][beta] * u[alpha] * u[beta];
                }
            }
        }

        // Return derivatives: (dx/dτ, dy/dτ, dz/dτ, du_x/dτ, du_y/dτ, du_z/dτ)
        return {ux, uy, uz, du[1], du[2], du[3]};
    };

    // RK4 integration
    auto k1 = computeDerivatives(p.x, p.y, p.z, p.u_x, p.u_y, p.u_z);

    auto k2 = computeDerivatives(
        p.x + 0.5f * dtau * k1[0],
        p.y + 0.5f * dtau * k1[1],
        p.z + 0.5f * dtau * k1[2],
        p.u_x + 0.5f * dtau * k1[3],
        p.u_y + 0.5f * dtau * k1[4],
        p.u_z + 0.5f * dtau * k1[5]
    );

    auto k3 = computeDerivatives(
        p.x + 0.5f * dtau * k2[0],
        p.y + 0.5f * dtau * k2[1],
        p.z + 0.5f * dtau * k2[2],
        p.u_x + 0.5f * dtau * k2[3],
        p.u_y + 0.5f * dtau * k2[4],
        p.u_z + 0.5f * dtau * k2[5]
    );

    auto k4 = computeDerivatives(
        p.x + dtau * k3[0],
        p.y + dtau * k3[1],
        p.z + dtau * k3[2],
        p.u_x + dtau * k3[3],
        p.u_y + dtau * k3[4],
        p.u_z + dtau * k3[5]
    );

    // Update state
    p.x += dtau * (k1[0] + 2.0f*k2[0] + 2.0f*k3[0] + k4[0]) / 6.0f;
    p.y += dtau * (k1[1] + 2.0f*k2[1] + 2.0f*k3[1] + k4[1]) / 6.0f;
    p.z += dtau * (k1[2] + 2.0f*k2[2] + 2.0f*k3[2] + k4[2]) / 6.0f;
    p.u_x += dtau * (k1[3] + 2.0f*k2[3] + 2.0f*k3[3] + k4[3]) / 6.0f;
    p.u_y += dtau * (k1[4] + 2.0f*k2[4] + 2.0f*k3[4] + k4[4]) / 6.0f;
    p.u_z += dtau * (k1[5] + 2.0f*k2[5] + 2.0f*k3[5] + k4[5]) / 6.0f;

    // Update u_t from normalization
    auto data = getRField(p.x, p.y, p.z);
    float u_spatial_sq = p.u_x*p.u_x + p.u_y*p.u_y + p.u_z*p.u_z;
    p.u_t = std::sqrt(1.0f / (data.R*data.R) + u_spatial_sq);

    p.tau += dtau;
}

/**
 * Test 1: Flat space (R=1) → Straight-line motion
 */
bool testFlatSpaceGeodesic() {
    std::cout << "\n=== Test 1: Flat Space Geodesic (R=1) ===\n";
    std::cout << "Expected: Straight-line motion (no curvature)\n\n";

    GeodesicParticle p;
    p.x = -5.0f; p.y = 0.0f; p.z = 0.0f;
    p.u_x = 0.5f; p.u_y = 0.0f; p.u_z = 0.0f;
    p.u_t = std::sqrt(1.0f + p.u_x*p.u_x);  // Normalization
    p.tau = 0.0f;

    const float dtau = 0.01f;
    const int num_steps = 1000;

    // Record initial state
    float x_initial = p.x;
    float E_initial = p.getEnergy(1.0f);

    // Evolve
    for (int step = 0; step < num_steps; ++step) {
        p.recordState();
        geodesicRK4Step(p, flatRField, dtau);
    }

    // Verify straight-line motion
    float x_expected = x_initial + p.u_x * num_steps * dtau;
    float x_error = std::abs(p.x - x_expected) / std::abs(x_expected);

    // Energy conservation
    float E_final = p.getEnergy(1.0f);
    float E_error = std::abs(E_final - E_initial) / E_initial;

    std::cout << "Initial x: " << x_initial << "\n";
    std::cout << "Final x: " << p.x << " (expected: " << x_expected << ")\n";
    std::cout << "Position error: " << (x_error * 100.0f) << "%\n";
    std::cout << "Energy drift: " << (E_error * 100.0f) << "%\n\n";

    bool pos_pass = x_error < 0.01f;  // 1%
    bool energy_pass = E_error < 0.001f;  // 0.1%

    std::cout << "Quality Gates:\n";
    std::cout << "  Position accuracy (< 1%): " << (pos_pass ? "PASS ✓" : "FAIL ✗") << "\n";
    std::cout << "  Energy conservation (< 0.1%): " << (energy_pass ? "PASS ✓" : "FAIL ✗") << "\n";

    return pos_pass && energy_pass;
}

/**
 * Test 2: Gaussian R-peak → Curved trajectory
 */
bool testCurvedSpaceGeodesic() {
    std::cout << "\n=== Test 2: Curved Space Geodesic (Gaussian R-peak) ===\n";
    std::cout << "Expected: Trajectory deflection due to curvature\n\n";

    const float A = 0.5f;  // Peak amplitude
    const float sigma = 2.0f;  // Peak width

    auto getRField = [=](float x, float y, float z) -> RFieldData {
        return gaussianRPeak(x, y, z, A, sigma);
    };

    GeodesicParticle p;
    p.x = -8.0f; p.y = 2.0f; p.z = 0.0f;  // Start with impact parameter
    p.u_x = 0.5f; p.u_y = 0.0f; p.u_z = 0.0f;  // Move in +x direction
    p.u_t = std::sqrt(1.0f + p.u_x*p.u_x);
    p.tau = 0.0f;

    const float dtau = 0.01f;
    const int num_steps = 3000;

    float y_initial = p.y;
    auto data_initial = getRField(p.x, p.y, p.z);
    float E_initial = p.getEnergy(data_initial.R);

    // Evolve
    for (int step = 0; step < num_steps; ++step) {
        if (step % 100 == 0) {
            p.recordState();
        }
        geodesicRK4Step(p, getRField, dtau);
    }

    // Measure deflection
    float y_final = p.y;
    float deflection = std::abs(y_final - y_initial);

    // Energy conservation
    auto data_final = getRField(p.x, p.y, p.z);
    float E_final = p.getEnergy(data_final.R);
    float E_error = std::abs(E_final - E_initial) / E_initial;

    std::cout << "Peak amplitude: " << A << ", width: " << sigma << "\n";
    std::cout << "Impact parameter: " << y_initial << "\n";
    std::cout << "Deflection: " << deflection << "\n";
    std::cout << "Energy drift: " << (E_error * 100.0f) << "%\n\n";

    // Quality gates
    bool deflection_exists = deflection > 0.01f;  // Non-zero deflection
    bool energy_pass = E_error < 0.001f;  // 0.1%

    std::cout << "Quality Gates:\n";
    std::cout << "  Deflection exists (> 0.01): " << (deflection_exists ? "PASS ✓" : "FAIL ✗") << "\n";
    std::cout << "  Energy conservation (< 0.1%): " << (energy_pass ? "PASS ✓" : "FAIL ✗") << "\n";

    return deflection_exists && energy_pass;
}

/**
 * Main test runner
 */
int runGeodesic3DTest() {
    std::cout << "========================================\n";
    std::cout << "  3D Geodesic Equation Validation\n";
    std::cout << "========================================\n";
    std::cout << "Metric: g_μν = R²(x,y,z)·η_μν\n";
    std::cout << "Equation: d²x^μ/dτ² = -Γ^μ_αβ u^α u^β\n\n";

    bool all_pass = true;

    all_pass &= testFlatSpaceGeodesic();
    all_pass &= testCurvedSpaceGeodesic();

    std::cout << "\n========================================\n";
    std::cout << "FINAL VERDICT: " << (all_pass ? "PASS ✓" : "FAIL ✗") << "\n";
    std::cout << "========================================\n";

    return all_pass ? 0 : 1;
}
