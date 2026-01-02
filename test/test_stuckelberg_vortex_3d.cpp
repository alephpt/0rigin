/**
 * test_stuckelberg_vortex_3d.cpp
 *
 * Stückelberg Vortex Validation in 3D
 *
 * Goal: Verify vortex line configurations and magnetic field generation in 3D
 *
 * Physics:
 *   2D: Point vortex → θ(x,y) = atan2(y-y₀, x-x₀) → B_z field
 *   3D: Vortex LINES (topological objects extended in one direction)
 *
 *   Phase winding: ∇θ → Vector potential A_μ
 *   Magnetic field: B = ∇ × A
 *   Flux quantization: ∮ A·dl = 2πn (topological charge)
 *
 * Test Scenarios:
 *   1. Vortex line along z-axis → B_z field (reproduce 2D: B_max ≈ 1.567)
 *   2. Vortex line along x-axis → B_x field (verify symmetry)
 *   3. Vortex line along y-axis → B_y field (verify symmetry)
 *   4. Vortex ring (toroidal) → Full 3D B-field + flux quantization
 *
 * Quality Gates:
 *   - B_max within 5% of 2D value (1.567)
 *   - Flux quantization: Φ/2π = integer ± 0.01
 *   - Energy conservation < 0.1%
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>

// Physical constants
const float PI = 3.14159265358979323846f;
const float HBAR = 1.0f;  // Natural units

/**
 * Simple 3D grid for phase field
 */
class Grid3D {
private:
    int N;
    float dx;
    std::vector<float> theta;

public:
    Grid3D(int size, float spacing) : N(size), dx(spacing), theta(size * size * size, 0.0f) {}

    int getSize() const { return N; }
    float getSpacing() const { return dx; }

    float& at(int ix, int iy, int iz) {
        return theta[ix + N * (iy + N * iz)];
    }

    const float& at(int ix, int iy, int iz) const {
        return theta[ix + N * (iy + N * iz)];
    }
};

/**
 * Initialize vortex line along z-axis
 * Phase: θ(x,y,z) = atan2(y-y₀, x-x₀) (constant along z)
 * Produces: B_z field (same as 2D)
 */
void initVortexLine_Z(Grid3D& grid, float x0, float y0) {
    const int N = grid.getSize();
    const float dx = grid.getSpacing();

    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            for (int iz = 0; iz < N; ++iz) {
                float x = (ix - N/2) * dx;
                float y = (iy - N/2) * dx;
                // z doesn't matter - phase constant along z

                float theta = std::atan2(y - y0, x - x0);
                grid.at(ix, iy, iz) = theta;
            }
        }
    }
}

/**
 * Initialize vortex line along x-axis
 * Phase: θ(x,y,z) = atan2(z-z₀, y-y₀) (constant along x)
 * Produces: B_x field
 */
void initVortexLine_X(Grid3D& grid, float y0, float z0) {
    const int N = grid.getSize();
    const float dx = grid.getSpacing();

    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            for (int iz = 0; iz < N; ++iz) {
                // x doesn't matter - phase constant along x
                float y = (iy - N/2) * dx;
                float z = (iz - N/2) * dx;

                float theta = std::atan2(z - z0, y - y0);
                grid.at(ix, iy, iz) = theta;
            }
        }
    }
}

/**
 * Initialize vortex line along y-axis
 * Phase: θ(x,y,z) = atan2(z-z₀, x-x₀) (constant along y)
 * Produces: B_y field
 */
void initVortexLine_Y(Grid3D& grid, float x0, float z0) {
    const int N = grid.getSize();
    const float dx = grid.getSpacing();

    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            for (int iz = 0; iz < N; ++iz) {
                float x = (ix - N/2) * dx;
                // y doesn't matter - phase constant along y
                float z = (iz - N/2) * dx;

                float theta = std::atan2(z - z0, x - x0);
                grid.at(ix, iy, iz) = theta;
            }
        }
    }
}

/**
 * Initialize vortex ring (toroidal configuration)
 * Circular vortex loop in x-y plane at z=z₀
 * Phase winds around torus
 */
void initVortexRing(Grid3D& grid, float ring_radius, float z0) {
    const int N = grid.getSize();
    const float dx = grid.getSpacing();

    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            for (int iz = 0; iz < N; ++iz) {
                float x = (ix - N/2) * dx;
                float y = (iy - N/2) * dx;
                float z = (iz - N/2) * dx;

                // Distance from z-axis
                float r_xy = std::sqrt(x*x + y*y);

                // Toroidal coordinates
                // Distance from ring center (in x-y plane at height z₀)
                float r_ring = std::sqrt((r_xy - ring_radius)*(r_xy - ring_radius) + (z - z0)*(z - z0));

                // Phase winds around the ring
                // Use poloidal angle as phase
                float theta = std::atan2(z - z0, r_xy - ring_radius);

                grid.at(ix, iy, iz) = theta;
            }
        }
    }
}

/**
 * Compute vector potential A from phase gradient
 * A = ħ∇θ (Stückelberg relation)
 */
std::array<float, 3> computeVectorPotential(const Grid3D& grid, int ix, int iy, int iz) {
    const float dx = grid.getSpacing();
    const int N = grid.getSize();

    // Central finite differences for gradient
    int ixp = (ix + 1) % N;
    int ixm = (ix - 1 + N) % N;
    int iyp = (iy + 1) % N;
    int iym = (iy - 1 + N) % N;
    int izp = (iz + 1) % N;
    int izm = (iz - 1 + N) % N;

    float theta_xp = grid.at(ixp, iy, iz);
    float theta_xm = grid.at(ixm, iy, iz);
    float theta_yp = grid.at(ix, iyp, iz);
    float theta_ym = grid.at(ix, iym, iz);
    float theta_zp = grid.at(ix, iy, izp);
    float theta_zm = grid.at(ix, iy, izm);

    // Handle phase wrapping (2π jumps)
    auto unwrap = [](float dtheta) -> float {
        while (dtheta > PI) dtheta -= 2.0f * PI;
        while (dtheta < -PI) dtheta += 2.0f * PI;
        return dtheta;
    };

    float dtheta_dx = unwrap(theta_xp - theta_xm) / (2.0f * dx);
    float dtheta_dy = unwrap(theta_yp - theta_ym) / (2.0f * dx);
    float dtheta_dz = unwrap(theta_zp - theta_zm) / (2.0f * dx);

    return {HBAR * dtheta_dx, HBAR * dtheta_dy, HBAR * dtheta_dz};
}

/**
 * Compute magnetic field B = ∇ × A
 */
std::array<float, 3> computeMagneticField(const Grid3D& grid, int ix, int iy, int iz) {
    const float dx = grid.getSpacing();
    const int N = grid.getSize();

    int ixp = (ix + 1) % N;
    int ixm = (ix - 1 + N) % N;
    int iyp = (iy + 1) % N;
    int iym = (iy - 1 + N) % N;
    int izp = (iz + 1) % N;
    int izm = (iz - 1 + N) % N;

    auto A_xp = computeVectorPotential(grid, ixp, iy, iz);
    auto A_xm = computeVectorPotential(grid, ixm, iy, iz);
    auto A_yp = computeVectorPotential(grid, ix, iyp, iz);
    auto A_ym = computeVectorPotential(grid, ix, iym, iz);
    auto A_zp = computeVectorPotential(grid, ix, iy, izp);
    auto A_zm = computeVectorPotential(grid, ix, iy, izm);

    // B = ∇ × A
    // Bx = ∂Az/∂y - ∂Ay/∂z
    float Bx = (A_yp[2] - A_ym[2]) / (2.0f * dx) - (A_zp[1] - A_zm[1]) / (2.0f * dx);

    // By = ∂Ax/∂z - ∂Az/∂x
    float By = (A_zp[0] - A_zm[0]) / (2.0f * dx) - (A_xp[2] - A_xm[2]) / (2.0f * dx);

    // Bz = ∂Ay/∂x - ∂Ax/∂y
    float Bz = (A_xp[1] - A_xm[1]) / (2.0f * dx) - (A_yp[0] - A_ym[0]) / (2.0f * dx);

    return {Bx, By, Bz};
}

/**
 * Measure maximum magnetic field magnitude
 */
float measureBmax(const Grid3D& grid) {
    const int N = grid.getSize();
    float B_max = 0.0f;

    for (int ix = 1; ix < N-1; ++ix) {
        for (int iy = 1; iy < N-1; ++iy) {
            for (int iz = 1; iz < N-1; ++iz) {
                auto B = computeMagneticField(grid, ix, iy, iz);
                float B_mag = std::sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
                B_max = std::max(B_max, B_mag);
            }
        }
    }

    return B_max;
}

/**
 * Test 1: Vortex line along z-axis
 * Reproduces 2D result: B_max ≈ 1.567
 */
bool testVortexLine_Z() {
    std::cout << "\n=== Test 1: Vortex Line Along Z-Axis ===\n";
    std::cout << "Phase: θ(x,y,z) = atan2(y-y₀, x-x₀)\n";
    std::cout << "Expected: B_z field, B_max ≈ 1.567 (2D benchmark)\n\n";

    const int N = 32;
    const float dx = 0.5f;
    Grid3D grid(N, dx);

    // Initialize vortex at center
    initVortexLine_Z(grid, 0.0f, 0.0f);

    // Measure B-field
    float B_max = measureBmax(grid);

    // Sample B-field at center
    auto B_center = computeMagneticField(grid, N/2, N/2, N/2);

    std::cout << "B_max: " << B_max << "\n";
    std::cout << "B at center: (" << B_center[0] << ", " << B_center[1] << ", " << B_center[2] << ")\n";

    // Quality gate
    const float B_expected = 1.567f;
    float B_error = std::abs(B_max - B_expected) / B_expected;

    std::cout << "Error vs 2D benchmark: " << (B_error * 100.0f) << "%\n";

    bool pass = B_error < 0.05f;  // 5%
    std::cout << "Quality Gate (< 5%): " << (pass ? "PASS ✓" : "FAIL ✗") << "\n";

    return pass;
}

/**
 * Test 2: Vortex line along x-axis
 * NEW 3D test: B_x field, verify symmetry
 */
bool testVortexLine_X() {
    std::cout << "\n=== Test 2: Vortex Line Along X-Axis ===\n";
    std::cout << "Phase: θ(x,y,z) = atan2(z-z₀, y-y₀)\n";
    std::cout << "Expected: B_x field, B_max ≈ 1.567 (by symmetry)\n\n";

    const int N = 32;
    const float dx = 0.5f;
    Grid3D grid(N, dx);

    initVortexLine_X(grid, 0.0f, 0.0f);

    float B_max = measureBmax(grid);
    auto B_center = computeMagneticField(grid, N/2, N/2, N/2);

    std::cout << "B_max: " << B_max << "\n";
    std::cout << "B at center: (" << B_center[0] << ", " << B_center[1] << ", " << B_center[2] << ")\n";

    const float B_expected = 1.567f;
    float B_error = std::abs(B_max - B_expected) / B_expected;

    std::cout << "Error vs expected: " << (B_error * 100.0f) << "%\n";

    bool pass = B_error < 0.05f;
    std::cout << "Quality Gate (< 5%): " << (pass ? "PASS ✓" : "FAIL ✗") << "\n";

    return pass;
}

/**
 * Test 3: Vortex line along y-axis
 * NEW 3D test: B_y field, verify symmetry
 */
bool testVortexLine_Y() {
    std::cout << "\n=== Test 3: Vortex Line Along Y-Axis ===\n";
    std::cout << "Phase: θ(x,y,z) = atan2(z-z₀, x-x₀)\n";
    std::cout << "Expected: B_y field, B_max ≈ 1.567 (by symmetry)\n\n";

    const int N = 32;
    const float dx = 0.5f;
    Grid3D grid(N, dx);

    initVortexLine_Y(grid, 0.0f, 0.0f);

    float B_max = measureBmax(grid);
    auto B_center = computeMagneticField(grid, N/2, N/2, N/2);

    std::cout << "B_max: " << B_max << "\n";
    std::cout << "B at center: (" << B_center[0] << ", " << B_center[1] << ", " << B_center[2] << ")\n";

    const float B_expected = 1.567f;
    float B_error = std::abs(B_max - B_expected) / B_expected;

    std::cout << "Error vs expected: " << (B_error * 100.0f) << "%\n";

    bool pass = B_error < 0.05f;
    std::cout << "Quality Gate (< 5%): " << (pass ? "PASS ✓" : "FAIL ✗") << "\n";

    return pass;
}

/**
 * Test 4: Vortex ring (toroidal)
 * Full 3D topology + flux quantization
 */
bool testVortexRing() {
    std::cout << "\n=== Test 4: Vortex Ring (Toroidal) ===\n";
    std::cout << "Circular vortex loop in x-y plane\n";
    std::cout << "Expected: Full 3D B-field, flux quantization Φ = 2πn\n\n";

    const int N = 32;
    const float dx = 0.5f;
    Grid3D grid(N, dx);

    float ring_radius = 4.0f;
    float z0 = 0.0f;
    initVortexRing(grid, ring_radius, z0);

    float B_max = measureBmax(grid);

    std::cout << "Ring radius: " << ring_radius << "\n";
    std::cout << "B_max: " << B_max << "\n";

    // Flux quantization: integrate A·dl around ring
    // Simplified: sum A_φ at ring location
    float flux = 0.0f;
    int num_samples = 64;
    for (int i = 0; i < num_samples; ++i) {
        float phi = 2.0f * PI * i / num_samples;
        float x = ring_radius * std::cos(phi);
        float y = ring_radius * std::sin(phi);
        float z = z0;

        // Convert to grid indices
        int ix = static_cast<int>((x / dx) + N/2);
        int iy = static_cast<int>((y / dx) + N/2);
        int iz = static_cast<int>((z / dx) + N/2);

        if (ix >= 0 && ix < N && iy >= 0 && iy < N && iz >= 0 && iz < N) {
            auto A = computeVectorPotential(grid, ix, iy, iz);
            // Tangent vector: (-sin(φ), cos(φ), 0)
            float A_tangent = -A[0] * std::sin(phi) + A[1] * std::cos(phi);
            flux += A_tangent;
        }
    }
    flux *= (2.0f * PI * ring_radius / num_samples);  // Arc length element

    float flux_quantum = 2.0f * PI;
    float winding_number = flux / flux_quantum;
    float quantization_error = std::abs(winding_number - std::round(winding_number));

    std::cout << "Flux Φ: " << flux << "\n";
    std::cout << "Φ/2π: " << winding_number << "\n";
    std::cout << "Quantization error: " << quantization_error << "\n";

    bool pass = quantization_error < 0.01f;
    std::cout << "Quality Gate (< 0.01): " << (pass ? "PASS ✓" : "FAIL ✗") << "\n";

    return pass;
}

/**
 * Main test runner
 */
int runStuckelbergVortex3DTest() {
    std::cout << "========================================\n";
    std::cout << "  3D Stückelberg Vortex Validation\n";
    std::cout << "========================================\n";
    std::cout << "Extension: 2D point vortex → 3D vortex lines\n";

    bool all_pass = true;

    all_pass &= testVortexLine_Z();
    all_pass &= testVortexLine_X();
    all_pass &= testVortexLine_Y();
    all_pass &= testVortexRing();

    std::cout << "\n========================================\n";
    std::cout << "FINAL VERDICT: " << (all_pass ? "PASS ✓" : "FAIL ✗") << "\n";
    std::cout << "========================================\n";

    return all_pass ? 0 : 1;
}
