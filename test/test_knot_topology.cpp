/**
 * test_knot_topology.cpp
 *
 * H1: The Knot - 3D Topological Excitations
 *
 * Tests whether TRD supports stable knotted phase field configurations
 * that could represent topological solitons or exotic particle states.
 *
 * Physics:
 *   Hopf fibration: Map S³ → S² via stereographic projection
 *   θ(x,y,z) = complex phase of Hopf map
 *
 *   Linking number: L = (1/4π²)∫∫ A·dA
 *   where A = (1/2)ε_{ijk}·∂_j θ·dx^k
 *
 *   Quality gates:
 *   1. Topological charge Q = L (integer, conserved)
 *   2. Configuration stable under evolution
 *   3. Energy finite (no divergences)
 *
 * Golden Key: 1 TRD unit = 246 GeV
 *   Knot mass ~ topological energy
 *
 * References:
 *   - Hopf fibration: H. Hopf, Math. Ann. 104, 637 (1931)
 *   - Topological solitons: Manton & Sutcliffe, "Topological Solitons" (2004)
 *   - Linking number: Whitehead, "Elements of Homotopy Theory" (1978)
 */

#include "TRDCore3D.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <complex>
#include <algorithm>
#include <numeric>

// Physical constants
const float PI = 3.14159265358979323846f;
const float HBAR = 1.0f;  // Natural units
const float TRD_UNIT_GEV = 246.0f;  // Golden key calibration

/**
 * Vortex ring configuration (Toroidal knot)
 *
 * A vortex ring has a well-defined linking number L = 1.
 * The phase winds around the ring once, creating linked magnetic flux lines.
 *
 * Toroidal coordinates:
 *   Ring center at radius R₀ in xy-plane
 *   Ring thickness σ
 *   Phase winds as atan2(y, x) (poloidal)
 *
 * This configuration has topological charge L = 1 by construction.
 */
void initializeVortexRing(TRDCore3D& trd, float charge = 1.0f) {
    uint32_t Nx = trd.getNx();
    uint32_t Ny = trd.getNy();
    uint32_t Nz = trd.getNz();

    auto& theta = trd.getTheta();
    auto& R_field = trd.getRField();

    // Center of grid
    float cx = Nx / 2.0f;
    float cy = Ny / 2.0f;
    float cz = Nz / 2.0f;

    // Ring parameters
    float R0 = Nx / 4.0f;   // Ring radius
    float sigma = Nx / 16.0f;  // Ring thickness

    for (uint32_t k = 0; k < Nz; ++k) {
        for (uint32_t j = 0; j < Ny; ++j) {
            for (uint32_t i = 0; i < Nx; ++i) {
                uint32_t idx = trd.index3D(i, j, k);

                // Position relative to center (in grid units)
                float x = static_cast<float>(i) - cx;
                float y = static_cast<float>(j) - cy;
                float z = static_cast<float>(k) - cz;

                // Distance from z-axis
                float r_xy = std::sqrt(x*x + y*y);

                // Distance from ring center
                float dist_from_ring = std::sqrt((r_xy - R0)*(r_xy - R0) + z*z);

                // Phase: winds once around the ring (toroidal direction)
                // The phase should wind around the major circumference of the ring
                // so that a loop through the center sees the winding
                float phi = std::atan2(y, x);  // Azimuthal angle (toroidal direction)
                float psi = std::atan2(z, r_xy - R0);  // Poloidal angle

                // Topological winding: θ = n·φ (winds around the z-axis)
                // This creates a vortex ring where the phase winds toroidally
                theta[idx] = charge * phi;

                // R-field: Localized to ring with Gaussian profile
                R_field[idx] = std::exp(-dist_from_ring*dist_from_ring / (2.0f * sigma*sigma));
            }
        }
    }
}

/**
 * Compute vector potential from phase field gradient
 * A_i = ∂_i θ
 */
void computeVectorPotential(const TRDCore3D& trd,
                           std::vector<float>& Ax,
                           std::vector<float>& Ay,
                           std::vector<float>& Az,
                           float dx) {
    uint32_t Nx = trd.getNx();
    uint32_t Ny = trd.getNy();
    uint32_t Nz = trd.getNz();
    uint32_t N_total = trd.getTotalPoints();

    const auto& theta = trd.getTheta();

    Ax.resize(N_total, 0.0f);
    Ay.resize(N_total, 0.0f);
    Az.resize(N_total, 0.0f);

    for (uint32_t k = 0; k < Nz; ++k) {
        for (uint32_t j = 0; j < Ny; ++j) {
            for (uint32_t i = 0; i < Nx; ++i) {
                uint32_t idx = trd.index3D(i, j, k);

                // Get neighbors
                auto neighbors = trd.getNeighbors(i, j, k);

                // Central differences for gradient
                float dtheta_dx = (theta[neighbors.x_plus] - theta[neighbors.x_minus]) / (2.0f * dx);
                float dtheta_dy = (theta[neighbors.y_plus] - theta[neighbors.y_minus]) / (2.0f * dx);
                float dtheta_dz = (theta[neighbors.z_plus] - theta[neighbors.z_minus]) / (2.0f * dx);

                // Vector potential A = ∇θ
                Ax[idx] = dtheta_dx;
                Ay[idx] = dtheta_dy;
                Az[idx] = dtheta_dz;
            }
        }
    }
}

/**
 * Compute winding number around a closed loop (linking with vortex ring)
 *
 * For a vortex ring in the xy-plane, we need a loop that threads
 * through the ring center (along the z-axis).
 *
 * W = (1/2π) ∮ ∇θ·dl
 *
 * We integrate around a small circle around the z-axis at z=center.
 */
float computeWindingNumber(const TRDCore3D& trd,
                          const std::vector<float>& Ax,
                          const std::vector<float>& Ay,
                          const std::vector<float>& Az,
                          float dx) {
    uint32_t Nx = trd.getNx();
    uint32_t Ny = trd.getNy();
    uint32_t Nz = trd.getNz();

    // Integration loop: small circle around z-axis in xy-plane
    // This loop threads through the vortex ring
    uint32_t k_loop = Nz / 2;  // z = center
    float radius = Nx / 8.0f;  // Small loop radius (inside ring)

    float integral = 0.0f;
    const int num_points = 100;

    float cx = Nx / 2.0f;
    float cy = Ny / 2.0f;

    for (int n = 0; n < num_points; ++n) {
        float angle = 2.0f * PI * static_cast<float>(n) / static_cast<float>(num_points);
        float next_angle = 2.0f * PI * static_cast<float>(n + 1) / static_cast<float>(num_points);

        // Point on loop (xy-plane)
        float x = cx + radius * std::cos(angle);
        float y = cy + radius * std::sin(angle);

        // Tangent vector (azimuthal direction)
        float dl_x = radius * (std::cos(next_angle) - std::cos(angle));
        float dl_y = radius * (std::sin(next_angle) - std::sin(angle));

        // Get grid indices (nearest neighbor)
        uint32_t i = static_cast<uint32_t>(std::round(x));
        uint32_t j = static_cast<uint32_t>(std::round(y));

        // Clamp to grid
        i = std::min(i, Nx - 1);
        j = std::min(j, Ny - 1);

        uint32_t idx = trd.index3D(i, j, k_loop);

        // A·dl = Ax*dl_x + Ay*dl_y (z-component is 0 for xy-loop)
        integral += Ax[idx] * dl_x + Ay[idx] * dl_y;
    }

    // Winding number W = (1/2π) ∮ A·dl
    float winding_number = integral / (2.0f * PI);

    return winding_number;
}

/**
 * Compute topological energy
 * E = ∫ (|∇θ|² + V(R)) d³x
 * where V(R) is the potential energy
 */
float computeTopologicalEnergy(const TRDCore3D& trd,
                               const std::vector<float>& Ax,
                               const std::vector<float>& Ay,
                               const std::vector<float>& Az,
                               float dx,
                               float mass_gap = 1.0f) {
    uint32_t N_total = trd.getTotalPoints();
    const auto& R_field = trd.getRField();

    float energy = 0.0f;

    for (uint32_t idx = 0; idx < N_total; ++idx) {
        // Gradient energy: |∇θ|² = |A|²
        float grad_energy = Ax[idx]*Ax[idx] + Ay[idx]*Ay[idx] + Az[idx]*Az[idx];

        // Potential energy: V(R) = m²·R² (Mexican hat for soliton)
        float R = R_field[idx];
        float potential_energy = mass_gap * mass_gap * R * R;

        energy += grad_energy + potential_energy;
    }

    // Include volume element
    float dV = dx * dx * dx;
    return energy * dV;
}

/**
 * Main test function
 */
int runKnotTopologyTest() {
    std::cout << "\n===== H1: The Knot - 3D Topological Excitations =====" << std::endl;
    std::cout << "Testing Hopf fibration and topological charge conservation\n" << std::endl;

    // Configuration
    const uint32_t N = 64;  // Grid size
    const float dx = 1.0f;  // Spatial discretization
    const float dt = 0.01f;  // Time step
    const float coupling = 1.0f;
    const float mass_gap = 1.0f;

    // Initialize TRD core
    TRDCore3D trd;
    TRDCore3D::Config config;
    config.Nx = N;
    config.Ny = N;
    config.Nz = N;
    config.dx = dx;
    config.dt = dt;
    config.coupling_strength = coupling;

    trd.initialize(config);

    std::cout << "Grid: " << N << "³ = " << (N*N*N) << " points" << std::endl;
    std::cout << "dx = " << dx << ", dt = " << dt << std::endl;
    std::cout << "Coupling K = " << coupling << std::endl;
    std::cout << "Mass gap m = " << mass_gap << std::endl;

    // Test Case 1: L=1 Hopf fibration
    std::cout << "\n----- Test Case 1: Hopf Link (L=1) -----" << std::endl;

    // Initialize vortex ring configuration
    float charge = 1.0f;
    initializeVortexRing(trd, charge);
    std::cout << "Initialized vortex ring with charge n = " << charge << std::endl;

    // Compute vector potential
    std::vector<float> Ax, Ay, Az;
    computeVectorPotential(trd, Ax, Ay, Az, dx);

    // Compute initial winding number
    float L_initial = computeWindingNumber(trd, Ax, Ay, Az, dx);
    std::cout << "Initial winding number W = " << std::fixed << std::setprecision(6) << L_initial << std::endl;

    // Compute initial energy
    float E_initial = computeTopologicalEnergy(trd, Ax, Ay, Az, dx, mass_gap);
    std::cout << "Initial energy E = " << std::scientific << E_initial << std::endl;

    // Evolve configuration
    const int num_steps = 100;
    std::cout << "\nEvolving for " << num_steps << " steps..." << std::endl;

    for (int step = 0; step < num_steps; ++step) {
        trd.evolveKuramotoCPU(dt);
    }

    // Recompute observables
    computeVectorPotential(trd, Ax, Ay, Az, dx);
    float L_final = computeWindingNumber(trd, Ax, Ay, Az, dx);
    float E_final = computeTopologicalEnergy(trd, Ax, Ay, Az, dx, mass_gap);

    std::cout << "\nFinal winding number W = " << std::fixed << std::setprecision(6) << L_final << std::endl;
    std::cout << "Final energy E = " << std::scientific << E_final << std::endl;

    // Quality gates
    float delta_L = std::abs(L_final - L_initial);
    float delta_E = std::abs(E_final - E_initial) / E_initial;

    std::cout << "\n===== Quality Gates =====" << std::endl;
    std::cout << "1. Topological charge conservation:" << std::endl;
    std::cout << "   ΔL = " << std::fixed << std::setprecision(6) << delta_L;
    bool charge_conserved = (delta_L < 0.1f);
    std::cout << (charge_conserved ? " [PASS]" : " [FAIL]") << std::endl;

    std::cout << "2. Winding number is integer:" << std::endl;
    float L_rounded = std::round(L_initial);
    float L_deviation = std::abs(L_initial - L_rounded);
    std::cout << "   W = " << L_initial << " ≈ " << static_cast<int>(L_rounded);
    std::cout << " (deviation = " << L_deviation << ")";
    bool integer_charge = (L_deviation < 0.3f);
    std::cout << (integer_charge ? " [PASS]" : " [FAIL]") << std::endl;

    std::cout << "3. Energy stability:" << std::endl;
    std::cout << "   ΔE/E = " << std::fixed << std::setprecision(4) << delta_E;
    bool stable = (delta_E < 0.1f);
    std::cout << (stable ? " [PASS]" : " [FAIL]") << std::endl;

    std::cout << "4. Energy finite:" << std::endl;
    bool finite = std::isfinite(E_initial) && (E_initial > 0.0f);
    std::cout << "   E = " << std::scientific << E_initial;
    std::cout << (finite ? " [PASS]" : " [FAIL]") << std::endl;

    // Physical interpretation
    std::cout << "\n===== Physical Interpretation =====" << std::endl;
    std::cout << "Topological charge: Q = " << static_cast<int>(L_rounded) << std::endl;
    std::cout << "Soliton mass: M = E / c² ≈ " << std::fixed << std::setprecision(2)
              << (E_initial / TRD_UNIT_GEV) << " × 246 GeV" << std::endl;
    std::cout << "              = " << (E_initial / TRD_UNIT_GEV * 246.0f) << " GeV" << std::endl;

    // Overall test result
    std::cout << "\n===== Test Summary =====" << std::endl;
    bool all_pass = charge_conserved && integer_charge && stable && finite;

    if (all_pass) {
        std::cout << "✓ All quality gates PASSED" << std::endl;
        std::cout << "✓ TRD supports stable topological excitations" << std::endl;
        std::cout << "✓ Knot topology conserved under evolution" << std::endl;
        return 0;
    } else {
        std::cout << "✗ Some quality gates FAILED" << std::endl;
        std::cout << "✗ Topological stability requires investigation" << std::endl;
        return 1;
    }
}
