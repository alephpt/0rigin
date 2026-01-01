/**
 * test/test_geodesic_verification.cpp
 *
 * Geodesic Equation Verification Test
 *
 * Goal: Verify that particles follow geodesics in SMFT curved spacetime
 * Metric: g_μν = R²(x,y) × diag(-(1-v²), 1, 1)
 * Equation: d²x^μ/dτ² + Γ^μ_νλ(dx^ν/dτ)(dx^λ/dτ) = 0
 *
 * Test Plan:
 * 1. Initialize particle as Dirac wavepacket Gaussian
 * 2. Evolve Dirac equation in SMFT curved spacetime
 * 3. Track wavepacket center of mass trajectory
 * 4. Independently solve geodesic equation with same initial conditions
 * 5. Compare trajectories: deviation should be <1%
 *
 * Quality Gate: max_deviation < 0.01 (1%)
 */

#include "GeodesicIntegrator.h"
#include "DiracEvolution.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <array>
#include <algorithm>

/**
 * Helper: Generate simple R-field for testing
 * R(x,y) = 1 + A*exp(-((x-x0)²+(y-y0)²)/σ²)
 */
std::vector<double> generateRField(int Nx, int Ny, double amplitude = 0.5, double sigma = 15.0) {
    std::vector<double> R(Nx * Ny);
    double x0 = Nx / 2.0;
    double y0 = Ny / 2.0;

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            double dx = i - x0;
            double dy = j - y0;
            double r2 = dx * dx + dy * dy;
            R[j * Nx + i] = 1.0 + amplitude * std::exp(-r2 / (sigma * sigma));
        }
    }
    return R;
}

/**
 * Helper: Generate velocity field from R-field
 * v(x,y) = sqrt(1 - 1/R²(x,y)), clamped to [0, 1)
 */
std::vector<double> generateVelocityField(const std::vector<double>& R_field) {
    std::vector<double> v = R_field;
    for (size_t i = 0; i < v.size(); ++i) {
        double R = v[i];
        if (R > 1.0) {
            v[i] = std::sqrt(1.0 - 1.0 / (R * R));
        } else {
            v[i] = 0.0;
        }
    }
    return v;
}

/**
 * Extract center of mass from Dirac wavepacket
 */
std::array<double, 2> computeWavepacketCenter(const DiracEvolution& dirac) {
    std::vector<float> density = dirac.getDensity();
    int Nx = dirac.getNx();
    int Ny = dirac.getNy();

    double x_sum = 0.0;
    double y_sum = 0.0;
    double density_sum = 0.0;

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            double rho = density[j * Nx + i];
            x_sum += i * rho;
            y_sum += j * rho;
            density_sum += rho;
        }
    }

    if (density_sum < 1e-10) {
        return {Nx / 2.0, Ny / 2.0};
    }

    return {x_sum / density_sum, y_sum / density_sum};
}

int main() {
    std::cout << "\n";
    std::cout << "========================================\n";
    std::cout << "    GEODESIC EQUATION VERIFICATION\n";
    std::cout << "========================================\n\n";

    // Configuration
    const int Nx = 128;
    const int Ny = 128;
    const double delta = 2.5;
    const double dt = 0.001;
    const int num_steps = 200;  // 0.2 time units

    std::cout << "Configuration:\n";
    std::cout << "  Grid: " << Nx << " × " << Ny << "\n";
    std::cout << "  Delta (vacuum potential): " << delta << "\n";
    std::cout << "  Time step dt: " << dt << "\n";
    std::cout << "  Total steps: " << num_steps << "\n";
    std::cout << "  Total time: " << (dt * num_steps) << "\n\n";

    // Step 1: Generate R-field and velocity field
    std::cout << "[Step 1] Generating curved spacetime fields...\n";
    std::vector<double> R_field = generateRField(Nx, Ny, 0.5, 15.0);
    std::vector<double> v_field = generateVelocityField(R_field);

    double R_min = *std::min_element(R_field.begin(), R_field.end());
    double R_max = *std::max_element(R_field.begin(), R_field.end());
    std::cout << "  R-field range: [" << R_min << ", " << R_max << "]\n";
    std::cout << "  Spacetime curvature created ✓\n\n";

    // Step 2: Initialize Dirac wavepacket
    std::cout << "[Step 2] Initializing Dirac wavepacket...\n";
    DiracEvolution dirac(Nx, Ny);
    double x0_init = Nx / 2.0;
    double y0_init = Ny / 2.0;
    dirac.initialize(x0_init, y0_init, 3.0);  // Gaussian with σ=3

    std::array<double, 2> center = computeWavepacketCenter(dirac);
    std::cout << "  Initial position: (" << center[0] << ", " << center[1] << ")\n";
    std::cout << "  Initial width: σ = 3.0\n";
    std::cout << "  Wavepacket initialized ✓\n\n";

    // Step 3: Extract initial velocity from Dirac evolution
    std::cout << "[Step 3] Computing initial velocity from Dirac state...\n";
    // For static initial state, velocity is approximately zero
    // but we'll perturb it slightly for testing
    double vx_init = 0.1;  // Small initial velocity
    double vy_init = 0.05;
    std::cout << "  Initial velocity: (" << vx_init << ", " << vy_init << ")\n";
    std::cout << "  (Small initial velocity for trajectory analysis)\n\n";

    // Step 4: Initialize geodesic integrator
    std::cout << "[Step 4] Initializing geodesic integrator...\n";
    GeodesicIntegrator integrator(Nx, Ny, delta);
    std::cout << "  Christoffel symbols computed dynamically\n";
    std::cout << "  RK4 geodesic solver ready\n\n";

    // Step 5: Solve geodesic equation
    std::cout << "[Step 5] Solving geodesic equation...\n";
    std::array<double, 2> geo_initial_pos = {center[0], center[1]};
    std::array<double, 2> geo_initial_vel = {vx_init, vy_init};
    std::vector<GeodesicIntegrator::GeodesicPoint> geodesic_trajectory =
        integrator.integrateGeodesic(geo_initial_pos, geo_initial_vel,
                                     R_field, v_field, dt, num_steps);
    std::cout << "  Geodesic trajectory computed: " << geodesic_trajectory.size() << " points\n";
    std::cout << "  Geodesic solved ✓\n\n";

    // Step 6: Evolve Dirac equation
    std::cout << "[Step 6] Evolving Dirac equation in curved spacetime...\n";
    std::vector<std::array<double, 3>> dirac_trajectory;  // {t, x, y}

    // Store initial state
    center = computeWavepacketCenter(dirac);
    dirac_trajectory.push_back({0.0, center[0], center[1]});

    // Evolve with R-field as mass field (curved spacetime effect)
    for (int step = 0; step < num_steps; ++step) {
        // Convert R-field to mass field: m(x) = Δ × R(x)
        std::vector<float> mass_field(Nx * Ny);
        for (size_t i = 0; i < R_field.size(); ++i) {
            mass_field[i] = static_cast<float>(delta * R_field[i]);
        }

        // Evolve Dirac equation
        dirac.step(mass_field, dt);

        // Extract center of mass
        center = computeWavepacketCenter(dirac);
        double t = (step + 1) * dt;
        dirac_trajectory.push_back({t, center[0], center[1]});
    }
    std::cout << "  Dirac evolution completed: " << dirac_trajectory.size() << " points\n";
    std::cout << "  Dirac equation solved ✓\n\n";

    // Step 7: Compare trajectories
    std::cout << "[Step 7] Comparing Dirac trajectory to geodesic prediction...\n";

    // Convert geodesic trajectory to expected format {t, x, y}
    std::vector<std::array<double, 3>> geodesic_traj;
    for (const auto& gp : geodesic_trajectory) {
        geodesic_traj.push_back({gp.t, gp.x, gp.y});
    }

    std::vector<double> deviations =
        integrator.compareTrajectories(dirac_trajectory, geodesic_traj);

    double max_deviation = 0.0;
    double avg_deviation = 0.0;
    if (!deviations.empty()) {
        max_deviation = *std::max_element(deviations.begin(), deviations.end());
        for (double dev : deviations) {
            avg_deviation += dev;
        }
        avg_deviation /= deviations.size();
    }

    std::cout << "  Trajectory comparison:\n";
    std::cout << "    Points compared: " << deviations.size() << "\n";
    std::cout << "    Average deviation: " << std::fixed << std::setprecision(6)
              << (avg_deviation * 100.0) << "%\n";
    std::cout << "    Maximum deviation: " << (max_deviation * 100.0) << "%\n\n";

    // Step 8: Validation
    std::cout << "[Step 8] Validation Result\n";
    const double DEVIATION_TOLERANCE = 0.01;  // 1%
    bool validation_pass = (max_deviation < DEVIATION_TOLERANCE);

    std::cout << "  Quality gate: max_deviation < " << (DEVIATION_TOLERANCE * 100.0) << "%\n";
    if (validation_pass) {
        std::cout << "  ✓ PASSED: Particles follow geodesics within tolerance\n";
    } else {
        std::cout << "  ✗ FAILED: Deviation exceeds tolerance\n";
        std::cout << "    Excess: " << ((max_deviation - DEVIATION_TOLERANCE) * 100.0) << "%\n";
    }
    std::cout << "\n";

    // Step 9: Generate output CSV
    std::cout << "[Step 9] Generating output...\n";
    std::ofstream csv("output/test/geodesic_verification.csv");
    csv << "time,dirac_x,dirac_y,geodesic_x,geodesic_y,deviation\n";

    for (size_t i = 0; i < std::min(dirac_trajectory.size(), geodesic_traj.size()); ++i) {
        csv << std::fixed << std::setprecision(6);
        csv << dirac_trajectory[i][0] << ","
            << dirac_trajectory[i][1] << ","
            << dirac_trajectory[i][2] << ","
            << geodesic_traj[i][1] << ","
            << geodesic_traj[i][2] << ","
            << (i < deviations.size() ? deviations[i] : 0.0) << "\n";
    }
    csv.close();
    std::cout << "  CSV written: output/test/geodesic_verification.csv\n\n";

    // Final summary
    std::cout << "========================================\n";
    std::cout << "TEST SUMMARY\n";
    std::cout << "========================================\n\n";

    if (validation_pass) {
        std::cout << "Result: ✓ PASSED\n\n";
        std::cout << "Conclusion:\n";
        std::cout << "  Particles follow geodesics in SMFT curved spacetime\n";
        std::cout << "  with high accuracy (deviation < 1%).\n";
        return 0;
    } else {
        std::cout << "Result: ✗ FAILED\n\n";
        std::cout << "Issue:\n";
        std::cout << "  Trajectory deviation exceeds 1% tolerance.\n";
        std::cout << "  Check Christoffel symbol computation or geodesic integrator.\n";
        return 1;
    }
}
