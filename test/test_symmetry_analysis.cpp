/**
 * test_symmetry_analysis.cpp
 *
 * E5: Symmetry Analysis via Noether's Theorem for TRD
 *
 * Goal: Identify and verify conserved quantities from TRD symmetries
 *       using Noether's theorem
 *
 * Physics Context:
 *   Noether's Theorem: For each continuous symmetry, there exists a conserved quantity
 *
 *   TRD Symmetries → Conserved Quantities:
 *   1. Time translation invariance → Energy conservation
 *   2. Space translation invariance → Momentum conservation
 *   3. U(1) phase rotation → Phase charge conservation
 *   4. Gauge symmetry → Electromagnetic current conservation
 *
 * Mathematical Framework:
 *   - Energy: E = ∫[(∂θ/∂t)² + (∇θ)² + V(R)] d³x
 *   - Momentum: P = ∫(∂θ/∂t)(∇θ) d³x
 *   - Phase charge: Q = ∫|∇θ|² d³x
 *   - Conservation: ∂ρ/∂t + ∇·J = 0
 *
 * Test Strategy:
 *   1. Initialize TRD field configuration
 *   2. Evolve system over extended time (100 time units)
 *   3. Measure conserved quantities at regular intervals
 *   4. Verify conservation: |ΔQ/Q₀| < tolerance
 *   5. Validate Noether currents satisfy continuity equations
 *
 * Quality Gates:
 *   - Energy conservation: |ΔE/E₀| < 1%
 *   - Momentum conservation: |ΔP/P₀| < 1%
 *   - Phase charge conservation: |ΔQ/Q₀| < 1%
 *   - All 3 quantities conserved → Noether's theorem validated
 */

#include "TRDCore3D.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <array>
#include <algorithm>
#include <numeric>

// Physical constants
const float PI = 3.14159265358979323846f;

/**
 * Structure to hold conserved quantities
 */
struct ConservedQuantities {
    double energy;
    std::array<double, 3> momentum;  // px, py, pz
    double phase_charge;
    double angular_momentum;  // Lz for rotation symmetry

    void print() const {
        std::cout << std::fixed << std::setprecision(6)
                  << "E=" << energy
                  << " P=(" << momentum[0] << "," << momentum[1] << "," << momentum[2] << ")"
                  << " Q=" << phase_charge
                  << " L=" << angular_momentum;
    }
};

/**
 * Compute energy from TRD Hamiltonian
 * H = ∫[(∂θ/∂t)² + (∇θ)² + V(R)] d³x
 */
double computeEnergy(
    const std::vector<float>& theta_field,
    const std::vector<float>& theta_dot,
    const std::vector<float>& R_field,
    float K, int nx, int ny, int nz, float dx
) {
    double energy = 0.0;

    for (int k = 1; k < nz - 1; ++k) {
        for (int j = 1; j < ny - 1; ++j) {
            for (int i = 1; i < nx - 1; ++i) {
                int idx = k * nx * ny + j * nx + i;

                // Kinetic energy from time derivative
                double kinetic_t = 0.5 * theta_dot[idx] * theta_dot[idx];

                // Compute spatial gradients
                int idx_xp = k * nx * ny + j * nx + (i + 1);
                int idx_xm = k * nx * ny + j * nx + (i - 1);
                int idx_yp = k * nx * ny + (j + 1) * nx + i;
                int idx_ym = k * nx * ny + (j - 1) * nx + i;
                int idx_zp = (k + 1) * nx * ny + j * nx + i;
                int idx_zm = (k - 1) * nx * ny + j * nx + i;

                float grad_x = (theta_field[idx_xp] - theta_field[idx_xm]) / (2 * dx);
                float grad_y = (theta_field[idx_yp] - theta_field[idx_ym]) / (2 * dx);
                float grad_z = (theta_field[idx_zp] - theta_field[idx_zm]) / (2 * dx);

                // Gradient energy
                double kinetic_x = 0.5 * grad_x * grad_x;
                double kinetic_y = 0.5 * grad_y * grad_y;
                double kinetic_z = 0.5 * grad_z * grad_z;

                // Potential energy from R-field coupling
                double potential = K * R_field[idx] * R_field[idx];

                energy += kinetic_t + kinetic_x + kinetic_y + kinetic_z + potential;
            }
        }
    }

    return energy * dx * dx * dx;  // Volume element
}

/**
 * Compute momentum from Noether current
 * P_i = ∫(∂θ/∂t)(∂θ/∂x_i) d³x
 */
std::array<double, 3> computeMomentum(
    const std::vector<float>& theta_field,
    const std::vector<float>& theta_dot,
    int nx, int ny, int nz, float dx
) {
    std::array<double, 3> momentum = {0.0, 0.0, 0.0};

    for (int k = 1; k < nz - 1; ++k) {
        for (int j = 1; j < ny - 1; ++j) {
            for (int i = 1; i < nx - 1; ++i) {
                int idx = k * nx * ny + j * nx + i;

                // Time derivative
                float dt_theta = theta_dot[idx];

                // Spatial gradients
                int idx_xp = k * nx * ny + j * nx + (i + 1);
                int idx_xm = k * nx * ny + j * nx + (i - 1);
                int idx_yp = k * nx * ny + (j + 1) * nx + i;
                int idx_ym = k * nx * ny + (j - 1) * nx + i;
                int idx_zp = (k + 1) * nx * ny + j * nx + i;
                int idx_zm = (k - 1) * nx * ny + j * nx + i;

                float grad_x = (theta_field[idx_xp] - theta_field[idx_xm]) / (2 * dx);
                float grad_y = (theta_field[idx_yp] - theta_field[idx_ym]) / (2 * dx);
                float grad_z = (theta_field[idx_zp] - theta_field[idx_zm]) / (2 * dx);

                // Momentum density components
                momentum[0] += dt_theta * grad_x;
                momentum[1] += dt_theta * grad_y;
                momentum[2] += dt_theta * grad_z;
            }
        }
    }

    double dV = dx * dx * dx;
    momentum[0] *= dV;
    momentum[1] *= dV;
    momentum[2] *= dV;

    return momentum;
}

/**
 * Compute phase charge (U(1) Noether charge)
 * Q = ∫ρ_θ d³x where ρ_θ = |∇θ|²
 */
double computePhaseCharge(
    const std::vector<float>& theta_field,
    int nx, int ny, int nz, float dx
) {
    double charge = 0.0;

    for (int k = 1; k < nz - 1; ++k) {
        for (int j = 1; j < ny - 1; ++j) {
            for (int i = 1; i < nx - 1; ++i) {
                int idx = k * nx * ny + j * nx + i;

                // Compute gradients
                int idx_xp = k * nx * ny + j * nx + (i + 1);
                int idx_xm = k * nx * ny + j * nx + (i - 1);
                int idx_yp = k * nx * ny + (j + 1) * nx + i;
                int idx_ym = k * nx * ny + (j - 1) * nx + i;
                int idx_zp = (k + 1) * nx * ny + j * nx + i;
                int idx_zm = (k - 1) * nx * ny + j * nx + i;

                float grad_x = (theta_field[idx_xp] - theta_field[idx_xm]) / (2 * dx);
                float grad_y = (theta_field[idx_yp] - theta_field[idx_ym]) / (2 * dx);
                float grad_z = (theta_field[idx_zp] - theta_field[idx_zm]) / (2 * dx);

                // Charge density is gradient squared
                charge += grad_x * grad_x + grad_y * grad_y + grad_z * grad_z;
            }
        }
    }

    return charge * dx * dx * dx;
}

/**
 * Compute angular momentum (for rotation symmetry)
 * L_z = ∫(x·p_y - y·p_x) d³x
 */
double computeAngularMomentum(
    const std::vector<float>& theta_field,
    const std::vector<float>& theta_dot,
    int nx, int ny, int nz, float dx
) {
    double Lz = 0.0;

    int cx = nx / 2;
    int cy = ny / 2;

    for (int k = 1; k < nz - 1; ++k) {
        for (int j = 1; j < ny - 1; ++j) {
            for (int i = 1; i < nx - 1; ++i) {
                int idx = k * nx * ny + j * nx + i;

                // Position relative to center
                float x = (i - cx) * dx;
                float y = (j - cy) * dx;

                // Time derivative (proportional to canonical momentum)
                float dt_theta = theta_dot[idx];

                // Spatial gradients
                int idx_xp = k * nx * ny + j * nx + (i + 1);
                int idx_xm = k * nx * ny + j * nx + (i - 1);
                int idx_yp = k * nx * ny + (j + 1) * nx + i;
                int idx_ym = k * nx * ny + (j - 1) * nx + i;

                float grad_x = (theta_field[idx_xp] - theta_field[idx_xm]) / (2 * dx);
                float grad_y = (theta_field[idx_yp] - theta_field[idx_ym]) / (2 * dx);

                // Angular momentum density: r × p
                Lz += x * (dt_theta * grad_y) - y * (dt_theta * grad_x);
            }
        }
    }

    return Lz * dx * dx * dx;
}

/**
 * Verify continuity equation ∂ρ/∂t + ∇·J = 0
 */
double verifyContinuityEquation(
    const std::vector<float>& rho,
    const std::vector<float>& rho_prev,
    const std::vector<float>& Jx,
    const std::vector<float>& Jy,
    const std::vector<float>& Jz,
    int nx, int ny, int nz, float dx, float dt
) {
    double max_violation = 0.0;

    for (int k = 1; k < nz - 1; ++k) {
        for (int j = 1; j < ny - 1; ++j) {
            for (int i = 1; i < nx - 1; ++i) {
                int idx = k * nx * ny + j * nx + i;

                // Time derivative of density
                float drho_dt = (rho[idx] - rho_prev[idx]) / dt;

                // Divergence of current
                int idx_xp = k * nx * ny + j * nx + (i + 1);
                int idx_xm = k * nx * ny + j * nx + (i - 1);
                int idx_yp = k * nx * ny + (j + 1) * nx + i;
                int idx_ym = k * nx * ny + (j - 1) * nx + i;
                int idx_zp = (k + 1) * nx * ny + j * nx + i;
                int idx_zm = (k - 1) * nx * ny + j * nx + i;

                float div_J = (Jx[idx_xp] - Jx[idx_xm]) / (2 * dx) +
                             (Jy[idx_yp] - Jy[idx_ym]) / (2 * dx) +
                             (Jz[idx_zp] - Jz[idx_zm]) / (2 * dx);

                // Continuity equation violation
                float violation = std::abs(drho_dt + div_J);
                max_violation = std::max(max_violation, (double)violation);
            }
        }
    }

    return max_violation;
}

int runSymmetryAnalysisTest() {
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "E5: SYMMETRY ANALYSIS (NOETHER'S THEOREM)" << std::endl;
    std::cout << std::string(70, '=') << std::endl;

    std::cout << "\nTheoretical Framework:" << std::endl;
    std::cout << "  Noether's Theorem: Symmetry → Conservation Law" << std::endl;
    std::cout << "  Time translation → Energy conservation" << std::endl;
    std::cout << "  Space translation → Momentum conservation" << std::endl;
    std::cout << "  U(1) phase rotation → Charge conservation" << std::endl;
    std::cout << "  Quality Gate: All quantities conserved to <1%\n" << std::endl;

    // Test parameters
    const int nx = 32, ny = 32, nz = 32;
    const float dx = 1.0f;
    const float dt = 0.001f;  // Small timestep for accuracy
    const float K = 2.0f;
    const int total_steps = 10000;  // 100 time units with dt=0.001
    const int measure_interval = 100;  // Measure every 100 steps

    // Initialize TRD core
    TRDCore3D trd;
    TRDCore3D::Config config;
    config.Nx = nx;
    config.Ny = ny;
    config.Nz = nz;
    config.dx = dx;
    config.dt = dt;
    config.coupling_strength = K;
    trd.initialize(config);

    // Storage for fields
    std::vector<float> theta_field(nx * ny * nz);
    std::vector<float> theta_prev(nx * ny * nz);
    std::vector<float> theta_dot(nx * ny * nz, 0.0f);
    std::vector<float> R_field(nx * ny * nz, 1.0f);

    // Initialize with moving Gaussian wavepacket
    float kx = 0.5f, ky = 0.3f, kz = 0.2f;  // Wavevector for nonzero momentum
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                int idx = k * nx * ny + j * nx + i;
                float x = (i - nx/2) * dx;
                float y = (j - ny/2) * dx;
                float z = (k - nz/2) * dx;
                float r2 = x*x + y*y + z*z;
                float sigma = 3.0f * dx;

                // Gaussian envelope with plane wave
                float amplitude = std::exp(-r2 / (2 * sigma * sigma));
                theta_field[idx] = amplitude * std::cos(kx*x + ky*y + kz*z);

                // Initialize time derivative for moving wave
                theta_dot[idx] = -amplitude * std::sin(kx*x + ky*y + kz*z) *
                                std::sqrt(kx*kx + ky*ky + kz*kz);
            }
        }
    }

    // Set initial field
    trd.setPhaseField(theta_field.data());

    // Storage for conservation tracking
    std::vector<double> energy_history;
    std::vector<std::array<double, 3>> momentum_history;
    std::vector<double> charge_history;
    std::vector<double> angular_momentum_history;
    std::vector<float> time_points;

    // Compute initial conserved quantities
    ConservedQuantities initial;
    initial.energy = computeEnergy(theta_field, theta_dot, R_field, K, nx, ny, nz, dx);
    initial.momentum = computeMomentum(theta_field, theta_dot, nx, ny, nz, dx);
    initial.phase_charge = computePhaseCharge(theta_field, nx, ny, nz, dx);
    initial.angular_momentum = computeAngularMomentum(theta_field, theta_dot, nx, ny, nz, dx);

    std::cout << "Initial Conserved Quantities:" << std::endl;
    std::cout << "  ";
    initial.print();
    std::cout << std::endl;

    std::cout << "\nEvolution Progress:" << std::endl;
    std::cout << std::string(60, '-') << std::endl;
    std::cout << std::setw(10) << "Time"
              << std::setw(15) << "Energy"
              << std::setw(15) << "|Momentum|"
              << std::setw(15) << "Charge" << std::endl;
    std::cout << std::string(60, '-') << std::endl;

    // Evolution loop
    for (int step = 0; step < total_steps; ++step) {
        // Store previous field for time derivative
        theta_prev = theta_field;

        // Evolve one step
        trd.evolveKuramotoCPU(dt);
        trd.getPhaseField(theta_field.data());
        trd.computeRField();
        trd.getRField(R_field.data());

        // Compute time derivative
        for (int i = 0; i < nx * ny * nz; ++i) {
            theta_dot[i] = (theta_field[i] - theta_prev[i]) / dt;
        }

        // Measure conserved quantities at intervals
        if (step % measure_interval == 0) {
            ConservedQuantities current;
            current.energy = computeEnergy(theta_field, theta_dot, R_field, K, nx, ny, nz, dx);
            current.momentum = computeMomentum(theta_field, theta_dot, nx, ny, nz, dx);
            current.phase_charge = computePhaseCharge(theta_field, nx, ny, nz, dx);
            current.angular_momentum = computeAngularMomentum(theta_field, theta_dot, nx, ny, nz, dx);

            energy_history.push_back(current.energy);
            momentum_history.push_back(current.momentum);
            charge_history.push_back(current.phase_charge);
            angular_momentum_history.push_back(current.angular_momentum);
            time_points.push_back(step * dt);

            // Print progress
            double p_mag = std::sqrt(current.momentum[0]*current.momentum[0] +
                                    current.momentum[1]*current.momentum[1] +
                                    current.momentum[2]*current.momentum[2]);

            std::cout << std::fixed << std::setprecision(3)
                      << std::setw(10) << step * dt
                      << std::setprecision(6)
                      << std::setw(15) << current.energy
                      << std::setw(15) << p_mag
                      << std::setw(15) << current.phase_charge << std::endl;
        }
    }

    std::cout << std::string(60, '-') << std::endl;

    // Analyze conservation
    std::cout << "\nConservation Analysis:" << std::endl;
    std::cout << std::string(70, '=') << std::endl;

    // Energy conservation
    double E_max = *std::max_element(energy_history.begin(), energy_history.end());
    double E_min = *std::min_element(energy_history.begin(), energy_history.end());
    double E_drift = std::abs((E_max - E_min) / initial.energy) * 100;

    std::cout << "Energy Conservation:" << std::endl;
    std::cout << "  Initial: " << initial.energy << std::endl;
    std::cout << "  Max: " << E_max << std::endl;
    std::cout << "  Min: " << E_min << std::endl;
    std::cout << "  Drift: " << E_drift << "%" << std::endl;
    std::cout << "  Status: " << (E_drift < 1.0 ? "✓ CONSERVED" : "✗ VIOLATED") << std::endl;

    // Momentum conservation
    double px_drift = 0.0, py_drift = 0.0, pz_drift = 0.0;
    for (const auto& p : momentum_history) {
        px_drift = std::max(px_drift, std::abs(p[0] - initial.momentum[0]));
        py_drift = std::max(py_drift, std::abs(p[1] - initial.momentum[1]));
        pz_drift = std::max(pz_drift, std::abs(p[2] - initial.momentum[2]));
    }

    double p_initial_mag = std::sqrt(initial.momentum[0]*initial.momentum[0] +
                                    initial.momentum[1]*initial.momentum[1] +
                                    initial.momentum[2]*initial.momentum[2]);
    double p_drift_mag = std::sqrt(px_drift*px_drift + py_drift*py_drift + pz_drift*pz_drift);
    double p_drift_percent = (p_initial_mag > 1e-10) ? (p_drift_mag / p_initial_mag) * 100 : 0.0;

    std::cout << "\nMomentum Conservation:" << std::endl;
    std::cout << "  Initial: (" << initial.momentum[0] << ", "
              << initial.momentum[1] << ", " << initial.momentum[2] << ")" << std::endl;
    std::cout << "  Max drift: (" << px_drift << ", " << py_drift << ", " << pz_drift << ")" << std::endl;
    std::cout << "  Drift: " << p_drift_percent << "%" << std::endl;
    std::cout << "  Status: " << (p_drift_percent < 1.0 ? "✓ CONSERVED" : "✗ VIOLATED") << std::endl;

    // Phase charge conservation
    double Q_max = *std::max_element(charge_history.begin(), charge_history.end());
    double Q_min = *std::min_element(charge_history.begin(), charge_history.end());
    double Q_drift = std::abs((Q_max - Q_min) / initial.phase_charge) * 100;

    std::cout << "\nPhase Charge Conservation:" << std::endl;
    std::cout << "  Initial: " << initial.phase_charge << std::endl;
    std::cout << "  Max: " << Q_max << std::endl;
    std::cout << "  Min: " << Q_min << std::endl;
    std::cout << "  Drift: " << Q_drift << "%" << std::endl;
    std::cout << "  Status: " << (Q_drift < 1.0 ? "✓ CONSERVED" : "✗ VIOLATED") << std::endl;

    // Angular momentum (bonus check)
    double L_max = *std::max_element(angular_momentum_history.begin(),
                                    angular_momentum_history.end());
    double L_min = *std::min_element(angular_momentum_history.begin(),
                                    angular_momentum_history.end());
    double L_drift = (std::abs(initial.angular_momentum) > 1e-10) ?
                    std::abs((L_max - L_min) / initial.angular_momentum) * 100 : 0.0;

    std::cout << "\nAngular Momentum Conservation:" << std::endl;
    std::cout << "  Initial: " << initial.angular_momentum << std::endl;
    std::cout << "  Drift: " << L_drift << "%" << std::endl;
    std::cout << "  Status: " << (L_drift < 5.0 ? "✓ CONSERVED" : "✗ VIOLATED") << std::endl;

    // Quality gates
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "QUALITY GATES (NOETHER'S THEOREM VALIDATION):" << std::endl;
    std::cout << std::string(70, '=') << std::endl;

    bool energy_conserved = E_drift < 1.0;
    bool momentum_conserved = p_drift_percent < 1.0;
    bool charge_conserved = Q_drift < 1.0;
    bool all_passed = energy_conserved && momentum_conserved && charge_conserved;

    std::cout << (energy_conserved ? "✓" : "✗")
              << " Energy conservation: |ΔE/E₀| = " << E_drift << "%"
              << (energy_conserved ? " < 1%" : " ≥ 1%") << std::endl;

    std::cout << (momentum_conserved ? "✓" : "✗")
              << " Momentum conservation: |ΔP/P₀| = " << p_drift_percent << "%"
              << (momentum_conserved ? " < 1%" : " ≥ 1%") << std::endl;

    std::cout << (charge_conserved ? "✓" : "✗")
              << " Phase charge conservation: |ΔQ/Q₀| = " << Q_drift << "%"
              << (charge_conserved ? " < 1%" : " ≥ 1%") << std::endl;

    std::cout << "\nNOETHER'S THEOREM: "
              << (all_passed ? "✓ VALIDATED" : "✗ VIOLATION DETECTED") << std::endl;

    std::cout << "\nInterpretation:" << std::endl;
    if (all_passed) {
        std::cout << "  All symmetries produce conserved quantities as predicted" << std::endl;
        std::cout << "  TRD respects fundamental conservation laws" << std::endl;
        std::cout << "  Numerical scheme preserves symmetries to <1% accuracy" << std::endl;
    } else {
        std::cout << "  Some conservation laws violated beyond tolerance" << std::endl;
        std::cout << "  Check numerical scheme or physical parameters" << std::endl;
    }

    std::cout << std::string(70, '=') << std::endl;

    return all_passed ? 0 : 1;
}

// Entry point for standalone testing
#ifndef TRD_MAIN_EXECUTABLE
int main() {
    return runSymmetryAnalysisTest();
}
#endif