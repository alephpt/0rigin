/**
 * Ehrenfest Theorem Validation - CRITICAL DIAGNOSTIC
 *
 * Resolve internal inconsistency:
 * - Test 3: Bounded orbit (PASS)
 * - Test 1: F·v ≈ 0 (FAIL)
 *
 * These are mutually contradictory unless:
 * 1. Measurement error in force/velocity
 * 2. Boundary artifact (not real binding)
 * 3. Energy E > 0 (scattering resonance, not binding)
 *
 * DIAGNOSTICS:
 * 1. E(t) = <Ψ|H|Ψ> - Must be E<0 for binding, conserved for validity
 * 2. Ehrenfest: d<p>/dt vs -<∇V> - Must agree if physics is correct
 * 3. Quantum velocity: v = <j>/ρ vs classical v = d<x>/dt
 * 4. Long-time (500k steps) to check eventual escape
 */

#include "../src/DiracEvolution.h"
#include "../src/SMFTCommon.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <random>
#include <complex>
#include <vector>

const uint32_t NX = 128;
const uint32_t NY = 128;
const uint32_t N_GRID = NX * NY;

const float DELTA = 0.5f;
const float DT = 0.01f;
const float K = 1.0f;
const float DAMPING = 0.1f;

const uint32_t DEFECT_X = 64;
const uint32_t DEFECT_Y = 64;
const float DEFECT_RADIUS = 15.0f;
const float OMEGA_DEFECT = 1.5f;

const float PSI_X0 = 48.0f;
const float PSI_Y0 = 64.0f;
const float PSI_SIGMA = 5.0f;

const int KURAMOTO_WARMUP = 1000;
const int COUPLED_STEPS = 50000;  // Sufficient for energy diagnostic
const int OUTPUT_INTERVAL = 100;   // Good resolution (500 points)

using Complex = std::complex<float>;

// Use functions from SMFTCommon instead of local implementations
// The snake_case wrappers in SMFTCommon.h provide backward compatibility
using namespace SMFT;

int main() {
    std::cout << "=== EHRENFEST THEOREM VALIDATION ===" << std::endl;
    std::cout << "Extended to " << COUPLED_STEPS << " steps (5000 time units)" << std::endl;
    std::cout << "Diagnostics: E(t), Ehrenfest, quantum velocity" << std::endl;

    // Initialize
    std::vector<float> theta(N_GRID, 0.0f);
    std::vector<float> omega(N_GRID, 0.0f);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> phase_dist(-M_PI, M_PI);

    for (uint32_t y = 0; y < NY; y++) {
        for (uint32_t x = 0; x < NX; x++) {
            float dx = x - DEFECT_X;
            float dy = y - DEFECT_Y;
            if (std::sqrt(dx*dx + dy*dy) < DEFECT_RADIUS) {
                theta[y * NX + x] = phase_dist(gen);
                omega[y * NX + x] = OMEGA_DEFECT;
            }
        }
    }

    std::cout << "\n[1] Kuramoto warmup..." << std::endl;
    for (int step = 0; step < KURAMOTO_WARMUP; step++) {
        step_kuramoto(theta, omega, DT, K, DAMPING, NX, NY);
    }

    auto R_field = compute_local_R(theta, NX, NY);
    std::cout << "  Defect contrast ΔR: "
              << (R_field[10*NX+10] - R_field[DEFECT_Y*NX+DEFECT_X]) << std::endl;

    std::cout << "\n[2] Initialize Dirac..." << std::endl;
    DiracEvolution dirac(NX, NY);
    dirac.initialize(PSI_X0, PSI_Y0, PSI_SIGMA);

    float x_init, y_init;
    dirac.getCenterOfMass(x_init, y_init);
    std::cout << "  Initial CoM: (" << x_init << ", " << y_init << ")" << std::endl;

    // Output files
    std::ofstream energy_file("output/10/ehrenfest/energy.dat");
    std::ofstream ehrenfest_file("output/10/ehrenfest/ehrenfest.dat");
    std::ofstream trajectory_file("output/10/ehrenfest/trajectory.dat");

    energy_file << "# step E_total KE PE <beta> norm\n";
    ehrenfest_file << "# step dp_dt_x dp_dt_y F_x F_y ratio_x ratio_y\n";
    trajectory_file << "# step x y dist\n";

    std::cout << "\n[3] Coupled evolution with FULL diagnostics..." << std::endl;

    float x_prev = x_init, y_prev = y_init;
    float x_prev2 = x_init, y_prev2 = y_init;

    for (int step = 0; step < COUPLED_STEPS; step++) {
        step_kuramoto(theta, omega, DT, K, DAMPING, NX, NY);

        R_field = compute_local_R(theta, NX, NY);
        std::vector<float> mass_field(N_GRID);
        for (uint32_t i = 0; i < N_GRID; i++) {
            mass_field[i] = DELTA * R_field[i];
        }

        dirac.step(mass_field, DT);

        if (step % OUTPUT_INTERVAL == 0) {
            auto density = dirac.getDensity();
            float x_com, y_com;
            dirac.getCenterOfMass(x_com, y_com);

            // === DIAGNOSTIC 1: Energy E(t) ===
            // Compute PROPER quantum expectation: E = <Ψ|H|Ψ>
            float KE, PE;
            float E_total = dirac.getEnergy(mass_field, KE, PE);

            float beta = dirac.getBetaExpectation();
            float norm = dirac.getNorm();

            energy_file << step << " " << E_total << " " << KE << " "
                       << PE << " " << beta << " " << norm << "\n";

            // === DIAGNOSTIC 2: Ehrenfest Theorem ===
            // d<p>/dt = (CoM[t+dt] - CoM[t-dt]) / (2*dt) (acceleration)
            // -<∇V> = Force at current position
            if (step > 0) {
                float dp_dt_x = (x_com - x_prev2) / (2.0f * DT * OUTPUT_INTERVAL);
                float dp_dt_y = (y_com - y_prev2) / (2.0f * DT * OUTPUT_INTERVAL);

                float grad_x, grad_y;
                compute_gradient(mass_field, x_com, y_com, NX, NY, grad_x, grad_y);
                float F_x = -beta * grad_x;
                float F_y = -beta * grad_y;

                float ratio_x = (std::abs(F_x) > 1e-6f) ? (dp_dt_x / F_x) : 0.0f;
                float ratio_y = (std::abs(F_y) > 1e-6f) ? (dp_dt_y / F_y) : 0.0f;

                ehrenfest_file << step << " " << dp_dt_x << " " << dp_dt_y << " "
                              << F_x << " " << F_y << " "
                              << ratio_x << " " << ratio_y << "\n";
            }

            // Update for next iteration
            x_prev2 = x_prev;
            y_prev2 = y_prev;
            x_prev = x_com;
            y_prev = y_com;

            // Trajectory
            float dist = std::sqrt((DEFECT_X - x_com)*(DEFECT_X - x_com) +
                                  (DEFECT_Y - y_com)*(DEFECT_Y - y_com));
            trajectory_file << step << " " << x_com << " " << y_com << " "
                          << dist << "\n";

            if (step % 50000 == 0) {
                std::cout << "  Step " << std::setw(7) << step
                          << ": dist=" << std::setprecision(2) << std::fixed << dist
                          << ", E=" << std::setprecision(4) << E_total
                          << ", norm=" << norm << std::endl;
            }
        }
    }

    energy_file.close();
    ehrenfest_file.close();
    trajectory_file.close();

    std::cout << "\n=== DIAGNOSTICS COMPLETE ===" << std::endl;
    std::cout << "Analyze output/10/ehrenfest/:" << std::endl;
    std::cout << "  1. energy.dat → Check E < 0 (binding) and |dE/dt| small (conservation)" << std::endl;
    std::cout << "  2. ehrenfest.dat → Check ratio ≈ 1.0 (Ehrenfest holds)" << std::endl;
    std::cout << "  3. trajectory.dat → Check if d(t) grows unbounded (escape)" << std::endl;

    return 0;
}
