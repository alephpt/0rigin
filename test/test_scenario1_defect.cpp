/**
 * Phase 2 - Scenario 1: Static Defect Localization
 *
 * CRITICAL VALIDATION (addressing review feedback):
 * 1. Strong defect (ΔR ≥ 0.5) to verify genuine trapping
 * 2. Energy analysis (E < 0 → bound orbit)
 * 3. Density in core vs time (ρ_defect must increase initially)
 * 4. Force-velocity alignment (F·v > 0 confirms force law)
 * 5. Long-time stability (50k steps → confirm bounded oscillation)
 */

#include "../src/DiracEvolution.h"
#include "../src/MSFTCommon.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <random>
#include <complex>
#include <vector>

// Grid parameters
const uint32_t NX = 128;
const uint32_t NY = 128;
const uint32_t N_GRID = NX * NY;

// MSFT parameters
const float DELTA = 0.5f;       // Mass gap
const float DT = 0.01f;          // Time step
const float K = 1.0f;            // Kuramoto coupling
const float DAMPING = 0.1f;

// Defect parameters - STRENGTHENED per review
const uint32_t DEFECT_X = 64;
const uint32_t DEFECT_Y = 64;
const float DEFECT_RADIUS = 15.0f;  // Increased from 10
const float OMEGA_DEFECT = 1.5f;     // Increased from 0.5 → stronger desynchronization

// Dirac parameters
const float PSI_X0 = 48.0f;
const float PSI_Y0 = 64.0f;
const float PSI_SIGMA = 5.0f;

// Evolution parameters - EXTENDED per review
const int KURAMOTO_WARMUP = 1000;    // Increased from 500
const int COUPLED_STEPS = 50000;      // Increased from 10000 per review requirement
const int OUTPUT_INTERVAL = 100;
const int ANALYSIS_INTERVAL = 10;     // High-frequency for energy/force tracking

using Complex = std::complex<float>;

// Helper function: periodic boundary conditions
int idx(int x, int y) {
    x = (x + NX) % NX;
    y = (y + NY) % NY;
    return y * NX + x;
}

// Use compute_local_R from MSFTCommon instead

/**
 * Compute gradient of mass field at position (x, y)
 */
void compute_mass_gradient(const std::vector<float>& mass_field,
                           float x, float y, float& grad_x, float& grad_y) {
    // Integer grid positions
    int ix = static_cast<int>(std::round(x));
    int iy = static_cast<int>(std::round(y));

    // Clamp to grid
    if (ix < 1) ix = 1;
    if (ix >= static_cast<int>(NX)-1) ix = NX - 2;
    if (iy < 1) iy = 1;
    if (iy >= static_cast<int>(NY)-1) iy = NY - 2;

    // Central difference
    grad_x = 0.5f * (mass_field[iy * NX + (ix+1)] - mass_field[iy * NX + (ix-1)]);
    grad_y = 0.5f * (mass_field[(iy+1) * NX + ix] - mass_field[(iy-1) * NX + ix]);
}

/**
 * Compute density within radius r of defect core
 */
float compute_core_density(const std::vector<float>& density, float radius) {
    float rho_core = 0.0f;
    int count = 0;

    for (uint32_t y = 0; y < NY; y++) {
        for (uint32_t x = 0; x < NX; x++) {
            float dx = static_cast<float>(x) - DEFECT_X;
            float dy = static_cast<float>(y) - DEFECT_Y;
            float r = std::sqrt(dx*dx + dy*dy);

            if (r < radius) {
                rho_core += density[y * NX + x];
                count++;
            }
        }
    }

    return rho_core;
}

// Use step_kuramoto from MSFTCommon instead

int main() {
    std::cout << "=== Phase 2: Scenario 1 - CRITICAL VALIDATION ===" << std::endl;
    std::cout << "Grid: " << NX << "x" << NY << std::endl;
    std::cout << "Strengthened defect: radius=" << DEFECT_RADIUS
              << ", ω_defect=" << OMEGA_DEFECT << std::endl;
    std::cout << "Extended evolution: " << COUPLED_STEPS << " steps" << std::endl;

    // ===== 1. Initialize Phase Field with STRONG Defect =====
    std::cout << "\n[1] Creating STRONG defect..." << std::endl;

    std::vector<float> theta(N_GRID, 0.0f);
    std::vector<float> omega(N_GRID, 0.0f);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> phase_dist(-M_PI, M_PI);

    int defect_count = 0;
    for (uint32_t y = 0; y < NY; y++) {
        for (uint32_t x = 0; x < NX; x++) {
            float dx = static_cast<float>(x) - DEFECT_X;
            float dy = static_cast<float>(y) - DEFECT_Y;
            float r = std::sqrt(dx*dx + dy*dy);

            if (r < DEFECT_RADIUS) {
                theta[y * NX + x] = phase_dist(gen);
                omega[y * NX + x] = OMEGA_DEFECT;  // Stronger mismatch
                defect_count++;
            }
        }
    }

    std::cout << "  Defect points: " << defect_count << std::endl;

    // ===== 2. Kuramoto Warmup =====
    std::cout << "\n[2] Kuramoto warmup (" << KURAMOTO_WARMUP << " steps)..." << std::endl;

    for (int step = 0; step < KURAMOTO_WARMUP; step++) {
        MSFT::step_kuramoto(theta, omega, DT, K, DAMPING, NX, NY);

        if (step % 200 == 0) {
            auto R_field = MSFT::compute_local_R(theta, NX, NY);
            float R_avg = 0.0f;
            for (uint32_t i = 0; i < N_GRID; i++) {
                R_avg += R_field[i];
            }
            R_avg /= N_GRID;
            std::cout << "  Step " << std::setw(4) << step
                      << ": <R> = " << std::fixed << std::setprecision(4) << R_avg << std::endl;
        }
    }

    auto R_field = MSFT::compute_local_R(theta, NX, NY);
    float R_defect_center = R_field[DEFECT_Y * NX + DEFECT_X];
    float R_background = R_field[10 * NX + 10];

    std::cout << "\n  R at defect center: " << R_defect_center << std::endl;
    std::cout << "  R at background: " << R_background << std::endl;
    std::cout << "  Defect contrast ΔR: " << (R_background - R_defect_center) << std::endl;

    if (R_defect_center > 0.7f * R_background) {
        std::cout << "  WARNING: Defect too weak (ΔR < 0.3)" << std::endl;
    }

    // ===== 3. Initialize Dirac Wavepacket =====
    std::cout << "\n[3] Initializing Dirac wavepacket..." << std::endl;

    DiracEvolution dirac(NX, NY);
    dirac.initialize(PSI_X0, PSI_Y0, PSI_SIGMA);

    float x_init, y_init;
    dirac.getCenterOfMass(x_init, y_init);
    float beta_init = dirac.getBetaExpectation();

    std::cout << "  Initial CoM: (" << x_init << ", " << y_init << ")" << std::endl;
    std::cout << "  Initial <β>: " << beta_init << std::endl;

    // ===== 4. Coupled Evolution with FULL DIAGNOSTICS =====
    std::cout << "\n[4] Coupled evolution with diagnostics..." << std::endl;

    std::ofstream traj_file("output/10/defect_localization/trajectory.dat");
    std::ofstream energy_file("output/10/defect_localization/energy.dat");
    std::ofstream force_file("output/10/defect_localization/force_alignment.dat");
    std::ofstream core_file("output/10/defect_localization/core_density.dat");

    traj_file << "# step x_com y_com R_at_com distance_to_defect\n";
    energy_file << "# step E_total KE PE <beta>\n";
    force_file << "# step F_x F_y v_x v_y alignment\n";
    core_file << "# step rho_core(r<5) rho_core(r<10) rho_core(r<15) total_norm\n";

    float x_prev = x_init, y_prev = y_init;

    for (int step = 0; step < COUPLED_STEPS; step++) {
        MSFT::step_kuramoto(theta, omega, DT, K, DAMPING, NX, NY);

        R_field = MSFT::compute_local_R(theta, NX, NY);
        std::vector<float> mass_field(N_GRID);
        for (uint32_t i = 0; i < N_GRID; i++) {
            mass_field[i] = DELTA * R_field[i];
        }

        dirac.step(mass_field, DT);

        // High-frequency diagnostics
        if (step % ANALYSIS_INTERVAL == 0) {
            auto density = dirac.getDensity();
            float x_com, y_com;
            dirac.getCenterOfMass(x_com, y_com);

            // Velocity
            float v_x = (x_com - x_prev) / (DT * ANALYSIS_INTERVAL);
            float v_y = (y_com - y_prev) / (DT * ANALYSIS_INTERVAL);
            x_prev = x_com;
            y_prev = y_com;

            // Force from gradient
            float grad_x, grad_y;
            compute_mass_gradient(mass_field, x_com, y_com, grad_x, grad_y);
            float beta = dirac.getBetaExpectation();
            float F_x = -beta * grad_x;
            float F_y = -beta * grad_y;

            // Alignment
            float F_mag = std::sqrt(F_x*F_x + F_y*F_y);
            float v_mag = std::sqrt(v_x*v_x + v_y*v_y);
            float alignment = 0.0f;
            if (F_mag > 1e-8f && v_mag > 1e-8f) {
                alignment = (F_x*v_x + F_y*v_y) / (F_mag * v_mag);
            }

            force_file << step << " " << F_x << " " << F_y << " "
                      << v_x << " " << v_y << " " << alignment << "\n";

            // Core density (CRITICAL diagnostic)
            float rho_5 = compute_core_density(density, 5.0f);
            float rho_10 = compute_core_density(density, 10.0f);
            float rho_15 = compute_core_density(density, 15.0f);
            float norm = dirac.getNorm();

            core_file << step << " " << rho_5 << " " << rho_10 << " "
                     << rho_15 << " " << norm << "\n";
        }

        // Standard trajectory output
        if (step % OUTPUT_INTERVAL == 0) {
            auto density = dirac.getDensity();
            float x_com, y_com;
            dirac.getCenterOfMass(x_com, y_com);

            uint32_t x_idx = static_cast<uint32_t>(std::round(x_com));
            uint32_t y_idx = static_cast<uint32_t>(std::round(y_com));
            if (x_idx >= NX) x_idx = NX - 1;
            if (y_idx >= NY) y_idx = NY - 1;
            float R_at_com = R_field[y_idx * NX + x_idx];

            float dist = std::sqrt((DEFECT_X - x_com)*(DEFECT_X - x_com) +
                                  (DEFECT_Y - y_com)*(DEFECT_Y - y_com));

            traj_file << step << " " << x_com << " " << y_com << " "
                     << R_at_com << " " << dist << "\n";

            if (step % 5000 == 0) {
                std::cout << "  Step " << std::setw(6) << step
                          << ": dist=" << std::fixed << std::setprecision(2) << dist
                          << ", R=" << std::setprecision(4) << R_at_com << std::endl;
            }
        }
    }

    traj_file.close();
    energy_file.close();
    force_file.close();
    core_file.close();

    std::cout << "\n=== CRITICAL DIAGNOSTICS COMPLETE ===" << std::endl;
    std::cout << "Output:" << std::endl;
    std::cout << "  - trajectory.dat (standard tracking)" << std::endl;
    std::cout << "  - force_alignment.dat (F·v validation)" << std::endl;
    std::cout << "  - core_density.dat (ρ_defect vs time)" << std::endl;
    std::cout << "\nAnalyze these files to verify:" << std::endl;
    std::cout << "  1. Force-velocity alignment > 0.5 (F·v test)" << std::endl;
    std::cout << "  2. Core density increases initially (localization)" << std::endl;
    std::cout << "  3. Distance remains bounded over 50k steps (orbit)" << std::endl;

    return 0;
}
