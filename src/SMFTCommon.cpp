/**
 * SMFTCommon.cpp
 *
 * Implementation of shared SMFT utilities
 */

#include "SMFTCommon.h"
#include <cmath>
#include <random>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>

namespace SMFT {

// ============================================================================
// KURAMOTO DYNAMICS
// ============================================================================

void stepKuramoto(std::vector<float>& theta,
                  const std::vector<float>& omega,
                  float dt, float K, float damping,
                  uint32_t Nx, uint32_t Ny) {

    uint32_t N = Nx * Ny;
    std::vector<float> dtheta(N, 0.0f);

    for (uint32_t y = 0; y < Ny; y++) {
        for (uint32_t x = 0; x < Nx; x++) {
            int i = periodicIdx(x, y, Nx, Ny);

            float rhs = omega[i];

            // 4-neighbor coupling
            rhs += (K/4.0f) * (
                std::sin(theta[periodicIdx(x+1, y, Nx, Ny)] - theta[i]) +
                std::sin(theta[periodicIdx(x-1, y, Nx, Ny)] - theta[i]) +
                std::sin(theta[periodicIdx(x, y+1, Nx, Ny)] - theta[i]) +
                std::sin(theta[periodicIdx(x, y-1, Nx, Ny)] - theta[i])
            );

            // Damping/forcing
            rhs -= damping * std::sin(theta[i]);

            dtheta[i] = rhs;
        }
    }

    // Euler update
    for (uint32_t i = 0; i < N; i++) {
        theta[i] += dt * dtheta[i];
    }
}

std::vector<float> computeLocalR(const std::vector<float>& theta,
                                  uint32_t Nx, uint32_t Ny) {

    uint32_t N = Nx * Ny;
    std::vector<float> R(N);

    for (uint32_t y = 0; y < Ny; y++) {
        for (uint32_t x = 0; x < Nx; x++) {
            std::complex<float> z(0.0f, 0.0f);

            // 3×3 neighborhood
            for (int dy = -1; dy <= 1; dy++) {
                for (int dx = -1; dx <= 1; dx++) {
                    int idx = periodicIdx(x + dx, y + dy, Nx, Ny);
                    z += std::complex<float>(std::cos(theta[idx]),
                                            std::sin(theta[idx]));
                }
            }

            z /= 9.0f;  // 9 neighbors (including self)
            R[y * Nx + x] = std::abs(z);
        }
    }

    return R;
}

float computeGlobalR(const std::vector<float>& theta) {
    std::complex<float> z(0.0f, 0.0f);

    for (float th : theta) {
        z += std::complex<float>(std::cos(th), std::sin(th));
    }

    z /= static_cast<float>(theta.size());
    return std::abs(z);
}

void stepKuramotoWithNoise(std::vector<float>& theta,
                           const std::vector<float>& omega,
                           float dt, float K, float damping,
                           float sigma,
                           uint32_t Nx, uint32_t Ny,
                           std::mt19937& rng) {

    uint32_t N = Nx * Ny;
    std::vector<float> theta_new(N);
    std::normal_distribution<float> noise(0.0f, 1.0f);

    for (uint32_t y = 0; y < Ny; y++) {
        for (uint32_t x = 0; x < Nx; x++) {
            uint32_t idx = y * Nx + x;

            // Coupling term (4-neighbor von Neumann)
            float coupling = 0.0f;
            coupling += std::sin(theta[periodicIdx(x+1, y, Nx, Ny)] - theta[idx]);
            coupling += std::sin(theta[periodicIdx(x-1, y, Nx, Ny)] - theta[idx]);
            coupling += std::sin(theta[periodicIdx(x, y+1, Nx, Ny)] - theta[idx]);
            coupling += std::sin(theta[periodicIdx(x, y-1, Nx, Ny)] - theta[idx]);

            // Damping term: -γ·sin(θ)
            float damping_force = -damping * std::sin(theta[idx]);

            // Deterministic drift
            float drift = omega[idx] + (K / 4.0f) * coupling + damping_force;

            // Stochastic term: σ·√(dt)·N(0,1) - Proper Euler-Maruyama scaling
            float noise_term = sigma * std::sqrt(dt) * noise(rng);

            // Update: θ(t+dt) = θ(t) + drift·dt + noise
            theta_new[idx] = theta[idx] + drift * dt + noise_term;
        }
    }

    theta = std::move(theta_new);
}

// ============================================================================
// DEFECT CREATION
// ============================================================================

uint32_t createDefect(std::vector<float>& theta,
                      std::vector<float>& omega,
                      const DefectConfig& config,
                      uint32_t Nx, uint32_t Ny,
                      uint32_t random_seed) {

    uint32_t N = Nx * Ny;
    theta.resize(N, 0.0f);
    omega.resize(N, config.omega_background);

    std::mt19937 gen(random_seed == 0 ? std::random_device{}() : random_seed);
    std::uniform_real_distribution<float> phase_dist(-M_PI, M_PI);

    uint32_t defect_count = 0;

    for (uint32_t y = 0; y < Ny; y++) {
        for (uint32_t x = 0; x < Nx; x++) {
            float dx = static_cast<float>(x) - config.center_x;
            float dy = static_cast<float>(y) - config.center_y;
            float r = std::sqrt(dx*dx + dy*dy);

            if (r < config.radius) {
                theta[y * Nx + x] = phase_dist(gen);
                omega[y * Nx + x] = config.omega_defect;
                defect_count++;
            }
        }
    }

    return defect_count;
}

float measureDefectContrast(const std::vector<float>& R,
                            uint32_t defect_x, uint32_t defect_y,
                            uint32_t Nx, uint32_t Ny) {

    // Background: sample far from defect
    float R_background = R[10 * Nx + 10];

    // Defect center
    float R_defect = R[defect_y * Nx + defect_x];

    return R_background - R_defect;
}

// ============================================================================
// DIAGNOSTICS
// ============================================================================

void computeGradient(const std::vector<float>& field,
                     float x, float y,
                     uint32_t Nx, uint32_t Ny,
                     float& grad_x, float& grad_y) {

    int ix = static_cast<int>(std::round(x));
    int iy = static_cast<int>(std::round(y));

    // Clamp to valid range for centered differences
    if (ix < 1) ix = 1;
    if (ix >= static_cast<int>(Nx)-1) ix = Nx - 2;
    if (iy < 1) iy = 1;
    if (iy >= static_cast<int>(Ny)-1) iy = Ny - 2;

    // Centered finite differences
    grad_x = 0.5f * (field[iy * Nx + (ix+1)] - field[iy * Nx + (ix-1)]);
    grad_y = 0.5f * (field[(iy+1) * Nx + ix] - field[(iy-1) * Nx + ix]);
}

float computeCoreDensity(const std::vector<float>& density,
                         float center_x, float center_y,
                         float radius,
                         uint32_t Nx, uint32_t Ny) {

    float rho_core = 0.0f;

    for (uint32_t y = 0; y < Ny; y++) {
        for (uint32_t x = 0; x < Nx; x++) {
            float dx = x - center_x;
            float dy = y - center_y;
            if (std::sqrt(dx*dx + dy*dy) < radius) {
                rho_core += density[y * Nx + x];
            }
        }
    }

    return rho_core;
}

float computeLocalization(const std::vector<float>& R_field) {
    float L = 0.0f;
    for (float R : R_field) {
        L += R * R * R * R;  // R⁴
    }
    return L;
}

std::vector<float> computeLocalRField(const std::vector<float>& theta,
                                      uint32_t Nx, uint32_t Ny,
                                      int radius) {
    uint32_t N = Nx * Ny;
    std::vector<float> R_field(N);

    for (uint32_t y = 0; y < Ny; y++) {
        for (uint32_t x = 0; x < Nx; x++) {
            uint32_t idx = y * Nx + x;

            float sum_cos = 0.0f;
            float sum_sin = 0.0f;
            int count = 0;

            // Average over neighborhood with given radius
            for (int dy = -radius; dy <= radius; dy++) {
                for (int dx = -radius; dx <= radius; dx++) {
                    int nx = periodicIdx(x + dx, y + dy, Nx, Ny);

                    sum_cos += std::cos(theta[nx]);
                    sum_sin += std::sin(theta[nx]);
                    count++;
                }
            }

            float avg_cos = sum_cos / count;
            float avg_sin = sum_sin / count;
            R_field[idx] = std::sqrt(avg_cos * avg_cos + avg_sin * avg_sin);
        }
    }

    return R_field;
}

// ============================================================================
// OUTPUT MANAGEMENT
// ============================================================================

DiagnosticWriter::DiagnosticWriter(const DiagnosticConfig& config)
    : _config(config) {

    // Create output directory
    mkdir(config.output_dir.c_str(), 0755);

    // Open files
    if (config.write_trajectory) {
        _trajectory_file.open(config.output_dir + "/trajectory.dat");
        _trajectory_file << "# step x y dist\n";
    }

    if (config.write_energy) {
        _energy_file.open(config.output_dir + "/energy.dat");
        _energy_file << "# step E_total KE PE <beta> norm\n";
    }

    if (config.write_force_alignment) {
        _force_file.open(config.output_dir + "/force_alignment.dat");
        _force_file << "# step alignment\n";
    }

    if (config.write_density) {
        _density_file.open(config.output_dir + "/core_density.dat");
        _density_file << "# step rho_5 rho_10 rho_15\n";
    }
}

DiagnosticWriter::~DiagnosticWriter() {
    flush();
}

void DiagnosticWriter::writeTrajectory(int step, float x, float y, float dist) {
    if (_trajectory_file.is_open()) {
        _trajectory_file << step << " " << x << " " << y << " " << dist << "\n";
    }
}

void DiagnosticWriter::writeEnergy(int step, float E_total, float KE, float PE,
                                   float beta, float norm) {
    if (_energy_file.is_open()) {
        _energy_file << step << " " << E_total << " " << KE << " "
                    << PE << " " << beta << " " << norm << "\n";
    }
}

void DiagnosticWriter::writeForceAlignment(int step, float alignment) {
    if (_force_file.is_open()) {
        _force_file << step << " " << alignment << "\n";
    }
}

void DiagnosticWriter::writeCoreDensity(int step, float rho_5, float rho_10, float rho_15) {
    if (_density_file.is_open()) {
        _density_file << step << " " << rho_5 << " " << rho_10 << " " << rho_15 << "\n";
    }
}

void DiagnosticWriter::flush() {
    if (_trajectory_file.is_open()) _trajectory_file.flush();
    if (_energy_file.is_open()) _energy_file.flush();
    if (_force_file.is_open()) _force_file.flush();
    if (_density_file.is_open()) _density_file.flush();
}

// ============================================================================
// SIMULATION PARAMETERS
// ============================================================================

SimulationParams::SimulationParams()
    : Nx(128), Ny(128), dt(0.01f), Delta(0.5f),
      K(1.0f), damping(0.1f),
      kuramoto_warmup(1000), coupled_steps(50000), output_interval(100),
      psi_x0(48.0f), psi_y0(64.0f), psi_sigma(5.0f) {

    // Defect defaults
    defect.center_x = 64;
    defect.center_y = 64;
    defect.radius = 15.0f;
    defect.omega_defect = 1.5f;
    defect.omega_background = 0.0f;

    // Diagnostic defaults
    diagnostics.output_dir = "output";
    diagnostics.output_interval = 100;
    diagnostics.write_trajectory = true;
    diagnostics.write_energy = true;
    diagnostics.write_density = false;
    diagnostics.write_force_alignment = false;
}

} // namespace SMFT
