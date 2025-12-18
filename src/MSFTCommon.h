/**
 * MSFTCommon.h
 *
 * Shared utilities for MSFT (Mass Synchronization Field Theory) simulations
 * Eliminates code duplication across test files
 */

#pragma once

#include <vector>
#include <complex>
#include <cstdint>
#include <string>
#include <fstream>
#include <random>

namespace MSFT {

// ============================================================================
// KURAMOTO DYNAMICS (CPU)
// ============================================================================

/**
 * Periodic boundary index wrapping
 */
inline int periodicIdx(int x, int y, uint32_t Nx, uint32_t Ny) {
    x = (x + Nx) % Nx;
    y = (y + Ny) % Ny;
    return y * Nx + x;
}

/**
 * Single Kuramoto evolution step (4-neighbor coupling)
 *
 * dθ/dt = ω + (K/4)·Σ sin(θ_neighbor - θ) - γ·sin(θ)
 */
void stepKuramoto(std::vector<float>& theta,
                  const std::vector<float>& omega,
                  float dt, float K, float damping,
                  uint32_t Nx, uint32_t Ny);

/**
 * Compute local order parameter R(x,y) using 3×3 neighborhood
 *
 * R = |⟨exp(iθ)⟩| = |Σ_{neighbors} e^{iθ_j}| / N_neighbors
 */
std::vector<float> computeLocalR(const std::vector<float>& theta,
                                  uint32_t Nx, uint32_t Ny);

/**
 * Compute global order parameter R_global
 */
float computeGlobalR(const std::vector<float>& theta);

/**
 * Kuramoto evolution step with noise (Euler-Maruyama)
 *
 * dθ/dt = ω + (K/4)·Σ sin(θ_neighbor - θ) - γ·sin(θ) + σ·dW
 *
 * Uses proper Euler-Maruyama: θ(t+dt) = θ(t) + drift·dt + σ·√(dt)·N(0,1)
 */
void stepKuramotoWithNoise(std::vector<float>& theta,
                           const std::vector<float>& omega,
                           float dt, float K, float damping,
                           float sigma,
                           uint32_t Nx, uint32_t Ny,
                           std::mt19937& rng);

// ============================================================================
// DEFECT CREATION
// ============================================================================

struct DefectConfig {
    uint32_t center_x;
    uint32_t center_y;
    float radius;
    float omega_defect;  // Natural frequency inside defect
    float omega_background;  // Natural frequency outside defect (default 0)

    DefectConfig()
        : center_x(64), center_y(64), radius(15.0f),
          omega_defect(1.5f), omega_background(0.0f) {}
};

/**
 * Initialize defect via natural frequency heterogeneity
 * Returns number of defect points
 */
uint32_t createDefect(std::vector<float>& theta,
                      std::vector<float>& omega,
                      const DefectConfig& config,
                      uint32_t Nx, uint32_t Ny,
                      uint32_t random_seed = 0);

/**
 * Measure defect contrast ΔR = R_background - R_defect
 */
float measureDefectContrast(const std::vector<float>& R,
                            uint32_t defect_x, uint32_t defect_y,
                            uint32_t Nx, uint32_t Ny);

// ============================================================================
// DIAGNOSTICS
// ============================================================================

struct DiagnosticConfig {
    std::string output_dir;
    int output_interval;
    bool write_trajectory;
    bool write_energy;
    bool write_density;
    bool write_force_alignment;

    DiagnosticConfig()
        : output_dir("output"), output_interval(100),
          write_trajectory(true), write_energy(true),
          write_density(false), write_force_alignment(false) {}
};

/**
 * Compute gradient of a scalar field at position (x,y)
 * Uses centered finite differences
 */
void computeGradient(const std::vector<float>& field,
                     float x, float y,
                     uint32_t Nx, uint32_t Ny,
                     float& grad_x, float& grad_y);

/**
 * Compute integrated density within radius r of center
 */
float computeCoreDensity(const std::vector<float>& density,
                         float center_x, float center_y,
                         float radius,
                         uint32_t Nx, uint32_t Ny);

/**
 * Compute localization metric L = ∫ R⁴ dA
 * Higher L indicates more localized synchronization patterns
 */
float computeLocalization(const std::vector<float>& R_field);

/**
 * Compute local R field with configurable neighborhood radius
 * Similar to computeLocalR but with adjustable radius (default=1)
 */
std::vector<float> computeLocalRField(const std::vector<float>& theta,
                                      uint32_t Nx, uint32_t Ny,
                                      int radius = 1);

// ============================================================================
// OUTPUT MANAGEMENT
// ============================================================================

/**
 * File handle manager for diagnostics
 */
class DiagnosticWriter {
public:
    DiagnosticWriter(const DiagnosticConfig& config);
    ~DiagnosticWriter();

    void writeTrajectory(int step, float x, float y, float dist);
    void writeEnergy(int step, float E_total, float KE, float PE,
                     float beta, float norm);
    void writeForceAlignment(int step, float alignment);
    void writeCoreDensity(int step, float rho_5, float rho_10, float rho_15);

    void flush();

private:
    DiagnosticConfig _config;
    std::ofstream _trajectory_file;
    std::ofstream _energy_file;
    std::ofstream _force_file;
    std::ofstream _density_file;
};

// ============================================================================
// SIMULATION PARAMETERS
// ============================================================================

struct SimulationParams {
    // Grid
    uint32_t Nx;
    uint32_t Ny;

    // Timestep
    float dt;

    // MSFT coupling
    float Delta;  // Coupling strength: m = Δ·R

    // Kuramoto parameters
    float K;       // Coupling strength
    float damping; // Damping/forcing: γ

    // Evolution
    int kuramoto_warmup;
    int coupled_steps;
    int output_interval;

    // Dirac initialization
    float psi_x0;
    float psi_y0;
    float psi_sigma;

    // Defect
    DefectConfig defect;

    // Output
    DiagnosticConfig diagnostics;

    // Default constructor with standard parameters
    SimulationParams();
};

// ============================================================================
// BACKWARD COMPATIBILITY: snake_case wrappers
// ============================================================================
// These inline wrappers allow test files to use snake_case naming
// while the main API uses camelCase

/**
 * Wrapper for computeLocalR - compute local order parameter with 3×3 neighborhood
 */
inline std::vector<float> compute_local_R(const std::vector<float>& theta,
                                          uint32_t Nx, uint32_t Ny) {
    return computeLocalR(theta, Nx, Ny);
}

/**
 * Wrapper for computeGlobalR - compute global order parameter
 */
inline float compute_global_R(const std::vector<float>& theta) {
    return computeGlobalR(theta);
}

/**
 * Wrapper for stepKuramoto - single evolution step with 4-neighbor coupling
 */
inline void step_kuramoto(std::vector<float>& theta,
                          const std::vector<float>& omega,
                          float dt, float K, float damping,
                          uint32_t Nx, uint32_t Ny) {
    stepKuramoto(theta, omega, dt, K, damping, Nx, Ny);
}

/**
 * Wrapper for computeGradient - centered finite differences
 */
inline void compute_gradient(const std::vector<float>& field,
                            float x, float y,
                            uint32_t Nx, uint32_t Ny,
                            float& grad_x, float& grad_y) {
    computeGradient(field, x, y, Nx, Ny, grad_x, grad_y);
}

/**
 * Wrapper for computeCoreDensity - integrated density within radius
 */
inline float compute_core_density(const std::vector<float>& density,
                                  float center_x, float center_y,
                                  float radius,
                                  uint32_t Nx, uint32_t Ny) {
    return computeCoreDensity(density, center_x, center_y, radius, Nx, Ny);
}

/**
 * Wrapper for createDefect - initialize defect via natural frequency heterogeneity
 */
inline uint32_t create_defect(std::vector<float>& theta,
                              std::vector<float>& omega,
                              const DefectConfig& config,
                              uint32_t Nx, uint32_t Ny,
                              uint32_t random_seed = 0) {
    return createDefect(theta, omega, config, Nx, Ny, random_seed);
}

/**
 * Wrapper for measureDefectContrast - measure ΔR = R_background - R_defect
 */
inline float measure_defect_contrast(const std::vector<float>& R,
                                     uint32_t defect_x, uint32_t defect_y,
                                     uint32_t Nx, uint32_t Ny) {
    return measureDefectContrast(R, defect_x, defect_y, Nx, Ny);
}

/**
 * Wrapper for stepKuramotoWithNoise - stochastic evolution step
 */
inline void step_kuramoto_with_noise(std::vector<float>& theta,
                                     const std::vector<float>& omega,
                                     float dt, float K, float damping,
                                     float sigma,
                                     uint32_t Nx, uint32_t Ny,
                                     std::mt19937& rng) {
    stepKuramotoWithNoise(theta, omega, dt, K, damping, sigma, Nx, Ny, rng);
}

/**
 * Wrapper for computeLocalization - localization metric
 */
inline float compute_localization(const std::vector<float>& R_field) {
    return computeLocalization(R_field);
}

/**
 * Wrapper for computeLocalRField - local R field with radius
 */
inline std::vector<float> compute_local_R_field(const std::vector<float>& theta,
                                                uint32_t Nx, uint32_t Ny,
                                                int radius = 1) {
    return computeLocalRField(theta, Nx, Ny, radius);
}

} // namespace MSFT
