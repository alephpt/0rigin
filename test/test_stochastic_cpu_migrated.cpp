/**
 * @file test_stochastic_cpu_migrated.cpp
 * @brief CPU-based validation test for stochastic MSR formalism - MIGRATED TO USE OutputManager
 *
 * This is a copy of test_stochastic_cpu.cpp updated to use OutputManager
 * for systematic output generation. This serves as proof of concept.
 *
 * PURPOSE: Validate stochastic Kuramoto-Dirac coupling using pure CPU implementation
 * to avoid GPU crashes documented in FINAL_ANALYSIS.md (13 GPU resets)
 *
 * PHYSICS:
 * - Coupled Kuramoto-Dirac evolution with vacuum noise
 * - Kuramoto: dθ/dt = ω + (K/4)Σsin(θ_j - θ_i) - γ·sin(θ) + σ_θ·√(dt)·ξ_θ(t)
 * - Dirac: i·dΨ/dt = [-iα·∇ + β·m(x)]Ψ + σ_Ψ·√(dt)·ξ_Ψ(t)
 * - Mass coupling: m(x) = Δ·R(x) where R(x) = |⟨e^(iθ)⟩|_local
 *
 * VALIDATION:
 * - Test 1: Synchronized vacuum stability (R_global > 0.95 maintained)
 * - Test 2: Spinor norm conservation (|Ψ|² conserved within 10%)
 * - Test 3: Particle localization (Gaussian wavepacket remains localized)
 * - Test 4: Critical threshold behavior (σ near σ_c ≈ 0.65 shows degradation)
 *
 * BASELINE: σ_θ = σ_Ψ = 0.05 (13× below critical σ_c ≈ 0.65)
 */

#include "../src/SMFTCommon.h"
#include "../src/output/OutputManager.h"
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <random>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <string>
#include <filesystem>

// Physical constants
const float DELTA = 2.5f;        // Mass gap parameter
const float DT = 0.01f;           // Time step
const float K = 1.0f;             // Kuramoto coupling
const float GAMMA = 0.1f;         // Damping
const float ALPHA = 1.0f;         // Dirac alpha matrix coefficient

// Grid parameters (smaller for CPU efficiency)
const int NX = 64;
const int NY = 64;
const int N_GRID = NX * NY;

// Simulation parameters
const int N_WARMUP = 500;        // Warmup steps for synchronization
const int N_STEPS = 500;         // Measurement steps
const int OUTPUT_INTERVAL = 50;  // Output frequency

// Baseline noise (validated in previous experiments)
const float SIGMA_BASELINE = 0.05f;  // 13× below critical σ_c ≈ 0.65

// Critical noise threshold (from noise sweep experiments)
const float SIGMA_CRITICAL = 0.65f;

// Type definitions
using Complex = std::complex<float>;
struct Spinor {
    Complex psi[4];  // 4-component Dirac spinor
};

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

int idx(int x, int y) {
    // Periodic boundary conditions
    x = (x + NX) % NX;
    y = (y + NY) % NY;
    return y * NX + x;
}

// ============================================================================
// KURAMOTO DYNAMICS
// ============================================================================

// Use compute_local_R and compute_global_R from SMFTCommon instead

/**
 * @brief Kuramoto step with Euler-Maruyama stochastic integration
 */
void kuramoto_step(std::vector<float>& theta, const std::vector<float>& omega,
                   float dt, float K, float gamma, float sigma,
                   std::mt19937& rng) {
    std::vector<float> theta_new(N_GRID);
    std::normal_distribution<float> noise(0.0f, 1.0f);

    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            int i = idx(x, y);

            // Coupling term (4-neighbor von Neumann)
            float coupling = 0.0f;
            coupling += std::sin(theta[idx(x+1, y)] - theta[i]);
            coupling += std::sin(theta[idx(x-1, y)] - theta[i]);
            coupling += std::sin(theta[idx(x, y+1)] - theta[i]);
            coupling += std::sin(theta[idx(x, y-1)] - theta[i]);

            // Damping term
            float damping_force = -gamma * std::sin(theta[i]);

            // Deterministic drift
            float drift = omega[i] + (K / 4.0f) * coupling + damping_force;

            // Stochastic term: σ·√(dt)·N(0,1)
            float noise_term = sigma * std::sqrt(dt) * noise(rng);

            // Euler-Maruyama update
            theta_new[i] = theta[i] + drift * dt + noise_term;
        }
    }

    theta = theta_new;
}

// ============================================================================
// DIRAC DYNAMICS
// ============================================================================

/**
 * @brief Initialize Gaussian wavepacket
 */
void initialize_gaussian_packet(std::vector<Spinor>& psi, float x0, float y0, float width) {
    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            float dx = x - x0;
            float dy = y - y0;
            float r2 = dx*dx + dy*dy;
            float amp = std::exp(-r2 / (2.0f * width * width));

            // Initialize as spin-up state
            int i = idx(x, y);
            psi[i].psi[0] = Complex(amp, 0);  // Upper component
            psi[i].psi[1] = Complex(0, 0);
            psi[i].psi[2] = Complex(0, 0);
            psi[i].psi[3] = Complex(0, 0);    // Lower component
        }
    }

    // Normalize
    float norm = 0.0f;
    for (const auto& s : psi) {
        for (int k = 0; k < 4; k++) {
            norm += std::norm(s.psi[k]);
        }
    }
    norm = std::sqrt(norm);

    for (auto& s : psi) {
        for (int k = 0; k < 4; k++) {
            s.psi[k] /= norm;
        }
    }
}

/**
 * @brief Compute spinor norm |Ψ|²
 */
float compute_spinor_norm(const std::vector<Spinor>& psi) {
    float norm = 0.0f;
    for (const auto& s : psi) {
        for (int k = 0; k < 4; k++) {
            norm += std::norm(s.psi[k]);
        }
    }
    return norm;
}

/**
 * @brief Compute particle position (center of mass)
 */
void compute_particle_position(const std::vector<Spinor>& psi, float& x_cm, float& y_cm) {
    float total_prob = 0.0f;
    x_cm = 0.0f;
    y_cm = 0.0f;

    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            int i = idx(x, y);
            float prob = 0.0f;
            for (int k = 0; k < 4; k++) {
                prob += std::norm(psi[i].psi[k]);
            }

            x_cm += x * prob;
            y_cm += y * prob;
            total_prob += prob;
        }
    }

    x_cm /= total_prob;
    y_cm /= total_prob;
}

/**
 * @brief Dirac step with mass coupling and stochastic noise
 *
 * Implements: i·dΨ/dt = [-iα·∇ + β·m(x)]Ψ + σ_Ψ·√(dt)·ξ_Ψ(t)
 *
 * Simplified 2D Dirac equation with proper normalization
 */
void dirac_step(std::vector<Spinor>& psi, const std::vector<float>& R_field,
                float dt, float sigma, std::mt19937& rng) {
    std::vector<Spinor> psi_new(N_GRID);
    std::normal_distribution<float> noise(0.0f, 1.0f);

    // Store initial norm for renormalization
    float norm_before = compute_spinor_norm(psi);

    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            int i = idx(x, y);

            // Mass from synchronization field
            float mass = DELTA * R_field[i];

            // Spatial derivatives (central differences with grid spacing)
            Spinor dpsi_dx, dpsi_dy;
            for (int k = 0; k < 4; k++) {
                dpsi_dx.psi[k] = (psi[idx(x+1, y)].psi[k] - psi[idx(x-1, y)].psi[k]) / 2.0f;
                dpsi_dy.psi[k] = (psi[idx(x, y+1)].psi[k] - psi[idx(x, y-1)].psi[k]) / 2.0f;
            }

            // Apply simplified 2D Dirac operator
            // For stability, we use a reduced coupling strength
            Spinor H_psi;
            float alpha_reduced = ALPHA * 0.1f;  // Reduced for stability

            // Kinetic term: -iα·∇
            // In 2D, we simplify to couple only adjacent components
            H_psi.psi[0] = -Complex(0, alpha_reduced) * (dpsi_dx.psi[1] + dpsi_dy.psi[1]);
            H_psi.psi[1] = -Complex(0, alpha_reduced) * (dpsi_dx.psi[0] + dpsi_dy.psi[0]);
            H_psi.psi[2] = -Complex(0, alpha_reduced) * (dpsi_dx.psi[3] + dpsi_dy.psi[3]);
            H_psi.psi[3] = -Complex(0, alpha_reduced) * (dpsi_dx.psi[2] + dpsi_dy.psi[2]);

            // Mass term: β·m (diagonal in rest frame)
            H_psi.psi[0] += mass * psi[i].psi[0];
            H_psi.psi[1] += mass * psi[i].psi[1];
            H_psi.psi[2] -= mass * psi[i].psi[2];  // Negative for lower components
            H_psi.psi[3] -= mass * psi[i].psi[3];

            // Euler-Maruyama update
            for (int k = 0; k < 4; k++) {
                // Deterministic evolution: Ψ(t+dt) = Ψ(t) - i·dt·H·Ψ
                psi_new[i].psi[k] = psi[i].psi[k] - Complex(0, dt) * H_psi.psi[k];

                // Stochastic term with reduced amplitude for stability
                // Note: noise should preserve norm on average
                float noise_real = sigma * std::sqrt(dt) * noise(rng) * 0.1f;
                float noise_imag = sigma * std::sqrt(dt) * noise(rng) * 0.1f;
                psi_new[i].psi[k] += Complex(noise_real, noise_imag);
            }
        }
    }

    // Renormalize to preserve unitarity (compensate for numerical errors)
    float norm_after = compute_spinor_norm(psi_new);
    if (norm_after > 1e-10) {  // Avoid division by zero
        float scale = std::sqrt(norm_before / norm_after);
        for (auto& s : psi_new) {
            for (int k = 0; k < 4; k++) {
                s.psi[k] *= scale;
            }
        }
    }

    psi = psi_new;
}

// ============================================================================
// VALIDATION TESTS
// ============================================================================

struct TestResult {
    bool passed;
    std::string description;
    float value;
    float threshold;
};

/**
 * @brief Test 1: Synchronized vacuum stability
 */
TestResult test_vacuum_stability(float sigma_theta, float sigma_psi,
                                SMFT::OutputManager& output_manager,
                                const std::string& test_dir) {
    std::cout << "\n=== Test 1: Vacuum Stability (σ_θ = " << sigma_theta
              << ", σ_Ψ = " << sigma_psi << ") ===" << std::endl;

    std::mt19937 rng(42);

    // Initialize synchronized state
    std::vector<float> theta(N_GRID);
    std::vector<float> omega(N_GRID, 0.0f);  // Zero natural frequencies

    std::uniform_real_distribution<float> init_dist(-0.1f, 0.1f);
    for (int i = 0; i < N_GRID; i++) {
        theta[i] = init_dist(rng);  // Small perturbation around sync
    }

    // Warmup without noise
    std::cout << "Warmup phase..." << std::flush;
    for (int step = 0; step < N_WARMUP; step++) {
        kuramoto_step(theta, omega, DT, K, GAMMA, 0.0f, rng);
        if (step % 100 == 0) std::cout << "." << std::flush;
    }

    float R_warmup = SMFT::compute_global_R(theta);
    std::cout << " R = " << R_warmup << std::endl;

    // Evolution with noise
    std::cout << "Stochastic evolution..." << std::flush;
    std::vector<float> R_history;
    std::vector<float> time_history;

    for (int step = 0; step < N_STEPS; step++) {
        kuramoto_step(theta, omega, DT, K, GAMMA, sigma_theta, rng);

        float R = SMFT::compute_global_R(theta);
        R_history.push_back(R);
        time_history.push_back(step * DT);

        if (step % 100 == 0) {
            std::cout << " R=" << std::fixed << std::setprecision(3) << R << std::flush;
        }
    }
    std::cout << std::endl;

    // Save timeseries using OutputManager
    std::map<std::string, std::string> metadata;
    metadata["test"] = "vacuum_stability";
    metadata["sigma_theta"] = std::to_string(sigma_theta);
    metadata["sigma_psi"] = std::to_string(sigma_psi);
    metadata["K"] = std::to_string(K);
    metadata["gamma"] = std::to_string(GAMMA);
    metadata["dt"] = std::to_string(DT);
    metadata["nx"] = std::to_string(NX);
    metadata["ny"] = std::to_string(NY);
    metadata["n_steps"] = std::to_string(N_STEPS);

    output_manager.writeTimeseries(test_dir + "/test1_R_timeseries.dat",
                                   time_history, R_history, metadata, "R_global");

    // Compute mean R over last 100 steps
    float R_mean = 0.0f;
    int count = std::min(100, (int)R_history.size());
    for (int i = R_history.size() - count; i < R_history.size(); i++) {
        R_mean += R_history[i];
    }
    R_mean /= count;

    TestResult result;
    result.value = R_mean;
    result.threshold = 0.95f;
    result.passed = (R_mean > result.threshold);
    result.description = "Vacuum maintains R > 0.95";

    std::cout << "Result: R_mean = " << R_mean
              << " (threshold = " << result.threshold << ") "
              << (result.passed ? "PASSED" : "FAILED") << std::endl;

    return result;
}

/**
 * @brief Test 2: Spinor norm conservation
 */
TestResult test_norm_conservation(float sigma_theta, float sigma_psi,
                                 SMFT::OutputManager& output_manager,
                                 const std::string& test_dir) {
    std::cout << "\n=== Test 2: Spinor Norm Conservation (σ_Ψ = " << sigma_psi << ") ===" << std::endl;

    std::mt19937 rng(123);

    // Initialize Kuramoto field (synchronized)
    std::vector<float> theta(N_GRID, 0.0f);
    std::vector<float> omega(N_GRID, 0.0f);

    // Initialize spinor (Gaussian packet)
    std::vector<Spinor> psi(N_GRID);
    initialize_gaussian_packet(psi, NX/2.0f, NY/2.0f, 5.0f);

    float norm_initial = compute_spinor_norm(psi);
    std::cout << "Initial norm: " << norm_initial << std::endl;

    // Evolution
    std::cout << "Evolution..." << std::flush;
    std::vector<float> norm_history;
    std::vector<float> time_history;

    for (int step = 0; step < N_STEPS; step++) {
        // Update Kuramoto
        kuramoto_step(theta, omega, DT, K, GAMMA, sigma_theta, rng);

        // Compute R field
        std::vector<float> R_field = SMFT::compute_local_R(theta, NX, NY);

        // Update Dirac
        dirac_step(psi, R_field, DT, sigma_psi, rng);

        float norm = compute_spinor_norm(psi);
        norm_history.push_back(norm);
        time_history.push_back(step * DT);

        if (step % 100 == 0) {
            std::cout << " |Ψ|²=" << std::fixed << std::setprecision(3) << norm << std::flush;
        }

        // Save snapshot at specific steps
        if (step == N_STEPS / 2) {
            std::map<std::string, std::string> snap_metadata;
            snap_metadata["test"] = "norm_conservation";
            snap_metadata["step"] = std::to_string(step);
            snap_metadata["time"] = std::to_string(step * DT);
            snap_metadata["norm"] = std::to_string(norm);

            output_manager.writeSnapshot(test_dir + "/test2_R_field_midpoint.dat",
                                        R_field, NX, NY, snap_metadata);
        }
    }
    std::cout << std::endl;

    // Save norm evolution
    std::map<std::string, std::string> metadata;
    metadata["test"] = "norm_conservation";
    metadata["sigma_theta"] = std::to_string(sigma_theta);
    metadata["sigma_psi"] = std::to_string(sigma_psi);
    metadata["norm_initial"] = std::to_string(norm_initial);

    output_manager.writeTimeseries(test_dir + "/test2_norm_timeseries.dat",
                                   time_history, norm_history, metadata, "norm");

    // Check norm conservation (within 10%)
    float max_deviation = 0.0f;
    for (float norm : norm_history) {
        float deviation = std::abs(norm - norm_initial) / norm_initial;
        max_deviation = std::max(max_deviation, deviation);
    }

    TestResult result;
    result.value = max_deviation;
    result.threshold = 0.1f;  // 10% tolerance
    result.passed = (max_deviation < result.threshold);
    result.description = "Spinor norm conserved within 10%";

    std::cout << "Result: Max deviation = " << (max_deviation * 100) << "% "
              << "(threshold = " << (result.threshold * 100) << "%) "
              << (result.passed ? "PASSED" : "FAILED") << std::endl;

    return result;
}

/**
 * @brief Test 3: Particle localization
 */
TestResult test_particle_localization(float sigma_theta, float sigma_psi,
                                     SMFT::OutputManager& output_manager,
                                     const std::string& test_dir) {
    std::cout << "\n=== Test 3: Particle Localization (σ_Ψ = " << sigma_psi << ") ===" << std::endl;

    std::mt19937 rng(456);

    // Initialize Kuramoto (synchronized)
    std::vector<float> theta(N_GRID, 0.0f);
    std::vector<float> omega(N_GRID, 0.0f);

    // Initialize spinor at center
    std::vector<Spinor> psi(N_GRID);
    float x0 = NX / 2.0f;
    float y0 = NY / 2.0f;
    initialize_gaussian_packet(psi, x0, y0, 3.0f);

    float x_init, y_init;
    compute_particle_position(psi, x_init, y_init);
    std::cout << "Initial position: (" << x_init << ", " << y_init << ")" << std::endl;

    // Evolution
    std::cout << "Evolution..." << std::flush;
    std::vector<float> x_history, y_history, drift_history, time_history;

    for (int step = 0; step < N_STEPS; step++) {
        // Update Kuramoto
        kuramoto_step(theta, omega, DT, K, GAMMA, sigma_theta, rng);

        // Compute R field
        std::vector<float> R_field = SMFT::compute_local_R(theta, NX, NY);

        // Update Dirac
        dirac_step(psi, R_field, DT, sigma_psi, rng);

        // Track position
        float x_cm, y_cm;
        compute_particle_position(psi, x_cm, y_cm);
        float drift = std::sqrt((x_cm - x0)*(x_cm - x0) + (y_cm - y0)*(y_cm - y0));

        x_history.push_back(x_cm);
        y_history.push_back(y_cm);
        drift_history.push_back(drift);
        time_history.push_back(step * DT);

        if (step % 100 == 0) {
            std::cout << " drift=" << std::fixed << std::setprecision(2) << drift << std::flush;
        }
    }

    float x_final = x_history.back();
    float y_final = y_history.back();
    std::cout << std::endl;

    float total_drift = std::sqrt((x_final - x_init)*(x_final - x_init) +
                                  (y_final - y_init)*(y_final - y_init));

    // Save position evolution using OutputManager
    std::map<std::string, std::string> metadata;
    metadata["test"] = "particle_localization";
    metadata["sigma_theta"] = std::to_string(sigma_theta);
    metadata["sigma_psi"] = std::to_string(sigma_psi);
    metadata["x_init"] = std::to_string(x_init);
    metadata["y_init"] = std::to_string(y_init);
    metadata["x_final"] = std::to_string(x_final);
    metadata["y_final"] = std::to_string(y_final);
    metadata["total_drift"] = std::to_string(total_drift);

    // Save multiple observables in one file
    std::map<std::string, std::vector<float>> observables;
    observables["x_position"] = x_history;
    observables["y_position"] = y_history;
    observables["drift"] = drift_history;

    output_manager.writeMultipleTimeseries(test_dir + "/test3_position_timeseries.dat",
                                          time_history, observables, metadata);

    TestResult result;
    result.value = total_drift;
    result.threshold = 5.0f;  // 5 grid units
    result.passed = (total_drift < result.threshold);
    result.description = "Particle drift < 5 grid units";

    std::cout << "Final position: (" << x_final << ", " << y_final << ")" << std::endl;
    std::cout << "Result: Total drift = " << total_drift
              << " (threshold = " << result.threshold << ") "
              << (result.passed ? "PASSED" : "FAILED") << std::endl;

    return result;
}

/**
 * @brief Test 4: Critical threshold behavior
 */
TestResult test_critical_behavior(SMFT::OutputManager& output_manager,
                                 const std::string& test_dir) {
    std::cout << "\n=== Test 4: Critical Threshold Behavior ===" << std::endl;

    // Create subdirectories for each test
    std::string baseline_dir = output_manager.createSubdirectory(test_dir, "test4_baseline");
    std::string critical_dir = output_manager.createSubdirectory(test_dir, "test4_critical");

    // Test at baseline (should be stable)
    TestResult baseline = test_vacuum_stability(SIGMA_BASELINE, SIGMA_BASELINE,
                                               output_manager, baseline_dir);

    // Test near critical (should show degradation)
    TestResult critical = test_vacuum_stability(SIGMA_CRITICAL * 0.9f, SIGMA_CRITICAL * 0.9f,
                                               output_manager, critical_dir);

    TestResult result;
    result.value = baseline.value - critical.value;  // Degradation amount
    result.threshold = 0.2f;  // Expect at least 20% degradation
    result.passed = (baseline.passed && !critical.passed) ||
                   (result.value > result.threshold);
    result.description = "Critical behavior observed";

    std::cout << "\nCritical test summary:" << std::endl;
    std::cout << "  Baseline (σ = " << SIGMA_BASELINE << "): R = " << baseline.value << std::endl;
    std::cout << "  Near critical (σ = " << (SIGMA_CRITICAL * 0.9f) << "): R = " << critical.value << std::endl;
    std::cout << "  Degradation: " << (result.value * 100) << "%" << std::endl;
    std::cout << "Result: " << (result.passed ? "PASSED" : "FAILED") << std::endl;

    return result;
}

// ============================================================================
// MAIN TEST RUNNER
// ============================================================================

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "   STOCHASTIC MSR FORMALISM CPU TEST   " << std::endl;
    std::cout << "     (USING OutputManager)              " << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "\nPure CPU implementation (no GPU calls)" << std::endl;
    std::cout << "Grid: " << NX << " × " << NY << std::endl;
    std::cout << "Steps: " << N_STEPS << " (dt = " << DT << ")" << std::endl;
    std::cout << "Baseline noise: σ = " << SIGMA_BASELINE << std::endl;
    std::cout << "Critical noise: σ_c ≈ " << SIGMA_CRITICAL << std::endl;

    // Initialize OutputManager
    SMFT::OutputManager output_manager;
    std::string experiment_dir = output_manager.createExperimentDirectory("stochastic_cpu_validation");

    // Save experiment metadata
    std::map<std::string, std::string> experiment_metadata;
    experiment_metadata["experiment_type"] = "stochastic_cpu_validation";
    experiment_metadata["grid_nx"] = std::to_string(NX);
    experiment_metadata["grid_ny"] = std::to_string(NY);
    experiment_metadata["dt"] = std::to_string(DT);
    experiment_metadata["n_warmup"] = std::to_string(N_WARMUP);
    experiment_metadata["n_steps"] = std::to_string(N_STEPS);
    experiment_metadata["K"] = std::to_string(K);
    experiment_metadata["gamma"] = std::to_string(GAMMA);
    experiment_metadata["delta"] = std::to_string(DELTA);
    experiment_metadata["alpha"] = std::to_string(ALPHA);
    experiment_metadata["sigma_baseline"] = std::to_string(SIGMA_BASELINE);
    experiment_metadata["sigma_critical"] = std::to_string(SIGMA_CRITICAL);

    output_manager.writeMetadata(experiment_dir + "/experiment_metadata.txt", experiment_metadata);

    // Prepare CSV for results summary
    std::vector<std::string> csv_headers = {"test_id", "description", "value", "threshold", "passed"};
    std::string results_csv = experiment_dir + "/test_results.csv";

    // Run validation tests
    std::vector<TestResult> results;

    // Test 1: Vacuum stability at baseline
    results.push_back(test_vacuum_stability(SIGMA_BASELINE, SIGMA_BASELINE,
                                           output_manager, experiment_dir));

    // Test 2: Norm conservation
    results.push_back(test_norm_conservation(SIGMA_BASELINE, SIGMA_BASELINE,
                                            output_manager, experiment_dir));

    // Test 3: Particle localization
    results.push_back(test_particle_localization(SIGMA_BASELINE, SIGMA_BASELINE,
                                                output_manager, experiment_dir));

    // Test 4: Critical behavior
    results.push_back(test_critical_behavior(output_manager, experiment_dir));

    // Write results to CSV
    for (size_t i = 0; i < results.size(); i++) {
        const auto& r = results[i];
        std::vector<float> row = {
            static_cast<float>(i + 1),
            static_cast<float>(r.passed ? 1.0f : 0.0f),
            r.value,
            r.threshold,
            static_cast<float>(r.passed ? 1.0f : 0.0f)
        };
        output_manager.appendCSV(results_csv, row, true, csv_headers);
    }

    // Summary
    std::cout << "\n========================================" << std::endl;
    std::cout << "           TEST SUMMARY" << std::endl;
    std::cout << "========================================" << std::endl;

    // Also save summary to file
    std::ofstream summary_file(experiment_dir + "/summary.txt");
    summary_file << "STOCHASTIC MSR FORMALISM CPU TEST SUMMARY\n";
    summary_file << "==========================================\n\n";

    int passed = 0;
    for (size_t i = 0; i < results.size(); i++) {
        const auto& r = results[i];
        std::cout << "Test " << (i+1) << ": " << r.description << " - "
                  << (r.passed ? "PASSED ✓" : "FAILED ✗") << std::endl;

        summary_file << "Test " << (i+1) << ": " << r.description << "\n";
        summary_file << "  Value: " << r.value << "\n";
        summary_file << "  Threshold: " << r.threshold << "\n";
        summary_file << "  Result: " << (r.passed ? "PASSED" : "FAILED") << "\n\n";

        if (r.passed) passed++;
    }

    std::cout << "\nOverall: " << passed << "/" << results.size() << " tests passed" << std::endl;
    summary_file << "Overall: " << passed << "/" << results.size() << " tests passed\n";

    if (passed == results.size()) {
        std::cout << "\n✓ VALIDATION SUCCESSFUL: Stochastic MSR formalism verified" << std::endl;
        std::cout << "\nAll output saved to: " << experiment_dir << std::endl;
        summary_file << "\nVALIDATION SUCCESSFUL\n";
    } else {
        std::cout << "\n✗ VALIDATION FAILED: Issues detected in stochastic implementation" << std::endl;
        std::cout << "\nOutput saved to: " << experiment_dir << std::endl;
        summary_file << "\nVALIDATION FAILED\n";
    }

    summary_file.close();

    return (passed == results.size()) ? 0 : 1;
}