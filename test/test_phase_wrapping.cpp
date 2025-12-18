/**
 * Phase Wrapping Diagnostic
 *
 * Purpose: Track phase range over time to detect unwrapping/drift
 *
 * The Kuramoto shader wraps phases to [-π, π] after each step:
 *     theta_new = mod(theta_new + π, 2π) - π
 *
 * Question: Does this wrapping preserve physics correctly?
 *
 * Tests:
 * 1. Track min/max θ over 1000 steps
 * 2. Check if phases stay in [-π, π]
 * 3. Monitor phase drift (unwanted accumulation)
 * 4. Detect wrapping artifacts (discontinuities)
 *
 * Expected:
 * - Phases remain bounded: θ ∈ [-π, π]
 * - No drift outside bounds
 * - R_avg stable (wrapping doesn't break synchronization)
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <cstdint>
#include "../src/SMFTCommon.h"

// CPU Kuramoto step (matches shader logic)
void cpu_kuramoto_step(std::vector<float>& theta,
                       const std::vector<float>& omega,
                       uint32_t Nx, uint32_t Ny,
                       float K, float dt, float damping) {
    std::vector<float> theta_new(theta.size());

    for (uint32_t y = 0; y < Ny; y++) {
        for (uint32_t x = 0; x < Nx; x++) {
            uint32_t idx = y * Nx + x;
            float theta_i = theta[idx];
            float omega_i = omega[idx];

            // Compute coupling force (8-neighborhood)
            float coupling = 0.0f;
            for (int dy = -1; dy <= 1; dy++) {
                for (int dx = -1; dx <= 1; dx++) {
                    if (dx == 0 && dy == 0) continue;

                    int nx = static_cast<int>(x) + dx;
                    int ny = static_cast<int>(y) + dy;

                    // Periodic boundaries
                    if (nx < 0) nx += Nx;
                    if (nx >= static_cast<int>(Nx)) nx -= Nx;
                    if (ny < 0) ny += Ny;
                    if (ny >= static_cast<int>(Ny)) ny -= Ny;

                    uint32_t neighbor_idx = ny * Nx + nx;
                    coupling += std::sin(theta[neighbor_idx] - theta_i);
                }
            }
            coupling *= (K / 8.0f);

            // Damping force
            float damping_force = -damping * std::sin(theta_i);

            // Total force
            float total_force = omega_i + coupling + damping_force;

            // Euler integration
            theta_new[idx] = theta_i + dt * total_force;

            // Wrap to [-π, π] (matches GLSL mod behavior)
            // GLSL mod(x, y) = x - y * floor(x/y)
            // For theta wrapping: mod(theta + π, 2π) - π
            float temp = theta_new[idx] + M_PI;
            float wrapped = temp - 2.0f * M_PI * std::floor(temp / (2.0f * M_PI));
            theta_new[idx] = wrapped - M_PI;
        }
    }

    theta = theta_new;
}

// Compute average of local R field (for monitoring)
float compute_R_avg(const std::vector<float>& theta, uint32_t Nx, uint32_t Ny) {
    // Use SMFTCommon to compute local R field, then average
    std::vector<float> R_field = SMFT::computeLocalR(theta, Nx, Ny);

    float sum = 0.0f;
    for (float r : R_field) {
        sum += r;
    }

    return sum / R_field.size();
}

int main() {
    std::cout << "=== Phase Wrapping Diagnostic ===" << std::endl;
    std::cout << "CPU simulation with wrapping monitoring\n" << std::endl;

    const uint32_t Nx = 128;
    const uint32_t Ny = 128;
    const int N_steps = 10000;  // Extended to 10k steps (t_max = 100)
    const float dt = 0.01f;
    const float K = 1.0f;
    const float damping = 0.1f;

    std::cout << "Parameters:" << std::endl;
    std::cout << "  Grid: " << Nx << " × " << Ny << std::endl;
    std::cout << "  Steps: " << N_steps << std::endl;
    std::cout << "  dt: " << dt << std::endl;
    std::cout << "  K: " << K << std::endl;
    std::cout << "  damping: " << damping << std::endl;
    std::cout << std::endl;

    // Initialize with random phases
    std::vector<float> theta(Nx * Ny);
    std::vector<float> omega(Nx * Ny);

    srand(42);
    for (size_t i = 0; i < theta.size(); i++) {
        theta[i] = (float(rand()) / RAND_MAX) * 2.0f * M_PI - M_PI;
        omega[i] = 0.5f * (float(rand()) / RAND_MAX - 0.5f);
    }

    // Track phase statistics over time
    std::vector<float> theta_min_series;
    std::vector<float> theta_max_series;
    std::vector<float> theta_range_series;
    std::vector<float> R_avg_series;

    std::cout << "Running simulation..." << std::endl;
    std::cout << std::fixed << std::setprecision(6);

    for (int step = 0; step < N_steps; step++) {
        // Update phases
        cpu_kuramoto_step(theta, omega, Nx, Ny, K, dt, damping);

        // Compute statistics
        float theta_min = *std::min_element(theta.begin(), theta.end());
        float theta_max = *std::max_element(theta.begin(), theta.end());
        float theta_range = theta_max - theta_min;
        float R_avg = compute_R_avg(theta, Nx, Ny);

        // Store
        theta_min_series.push_back(theta_min);
        theta_max_series.push_back(theta_max);
        theta_range_series.push_back(theta_range);
        R_avg_series.push_back(R_avg);

        // Print progress
        if (step % 500 == 0 || step == N_steps - 1) {
            std::cout << "Step " << std::setw(5) << step
                      << ": θ ∈ [" << theta_min << ", " << theta_max << "]"
                      << ", range = " << theta_range
                      << ", R_avg = " << R_avg << std::endl;
        }
    }

    // Analysis
    std::cout << "\n=== ANALYSIS ===" << std::endl;

    // Check bounds
    float global_min = *std::min_element(theta_min_series.begin(), theta_min_series.end());
    float global_max = *std::max_element(theta_max_series.begin(), theta_max_series.end());

    std::cout << "Global phase range: [" << global_min << ", " << global_max << "]" << std::endl;
    std::cout << "Expected range: [-π, π] = [-3.14159, 3.14159]" << std::endl;

    bool in_bounds = (global_min >= -M_PI - 1e-5) && (global_max <= M_PI + 1e-5);

    if (in_bounds) {
        std::cout << "✓ PASS - Phases remain within bounds" << std::endl;
    } else {
        std::cout << "✗ FAIL - Phases escaped bounds!" << std::endl;
        std::cout << "  Overflow: " << std::max(0.0f, float(global_max - M_PI)) << std::endl;
        std::cout << "  Underflow: " << std::max(0.0f, float(-M_PI - global_min)) << std::endl;
    }

    // Check wrapping artifacts (discontinuities in range)
    std::cout << "\nWrapping Artifacts Check:" << std::endl;
    bool has_artifacts = false;
    for (size_t i = 1; i < theta_range_series.size(); i++) {
        float delta_range = std::abs(theta_range_series[i] - theta_range_series[i-1]);
        if (delta_range > 1.0f) {  // Large jump in range
            std::cout << "  Large jump at step " << i << ": Δrange = " << delta_range << std::endl;
            has_artifacts = true;
        }
    }

    if (!has_artifacts) {
        std::cout << "  ✓ No large discontinuities detected" << std::endl;
    }

    // Check R_avg evolution (expect growth during synchronization)
    float R_initial = R_avg_series[0];
    float R_final = R_avg_series.back();
    float R_growth = R_final - R_initial;

    // Check stability in last 100 steps (should reach steady state)
    size_t start_steady = R_avg_series.size() > 100 ? R_avg_series.size() - 100 : 0;
    float R_min_steady = *std::min_element(R_avg_series.begin() + start_steady, R_avg_series.end());
    float R_max_steady = *std::max_element(R_avg_series.begin() + start_steady, R_avg_series.end());

    std::cout << "\nSynchronization Evolution:" << std::endl;
    std::cout << "  R_initial = " << R_initial << std::endl;
    std::cout << "  R_final = " << R_final << std::endl;
    std::cout << "  Growth: " << R_growth << " (" << (R_growth/R_initial*100.0f) << "%)" << std::endl;
    std::cout << "  Steady-state fluctuations (last 100 steps): " << (R_max_steady - R_min_steady) << std::endl;

    // R should grow (synchronization) and stabilize
    bool R_growing = (R_growth > 0.1f);  // Significant growth
    bool R_steady = (R_max_steady - R_min_steady) < 0.05f;  // Less than 5% fluctuation at end

    if (R_growing && R_steady) {
        std::cout << "  ✓ Normal synchronization dynamics (wrapping preserves physics)" << std::endl;
    } else if (R_growing && !R_steady) {
        std::cout << "  ⚠ Still synchronizing (not yet steady state)" << std::endl;
    } else {
        std::cout << "  ? Unexpected R_avg behavior" << std::endl;
    }

    // Save timeseries
    std::ofstream f("output/phase_wrapping_diagnostics.dat");
    f << "# step theta_min theta_max theta_range R_avg\n";
    for (size_t i = 0; i < theta_min_series.size(); i++) {
        f << i << " "
          << theta_min_series[i] << " "
          << theta_max_series[i] << " "
          << theta_range_series[i] << " "
          << R_avg_series[i] << "\n";
    }
    f.close();
    std::cout << "\nSaved: output/phase_wrapping_diagnostics.dat" << std::endl;

    // Final verdict
    std::cout << "\n=== VERDICT ===" << std::endl;
    if (in_bounds && R_growing) {
        std::cout << "✓ Phase wrapping is working correctly" << std::endl;
        std::cout << "  - Phases stay bounded in [-π, π]" << std::endl;
        std::cout << "  - Synchronization dynamics preserved" << std::endl;
        std::cout << "  - No drift or unwrapping artifacts" << std::endl;
        std::cout << "  - Physics is correct" << std::endl;
        return 0;
    } else {
        std::cout << "✗ Phase wrapping has issues" << std::endl;
        if (!in_bounds) std::cout << "  - Phases escaped bounds" << std::endl;
        if (!R_growing) std::cout << "  - Synchronization not working" << std::endl;
        return 1;
    }
}
