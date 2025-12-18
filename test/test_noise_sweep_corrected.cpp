/**
 * CORRECTED CPU Noise Sweep - Proper Noise Range
 * Based on diagnostic findings: σ_c ≈ 0.8
 * This sweep covers the actual transition region
 */

#include "../src/MSFTCommon.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <iomanip>
#include <filesystem>

using namespace MSFT;

float computeCircularVariance(const std::vector<float>& theta) {
    return 1.0f - computeGlobalR(theta);
}

int main() {
    std::cout << "=== CORRECTED CPU Noise Sweep ===" << std::endl;
    std::cout << "Covering actual transition region (σ = 10⁻⁵ to 3.0)" << std::endl << std::endl;

    const uint32_t Nx = 128;
    const uint32_t Ny = 128;
    const float dt = 0.01f;
    const float K = 1.0f;
    const float damping = 0.1f;
    const int warmup_steps = 2000;
    const int measurement_steps = 5000;

    std::mt19937 rng(42);

    // Sweep covering both regimes
    std::vector<float> sigma_values = {
        // Sub-threshold (for reference)
        1e-5f, 1e-4f, 1e-3f, 1e-2f,
        // Transition region
        0.05f, 0.1f, 0.15f, 0.2f, 0.25f, 0.3f, 0.35f, 0.4f,
        0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f,
        // Above threshold
        1.5f, 2.0f, 3.0f
    };

    std::filesystem::create_directories("/home/persist/neotec/0rigin/output/corrected_sweep");

    std::ofstream summary("/home/persist/neotec/0rigin/output/corrected_sweep/summary.dat");
    summary << "# CORRECTED CPU Noise Sweep Results\n";
    summary << "# Grid: " << Nx << "×" << Ny << ", dt = " << dt << ", K = " << K << ", γ = " << damping << "\n";
    summary << "# Warmup: " << warmup_steps << " steps, Measurement: " << measurement_steps << " steps\n";
    summary << "#\n";
    summary << "# sigma  R_initial  R_final  R_mean  R_std  circular_var\n";
    summary << std::fixed << std::setprecision(6);

    for (float sigma : sigma_values) {
        std::cout << "\nσ = " << std::scientific << sigma << std::flush;

        std::vector<float> theta(Nx * Ny);
        std::vector<float> omega(Nx * Ny, 0.0f);

        // Initialize near-synchronized
        std::uniform_real_distribution<float> init_dist(-0.1f, 0.1f);
        for (uint32_t i = 0; i < Nx * Ny; i++) {
            theta[i] = init_dist(rng);
        }

        float R_initial = computeGlobalR(theta);

        // Warmup
        for (int step = 0; step < warmup_steps; step++) {
            stepKuramotoWithNoise(theta, omega, dt, K, damping, 0.0f, Nx, Ny, rng);
        }
        float R_warmup = computeGlobalR(theta);

        // Measurement with noise
        std::vector<float> R_series(measurement_steps);
        for (int step = 0; step < measurement_steps; step++) {
            stepKuramotoWithNoise(theta, omega, dt, K, damping, sigma, Nx, Ny, rng);
            R_series[step] = computeGlobalR(theta);
        }

        // Statistics (last 50%)
        int steady_start = measurement_steps / 2;
        double R_sum = 0.0, R_sum_sq = 0.0;
        for (int i = steady_start; i < measurement_steps; i++) {
            R_sum += R_series[i];
            R_sum_sq += R_series[i] * R_series[i];
        }

        int N_samples = measurement_steps - steady_start;
        float R_mean = R_sum / N_samples;
        float R_std = std::sqrt(R_sum_sq / N_samples - R_mean * R_mean);
        float R_final = R_series.back();
        float circular_var = computeCircularVariance(theta);

        std::cout << " → R = " << std::fixed << std::setprecision(4) << R_final;

        summary << std::scientific << sigma << "  "
                << std::fixed << R_warmup << "  " << R_final << "  "
                << R_mean << "  " << R_std << "  " << circular_var << "\n";
        summary.flush();

        // Save timeseries for key points
        if (sigma == 1e-5f || sigma == 1e-4f || sigma == 0.1f || sigma == 0.5f ||
            sigma == 0.8f || sigma == 1.0f || sigma == 2.0f) {

            std::string filename = "/home/persist/neotec/0rigin/output/corrected_sweep/timeseries_" +
                                   std::to_string(sigma) + ".dat";
            std::ofstream ts(filename);
            ts << "# Timeseries for σ = " << sigma << "\n";
            ts << "# step  R_global\n";
            ts << std::fixed << std::setprecision(6);
            for (int i = 0; i < measurement_steps; i++) {
                ts << i << "  " << R_series[i] << "\n";
            }
            ts.close();
        }
    }

    summary.close();
    std::cout << "\n\n=== COMPLETE ===" << std::endl;
    std::cout << "Output: /home/persist/neotec/0rigin/output/corrected_sweep/" << std::endl;

    return 0;
}
