/**
 * test_dirac_norm_scaling.cpp
 *
 * Test how norm drift scales with timestep size
 * to verify symplectic/unitary properties
 */

#include "Dirac3D.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

int main() {
    std::cout << "=== Dirac Norm Drift Scaling Test ===" << std::endl;
    std::cout << "Testing how norm drift scales with dt\n" << std::endl;

    // Small grid for fast testing
    const uint32_t N = 16;
    const float Delta = 0.5f;
    const float total_time = 10.0f;

    // Test different timesteps
    std::vector<float> dt_values = {0.01f, 0.005f, 0.002f, 0.001f};
    std::vector<float> norm_drifts;
    std::vector<int> step_counts;

    std::cout << std::setw(10) << "dt"
              << std::setw(10) << "steps"
              << std::setw(15) << "norm drift"
              << std::setw(15) << "drift/step"
              << std::setw(20) << "drift/(dt²)"
              << std::endl;
    std::cout << std::string(70, '-') << std::endl;

    for (float dt : dt_values) {
        int num_steps = static_cast<int>(total_time / dt);

        // Create fresh solver
        Dirac3D dirac(N, N, N);
        dirac.initializeGaussian(0.0f, 0.0f, 0.0f, 2.0f);

        // Create smooth vacuum fields
        std::vector<float> R_field(N*N*N);
        std::vector<float> theta_field(N*N*N);

        for (uint32_t k = 0; k < N; ++k) {
            for (uint32_t j = 0; j < N; ++j) {
                for (uint32_t i = 0; i < N; ++i) {
                    uint32_t idx = k*N*N + j*N + i;

                    float x = (float)i - N/2.0f;
                    float y = (float)j - N/2.0f;
                    float z = (float)k - N/2.0f;
                    float r2 = x*x + y*y + z*z;

                    R_field[idx] = 0.8f * std::exp(-r2/8.0f);
                    theta_field[idx] = 0.5f * std::atan2(y, x);
                }
            }
        }

        float norm_initial = dirac.getNorm();

        // Evolve
        for (int step = 0; step < num_steps; ++step) {
            dirac.stepWithChiralMass(R_field, theta_field, Delta, dt);
        }

        float norm_final = dirac.getNorm();
        float drift = std::abs(norm_final - norm_initial) / norm_initial;
        float drift_per_step = drift / num_steps;
        float drift_per_dt2 = drift / (dt * dt);

        norm_drifts.push_back(drift);
        step_counts.push_back(num_steps);

        std::cout << std::fixed << std::setprecision(4) << std::setw(10) << dt
                  << std::setw(10) << num_steps
                  << std::scientific << std::setprecision(4)
                  << std::setw(15) << drift
                  << std::setw(15) << drift_per_step
                  << std::setw(20) << drift_per_dt2
                  << std::endl;
    }

    std::cout << "\n=== Analysis ===" << std::endl;

    // Check scaling
    if (norm_drifts.size() >= 2) {
        float ratio = norm_drifts[0] / norm_drifts[1];
        float dt_ratio = dt_values[0] / dt_values[1];
        float expected_ratio = dt_ratio * dt_ratio; // For O(dt²) scaling

        std::cout << "Drift ratio (dt=" << dt_values[0] << " vs dt=" << dt_values[1] << "): "
                  << ratio << std::endl;
        std::cout << "Expected for O(dt²): " << expected_ratio << std::endl;
        std::cout << "Expected for O(dt): " << dt_ratio << std::endl;

        if (std::abs(ratio - expected_ratio) < 0.5 * expected_ratio) {
            std::cout << "✓ Scaling is approximately O(dt²) - good symplectic behavior" << std::endl;
        } else if (std::abs(ratio - dt_ratio) < 0.5 * dt_ratio) {
            std::cout << "⚠ Scaling is approximately O(dt) - first-order error" << std::endl;
        } else {
            std::cout << "✗ Scaling doesn't match O(dt²) or O(dt) - problematic" << std::endl;
        }
    }

    // Check if drift accumulates with steps
    bool accumulates = true;
    for (size_t i = 1; i < norm_drifts.size(); ++i) {
        float drift_per_step_i = norm_drifts[i] / step_counts[i];
        float drift_per_step_0 = norm_drifts[0] / step_counts[0];
        if (drift_per_step_i < drift_per_step_0 * 0.5) {
            accumulates = false;
            break;
        }
    }

    if (accumulates) {
        std::cout << "✗ Drift accumulates linearly with steps - non-unitary evolution!" << std::endl;
    } else {
        std::cout << "✓ Drift scales properly with dt - unitary evolution preserved" << std::endl;
    }

    // Overall assessment
    std::cout << "\n=== Verdict ===" << std::endl;
    if (norm_drifts[0] < 1e-3) { // Best case drift < 0.1%
        std::cout << "✓ PASS: Norm conservation meets TRD standards (<0.1%)" << std::endl;
    } else if (norm_drifts[0] < 1e-2) { // Drift < 1%
        std::cout << "⚠ MARGINAL: Norm conservation acceptable but could be improved" << std::endl;
    } else {
        std::cout << "✗ FAIL: Excessive norm drift (>1%)" << std::endl;
    }

    return 0;
}