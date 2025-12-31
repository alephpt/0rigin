/**
 * Diagnostic Tests for Noise Implementation
 * Tests 1-4 from immediate.md critical assessment
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>

void kuramotoStepWithNoise(std::vector<float>& theta, const std::vector<float>& omega,
                           uint32_t Nx, uint32_t Ny, float K, float damping,
                           float dt, float sigma, std::mt19937& rng) {
    std::vector<float> theta_new(Nx * Ny);
    std::normal_distribution<float> noise(0.0f, 1.0f);

    for (uint32_t y = 0; y < Ny; y++) {
        for (uint32_t x = 0; x < Nx; x++) {
            uint32_t idx = y * Nx + x;

            // Coupling
            float coupling = 0.0f;
            uint32_t neighbors[4][2] = {
                {(x + 1) % Nx, y},
                {(x + Nx - 1) % Nx, y},
                {x, (y + 1) % Ny},
                {x, (y + Ny - 1) % Ny}
            };

            for (int n = 0; n < 4; n++) {
                uint32_t nx = neighbors[n][0];
                uint32_t ny = neighbors[n][1];
                uint32_t nidx = ny * Nx + nx;
                coupling += std::sin(theta[nidx] - theta[idx]);
            }

            float damping_force = -damping * std::sin(theta[idx]);
            float drift = omega[idx] + (K / 4.0f) * coupling + damping_force;

            // CRITICAL: Noise term
            float noise_term = sigma * std::sqrt(dt) * noise(rng);

            theta_new[idx] = theta[idx] + drift * dt + noise_term;
        }
    }

    theta = theta_new;
}

float computeGlobalR(const std::vector<float>& theta) {
    double sum_real = 0.0;
    double sum_imag = 0.0;

    for (float t : theta) {
        sum_real += std::cos(t);
        sum_imag += std::sin(t);
    }

    sum_real /= theta.size();
    sum_imag /= theta.size();

    return std::sqrt(sum_real * sum_real + sum_imag * sum_imag);
}

float computePhaseVariance(const std::vector<float>& theta) {
    // Compute circular variance: 1 - |<e^(iθ)>|
    double sum_real = 0.0;
    double sum_imag = 0.0;

    for (float t : theta) {
        sum_real += std::cos(t);
        sum_imag += std::sin(t);
    }

    sum_real /= theta.size();
    sum_imag /= theta.size();

    float R = std::sqrt(sum_real * sum_real + sum_imag * sum_imag);
    return 1.0f - R;  // Circular variance
}

int main() {
    std::cout << "=== Noise Implementation Diagnostic Tests ===" << std::endl << std::endl;

    const uint32_t Nx = 32;  // Smaller for quick test
    const uint32_t Ny = 32;
    const float dt = 0.01f;
    const float K = 1.0f;

    std::mt19937 rng(42);

    // ========================================
    // TEST 1: High-noise test (σ = 1.0)
    // ========================================
    std::cout << "TEST 1: High-noise test (σ = 1.0)" << std::endl;
    std::cout << "Expected: R < 0.1 (thermal gas)" << std::endl;
    std::cout << "If R = 1.0: Noise NOT being applied" << std::endl << std::endl;

    {
        std::vector<float> theta(Nx * Ny, 0.0f);  // Start synchronized
        std::vector<float> omega(Nx * Ny, 0.0f);
        float damping = 0.1f;
        float sigma = 1.0f;  // HUGE noise

        std::cout << "Running 1000 steps with σ = 1.0..." << std::endl;
        for (int step = 0; step < 1000; step++) {
            kuramotoStepWithNoise(theta, omega, Nx, Ny, K, damping, dt, sigma, rng);

            if (step % 100 == 0) {
                float R = computeGlobalR(theta);
                float var = computePhaseVariance(theta);
                std::cout << "  Step " << step << ": R = " << std::fixed << std::setprecision(4)
                          << R << ", variance = " << std::scientific << var << std::endl;
            }
        }

        float R_final = computeGlobalR(theta);
        std::cout << "RESULT: R_final = " << std::fixed << std::setprecision(4) << R_final << std::endl;

        if (R_final > 0.5f) {
            std::cout << "❌ FAIL: R too high! Noise not being applied correctly." << std::endl;
        } else {
            std::cout << "✅ PASS: Noise destroys synchronization as expected." << std::endl;
        }
        std::cout << std::endl;
    }

    // ========================================
    // TEST 2: Phase variance tracking
    // ========================================
    std::cout << "TEST 2: Phase variance tracking" << std::endl;
    std::cout << "Expected: variance grows with time for σ > 0" << std::endl;
    std::cout << "If variance = 0: Evolution is frozen" << std::endl << std::endl;

    {
        std::vector<float> theta(Nx * Ny, 0.0f);
        std::vector<float> omega(Nx * Ny, 0.0f);
        float damping = 0.1f;
        float sigma = 1e-4f;

        std::cout << "Running 1000 steps with σ = 1e-4, tracking variance..." << std::endl;

        for (int step = 0; step <= 1000; step++) {
            if (step > 0) {
                kuramotoStepWithNoise(theta, omega, Nx, Ny, K, damping, dt, sigma, rng);
            }

            if (step % 200 == 0) {
                float var = computePhaseVariance(theta);
                float R = computeGlobalR(theta);

                // Also compute standard deviation of phases
                float mean_theta = 0.0f;
                for (float t : theta) mean_theta += t;
                mean_theta /= theta.size();

                float std_dev = 0.0f;
                for (float t : theta) {
                    float diff = t - mean_theta;
                    std_dev += diff * diff;
                }
                std_dev = std::sqrt(std_dev / theta.size());

                std::cout << "  Step " << step << ": R = " << std::fixed << std::setprecision(6) << R
                          << ", circular_var = " << std::scientific << std::setprecision(3) << var
                          << ", std_dev = " << std_dev << std::endl;
            }
        }

        float var_final = computePhaseVariance(theta);
        if (var_final < 1e-10f) {
            std::cout << "❌ FAIL: Variance is zero! Phases frozen or noise not applied." << std::endl;
        } else {
            std::cout << "✅ PASS: Variance non-zero, phases evolving." << std::endl;
        }
        std::cout << std::endl;
    }

    // ========================================
    // TEST 3: Manual phase inspection
    // ========================================
    std::cout << "TEST 3: Manual phase inspection" << std::endl;
    std::cout << "Check: Are phases identical or distributed?" << std::endl << std::endl;

    {
        std::vector<float> theta(Nx * Ny, 0.0f);
        std::vector<float> omega(Nx * Ny, 0.0f);
        float damping = 0.1f;
        float sigma = 1e-4f;

        for (int step = 0; step < 1000; step++) {
            kuramotoStepWithNoise(theta, omega, Nx, Ny, K, damping, dt, sigma, rng);
        }

        std::cout << "First 20 phases after 1000 steps (σ = 1e-4):" << std::endl;
        std::cout << std::fixed << std::setprecision(8);
        for (int i = 0; i < 20; i++) {
            std::cout << "  θ[" << std::setw(2) << i << "] = " << theta[i] << std::endl;
        }

        // Check if all identical
        bool all_identical = true;
        for (size_t i = 1; i < theta.size(); i++) {
            if (std::abs(theta[i] - theta[0]) > 1e-10f) {
                all_identical = false;
                break;
            }
        }

        if (all_identical) {
            std::cout << "❌ FAIL: All phases IDENTICAL! This is frozen state, not synchronization." << std::endl;
        } else {
            std::cout << "✅ PASS: Phases have variation (not frozen)." << std::endl;
        }
        std::cout << std::endl;
    }

    // ========================================
    // TEST 4: No damping test
    // ========================================
    std::cout << "TEST 4: No damping test (γ = 0)" << std::endl;
    std::cout << "Expected: R decreases over time (noise accumulates)" << std::endl;
    std::cout << "If R = 1.0: Noise implementation broken" << std::endl << std::endl;

    {
        std::vector<float> theta(Nx * Ny, 0.0f);
        std::vector<float> omega(Nx * Ny, 0.0f);
        float damping = 0.0f;  // NO DAMPING
        float sigma = 1e-4f;

        std::cout << "Running 2000 steps with σ = 1e-4, γ = 0..." << std::endl;

        float R_initial = computeGlobalR(theta);
        std::cout << "  Step 0: R = " << std::fixed << std::setprecision(6) << R_initial << std::endl;

        for (int step = 0; step < 2000; step++) {
            kuramotoStepWithNoise(theta, omega, Nx, Ny, K, damping, dt, sigma, rng);

            if ((step + 1) % 500 == 0) {
                float R = computeGlobalR(theta);
                std::cout << "  Step " << (step + 1) << ": R = " << R << std::endl;
            }
        }

        float R_final = computeGlobalR(theta);

        if (R_final > 0.999f) {
            std::cout << "❌ FAIL: R unchanged! Noise not accumulating." << std::endl;
        } else {
            std::cout << "✅ PASS: R decreased, noise is accumulating." << std::endl;
        }
        std::cout << std::endl;
    }

    std::cout << "=== Diagnostic Tests Complete ===" << std::endl;

    return 0;
}
