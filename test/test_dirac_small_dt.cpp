#include "Dirac3D.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <numeric>

float computeNorm(const Dirac3D& dirac) {
    auto density = dirac.getDensity();
    return std::accumulate(density.begin(), density.end(), 0.0f);
}

int main() {
    std::cout << "=== Testing Chiral Mass with dt=0.001 ===" << std::endl;

    const uint32_t N = 32;
    const float dt = 0.001f;  // 10x smaller timestep
    const float Delta = 2.5f;
    const int num_steps = 10000;  // Same total time

    Dirac3D dirac(N, N, N);
    dirac.initializeGaussian(0.0f, 0.0f, 0.0f, 3.0f);

    // Create vacuum fields (simple vortex)
    std::vector<float> R_field(N*N*N);
    std::vector<float> theta_field(N*N*N);

    for (uint32_t k = 0; k < N; ++k) {
        for (uint32_t j = 0; j < N; ++j) {
            for (uint32_t i = 0; i < N; ++i) {
                uint32_t idx = k * N * N + j * N + i;
                float x = (float)i - N/2.0f;
                float y = (float)j - N/2.0f;
                float z = (float)k - N/2.0f;
                float r_xy = std::sqrt(x*x + y*y) + 1e-6f;

                theta_field[idx] = std::atan2(y, x);
                R_field[idx] = std::tanh(r_xy / 5.0f);
            }
        }
    }

    float initial_norm = computeNorm(dirac);
    std::cout << "Initial norm: " << initial_norm << std::endl;

    // Evolve with smaller timestep
    for (int step = 0; step < num_steps; ++step) {
        dirac.stepWithChiralMass(R_field, theta_field, Delta, dt);

        if (step % 1000 == 999) {
            float norm = computeNorm(dirac);
            float drift = std::abs(norm - initial_norm) / initial_norm;
            std::cout << "Step " << (step+1) << ": norm drift = "
                      << drift * 100 << "%" << std::endl;
        }
    }

    float final_norm = computeNorm(dirac);
    float norm_drift = std::abs(final_norm - initial_norm) / initial_norm;

    std::cout << "\nFinal norm drift: " << norm_drift * 100 << "%\n";
    std::cout << "Target: < 0.01%\n";
    std::cout << "Result: " << (norm_drift < 0.0001f ? "PASS" : "FAIL") << std::endl;

    return norm_drift < 0.0001f ? 0 : 1;
}