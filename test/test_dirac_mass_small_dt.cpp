#include "Dirac3D.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <numeric>

float computeEnergy(const Dirac3D& dirac) {
    auto density = dirac.getDensity();
    return std::accumulate(density.begin(), density.end(), 0.0f);
}

float computeNorm(const Dirac3D& dirac) {
    auto density = dirac.getDensity();
    float norm = 0.0f;
    for (float d : density) norm += d;
    return norm;
}

int main() {
    std::cout << "=== Testing with Different Timesteps ===" << std::endl;

    const uint32_t N = 32;
    const float mass = 2.5f;
    const float total_time = 10.0f;

    std::vector<float> timesteps = {0.01f, 0.005f, 0.002f, 0.001f};

    for (float dt : timesteps) {
        int num_steps = static_cast<int>(total_time / dt);

        Dirac3D dirac(N, N, N);
        dirac.initializeGaussian(0.0f, 0.0f, 0.0f, 3.0f);

        std::vector<float> mass_field(N*N*N, mass);

        float initial_norm = computeNorm(dirac);

        // Evolve
        for (int step = 0; step < num_steps; ++step) {
            dirac.step(mass_field, dt);
        }

        float final_norm = computeNorm(dirac);
        float norm_drift = std::abs(final_norm - initial_norm) / initial_norm;

        std::cout << "dt = " << dt << ": "
                  << "steps = " << num_steps
                  << ", norm drift = " << norm_drift * 100 << "%"
                  << (norm_drift < 0.0001f ? " PASS" : " FAIL") << std::endl;
    }

    return 0;
}