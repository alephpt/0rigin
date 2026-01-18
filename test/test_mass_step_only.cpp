#include "Dirac3D.h"
#include <iostream>
#include <cmath>
#include <numeric>

float computeNorm(const Dirac3D& dirac) {
    auto density = dirac.getDensity();
    return std::accumulate(density.begin(), density.end(), 0.0f);
}

int main() {
    std::cout << "Testing applyMassStep stability..." << std::endl;

    const uint32_t N = 32;
    const float mass = 2.5f;
    const float dt = 0.01f;

    Dirac3D dirac(N, N, N);
    dirac.initializeGaussian(0.0f, 0.0f, 0.0f, 3.0f);

    std::vector<float> mass_field(N*N*N, mass);

    float initial = computeNorm(dirac);

    // Apply ONLY mass step (no kinetic)
    for (int i = 0; i < 1000; ++i) {
        dirac.applyMassStep(mass_field, dt);
    }

    float final = computeNorm(dirac);
    float drift = std::abs(final - initial) / initial;

    std::cout << "Initial norm: " << initial << std::endl;
    std::cout << "Final norm: " << final << std::endl;
    std::cout << "Drift: " << drift * 100 << "%" << std::endl;
    std::cout << "Expected: 0% (mass step alone should preserve norm exactly)" << std::endl;

    return drift < 1e-6 ? 0 : 1;
}