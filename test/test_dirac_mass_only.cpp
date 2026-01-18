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
    std::cout << "=== Testing Standard Mass Operator ===" << std::endl;

    const uint32_t N = 32;
    const float dt = 0.01f;
    const float mass = 2.5f;
    const int num_steps = 1000;

    Dirac3D dirac(N, N, N);
    dirac.initializeGaussian(0.0f, 0.0f, 0.0f, 3.0f);

    // Create uniform mass field
    std::vector<float> mass_field(N*N*N, mass);

    float initial_energy = computeEnergy(dirac);
    float initial_norm = computeNorm(dirac);

    std::cout << "Initial energy: " << initial_energy << std::endl;
    std::cout << "Initial norm: " << initial_norm << std::endl;

    // Evolve
    for (int step = 0; step < num_steps; ++step) {
        dirac.step(mass_field, dt);
    }

    float final_energy = computeEnergy(dirac);
    float final_norm = computeNorm(dirac);

    float energy_drift = std::abs(final_energy - initial_energy) / initial_energy;
    float norm_drift = std::abs(final_norm - initial_norm) / initial_norm;

    std::cout << "Final energy: " << final_energy << std::endl;
    std::cout << "Final norm: " << final_norm << std::endl;
    std::cout << "Energy drift: " << energy_drift * 100 << "%" << std::endl;
    std::cout << "Norm drift: " << norm_drift * 100 << "%" << std::endl;

    bool pass = (energy_drift < 0.0001f) && (norm_drift < 1e-6f);
    std::cout << "\nResult: " << (pass ? "PASS" : "FAIL") << std::endl;

    return pass ? 0 : 1;
}