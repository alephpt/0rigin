#include "Dirac3D.h"
#include <iostream>
#include <cmath>
#include <numeric>
#include <iomanip>

float computeNorm(const Dirac3D& dirac) {
    auto density = dirac.getDensity();
    float norm = 0.0f;
    for (size_t i = 0; i < density.size(); ++i) {
        norm += density[i];
    }
    return norm;
}

int main() {
    std::cout << "=== Testing Norm Conservation ===" << std::endl;
    std::cout << std::scientific << std::setprecision(6);

    const uint32_t N = 16;  // Smaller grid for testing
    const float dt = 0.001f;

    Dirac3D dirac(N, N, N);
    dirac.initializeGaussian(0.0f, 0.0f, 0.0f, 3.0f);

    // Test 1: Pure phase rotation should preserve norm EXACTLY
    std::cout << "\nTest 1: Zero mass field (should preserve norm exactly)" << std::endl;
    {
        std::vector<float> zero_mass(N*N*N, 0.0f);
        float norm0 = computeNorm(dirac);
        for (int i = 0; i < 100; ++i) {
            dirac.step(zero_mass, dt);
        }
        float norm1 = computeNorm(dirac);
        float drift = std::abs(norm1 - norm0);
        std::cout << "Initial: " << norm0 << ", Final: " << norm1
                  << ", Absolute drift: " << drift << std::endl;
    }

    // Test 2: Uniform mass field
    std::cout << "\nTest 2: Uniform mass field" << std::endl;
    {
        Dirac3D dirac2(N, N, N);
        dirac2.initializeGaussian(0.0f, 0.0f, 0.0f, 3.0f);

        std::vector<float> mass(N*N*N, 1.0f);
        float norm0 = computeNorm(dirac2);
        for (int i = 0; i < 100; ++i) {
            dirac2.step(mass, dt);
        }
        float norm1 = computeNorm(dirac2);
        float drift = std::abs(norm1 - norm0) / norm0;
        std::cout << "Initial: " << norm0 << ", Final: " << norm1
                  << ", Relative drift: " << drift * 100 << "%" << std::endl;
    }

    // Test 3: Chiral mass - simplified
    std::cout << "\nTest 3: Chiral mass (simplified)" << std::endl;
    {
        Dirac3D dirac3(N, N, N);
        dirac3.initializeGaussian(0.0f, 0.0f, 0.0f, 3.0f);

        std::vector<float> R_field(N*N*N, 0.5f);  // Uniform R
        std::vector<float> theta_field(N*N*N, 0.0f);  // Zero theta

        float norm0 = computeNorm(dirac3);
        for (int i = 0; i < 100; ++i) {
            dirac3.stepWithChiralMass(R_field, theta_field, 1.0f, dt);
        }
        float norm1 = computeNorm(dirac3);
        float drift = std::abs(norm1 - norm0) / norm0;
        std::cout << "Initial: " << norm0 << ", Final: " << norm1
                  << ", Relative drift: " << drift * 100 << "%" << std::endl;
    }

    return 0;
}