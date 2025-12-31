/**
 * Find actual critical noise threshold
 * Sweep from σ = 10⁻⁴ to σ = 10
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

int main() {
    std::cout << "=== Finding Critical Noise Threshold ===" << std::endl << std::endl;

    const uint32_t Nx = 64;
    const uint32_t Ny = 64;
    const float dt = 0.01f;
    const float K = 1.0f;
    const float damping = 0.1f;

    std::mt19937 rng(42);

    // Broader sweep
    std::vector<float> sigma_values = {
        1e-4f, 3e-4f,
        1e-3f, 3e-3f,
        1e-2f, 3e-2f,
        1e-1f, 3e-1f,
        1.0f, 3.0f
    };

    std::cout << "Grid: " << Nx << "×" << Ny << std::endl;
    std::cout << "K = " << K << ", γ = " << damping << ", dt = " << dt << std::endl;
    std::cout << "Simulation: 2000 steps per σ value" << std::endl << std::endl;

    std::cout << std::setw(12) << "sigma" << std::setw(12) << "R_final" << std::setw(15) << "interpretation" << std::endl;
    std::cout << std::string(40, '-') << std::endl;

    for (float sigma : sigma_values) {
        std::vector<float> theta(Nx * Ny, 0.0f);  // Start synchronized
        std::vector<float> omega(Nx * Ny, 0.0f);

        // Run simulation
        for (int step = 0; step < 2000; step++) {
            kuramotoStepWithNoise(theta, omega, Nx, Ny, K, damping, dt, sigma, rng);
        }

        float R_final = computeGlobalR(theta);

        std::cout << std::scientific << std::setprecision(2) << std::setw(12) << sigma
                  << std::fixed << std::setprecision(4) << std::setw(12) << R_final;

        if (R_final > 0.95f) {
            std::cout << "    synchronized";
        } else if (R_final > 0.5f) {
            std::cout << "    partial sync";
        } else if (R_final > 0.2f) {
            std::cout << "    weak sync";
        } else {
            std::cout << "    thermal gas";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
    std::cout << "Critical threshold σ_c is where R drops below 0.95" << std::endl;

    return 0;
}
