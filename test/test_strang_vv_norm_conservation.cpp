/**
 * test_strang_vv_norm_conservation.cpp
 *
 * Focused validation of Strang + Velocity Verlet for norm conservation
 * Norm conservation is the definitive test for unitary evolution
 */

#include "Dirac3D.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <complex>
#include <chrono>

// Test configuration
struct TestResult {
    std::string test_name;
    float initial_norm;
    float final_norm;
    float max_norm_drift;
    float norm_drift_percent;
    bool passed;
};

// Test with different configurations
TestResult runTest(const std::string& name,
                   uint32_t Nx, uint32_t Ny, uint32_t Nz,
                   float Delta, float dt, int num_steps,
                   const std::vector<float>& R_field,
                   const std::vector<float>& theta_field) {

    std::cout << "\n--- " << name << " ---" << std::endl;
    std::cout << "Grid: " << Nx << "×" << Ny << "×" << Nz << std::endl;
    std::cout << "dt = " << dt << ", steps = " << num_steps << std::endl;

    Dirac3D dirac(Nx, Ny, Nz);
    dirac.initializeGaussian(Nx/2.0f, Ny/2.0f, Nz/2.0f, 3.0f);

    float norm0 = dirac.getNorm();
    float max_drift = 0.0f;

    std::cout << "Initial norm: " << std::fixed << std::setprecision(6) << norm0 << std::endl;
    std::cout << "Evolving..." << std::endl;

    for (int step = 1; step <= num_steps; ++step) {
        dirac.stepWithChiralMass(R_field, theta_field, Delta, dt);

        float norm = dirac.getNorm();
        float drift = std::abs(norm - norm0);
        max_drift = std::max(max_drift, drift);

        if (step % 200 == 0) {
            std::cout << "  Step " << std::setw(4) << step
                     << ": norm = " << norm
                     << ", drift = " << std::scientific << drift << std::endl;
        }
    }

    float norm_final = dirac.getNorm();
    float drift_percent = std::abs(norm_final - norm0) / norm0 * 100.0f;

    TestResult result;
    result.test_name = name;
    result.initial_norm = norm0;
    result.final_norm = norm_final;
    result.max_norm_drift = max_drift;
    result.norm_drift_percent = drift_percent;
    result.passed = (max_drift < 1e-3);  // 0.1% threshold

    std::cout << std::fixed;
    std::cout << "Final norm: " << norm_final << std::endl;
    std::cout << "Max drift: " << std::scientific << max_drift << std::endl;
    std::cout << "Drift %: " << std::setprecision(4) << drift_percent << "%" << std::endl;

    if (result.passed) {
        std::cout << "✓ PASS" << std::endl;
    } else {
        std::cout << "✗ FAIL" << std::endl;
    }

    return result;
}

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "  STRANG-VV NORM CONSERVATION TEST" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "\nValidating unitary evolution via norm conservation" << std::endl;
    std::cout << "Target: Δ||Ψ||²/||Ψ||² < 0.1%" << std::endl;

    auto start_time = std::chrono::high_resolution_clock::now();

    std::vector<TestResult> results;

    // Test 1: Small grid, uniform scalar mass
    {
        const uint32_t N = 16;
        std::vector<float> R_field(N*N*N, 1.0f);
        std::vector<float> theta_field(N*N*N, 0.0f);

        results.push_back(runTest(
            "Small grid, uniform scalar (θ=0)",
            N, N, N, 2.5f, 0.01f, 500,
            R_field, theta_field
        ));
    }

    // Test 2: Small grid, uniform pseudoscalar
    {
        const uint32_t N = 16;
        std::vector<float> R_field(N*N*N, 1.0f);
        std::vector<float> theta_field(N*N*N, M_PI/2);

        results.push_back(runTest(
            "Small grid, pure pseudoscalar (θ=π/2)",
            N, N, N, 2.5f, 0.01f, 500,
            R_field, theta_field
        ));
    }

    // Test 3: Medium grid, spatially varying
    {
        const uint32_t N = 24;
        std::vector<float> R_field(N*N*N);
        std::vector<float> theta_field(N*N*N);

        // Create radially varying R field
        for (uint32_t k = 0; k < N; ++k) {
            for (uint32_t j = 0; j < N; ++j) {
                for (uint32_t i = 0; i < N; ++i) {
                    float x = (float)i - N/2;
                    float y = (float)j - N/2;
                    float z = (float)k - N/2;
                    float r = std::sqrt(x*x + y*y + z*z);

                    uint32_t idx = k * N * N + j * N + i;
                    R_field[idx] = 1.0f / (1.0f + 0.01f * r * r);
                    theta_field[idx] = std::atan2(y, x);
                }
            }
        }

        results.push_back(runTest(
            "Medium grid, vortex configuration",
            N, N, N, 2.5f, 0.01f, 500,
            R_field, theta_field
        ));
    }

    // Test 4: Longer evolution
    {
        const uint32_t N = 16;
        std::vector<float> R_field(N*N*N, 0.8f);
        std::vector<float> theta_field(N*N*N, M_PI/4);

        results.push_back(runTest(
            "Long evolution (1000 steps), θ=π/4",
            N, N, N, 2.5f, 0.01f, 1000,
            R_field, theta_field
        ));
    }

    // Test 5: Different time step
    {
        const uint32_t N = 16;
        std::vector<float> R_field(N*N*N, 1.0f);
        std::vector<float> theta_field(N*N*N, M_PI/6);

        results.push_back(runTest(
            "Smaller time step (dt=0.005)",
            N, N, N, 2.5f, 0.005f, 500,
            R_field, theta_field
        ));
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);

    // Summary
    std::cout << "\n========================================" << std::endl;
    std::cout << "           TEST SUMMARY" << std::endl;
    std::cout << "========================================" << std::endl;

    int passed = 0;
    int total = results.size();

    std::cout << std::fixed << std::setprecision(4);
    for (const auto& result : results) {
        std::cout << "\n" << result.test_name << ": ";
        if (result.passed) {
            std::cout << "✓ PASS";
            passed++;
        } else {
            std::cout << "✗ FAIL";
        }
        std::cout << " (drift = " << result.norm_drift_percent << "%)" << std::endl;
    }

    std::cout << "\n========================================" << std::endl;
    std::cout << "Tests passed: " << passed << " / " << total << std::endl;
    std::cout << "Test duration: " << duration.count() << " seconds" << std::endl;

    if (passed == total) {
        std::cout << "\n✓✓✓ ALL TESTS PASSED" << std::endl;
        std::cout << "Norm conservation validated: <0.1% drift" << std::endl;
        std::cout << "Unitary evolution confirmed" << std::endl;
        return 0;
    } else {
        std::cout << "\n✗✗✗ SOME TESTS FAILED" << std::endl;
        std::cout << "Norm conservation violated" << std::endl;
        return 1;
    }
}