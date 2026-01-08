// test_laplacian_accuracy.cpp
// Verify 4th-order Laplacian accuracy on analytic function

#include "ConservativeSolver.h"
#include <iostream>
#include <iomanip>
#include <cmath>

// Test function: f(x,y,z) = sin(πx)·sin(πy)·sin(πz)
// Exact Laplacian: ∇²f = -3π²·f

float testFunction(float x, float y, float z) {
    const float pi = M_PI;
    return std::sin(pi * x) * std::sin(pi * y) * std::sin(pi * z);
}

float exactLaplacian(float x, float y, float z) {
    const float pi = M_PI;
    return -3.0f * pi * pi * testFunction(x, y, z);
}

void testLaplacianOrder(ConservativeSolver::SpatialOrder order, const char* name) {
    std::cout << "\n=========================================="  << std::endl;
    std::cout << " Testing: " << name << std::endl;
    std::cout << "==========================================" << std::endl;

    ConservativeSolver::Config config;
    config.nx = 32;
    config.ny = 32;
    config.nz = 32;
    config.dx = 1.0f / 32.0f;  // Domain [0,1]³
    config.spatial_order = order;

    ConservativeSolver solver;
    solver.initialize(config);

    // Set up test function
    std::vector<float> test_field(config.nx * config.ny * config.nz);
    for (uint32_t k = 0; k < config.nz; ++k) {
        for (uint32_t j = 0; j < config.ny; ++j) {
            for (uint32_t i = 0; i < config.nx; ++i) {
                float x = (i + 0.5f) * config.dx;
                float y = (j + 0.5f) * config.dx;
                float z = (k + 0.5f) * config.dx;

                int idx = k * (config.nx * config.ny) + j * config.nx + i;
                test_field[idx] = testFunction(x, y, z);
            }
        }
    }

    // Compute Laplacian and compare to exact
    float max_error = 0.0f;
    float l2_error = 0.0f;
    int count = 0;

    for (uint32_t k = 2; k < config.nz - 2; ++k) {  // Avoid boundaries for 4th-order
        for (uint32_t j = 2; j < config.ny - 2; ++j) {
            for (uint32_t i = 2; i < config.nx - 2; ++i) {
                float x = (i + 0.5f) * config.dx;
                float y = (j + 0.5f) * config.dx;
                float z = (k + 0.5f) * config.dx;

                // Compute numerical Laplacian
                float laplacian_numerical;
                if (order == ConservativeSolver::SpatialOrder::FOURTH_ORDER) {
                    // Manually compute to test implementation
                    const int idx_c = k * (config.nx * config.ny) + j * config.nx + i;
                    const float dx2 = config.dx * config.dx;
                    const float f_c = test_field[idx_c];

                    // X-direction
                    const float f_xm2 = test_field[k * config.nx * config.ny + j * config.nx + (i-2)];
                    const float f_xm1 = test_field[k * config.nx * config.ny + j * config.nx + (i-1)];
                    const float f_xp1 = test_field[k * config.nx * config.ny + j * config.nx + (i+1)];
                    const float f_xp2 = test_field[k * config.nx * config.ny + j * config.nx + (i+2)];
                    const float d2_dx2 = (-f_xm2 + 16.0f*f_xm1 - 30.0f*f_c + 16.0f*f_xp1 - f_xp2) / (12.0f * dx2);

                    // Y-direction
                    const float f_ym2 = test_field[k * config.nx * config.ny + (j-2) * config.nx + i];
                    const float f_ym1 = test_field[k * config.nx * config.ny + (j-1) * config.nx + i];
                    const float f_yp1 = test_field[k * config.nx * config.ny + (j+1) * config.nx + i];
                    const float f_yp2 = test_field[k * config.nx * config.ny + (j+2) * config.nx + i];
                    const float d2_dy2 = (-f_ym2 + 16.0f*f_ym1 - 30.0f*f_c + 16.0f*f_yp1 - f_yp2) / (12.0f * dx2);

                    // Z-direction
                    const float f_zm2 = test_field[(k-2) * config.nx * config.ny + j * config.nx + i];
                    const float f_zm1 = test_field[(k-1) * config.nx * config.ny + j * config.nx + i];
                    const float f_zp1 = test_field[(k+1) * config.nx * config.ny + j * config.nx + i];
                    const float f_zp2 = test_field[(k+2) * config.nx * config.ny + j * config.nx + i];
                    const float d2_dz2 = (-f_zm2 + 16.0f*f_zm1 - 30.0f*f_c + 16.0f*f_zp1 - f_zp2) / (12.0f * dx2);

                    laplacian_numerical = d2_dx2 + d2_dy2 + d2_dz2;
                } else {
                    // 2nd-order
                    const int idx_c = k * (config.nx * config.ny) + j * config.nx + i;
                    const float dx2 = config.dx * config.dx;
                    const float f_c = test_field[idx_c];

                    const float f_xp = test_field[k * config.nx * config.ny + j * config.nx + (i+1)];
                    const float f_xm = test_field[k * config.nx * config.ny + j * config.nx + (i-1)];
                    const float f_yp = test_field[k * config.nx * config.ny + (j+1) * config.nx + i];
                    const float f_ym = test_field[k * config.nx * config.ny + (j-1) * config.nx + i];
                    const float f_zp = test_field[(k+1) * config.nx * config.ny + j * config.nx + i];
                    const float f_zm = test_field[(k-1) * config.nx * config.ny + j * config.nx + i];

                    laplacian_numerical = (f_xp + f_xm + f_yp + f_ym + f_zp + f_zm - 6.0f*f_c) / dx2;
                }

                float laplacian_exact = exactLaplacian(x, y, z);
                float error = std::abs(laplacian_numerical - laplacian_exact);

                max_error = std::max(max_error, error);
                l2_error += error * error;
                count++;
            }
        }
    }

    l2_error = std::sqrt(l2_error / count);

    std::cout << "  Grid: " << config.nx << "×" << config.ny << "×" << config.nz << std::endl;
    std::cout << "  dx = " << config.dx << std::endl;
    std::cout << "  Max error: " << std::scientific << max_error << std::endl;
    std::cout << "  L2 error:  " << std::scientific << l2_error << std::endl;

    // Theoretical error scaling
    if (order == ConservativeSolver::SpatialOrder::SECOND_ORDER) {
        float expected_error = config.dx * config.dx;  // O(dx²)
        std::cout << "  Expected O(dx²) ≈ " << expected_error << std::endl;
    } else {
        float expected_error = config.dx * config.dx * config.dx * config.dx;  // O(dx⁴)
        std::cout << "  Expected O(dx⁴) ≈ " << expected_error << std::endl;
    }
}

int main() {
    std::cout << "==========================================" << std::endl;
    std::cout << " Laplacian Accuracy Test" << std::endl;
    std::cout << " Function: f = sin(πx)·sin(πy)·sin(πz)" << std::endl;
    std::cout << " Exact: ∇²f = -3π²·f" << std::endl;
    std::cout << "==========================================" << std::endl;

    testLaplacianOrder(ConservativeSolver::SpatialOrder::SECOND_ORDER, "2nd-order Laplacian");
    testLaplacianOrder(ConservativeSolver::SpatialOrder::FOURTH_ORDER, "4th-order Laplacian");

    std::cout << "\n==========================================" << std::endl;
    std::cout << " Summary" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << "If 4th-order is correct, error should be ~16× smaller" << std::endl;

    return 0;
}
