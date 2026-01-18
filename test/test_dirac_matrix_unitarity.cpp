/**
 * test_dirac_matrix_unitarity.cpp
 *
 * Test that the kinetic evolution matrix is actually unitary
 */

#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <array>

// Test the matrix exponential formula used in Dirac3D
void testKineticUnitarity() {
    // Dirac alpha matrices (same as in Dirac3D.cpp)
    const std::array<std::complex<float>, 16> alpha_x = {{
        {0,0}, {0,0}, {0,0}, {1,0},
        {0,0}, {0,0}, {1,0}, {0,0},
        {0,0}, {1,0}, {0,0}, {0,0},
        {1,0}, {0,0}, {0,0}, {0,0}
    }};

    const std::array<std::complex<float>, 16> alpha_y = {{
        {0,0}, {0,0}, {0,0}, {0,-1},
        {0,0}, {0,0}, {0,1}, {0,0},
        {0,0}, {0,-1}, {0,0}, {0,0},
        {0,1}, {0,0}, {0,0}, {0,0}
    }};

    const std::array<std::complex<float>, 16> alpha_z = {{
        {0,0}, {0,0}, {1,0}, {0,0},
        {0,0}, {0,0}, {0,0}, {-1,0},
        {1,0}, {0,0}, {0,0}, {0,0},
        {0,0}, {-1,0}, {0,0}, {0,0}
    }};

    // Test momentum
    float kx = 0.5f, ky = 0.3f, kz = 0.2f;
    float dt = 0.01f;

    float k_mag = std::sqrt(kx*kx + ky*ky + kz*kz);
    float k_inv = 1.0f / k_mag;
    float phase = k_mag * dt;
    float cos_phase = std::cos(phase);
    float sin_phase = std::sin(phase);

    // Normalized k
    float kx_norm = kx * k_inv;
    float ky_norm = ky * k_inv;
    float kz_norm = kz * k_inv;

    // Construct α·k̂ matrix
    std::array<std::complex<float>, 16> alpha_k;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            alpha_k[i*4 + j] =
                alpha_x[i*4 + j] * kx_norm +
                alpha_y[i*4 + j] * ky_norm +
                alpha_z[i*4 + j] * kz_norm;
        }
    }

    // Construct evolution matrix: U = cos(phase)*I - i*sin(phase)*α·k̂
    std::array<std::complex<float>, 16> U;
    const std::complex<float> i_unit(0.0f, 1.0f);

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            U[i*4 + j] = -i_unit * sin_phase * alpha_k[i*4 + j];
            if (i == j) {
                U[i*4 + j] += cos_phase;
            }
        }
    }

    // Check unitarity: U† U = I
    std::array<std::complex<float>, 16> UdaggerU;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            UdaggerU[i*4 + j] = {0, 0};
            for (int k = 0; k < 4; ++k) {
                UdaggerU[i*4 + j] += std::conj(U[k*4 + i]) * U[k*4 + j];
            }
        }
    }

    // Check if identity
    float max_error = 0.0f;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            std::complex<float> expected = (i == j) ? std::complex<float>(1,0) : std::complex<float>(0,0);
            float error = std::abs(UdaggerU[i*4 + j] - expected);
            max_error = std::max(max_error, error);
        }
    }

    std::cout << "=== Kinetic Evolution Matrix Unitarity Test ===" << std::endl;
    std::cout << "k = (" << kx << ", " << ky << ", " << kz << ")" << std::endl;
    std::cout << "|k| = " << k_mag << std::endl;
    std::cout << "dt = " << dt << std::endl;
    std::cout << "phase = |k|*dt = " << phase << std::endl;
    std::cout << "\nMax deviation from unitarity (U†U - I): "
              << std::scientific << max_error << std::endl;

    if (max_error < 1e-6f) {
        std::cout << "✓ Matrix is unitary to machine precision" << std::endl;
    } else {
        std::cout << "✗ Matrix is NOT unitary!" << std::endl;

        // Print the problematic elements
        std::cout << "\nU†U matrix (should be identity):" << std::endl;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                auto val = UdaggerU[i*4 + j];
                std::cout << "(" << std::setw(8) << val.real()
                          << ", " << std::setw(8) << val.imag() << ") ";
            }
            std::cout << std::endl;
        }
    }

    // Also check norm preservation for a test vector
    std::complex<float> psi[4] = {
        {0.5f, 0.0f},
        {0.5f, 0.0f},
        {0.5f, 0.0f},
        {0.5f, 0.0f}
    };

    float norm_before = 0;
    for (int i = 0; i < 4; ++i) {
        norm_before += std::norm(psi[i]);
    }

    // Apply U
    std::complex<float> psi_new[4] = {{0,0}, {0,0}, {0,0}, {0,0}};
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            psi_new[i] += U[i*4 + j] * psi[j];
        }
    }

    float norm_after = 0;
    for (int i = 0; i < 4; ++i) {
        norm_after += std::norm(psi_new[i]);
    }

    std::cout << "\nTest vector norm preservation:" << std::endl;
    std::cout << "Norm before: " << norm_before << std::endl;
    std::cout << "Norm after: " << norm_after << std::endl;
    std::cout << "Relative change: " << (norm_after - norm_before)/norm_before << std::endl;
}

int main() {
    testKineticUnitarity();
    return 0;
}