/**
 * test_matrix_unitarity.cpp
 *
 * Verify that the Dirac kinetic matrix exponential is unitary: U†U = I
 */

#include <iostream>
#include <complex>
#include <cmath>
#include <iomanip>

using namespace std;

int main() {
    // Test for several k-vectors
    float test_k[][2] = {
        {0.1f, 0.2f},
        {1.0f, 0.5f},
        {2.0f, 3.0f},
        {0.5f, 0.5f}
    };
    
    float dt = 0.01f;
    
    for (int test = 0; test < 4; test++) {
        float kx = test_k[test][0];
        float ky = test_k[test][1];
        float k_mag = sqrt(kx*kx + ky*ky);
        
        float cos_term = cos(k_mag * dt);
        float sin_term = sin(k_mag * dt);
        float kx_norm = kx / k_mag;
        float ky_norm = ky / k_mag;
        
        // Build matrix U = exp(-i(σ·k)Δt)
        complex<float> U[2][2];
        U[0][0] = complex<float>(cos_term, 0.0f);
        U[0][1] = complex<float>(-sin_term * ky_norm, -sin_term * kx_norm);
        U[1][0] = complex<float>(sin_term * ky_norm, -sin_term * kx_norm);
        U[1][1] = complex<float>(cos_term, 0.0f);
        
        // Compute U† (complex conjugate transpose)
        complex<float> Udag[2][2];
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                Udag[i][j] = conj(U[j][i]);
            }
        }
        
        // Compute U†U
        complex<float> UdagU[2][2];
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                UdagU[i][j] = complex<float>(0.0f, 0.0f);
                for (int k = 0; k < 2; k++) {
                    UdagU[i][j] += Udag[i][k] * U[k][j];
                }
            }
        }
        
        // Check if U†U = I
        float max_error = 0.0f;
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                complex<float> expected = (i == j) ? complex<float>(1.0f, 0.0f) : complex<float>(0.0f, 0.0f);
                complex<float> diff = UdagU[i][j] - expected;
                float error = abs(diff);
                max_error = max(max_error, error);
            }
        }
        
        cout << "[TEST " << test << "] k=(" << kx << ", " << ky << "), |k|=" << k_mag;
        cout << ", max|U†U - I|=" << scientific << max_error;
        if (max_error < 1e-6f) {
            cout << " [PASS]" << endl;
        } else {
            cout << " [FAIL]" << endl;
            cout << "  U†U = [" << UdagU[0][0] << ", " << UdagU[0][1] << "]" << endl;
            cout << "        [" << UdagU[1][0] << ", " << UdagU[1][1] << "]" << endl;
        }
    }
    
    return 0;
}
