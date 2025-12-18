/**
 * test_dispersion_full.cpp
 * 
 * Full dispersion relation E(k) using momentum boost
 * Apply phase gradient to create plane wave with definite k
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <complex>
#include "fftw3.h"
#include <cstdint>

int main() {
    const uint32_t Nx = 64;
    const uint32_t Ny = 64;
    const uint32_t N = Nx * Ny;
    const float dt = 0.001f;
    const int steps = 5000;
    const float m0 = 0.5f;  // Rest mass
    
    std::cout << "=== Full Dispersion Relation E(k) ===" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    
    std::ofstream out("output/dispersion_full.dat");
    out << "# k_mag E_measured E_theory error_percent\n";
    
    // Test different k values by applying momentum boost
    std::vector<float> k_test = {0.0f, 0.2f, 0.4f, 0.6f, 0.8f, 1.0f, 1.5f, 2.0f};
    
    for (float k_mag : k_test) {
        // Initialize spinor with momentum k in +x direction
        std::vector<std::complex<float>> psi(N);
        
        // Create plane wave: ψ(x) = exp(i k·x) with Gaussian envelope
        float sigma = 10.0f;
        float norm_sum = 0.0f;
        for (uint32_t y = 0; y < Ny; y++) {
            for (uint32_t x = 0; x < Nx; x++) {
                uint32_t idx = y * Nx + x;
                float dx = x - Nx/2.0f;
                float dy = y - Ny/2.0f;
                float r2 = dx*dx + dy*dy;
                float envelope = std::exp(-r2 / (2.0f * sigma * sigma));
                float phase = k_mag * x;
                psi[idx] = envelope * std::exp(std::complex<float>(0.0f, phase));
                norm_sum += std::norm(psi[idx]);
            }
        }
        
        // Normalize
        float norm_factor = std::sqrt(norm_sum);
        for (uint32_t i = 0; i < N; i++) {
            psi[i] /= norm_factor;
        }
        
        // Setup FFTW
        std::vector<std::complex<float>> psi_k(N);
        fftwf_plan fft_forward = fftwf_plan_dft_2d(Ny, Nx,
            reinterpret_cast<fftwf_complex*>(psi.data()),
            reinterpret_cast<fftwf_complex*>(psi_k.data()),
            FFTW_FORWARD, FFTW_ESTIMATE);
        fftwf_plan fft_backward = fftwf_plan_dft_2d(Ny, Nx,
            reinterpret_cast<fftwf_complex*>(psi_k.data()),
            reinterpret_cast<fftwf_complex*>(psi.data()),
            FFTW_BACKWARD, FFTW_ESTIMATE);
        
        // Setup momentum grid
        std::vector<float> kx(Nx), ky(Ny);
        float dk_x = 2.0f * M_PI / Nx;
        float dk_y = 2.0f * M_PI / Ny;
        for (uint32_t i = 0; i < Nx; i++) {
            int freq = (i < Nx/2) ? i : i - Nx;
            kx[i] = freq * dk_x;
        }
        for (uint32_t j = 0; j < Ny; j++) {
            int freq = (j < Ny/2) ? j : j - Ny;
            ky[j] = freq * dk_y;
        }
        
        // Get initial global phase
        std::complex<float> psi_center_init = psi[Ny/2 * Nx + Nx/2];
        float phase_init = std::arg(psi_center_init);
        
        // Time evolution with split-operator
        for (int step = 0; step < steps; step++) {
            // K/2
            fftwf_execute(fft_forward);
            for (uint32_t j = 0; j < Ny; j++) {
                for (uint32_t i = 0; i < Nx; i++) {
                    uint32_t idx = j * Nx + i;
                    float k = std::sqrt(kx[i]*kx[i] + ky[j]*ky[j]);
                    // Dirac: exp(-i|k|Δt/2) for massless, exp(-i√(k²+m²)Δt/2) for massive
                    float E_k = std::sqrt(k*k + m0*m0);
                    psi_k[idx] *= std::exp(std::complex<float>(0.0f, -E_k * dt / 2.0f));
                }
            }
            fftwf_execute(fft_backward);
            for (uint32_t i = 0; i < N; i++) {
                psi[i] /= float(N);
            }
            
            // V (mass term - constant, so just global phase)
            for (uint32_t i = 0; i < N; i++) {
                psi[i] *= std::exp(std::complex<float>(0.0f, -m0 * dt));
            }
            
            // K/2
            fftwf_execute(fft_forward);
            for (uint32_t j = 0; j < Ny; j++) {
                for (uint32_t i = 0; i < Nx; i++) {
                    uint32_t idx = j * Nx + i;
                    float k = std::sqrt(kx[i]*kx[i] + ky[j]*ky[j]);
                    float E_k = std::sqrt(k*k + m0*m0);
                    psi_k[idx] *= std::exp(std::complex<float>(0.0f, -E_k * dt / 2.0f));
                }
            }
            fftwf_execute(fft_backward);
            for (uint32_t i = 0; i < N; i++) {
                psi[i] /= float(N);
            }
        }
        
        // Measure final phase
        std::complex<float> psi_center_final = psi[Ny/2 * Nx + Nx/2];
        float phase_final = std::arg(psi_center_final);
        
        // Extract frequency
        float dphase = phase_final - phase_init;
        while (dphase > M_PI) dphase -= 2.0f * M_PI;
        while (dphase < -M_PI) dphase += 2.0f * M_PI;
        
        float omega_measured = -dphase / (dt * steps);
        float E_theory = std::sqrt(k_mag * k_mag + m0 * m0);
        float error = std::abs(omega_measured - E_theory) / E_theory * 100.0f;
        
        out << k_mag << " " << omega_measured << " " << E_theory << " " << error << "\n";
        
        std::cout << "k=" << std::setw(5) << k_mag;
        std::cout << " → E_meas=" << std::setw(8) << omega_measured;
        std::cout << ", E_theo=" << std::setw(8) << E_theory;
        std::cout << ", error=" << std::setw(6) << error << "%";
        
        if (error < 2.0f) {
            std::cout << " [PASS]" << std::endl;
        } else if (error < 10.0f) {
            std::cout << " [WARN]" << std::endl;
        } else {
            std::cout << " [FAIL]" << std::endl;
        }
        
        fftwf_destroy_plan(fft_forward);
        fftwf_destroy_plan(fft_backward);
    }
    
    out.close();
    std::cout << "\nOutput: output/dispersion_full.dat" << std::endl;
    
    return 0;
}
