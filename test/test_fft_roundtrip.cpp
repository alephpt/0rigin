/**
 * test_fft_roundtrip.cpp
 *
 * Test FFT forward/backward roundtrip for norm preservation
 */

#include <fftw3.h>
#include <complex>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>

int main() {
    std::cout << "=== FFT Roundtrip Norm Preservation Test ===" << std::endl;

    const int N = 16;
    const int N_total = N * N * N;

    std::vector<std::complex<float>> data(N_total);
    std::vector<std::complex<float>> data_k(N_total);

    // Initialize with random data
    std::mt19937 gen(42);
    std::normal_distribution<float> dist(0.0f, 1.0f);

    for (int i = 0; i < N_total; ++i) {
        data[i] = std::complex<float>(dist(gen), dist(gen));
    }

    // Compute initial norm
    float norm_initial = 0;
    for (const auto& val : data) {
        norm_initial += std::norm(val);
    }

    std::cout << "Grid: " << N << "^3 = " << N_total << " points" << std::endl;
    std::cout << "Initial norm: " << std::scientific << norm_initial << std::endl;

    // Create FFTW plans
    fftwf_plan forward = fftwf_plan_dft_3d(
        N, N, N,
        reinterpret_cast<fftwf_complex*>(data.data()),
        reinterpret_cast<fftwf_complex*>(data_k.data()),
        FFTW_FORWARD,
        FFTW_ESTIMATE
    );

    fftwf_plan backward = fftwf_plan_dft_3d(
        N, N, N,
        reinterpret_cast<fftwf_complex*>(data_k.data()),
        reinterpret_cast<fftwf_complex*>(data.data()),
        FFTW_BACKWARD,
        FFTW_ESTIMATE
    );

    // Test multiple roundtrips
    std::cout << "\nRoundtrip test:" << std::endl;
    std::cout << std::setw(10) << "Iteration"
              << std::setw(15) << "Norm"
              << std::setw(15) << "Drift"
              << std::endl;
    std::cout << std::string(40, '-') << std::endl;

    for (int iter = 0; iter <= 10; ++iter) {
        if (iter > 0) {
            // Forward FFT
            fftwf_execute(forward);

            // Backward FFT
            fftwf_execute(backward);

            // Normalize
            for (auto& val : data) {
                val /= N_total;
            }
        }

        // Compute norm
        float norm = 0;
        for (const auto& val : data) {
            norm += std::norm(val);
        }

        float drift = (norm - norm_initial) / norm_initial;

        std::cout << std::setw(10) << iter
                  << std::scientific << std::setprecision(6)
                  << std::setw(15) << norm
                  << std::setw(15) << drift
                  << std::endl;
    }

    // Cleanup
    fftwf_destroy_plan(forward);
    fftwf_destroy_plan(backward);

    std::cout << "\n=== Result ===" << std::endl;
    float final_norm = 0;
    for (const auto& val : data) {
        final_norm += std::norm(val);
    }
    float total_drift = std::abs(final_norm - norm_initial) / norm_initial;

    if (total_drift < 1e-6) {
        std::cout << "✓ FFT roundtrip preserves norm perfectly (drift < 1e-6)" << std::endl;
    } else if (total_drift < 1e-4) {
        std::cout << "⚠ Small FFT roundtrip error (drift ~ " << total_drift << ")" << std::endl;
    } else {
        std::cout << "✗ Significant FFT roundtrip error!" << std::endl;
    }

    return 0;
}