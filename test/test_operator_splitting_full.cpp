/**
 * Full 10,000-step operator splitting simulation
 * CPU-based implementation with complete output
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <complex>
#include <iomanip>
#include <filesystem>
#include <chrono>
#include "../src/MSFTCommon.h"

namespace fs = std::filesystem;

// Simplified Kuramoto-Dirac coupled system (CPU-only)
class OperatorSplittingSimulation {
public:
    OperatorSplittingSimulation(int grid_size, int substep_ratio)
        : N(grid_size), N_points(grid_size * grid_size), substep_ratio(substep_ratio),
          theta(N_points), omega(N_points), R(N_points),
          theta_sum(N_points, 0.0f), R_sum(N_points, 0.0f),
          theta_avg(N_points), R_avg(N_points),
          psi(N_points, std::complex<float>(0.0f, 0.0f)),
          substep_count(0) {

        std::cout << "Initializing simulation: " << N << "x" << N << " grid" << std::endl;
        std::cout << "Substep ratio N = " << substep_ratio << std::endl;
    }

    void initialize(float K, float Delta) {
        this->K = K;
        this->Delta = Delta;

        // Random initial phases and frequencies
        std::random_device rd;
        std::mt19937 gen(42); // Fixed seed for reproducibility
        std::uniform_real_distribution<float> phase_dist(0.0f, 2.0f * M_PI);
        std::uniform_real_distribution<float> omega_dist(-0.1f, 0.1f);

        for (int i = 0; i < N_points; i++) {
            theta[i] = phase_dist(gen);
            omega[i] = omega_dist(gen);
        }

        // Initialize Dirac field as Gaussian wavepacket at center
        float x0 = N / 2.0f;
        float y0 = N / 2.0f;
        float sigma = 5.0f;
        float norm_sum = 0.0f;

        for (int y = 0; y < N; y++) {
            for (int x = 0; x < N; x++) {
                int idx = y * N + x;
                float dx = x - x0;
                float dy = y - y0;
                float r2 = dx*dx + dy*dy;
                float amplitude = std::exp(-r2 / (2.0f * sigma * sigma));
                psi[idx] = std::complex<float>(amplitude, 0.0f);
                norm_sum += amplitude * amplitude;
            }
        }

        // Normalize
        float norm_factor = std::sqrt(norm_sum);
        for (auto& p : psi) p /= norm_factor;

        std::cout << "✓ Initialized with Gaussian wavepacket at (" << x0 << ", " << y0 << ")" << std::endl;
    }

    void step(float dt) {
        // Fast subsystem: Kuramoto evolution
        stepKuramoto(dt);

        // Accumulate for time averaging
        for (int i = 0; i < N_points; i++) {
            theta_sum[i] += theta[i];
            R_sum[i] += R[i];
        }

        substep_count++;

        // Slow subsystem: Dirac evolution (every N steps)
        if (substep_count >= substep_ratio) {
            // Compute time averages
            for (int i = 0; i < N_points; i++) {
                theta_avg[i] = theta_sum[i] / substep_ratio;
                R_avg[i] = R_sum[i] / substep_ratio;
            }

            // Evolve Dirac with averaged fields
            stepDirac(dt * substep_ratio);

            // Reset accumulators
            std::fill(theta_sum.begin(), theta_sum.end(), 0.0f);
            std::fill(R_sum.begin(), R_sum.end(), 0.0f);
            substep_count = 0;
        }
    }

    void stepKuramoto(float dt) {
        std::vector<float> theta_new = theta;

        // Compute sync field R using MSFTCommon
        R = MSFT::computeLocalR(theta, N, N);

        // Evolve phases with Kuramoto dynamics + Dirac coupling
        for (int i = 0; i < N_points; i++) {
            int x = i % N;
            int y = i / N;

            float coupling = 0.0f;

            // Nearest neighbors
            for (int dy = -1; dy <= 1; dy++) {
                for (int dx = -1; dx <= 1; dx++) {
                    if (dx == 0 && dy == 0) continue;
                    int nx = (x + dx + N) % N;
                    int ny = (y + dy + N) % N;
                    int nidx = ny * N + nx;
                    coupling += std::sin(theta[nidx] - theta[i]);
                }
            }

            // Dirac density feedback
            float psi_density = std::norm(psi[i]);
            float mass_coupling = Delta * R[i] * psi_density;

            theta_new[i] = theta[i] + dt * (omega[i] + K * coupling / 8.0f + 0.1f * mass_coupling);
        }

        theta = theta_new;
    }

    // Compute Dirac Hamiltonian H*psi at given index
    std::complex<float> computeDiracH(const std::vector<std::complex<float>>& psi_state, int idx) {
        int x = idx % N;
        int y = idx / N;

        // Mass field from averaged sync field
        float mass = Delta * R_avg[idx];

        // Spatial derivatives (finite differences, periodic BC)
        int xp = (x + 1) % N;
        int xm = (x + N - 1) % N;
        int yp = (y + 1) % N;
        int ym = (y + N - 1) % N;

        std::complex<float> dpsi_dx = (psi_state[y * N + xp] - psi_state[y * N + xm]) * 0.5f;
        std::complex<float> dpsi_dy = (psi_state[yp * N + x] - psi_state[ym * N + x]) * 0.5f;

        // Simplified Dirac: H*Ψ = (-i∇ + m)Ψ = -i∇Ψ + m·Ψ
        return -std::complex<float>(0.0f, 1.0f) * (dpsi_dx + dpsi_dy)
               + mass * psi_state[idx];
    }

    void stepDiracEuler(float dt) {
        std::vector<std::complex<float>> psi_new = psi;

        for (int idx = 0; idx < N_points; idx++) {
            // dΨ/dt = -i*H*Ψ
            psi_new[idx] = psi[idx] - std::complex<float>(0.0f, dt) * computeDiracH(psi, idx);
        }

        psi = psi_new;

        // Normalize to prevent exponential growth
        float norm = 0.0f;
        for (int idx = 0; idx < N_points; idx++) {
            norm += std::norm(psi[idx]);
        }
        norm = std::sqrt(norm);

        for (int idx = 0; idx < N_points; idx++) {
            psi[idx] /= norm;
        }
    }

    void stepDiracRK4(float dt) {
        std::vector<std::complex<float>> k1(N_points), k2(N_points), k3(N_points), k4(N_points);
        std::vector<std::complex<float>> temp(N_points);

        // k1 = -i*H*psi
        for (int idx = 0; idx < N_points; idx++) {
            k1[idx] = -std::complex<float>(0.0f, 1.0f) * computeDiracH(psi, idx);
        }

        // k2 = -i*H*(psi + 0.5*dt*k1)
        for (int idx = 0; idx < N_points; idx++) {
            temp[idx] = psi[idx] + 0.5f * dt * k1[idx];
        }
        for (int idx = 0; idx < N_points; idx++) {
            k2[idx] = -std::complex<float>(0.0f, 1.0f) * computeDiracH(temp, idx);
        }

        // k3 = -i*H*(psi + 0.5*dt*k2)
        for (int idx = 0; idx < N_points; idx++) {
            temp[idx] = psi[idx] + 0.5f * dt * k2[idx];
        }
        for (int idx = 0; idx < N_points; idx++) {
            k3[idx] = -std::complex<float>(0.0f, 1.0f) * computeDiracH(temp, idx);
        }

        // k4 = -i*H*(psi + dt*k3)
        for (int idx = 0; idx < N_points; idx++) {
            temp[idx] = psi[idx] + dt * k3[idx];
        }
        for (int idx = 0; idx < N_points; idx++) {
            k4[idx] = -std::complex<float>(0.0f, 1.0f) * computeDiracH(temp, idx);
        }

        // psi_new = psi + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
        for (int idx = 0; idx < N_points; idx++) {
            psi[idx] += (dt / 6.0f) * (k1[idx] + 2.0f * k2[idx] + 2.0f * k3[idx] + k4[idx]);
        }

        // Normalize to ensure unitarity
        float norm = 0.0f;
        for (int idx = 0; idx < N_points; idx++) {
            norm += std::norm(psi[idx]);
        }
        norm = std::sqrt(norm);

        for (int idx = 0; idx < N_points; idx++) {
            psi[idx] /= norm;
        }
    }

    void stepDirac(float dt) {
        stepDiracRK4(dt);  // Use RK4 by default
    }

    void saveSnapshot(const std::string& prefix, int step) {
        std::ofstream file(prefix + "_step_" + std::to_string(step) + ".dat");
        file << "# Step " << step << "\n";
        file << "# x y theta R psi_density\n";

        for (int y = 0; y < N; y++) {
            for (int x = 0; x < N; x++) {
                int idx = y * N + x;
                float psi_density = std::norm(psi[idx]);
                file << x << " " << y << " "
                     << theta[idx] << " " << R[idx] << " "
                     << psi_density << "\n";
            }
        }
    }

    void saveTimeseries(const std::string& filename) {
        std::ofstream file(filename);
        file << "# Time series data\n";
        file << "# step mean_R std_R mean_psi_density max_psi_density\n";

        for (const auto& record : timeseries) {
            file << record.step << " "
                 << record.mean_R << " " << record.std_R << " "
                 << record.mean_psi << " " << record.max_psi << "\n";
        }
    }

    void recordTimeseries(int step) {
        float mean_R = 0.0f, mean_psi = 0.0f, max_psi = 0.0f;
        float sum_R2 = 0.0f;

        for (int i = 0; i < N_points; i++) {
            mean_R += R[i];
            sum_R2 += R[i] * R[i];
            float psi_density = std::norm(psi[i]);
            mean_psi += psi_density;
            max_psi = std::max(max_psi, psi_density);
        }

        mean_R /= N_points;
        mean_psi /= N_points;
        float std_R = std::sqrt(sum_R2 / N_points - mean_R * mean_R);

        timeseries.push_back({step, mean_R, std_R, mean_psi, max_psi});
    }

private:
    int N;
    int N_points;
    int substep_ratio;
    int substep_count;
    float K, Delta;

    std::vector<float> theta, omega, R;
    std::vector<float> theta_sum, R_sum, theta_avg, R_avg;
    std::vector<std::complex<float>> psi;

    struct TimeseriesRecord {
        int step;
        float mean_R, std_R, mean_psi, max_psi;
    };
    std::vector<TimeseriesRecord> timeseries;
};

int main(int argc, char* argv[]) {
    std::cout << "==================================================" << std::endl;
    std::cout << "Operator Splitting Full Simulation" << std::endl;
    std::cout << "==================================================" << std::endl;

    // Parse command line arguments
    int grid_size = 256;  // Default to 256×256
    int total_steps = 10000;
    int substep_ratio = 10;

    if (argc > 1) grid_size = std::atoi(argv[1]);
    if (argc > 2) total_steps = std::atoi(argv[2]);
    if (argc > 3) substep_ratio = std::atoi(argv[3]);

    std::cout << "\nCommand line options:" << std::endl;
    std::cout << "  Usage: " << argv[0] << " [grid_size] [total_steps] [substep_ratio]" << std::endl;
    std::cout << "  Current: grid_size=" << grid_size
              << ", total_steps=" << total_steps
              << ", N=" << substep_ratio << std::endl;

    // Validate grid size
    if (grid_size > 512) {
        std::cerr << "Warning: Grid size " << grid_size << " may cause memory issues" << std::endl;
        std::cerr << "Maximum recommended: 512×512 (262,144 points)" << std::endl;
        return 1;
    }

    // Create output directory
    std::string output_dir = "output/07_" + std::to_string(grid_size) + "x" + std::to_string(grid_size);
    fs::create_directories(output_dir);

    // Simulation parameters
    const float dt = 0.01f;
    const float K = 1.0f;
    const float Delta = 1.0f;

    OperatorSplittingSimulation sim(grid_size, substep_ratio);
    sim.initialize(K, Delta);

    std::cout << "\nRunning simulation..." << std::endl;
    std::cout << "  Total steps: " << total_steps << std::endl;
    std::cout << "  dt = " << dt << std::endl;
    std::cout << "  K = " << K << ", Δ = " << Delta << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    for (int step = 0; step < total_steps; step++) {
        sim.step(dt);
        sim.recordTimeseries(step);

        // Save snapshots at key points (0%, 10%, 50%, 100%)
        if (step == 0 || step == total_steps/10 || step == total_steps/2 || step == total_steps - 1) {
            sim.saveSnapshot(output_dir + "/snapshot", step);
            std::cout << "  ✓ Saved snapshot at step " << step << std::endl;
        }

        // Progress indicator
        if ((step + 1) % 1000 == 0) {
            std::cout << "  Progress: " << (step + 1) << " / " << total_steps
                     << " steps complete" << std::endl;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "\n✓ Simulation complete in " << duration.count() << " ms" << std::endl;

    // Save timeseries
    sim.saveTimeseries(output_dir + "/timeseries.dat");
    std::cout << "✓ Saved timeseries to " << output_dir << "/timeseries.dat" << std::endl;

    std::cout << "\n=== Output Files ===" << std::endl;
    std::cout << output_dir << "/snapshot_step_0.dat           - Initial state" << std::endl;
    std::cout << output_dir << "/snapshot_step_" << (total_steps/10) << ".dat     - 10% complete" << std::endl;
    std::cout << output_dir << "/snapshot_step_" << (total_steps/2) << ".dat    - 50% complete" << std::endl;
    std::cout << output_dir << "/snapshot_step_" << (total_steps-1) << ".dat    - Final state" << std::endl;
    std::cout << output_dir << "/timeseries.dat                - Full evolution" << std::endl;

    std::cout << "\n=== SUCCESS ===" << std::endl;
    std::cout << "Operator splitting simulation complete!" << std::endl;
    std::cout << "Results saved to " << output_dir << "/" << std::endl;

    return 0;
}
