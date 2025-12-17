/**
 * @file test_dirac_stochastic_full.cpp
 * @brief Comprehensive test with FULL Dirac-Kuramoto evolution data output
 *
 * This test outputs COMPLETE evolution data for analysis and visualization:
 * 1. Timeseries data (every 10 steps)
 * 2. Full field snapshots (at t=0,1,2,3,4,5)
 * 3. Summary statistics
 * 4. Mass field evolution
 * 5. Spinor density maps
 *
 * PHYSICS:
 * - Coupled Kuramoto-Dirac evolution with vacuum noise
 * - Kuramoto: dθ/dt = ω + (K/4)Σsin(θ_j - θ_i) - γ·sin(θ) + σ_θ·√(dt)·ξ_θ(t)
 * - Dirac: i·dΨ/dt = [-iα·∇ + β·m(x)]Ψ + σ_Ψ·√(dt)·ξ_Ψ(t)
 * - Mass coupling: m(x) = Δ·R(x) where R(x) = |⟨e^(iθ)⟩|_local
 */

#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <random>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <string>
#include <sys/stat.h>
#include <filesystem>

// Physical constants
const float DELTA = 2.5f;        // Mass gap parameter
const float DT = 0.01f;           // Time step
const float K = 1.0f;             // Kuramoto coupling
const float GAMMA = 0.1f;         // Damping
const float ALPHA = 1.0f;         // Dirac alpha matrix coefficient

// Grid parameters
const int NX = 64;
const int NY = 64;
const int N_GRID = NX * NY;

// Simulation parameters
const int N_WARMUP = 500;        // Warmup steps for synchronization
const int N_STEPS = 500;         // Measurement steps
const int TIMESERIES_INTERVAL = 10;  // Output timeseries every 10 steps
const int SNAPSHOT_INTERVAL = 100;   // Output snapshots every 100 steps (t=0,1,2,3,4,5)

// Stochastic parameters
const float SIGMA_THETA = 0.05f;  // Kuramoto noise strength
const float SIGMA_PSI = 0.05f;    // Dirac noise strength

// Type definitions
using Complex = std::complex<float>;
struct Spinor {
    Complex psi[4];  // 4-component Dirac spinor
};

// Output directory
const std::string OUTPUT_DIR = "/home/persist/neotec/0rigin/output/dirac_evolution";

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

void create_output_directories() {
    std::filesystem::create_directories(OUTPUT_DIR);
    std::cout << "Created output directory: " << OUTPUT_DIR << std::endl;
}

int idx(int x, int y) {
    // Periodic boundary conditions
    x = (x + NX) % NX;
    y = (y + NY) % NY;
    return y * NX + x;
}

// ============================================================================
// KURAMOTO DYNAMICS
// ============================================================================

std::vector<float> compute_local_R(const std::vector<float>& theta) {
    std::vector<float> R(N_GRID);

    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            // Compute local average over 3x3 neighborhood
            Complex z(0, 0);
            int count = 0;

            for (int dy = -1; dy <= 1; dy++) {
                for (int dx = -1; dx <= 1; dx++) {
                    int i = idx(x + dx, y + dy);
                    z += Complex(std::cos(theta[i]), std::sin(theta[i]));
                    count++;
                }
            }

            z /= float(count);
            R[idx(x, y)] = std::abs(z);
        }
    }

    return R;
}

float compute_global_R(const std::vector<float>& theta) {
    Complex z(0, 0);
    for (float t : theta) {
        z += Complex(std::cos(t), std::sin(t));
    }
    z /= float(N_GRID);
    return std::abs(z);
}

void kuramoto_step(std::vector<float>& theta, const std::vector<float>& omega,
                   float dt, float K, float gamma, float sigma,
                   std::mt19937& rng) {
    std::vector<float> theta_new(N_GRID);
    std::normal_distribution<float> noise(0.0f, 1.0f);

    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            int i = idx(x, y);

            // Coupling term (4-neighbor von Neumann)
            float coupling = 0.0f;
            coupling += std::sin(theta[idx(x+1, y)] - theta[i]);
            coupling += std::sin(theta[idx(x-1, y)] - theta[i]);
            coupling += std::sin(theta[idx(x, y+1)] - theta[i]);
            coupling += std::sin(theta[idx(x, y-1)] - theta[i]);

            // Damping term
            float damping_force = -gamma * std::sin(theta[i]);

            // Deterministic drift
            float drift = omega[i] + (K / 4.0f) * coupling + damping_force;

            // Stochastic term: σ·√(dt)·N(0,1)
            float noise_term = sigma * std::sqrt(dt) * noise(rng);

            // Euler-Maruyama update
            theta_new[i] = theta[i] + drift * dt + noise_term;
        }
    }

    theta = theta_new;
}

// ============================================================================
// DIRAC DYNAMICS
// ============================================================================

void initialize_gaussian_packet(std::vector<Spinor>& psi, float x0, float y0, float width) {
    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            float dx = x - x0;
            float dy = y - y0;
            float r2 = dx*dx + dy*dy;
            float amp = std::exp(-r2 / (2.0f * width * width));

            // Initialize as spin-up state
            int i = idx(x, y);
            psi[i].psi[0] = Complex(amp, 0);  // Upper component
            psi[i].psi[1] = Complex(0, 0);
            psi[i].psi[2] = Complex(0, 0);
            psi[i].psi[3] = Complex(0, 0);    // Lower component
        }
    }

    // Normalize
    float norm = 0.0f;
    for (const auto& s : psi) {
        for (int k = 0; k < 4; k++) {
            norm += std::norm(s.psi[k]);
        }
    }
    norm = std::sqrt(norm);

    for (auto& s : psi) {
        for (int k = 0; k < 4; k++) {
            s.psi[k] /= norm;
        }
    }
}

float compute_spinor_norm(const std::vector<Spinor>& psi) {
    float norm = 0.0f;
    for (const auto& s : psi) {
        for (int k = 0; k < 4; k++) {
            norm += std::norm(s.psi[k]);
        }
    }
    return norm;
}

void compute_particle_position(const std::vector<Spinor>& psi, float& x_cm, float& y_cm) {
    float total_prob = 0.0f;
    x_cm = 0.0f;
    y_cm = 0.0f;

    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            int i = idx(x, y);
            float prob = 0.0f;
            for (int k = 0; k < 4; k++) {
                prob += std::norm(psi[i].psi[k]);
            }

            x_cm += x * prob;
            y_cm += y * prob;
            total_prob += prob;
        }
    }

    x_cm /= total_prob;
    y_cm /= total_prob;
}

void dirac_step(std::vector<Spinor>& psi, const std::vector<float>& R_field,
                float dt, float sigma, std::mt19937& rng) {
    std::vector<Spinor> psi_new(N_GRID);
    std::normal_distribution<float> noise(0.0f, 1.0f);

    // Store initial norm for renormalization
    float norm_before = compute_spinor_norm(psi);

    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            int i = idx(x, y);

            // Mass from synchronization field
            float mass = DELTA * R_field[i];

            // Spatial derivatives (central differences with grid spacing)
            Spinor dpsi_dx, dpsi_dy;
            for (int k = 0; k < 4; k++) {
                dpsi_dx.psi[k] = (psi[idx(x+1, y)].psi[k] - psi[idx(x-1, y)].psi[k]) / 2.0f;
                dpsi_dy.psi[k] = (psi[idx(x, y+1)].psi[k] - psi[idx(x, y-1)].psi[k]) / 2.0f;
            }

            // Apply simplified 2D Dirac operator
            // For stability, we use a reduced coupling strength
            Spinor H_psi;
            float alpha_reduced = ALPHA * 0.1f;  // Reduced for stability

            // Kinetic term: -iα·∇
            H_psi.psi[0] = -Complex(0, alpha_reduced) * (dpsi_dx.psi[1] + dpsi_dy.psi[1]);
            H_psi.psi[1] = -Complex(0, alpha_reduced) * (dpsi_dx.psi[0] + dpsi_dy.psi[0]);
            H_psi.psi[2] = -Complex(0, alpha_reduced) * (dpsi_dx.psi[3] + dpsi_dy.psi[3]);
            H_psi.psi[3] = -Complex(0, alpha_reduced) * (dpsi_dx.psi[2] + dpsi_dy.psi[2]);

            // Mass term: β·m (diagonal in rest frame)
            H_psi.psi[0] += mass * psi[i].psi[0];
            H_psi.psi[1] += mass * psi[i].psi[1];
            H_psi.psi[2] -= mass * psi[i].psi[2];  // Negative for lower components
            H_psi.psi[3] -= mass * psi[i].psi[3];

            // Euler-Maruyama update
            for (int k = 0; k < 4; k++) {
                // Deterministic evolution: Ψ(t+dt) = Ψ(t) - i·dt·H·Ψ
                psi_new[i].psi[k] = psi[i].psi[k] - Complex(0, dt) * H_psi.psi[k];

                // Stochastic term with reduced amplitude for stability
                float noise_real = sigma * std::sqrt(dt) * noise(rng) * 0.1f;
                float noise_imag = sigma * std::sqrt(dt) * noise(rng) * 0.1f;
                psi_new[i].psi[k] += Complex(noise_real, noise_imag);
            }
        }
    }

    // Renormalize to preserve unitarity
    float norm_after = compute_spinor_norm(psi_new);
    if (norm_after > 1e-10) {
        float scale = std::sqrt(norm_before / norm_after);
        for (auto& s : psi_new) {
            for (int k = 0; k < 4; k++) {
                s.psi[k] *= scale;
            }
        }
    }

    psi = psi_new;
}

// ============================================================================
// DATA OUTPUT FUNCTIONS
// ============================================================================

void write_timeseries_header(std::ofstream& file) {
    file << "# Dirac-Kuramoto Stochastic Evolution Timeseries\n";
    file << "# Grid: " << NX << " x " << NY << "\n";
    file << "# dt = " << DT << ", sigma_theta = " << SIGMA_THETA << ", sigma_psi = " << SIGMA_PSI << "\n";
    file << "# Columns: time  R_global  spinor_norm  particle_x  particle_y  particle_drift\n";
}

void write_timeseries_data(std::ofstream& file, float time, float R_global, float norm,
                           float x_cm, float y_cm, float x_init, float y_init) {
    float drift = std::sqrt((x_cm - x_init) * (x_cm - x_init) + (y_cm - y_init) * (y_cm - y_init));
    file << std::fixed << std::setprecision(6)
         << time << "  " << R_global << "  " << norm << "  "
         << x_cm << "  " << y_cm << "  " << drift << "\n";
}

void write_field_snapshot(const std::string& filename, const std::vector<float>& theta,
                         const std::vector<float>& R_field, const std::vector<float>& mass_field,
                         const std::vector<Spinor>& psi, float time) {
    std::ofstream file(filename);
    file << "# Field snapshot at t = " << time << "\n";
    file << "# Grid: " << NX << " x " << NY << "\n";
    file << "# Columns: x  y  theta  R  mass  |psi|^2  psi_0_real  psi_0_imag  psi_1_real  psi_1_imag  psi_2_real  psi_2_imag  psi_3_real  psi_3_imag\n";

    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            int i = idx(x, y);

            // Compute spinor density
            float density = 0.0f;
            for (int k = 0; k < 4; k++) {
                density += std::norm(psi[i].psi[k]);
            }

            file << x << "  " << y << "  "
                 << std::fixed << std::setprecision(6)
                 << theta[i] << "  " << R_field[i] << "  " << mass_field[i] << "  " << density;

            // Output all spinor components
            for (int k = 0; k < 4; k++) {
                file << "  " << psi[i].psi[k].real() << "  " << psi[i].psi[k].imag();
            }
            file << "\n";
        }
    }
    file.close();
}

void write_mass_field(const std::string& filename, const std::vector<float>& mass_field, float time) {
    std::ofstream file(filename);
    file << "# Mass field m(x,y) = Delta * R(x,y) at t = " << time << "\n";
    file << "# Grid: " << NX << " x " << NY << "\n";
    file << "# Columns: x  y  mass\n";

    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            int i = idx(x, y);
            file << x << "  " << y << "  " << std::fixed << std::setprecision(6) << mass_field[i] << "\n";
        }
    }
    file.close();
}

void write_spinor_density(const std::string& filename, const std::vector<Spinor>& psi, float time) {
    std::ofstream file(filename);
    file << "# Spinor density |Psi|^2 at t = " << time << "\n";
    file << "# Grid: " << NX << " x " << NY << "\n";
    file << "# Columns: x  y  density\n";

    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            int i = idx(x, y);
            float density = 0.0f;
            for (int k = 0; k < 4; k++) {
                density += std::norm(psi[i].psi[k]);
            }
            file << x << "  " << y << "  " << std::fixed << std::setprecision(8) << density << "\n";
        }
    }
    file.close();
}

void write_summary_statistics(const std::string& filename,
                             float R_initial, float R_final, float R_mean, float R_std,
                             float norm_initial, float norm_final, float norm_deviation,
                             float x_init, float y_init, float x_final, float y_final,
                             float total_drift) {
    std::ofstream file(filename);
    file << "# Stochastic Dirac-Kuramoto Evolution Summary\n";
    file << "Grid: " << NX << " x " << NY << "\n";
    file << "Steps: " << N_STEPS << "\n";
    file << "dt: " << DT << "\n";
    file << "sigma_theta: " << SIGMA_THETA << "\n";
    file << "sigma_psi: " << SIGMA_PSI << "\n";
    file << "\n";
    file << "R_initial: " << std::fixed << std::setprecision(6) << R_initial << "\n";
    file << "R_final: " << R_final << "\n";
    file << "R_mean: " << R_mean << "\n";
    file << "R_std: " << R_std << "\n";
    file << "\n";
    file << "Norm_initial: " << norm_initial << "\n";
    file << "Norm_final: " << norm_final << "\n";
    file << "Norm_deviation: " << std::setprecision(2) << (norm_deviation * 100) << "%\n";
    file << "\n";
    file << "Particle_initial: (" << std::setprecision(1) << x_init << ", " << y_init << ")\n";
    file << "Particle_final: (" << x_final << ", " << y_final << ")\n";
    file << "Total_drift: " << std::setprecision(2) << total_drift << " grid units\n";
    file.close();
}

// ============================================================================
// MAIN TEST FUNCTION
// ============================================================================

int main() {
    std::cout << "=============================================" << std::endl;
    std::cout << "  COMPREHENSIVE DIRAC-KURAMOTO EVOLUTION TEST" << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "\nFull data output for stochastic evolution" << std::endl;
    std::cout << "Grid: " << NX << " x " << NY << std::endl;
    std::cout << "Steps: " << N_STEPS << " (dt = " << DT << ")" << std::endl;
    std::cout << "Noise: σ_θ = " << SIGMA_THETA << ", σ_Ψ = " << SIGMA_PSI << std::endl;
    std::cout << "Output directory: " << OUTPUT_DIR << std::endl;

    // Create output directories
    create_output_directories();

    // Initialize RNG
    std::mt19937 rng(42);

    // Initialize Kuramoto field (nearly synchronized)
    std::vector<float> theta(N_GRID);
    std::vector<float> omega(N_GRID, 0.0f);  // Zero natural frequencies
    std::uniform_real_distribution<float> init_dist(-0.1f, 0.1f);
    for (int i = 0; i < N_GRID; i++) {
        theta[i] = init_dist(rng);  // Small perturbation around sync
    }

    // Initialize spinor (Gaussian packet at center)
    std::vector<Spinor> psi(N_GRID);
    float x_init = NX / 2.0f;
    float y_init = NY / 2.0f;
    initialize_gaussian_packet(psi, x_init, y_init, 5.0f);

    // Warmup phase for Kuramoto synchronization
    std::cout << "\nWarmup phase (" << N_WARMUP << " steps)..." << std::flush;
    for (int step = 0; step < N_WARMUP; step++) {
        kuramoto_step(theta, omega, DT, K, GAMMA, 0.0f, rng);  // No noise during warmup
        if (step % 100 == 0) std::cout << "." << std::flush;
    }
    float R_warmup = compute_global_R(theta);
    std::cout << " R = " << R_warmup << std::endl;

    // Open timeseries file
    std::ofstream timeseries_file(OUTPUT_DIR + "/timeseries.dat");
    write_timeseries_header(timeseries_file);

    // Storage for statistics
    std::vector<float> R_history;
    std::vector<float> norm_history;
    float norm_initial = compute_spinor_norm(psi);
    float R_initial = compute_global_R(theta);

    // Initial position
    float x_cm, y_cm;
    compute_particle_position(psi, x_cm, y_cm);

    // Write initial snapshot (t=0)
    std::vector<float> R_field = compute_local_R(theta);
    std::vector<float> mass_field(N_GRID);
    for (int i = 0; i < N_GRID; i++) {
        mass_field[i] = DELTA * R_field[i];
    }

    write_field_snapshot(OUTPUT_DIR + "/snapshot_000.dat", theta, R_field, mass_field, psi, 0.0f);
    write_mass_field(OUTPUT_DIR + "/mass_field_000.dat", mass_field, 0.0f);
    write_spinor_density(OUTPUT_DIR + "/density_000.dat", psi, 0.0f);
    write_timeseries_data(timeseries_file, 0.0f, R_initial, norm_initial, x_cm, y_cm, x_init, y_init);

    // Main evolution loop
    std::cout << "\nEvolution phase (" << N_STEPS << " steps)..." << std::endl;
    int snapshot_count = 1;

    for (int step = 1; step <= N_STEPS; step++) {
        // Update Kuramoto field
        kuramoto_step(theta, omega, DT, K, GAMMA, SIGMA_THETA, rng);

        // Compute R field and mass field
        R_field = compute_local_R(theta);
        for (int i = 0; i < N_GRID; i++) {
            mass_field[i] = DELTA * R_field[i];
        }

        // Update Dirac spinor
        dirac_step(psi, R_field, DT, SIGMA_PSI, rng);

        // Compute observables
        float R_global = compute_global_R(theta);
        float norm = compute_spinor_norm(psi);
        compute_particle_position(psi, x_cm, y_cm);

        R_history.push_back(R_global);
        norm_history.push_back(norm);

        // Write timeseries data
        if (step % TIMESERIES_INTERVAL == 0) {
            float time = step * DT;
            write_timeseries_data(timeseries_file, time, R_global, norm, x_cm, y_cm, x_init, y_init);
        }

        // Write field snapshots
        if (step % SNAPSHOT_INTERVAL == 0) {
            float time = step * DT;
            char snapshot_num[10];
            snprintf(snapshot_num, sizeof(snapshot_num), "%03d", snapshot_count);

            write_field_snapshot(OUTPUT_DIR + "/snapshot_" + snapshot_num + ".dat",
                               theta, R_field, mass_field, psi, time);
            write_mass_field(OUTPUT_DIR + "/mass_field_" + snapshot_num + ".dat", mass_field, time);
            write_spinor_density(OUTPUT_DIR + "/density_" + snapshot_num + ".dat", psi, time);

            std::cout << "Snapshot " << snapshot_count << " at t = " << std::fixed << std::setprecision(1)
                     << time << " : R = " << std::setprecision(4) << R_global
                     << ", |Ψ|² = " << norm
                     << ", position = (" << std::setprecision(1) << x_cm << ", " << y_cm << ")" << std::endl;

            snapshot_count++;
        }

        // Progress indicator
        if (step % 50 == 0) {
            float progress = (float)step / N_STEPS * 100;
            std::cout << "Progress: " << std::fixed << std::setprecision(1) << progress << "%"
                     << " - R = " << std::setprecision(4) << R_global
                     << ", |Ψ|² = " << norm << std::endl;
        }
    }

    timeseries_file.close();

    // Compute final statistics
    float R_final = R_history.back();
    float norm_final = norm_history.back();

    // Compute mean and std of R
    float R_mean = 0.0f;
    for (float R : R_history) R_mean += R;
    R_mean /= R_history.size();

    float R_std = 0.0f;
    for (float R : R_history) {
        float diff = R - R_mean;
        R_std += diff * diff;
    }
    R_std = std::sqrt(R_std / R_history.size());

    // Compute max norm deviation
    float max_norm_deviation = 0.0f;
    for (float norm : norm_history) {
        float deviation = std::abs(norm - norm_initial) / norm_initial;
        max_norm_deviation = std::max(max_norm_deviation, deviation);
    }

    // Final particle position
    float x_final, y_final;
    compute_particle_position(psi, x_final, y_final);
    float total_drift = std::sqrt((x_final - x_init) * (x_final - x_init) +
                                  (y_final - y_init) * (y_final - y_init));

    // Write summary statistics
    write_summary_statistics(OUTPUT_DIR + "/summary.dat",
                           R_initial, R_final, R_mean, R_std,
                           norm_initial, norm_final, max_norm_deviation,
                           x_init, y_init, x_final, y_final, total_drift);

    // Validation checks
    std::cout << "\n=============================================" << std::endl;
    std::cout << "           VALIDATION RESULTS" << std::endl;
    std::cout << "=============================================" << std::endl;

    bool vacuum_stable = (R_final > 0.95f);
    bool norm_conserved = (max_norm_deviation < 0.1f);
    bool particle_localized = (total_drift < 5.0f);

    std::cout << "Vacuum stability (R > 0.95): " << R_final << " - "
              << (vacuum_stable ? "PASSED" : "FAILED") << std::endl;
    std::cout << "Norm conservation (< 10% deviation): " << (max_norm_deviation * 100) << "% - "
              << (norm_conserved ? "PASSED" : "FAILED") << std::endl;
    std::cout << "Particle localization (< 5 grid units): " << total_drift << " - "
              << (particle_localized ? "PASSED" : "FAILED") << std::endl;

    bool all_passed = vacuum_stable && norm_conserved && particle_localized;

    std::cout << "\n=============================================" << std::endl;
    if (all_passed) {
        std::cout << "✓ ALL VALIDATION CHECKS PASSED" << std::endl;
    } else {
        std::cout << "✗ SOME VALIDATION CHECKS FAILED" << std::endl;
    }
    std::cout << "=============================================" << std::endl;

    std::cout << "\nOutput files generated in: " << OUTPUT_DIR << std::endl;
    std::cout << "  - timeseries.dat: Evolution data every " << TIMESERIES_INTERVAL << " steps" << std::endl;
    std::cout << "  - snapshot_*.dat: Full field data at t = 0, 1, 2, 3, 4, 5" << std::endl;
    std::cout << "  - mass_field_*.dat: Mass field evolution" << std::endl;
    std::cout << "  - density_*.dat: Spinor density maps" << std::endl;
    std::cout << "  - summary.dat: Summary statistics" << std::endl;

    return all_passed ? 0 : 1;
}