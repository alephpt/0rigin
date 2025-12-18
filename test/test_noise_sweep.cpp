/**
 * @file test_noise_sweep.cpp
 * @brief Path B Experimental Protocol: Noise Sweep (σ ∈ [0, 1])
 *
 * From Directive.md:
 * - Phase 1: Validation (Week 1) - PRNG verification
 * - Phase 2: Noise Sweep (Week 2) - σ_c measurement
 * - Phase 3: Interpretation (Week 3) - Path A vs Path B decision
 *
 * Success Criteria:
 * - If σ_c > 10^-5: Path B (stochastic) validated ✓
 * - If σ_c < 10^-5: Path A (deterministic) required ✗
 *
 * Outputs: build/output/noise_sweep/
 */

#include "../lib/Nova/Nova.h"
#include "../src/SMFTEngine.h"
#include "../src/SMFTCommon.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sys/stat.h>
#include <sys/types.h>

void create_output_dirs() {
    mkdir("/home/persist/neotec/0rigin/output", 0755);
    mkdir("/home/persist/neotec/0rigin/output/noise_sweep", 0755);
}

void save_field(const std::vector<float>& field, int Nx, int Ny,
                const std::string& dir, const std::string& name) {
    std::string filename = dir + "/" + name + ".dat";
    std::ofstream f(filename);
    if (!f) {
        std::cerr << "Failed to open " << filename << std::endl;
        return;
    }

    for (int y = 0; y < Ny; y++) {
        for (int x = 0; x < Nx; x++) {
            f << field[y * Nx + x];
            if (x < Nx - 1) f << " ";
        }
        f << "\n";
    }
    f.close();
}

void save_timeseries(const std::vector<float>& data, const std::string& filename) {
    std::ofstream f(filename);
    if (!f) {
        std::cerr << "Failed to open " << filename << std::endl;
        return;
    }

    for (size_t i = 0; i < data.size(); i++) {
        f << i << " " << data[i] << "\n";
    }
    f.close();
}

struct Statistics {
    float mean;
    float min;
    float max;
    float variance;
};

Statistics compute_stats(const std::vector<float>& field) {
    Statistics s = {0, 1e10, -1e10, 0};

    for (float val : field) {
        s.mean += val;
        if (val < s.min) s.min = val;
        if (val > s.max) s.max = val;
    }
    s.mean /= field.size();

    for (float val : field) {
        float diff = val - s.mean;
        s.variance += diff * diff;
    }
    s.variance /= field.size();

    return s;
}

float compute_localization(const std::vector<float>& R_field) {
    // L = ∫|R|⁴ dx (measure of localization)
    float L = 0.0f;
    for (float r : R_field) {
        L += r * r * r * r;
    }
    return L;
}

// compute_global_R now provided by SMFTCommon.h

int main() {
    std::cout << "=== SMFT Noise Sweep Experiment ===" << std::endl;
    std::cout << "Directive.md Phase 2: Measuring σ_c\n" << std::endl;

    create_output_dirs();

    // Noise sweep values (logarithmic spacing)
    std::vector<float> sigma_values = {
        0.0,      // Control: deterministic
        1e-6,     // Very small
        1e-5,     // Threshold for Path B validation
        1e-4,
        1e-3,
        1e-2,
        1e-1,
        1.0       // Order unity
    };

    // Simulation parameters (from previous successful runs)
    const int Nx = 256;
    const int Ny = 256;
    const int N_warmup = 5000;    // Steps to reach equilibrium
    const int N_measure = 1000;   // Steps to measure statistics
    const float dt = 0.01f;
    const float K = 27.21f;       // Known stable coupling
    const float omega_mean = 0.0f;

    std::cout << "Grid: " << Nx << " x " << Ny << std::endl;
    std::cout << "Warmup steps: " << N_warmup << std::endl;
    std::cout << "Measurement steps: " << N_measure << std::endl;
    std::cout << "K: " << K << std::endl;
    std::cout << "σ values: " << sigma_values.size() << " points\n" << std::endl;

    // Results storage
    std::ofstream results("/home/persist/neotec/0rigin/output/noise_sweep/results.csv");
    results << "sigma,R_global_mean,R_global_std,R_local_mean,R_local_std,L_mean,L_std,phase_variance\n";

    // Run sweep
    for (float sigma : sigma_values) {
        std::cout << "\n=== Running σ = " << sigma << " ===" << std::endl;

        // Create output directory for this sigma
        char dirname[512];
        snprintf(dirname, sizeof(dirname), "/home/persist/neotec/0rigin/output/noise_sweep/sigma_%.2e", sigma);
        mkdir(dirname, 0755);

        // Initialize Nova (compute-only mode)
        NovaConfig config = {
            .name = "SMFT Noise Sweep",
            .screen = {800, 600},
            .debug_level = "info",
            .dimensions = "2D",
            .camera_type = "orthographic",
            .compute = true,
        };
        Nova nova(config);
        nova.initialized = true;

        // Initialize SMFT engine
        SMFTEngine engine(&nova);
        engine.initialize(Nx, Ny, 2.5f, 0.0f);

        // Set initial conditions (random phases)
        std::vector<float> theta_init(Nx * Ny);
        srand(42);  // Fixed seed for reproducibility
        for (int i = 0; i < Nx * Ny; i++) {
            theta_init[i] = (float(rand()) / RAND_MAX) * 2.0f * M_PI - M_PI;
        }
        engine.setInitialPhases(theta_init);

        // Warmup phase (no noise) - let system equilibrate
        std::cout << "Warmup phase (σ=0)..." << std::flush;
        for (int step = 0; step < N_warmup; step++) {
            engine.step(dt, K, 0.0f);  // No noise during warmup

            if (step % 1000 == 0) {
                std::cout << "." << std::flush;
            }
        }
        std::cout << " done" << std::endl;

        // Check warmup success (compute R_global)
        std::vector<float> theta_after_warmup = engine.getPhaseField();
        float R_global_warmup = SMFT::compute_global_R(theta_after_warmup);
        std::cout << "\n  R_global after warmup: " << R_global_warmup << std::endl;

        if (R_global_warmup < 0.5) {
            std::cerr << "  WARNING: Warmup failed to synchronize (R_global < 0.5)" << std::endl;
            std::cerr << "  This invalidates the noise stability test!" << std::endl;
        }

        // Measurement phase (with noise)
        std::cout << "Measurement phase (σ=" << sigma << ")..." << std::flush;

        std::vector<float> L_series;
        std::vector<float> R_global_series;  // CRITICAL: Global order parameter
        std::vector<float> R_local_mean_series;  // Spatial average of local R
        std::vector<float> phase_var_series;

        for (int step = 0; step < N_measure; step++) {
            // Step with stochastic noise (damping=0, sigma_theta=sigma, sigma_psi=0)
            engine.stepStochastic(dt, K, 0.0f, sigma, 0.0f);

            // Get fields
            std::vector<float> R_field = engine.getSyncField();
            std::vector<float> theta_field = engine.getPhaseField();

            // Compute observables
            float R_global = SMFT::compute_global_R(theta_field);  // CRITICAL: True order parameter
            float L = compute_localization(R_field);
            Statistics R_local_stats = compute_stats(R_field);  // Local R average (for SMFT spatial structure)
            Statistics theta_stats = compute_stats(theta_field);

            R_global_series.push_back(R_global);  // What we need for falsification test
            R_local_mean_series.push_back(R_local_stats.mean);  // For SMFT defect analysis
            L_series.push_back(L);
            phase_var_series.push_back(theta_stats.variance);

            if (step % 200 == 0) {
                std::cout << "." << std::flush;
            }

            // Save snapshots at key points
            if (step == 0 || step == N_measure / 2 || step == N_measure - 1) {
                char stepdir[512];
                snprintf(stepdir, sizeof(stepdir), "%s/step_%04d", dirname, step);
                mkdir(stepdir, 0755);

                save_field(R_field, Nx, Ny, stepdir, "R_field");
                save_field(theta_field, Nx, Ny, stepdir, "theta_field");
            }
        }
        std::cout << " done" << std::endl;

        // Compute statistics over measurement window
        Statistics R_global_stats = compute_stats(R_global_series);
        Statistics R_local_stats = compute_stats(R_local_mean_series);
        Statistics L_stats = compute_stats(L_series);
        float phase_var_mean = 0.0f;
        for (float v : phase_var_series) phase_var_mean += v;
        phase_var_mean /= phase_var_series.size();

        std::cout << "Results:" << std::endl;
        std::cout << "  R_global = " << R_global_stats.mean << " ± " << sqrt(R_global_stats.variance) << " (Kuramoto order parameter)" << std::endl;
        std::cout << "  R_local  = " << R_local_stats.mean << " ± " << sqrt(R_local_stats.variance) << " (Spatial avg of local sync)" << std::endl;
        std::cout << "  L = " << L_stats.mean << " ± " << sqrt(L_stats.variance) << std::endl;
        std::cout << "  Phase variance = " << phase_var_mean << std::endl;

        // Save time series
        save_timeseries(R_global_series, std::string(dirname) + "/R_global_timeseries.dat");
        save_timeseries(R_local_mean_series, std::string(dirname) + "/R_local_timeseries.dat");
        save_timeseries(L_series, std::string(dirname) + "/L_timeseries.dat");
        save_timeseries(phase_var_series, std::string(dirname) + "/phase_var_timeseries.dat");

        // Append to results CSV
        results << sigma << ","
                << R_global_stats.mean << "," << sqrt(R_global_stats.variance) << ","
                << R_local_stats.mean << "," << sqrt(R_local_stats.variance) << ","
                << L_stats.mean << "," << sqrt(L_stats.variance) << ","
                << phase_var_mean << "\n";
        results.flush();
    }

    results.close();

    std::cout << "\n=== Noise Sweep Complete ===" << std::endl;
    std::cout << "Results saved to: build/output/noise_sweep/" << std::endl;
    std::cout << "\nNext steps:" << std::endl;
    std::cout << "1. Plot L(σ) to identify σ_c" << std::endl;
    std::cout << "2. Fit exponential: L(σ) = L₀·exp(-σ/σ_c)" << std::endl;
    std::cout << "3. Decision:" << std::endl;
    std::cout << "   - If σ_c > 10⁻⁵: Path B validated ✓" << std::endl;
    std::cout << "   - If σ_c < 10⁻⁵: Path A required ✗" << std::endl;

    return 0;
}
