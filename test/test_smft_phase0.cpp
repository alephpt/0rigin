/**
 * Phase 0: Deterministic Characterization
 *
 * Run GPU SMFT for 1000 timesteps and save:
 * - theta(x,y,t) - phase field evolution
 * - R(x,y,t) - synchronization field evolution
 * - Statistics: R_avg(t), R_min(t), R_max(t)
 *
 * For offline analysis:
 * - Lyapunov exponent
 * - Power spectrum
 * - Autocorrelation
 */

#include "../lib/Nova/Nova.h"
#include "../src/SMFTEngine.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sys/stat.h>

void create_output_dir() {
    mkdir("output", 0755);
}

void save_field(const std::vector<float>& field, int Nx, int Ny, int step, const std::string& name) {
    char filename[256];
    snprintf(filename, sizeof(filename), "output/%s_step%04d.dat", name.c_str(), step);

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

void save_timeseries(const std::vector<float>& data, const std::string& name) {
    std::ofstream f("output/" + name + "_timeseries.dat");
    if (!f) {
        std::cerr << "Failed to open output/" << name << "_timeseries.dat" << std::endl;
        return;
    }

    for (size_t i = 0; i < data.size(); i++) {
        f << i << " " << data[i] << "\n";
    }
    f.close();
}

int main() {
    std::cout << "=== SMFT Phase 0: Deterministic Characterization ===" << std::endl;
    std::cout << "GPU-accelerated, multi-timestep batch mode\n" << std::endl;

    create_output_dir();

    // Initialize Nova with minimal config
    NovaConfig config = {
        .name = "SMFT Phase 0 - Batch Processing",
        .screen = {800, 600},
        .debug_level = "info",
        .dimensions = "2D",
        .camera_type = "orthographic",
        .compute = true,
    };
    Nova nova(config);
    nova.initialized = true;

    // Simulation parameters
    const int Nx = 128;
    const int Ny = 128;
    const int N_steps = 10000;  // Extended to 10k steps (t_max = 100)
    const float dt = 0.01f;
    const float K = 1.0f;
    const float Delta = 2.5f;
    const float chiral_angle = 0.0f;
    const int save_every = 100;  // Save field every 100 steps (10x less frequent)

    std::cout << "Grid: " << Nx << " x " << Ny << std::endl;
    std::cout << "Steps: " << N_steps << std::endl;
    std::cout << "dt: " << dt << std::endl;
    std::cout << "K: " << K << std::endl;
    std::cout << "Delta: " << Delta << std::endl;
    std::cout << "Save interval: " << save_every << " steps\n" << std::endl;

    // Initialize SMFT engine
    SMFTEngine engine(&nova);
    engine.initialize(Nx, Ny, Delta, chiral_angle);

    // Set initial conditions (random phases)
    std::vector<float> theta_init(Nx * Ny);
    srand(42);  // Fixed seed for reproducibility
    for (int i = 0; i < Nx * Ny; i++) {
        theta_init[i] = (float(rand()) / RAND_MAX) * 2.0f * M_PI - M_PI;
    }
    engine.setInitialPhases(theta_init);

    // Time series storage
    std::vector<float> R_avg_series;
    std::vector<float> R_min_series;
    std::vector<float> R_max_series;

    std::cout << "Running simulation..." << std::endl;

    auto start_time = std::chrono::high_resolution_clock::now();

    for (int step = 0; step < N_steps; step++) {
        // Run physics step
        engine.step(dt, K, 0.1f);

        // Get R field
        std::vector<float> R_field = engine.getSyncField();

        // Compute statistics
        float R_sum = 0.0f;
        float R_min = 1e10f;
        float R_max = -1e10f;

        for (float r : R_field) {
            R_sum += r;
            if (r < R_min) R_min = r;
            if (r > R_max) R_max = r;
        }

        float R_avg = R_sum / R_field.size();

        R_avg_series.push_back(R_avg);
        R_min_series.push_back(R_min);
        R_max_series.push_back(R_max);

        // Progress output
        if (step % 500 == 0) {
            std::cout << "Step " << step << ": R_avg = " << R_avg
                      << ", R_min = " << R_min << ", R_max = " << R_max << std::endl;
        }

        // Save field snapshots
        if (step % save_every == 0) {
            save_field(R_field, Nx, Ny, step, "R_field");

            // Also save theta field for Lyapunov analysis
            std::vector<float> theta_field = engine.getPhaseField();
            save_field(theta_field, Nx, Ny, step, "theta");
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    std::cout << "\nSimulation completed in " << duration.count() / 1000.0 << " seconds" << std::endl;
    std::cout << "Steps per second: " << N_steps * 1000.0 / duration.count() << std::endl;

    // Save time series
    std::cout << "\nSaving time series data..." << std::endl;
    save_timeseries(R_avg_series, "R_avg");
    save_timeseries(R_min_series, "R_min");
    save_timeseries(R_max_series, "R_max");

    // Save final state
    std::vector<float> R_final = engine.getSyncField();
    std::vector<float> theta_final = engine.getPhaseField();
    save_field(R_final, Nx, Ny, N_steps, "R_field");
    save_field(theta_final, Nx, Ny, N_steps, "theta");

    std::cout << "\nOutput saved to output/ directory" << std::endl;
    std::cout << "Ready for Phase 0 analysis:" << std::endl;
    std::cout << "  - Lyapunov exponent (from theta evolution)" << std::endl;
    std::cout << "  - Power spectrum (from R_avg timeseries)" << std::endl;
    std::cout << "  - Autocorrelation (from R_avg timeseries)" << std::endl;

    return 0;
}
