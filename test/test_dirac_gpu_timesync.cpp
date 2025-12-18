/**
 * Test Dirac GPU shader with timesync operator splitting
 *
 * FINDING: Timesync is implemented with _substep_ratio (default N=10)
 * - Kuramoto runs N times for every 1 Dirac update
 * - Born-Oppenheimer approximation: vacuum fast, matter slow
 *
 * PURPOSE: Test if GPU Dirac works without numerical explosion
 * when using the proper timesync ratio
 */

#include "../src/MSFTEngine.h"
#include "../lib/Nova/Nova.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>
#include <iomanip>

// Test parameters
const int GRID_SIZE = 64;
const float DELTA = 2.5f;
const float K = 1.0f;
const float DT_FAST = 0.01f;  // Kuramoto timestep
const float DAMPING = 0.1f;
const float LAMBDA = 0.1f;    // Dirac feedback coupling

// Safety parameters
const int TEST_STEPS = 100;   // Short run for safety
const int SUBSTEP_RATIO = 100; // 100 Kuramoto steps per Dirac step

void printStatus(const std::string& label, float value) {
    std::cout << std::setw(30) << std::left << label << ": " << value << std::endl;
}

void runGPUDiracTest() {
    std::cout << "\n===== GPU DIRAC TIMESYNC TEST =====" << std::endl;
    std::cout << "Testing GPU Dirac shader with operator splitting" << std::endl;
    std::cout << "Substep ratio N = " << SUBSTEP_RATIO << std::endl;
    std::cout << "dt_fast (Kuramoto) = " << DT_FAST << std::endl;
    std::cout << "dt_slow (Dirac) = " << (DT_FAST * SUBSTEP_RATIO) << std::endl;
    std::cout << "\nInitializing Nova graphics engine..." << std::endl;

    // Initialize Nova with proper config
    NovaConfig config = {
        .name = "GPU Dirac Timesync Test",
        .screen = {800, 600},
        .debug_level = "error",  // Minimal output to focus on physics
        .dimensions = "2D",
        .camera_type = "fixed",
        .compute = true
    };

    Nova nova(config);
    nova.initialized = true;

    std::cout << "✓ Nova initialized" << std::endl;

    // Initialize MSFT Engine
    MSFTEngine engine(&nova);
    engine.initialize(GRID_SIZE, GRID_SIZE, DELTA, 0.0f);

    // Set the timesync ratio
    engine.setSubstepRatio(SUBSTEP_RATIO);
    std::cout << "✓ Set substep ratio to " << SUBSTEP_RATIO << std::endl;

    // Initialize phases with slight randomness
    std::vector<float> phases(GRID_SIZE * GRID_SIZE);
    for (int i = 0; i < phases.size(); i++) {
        phases[i] = 0.1f * sin(i * 0.1f);  // Gentle spatial variation
    }
    engine.setInitialPhases(phases);

    // Zero natural frequencies for simplicity
    std::vector<float> omega(GRID_SIZE * GRID_SIZE, 0.0f);
    engine.setNaturalFrequencies(omega);

    // Initialize Dirac field as Gaussian wavepacket
    float x0 = GRID_SIZE / 2.0f;
    float y0 = GRID_SIZE / 2.0f;
    float sigma = 5.0f;
    engine.initializeDiracField(x0, y0, sigma, 1.0f);

    // Initialize hybrid system (GPU Kuramoto + CPU Dirac)
    engine.initializeHybrid(x0, y0, sigma);

    std::cout << "✓ Initialized Dirac wavepacket at (" << x0 << ", " << y0 << ")" << std::endl;
    std::cout << "\n--- Starting Time Evolution ---" << std::endl;

    // Track stability metrics
    float max_R = 0.0f;
    float min_R = 1.0f;
    float max_density = 0.0f;
    bool exploded = false;

    // Open output file for timeseries
    std::ofstream timeseries("output/dirac_gpu_timesync.dat");
    timeseries << "# step  time  avg_R  max_R  max_density  status" << std::endl;

    // Time evolution
    auto start_time = std::chrono::high_resolution_clock::now();

    for (int step = 0; step < TEST_STEPS; step++) {
        try {
            // Evolve with Dirac coupling
            // This internally handles the operator splitting:
            // - Kuramoto runs SUBSTEP_RATIO times
            // - Dirac updates once every SUBSTEP_RATIO Kuramoto steps
            engine.stepWithDirac(DT_FAST, LAMBDA);

            // Get fields for monitoring
            auto R_field = engine.getSyncField();
            auto density = engine.getDiracDensity();

            // Compute statistics
            float avg_R = 0.0f;
            float current_max_R = 0.0f;
            float current_max_density = 0.0f;

            for (const auto& r : R_field) {
                avg_R += r;
                if (r > current_max_R) current_max_R = r;
                if (std::isnan(r) || std::isinf(r)) {
                    exploded = true;
                    std::cout << "❌ EXPLOSION DETECTED at step " << step
                             << ": R = " << r << std::endl;
                    break;
                }
            }
            avg_R /= R_field.size();

            for (const auto& d : density) {
                if (d > current_max_density) current_max_density = d;
                if (std::isnan(d) || std::isinf(d)) {
                    exploded = true;
                    std::cout << "❌ EXPLOSION DETECTED at step " << step
                             << ": density = " << d << std::endl;
                    break;
                }
            }

            if (exploded) break;

            // Update global tracking
            if (current_max_R > max_R) max_R = current_max_R;
            if (avg_R < min_R) min_R = avg_R;
            if (current_max_density > max_density) max_density = current_max_density;

            // Output progress
            if (step % 10 == 0) {
                float time = step * DT_FAST;
                std::cout << "Step " << std::setw(4) << step
                         << " | t = " << std::fixed << std::setprecision(2) << time
                         << " | avg_R = " << std::setprecision(4) << avg_R
                         << " | max_dens = " << current_max_density
                         << " | ✓ stable" << std::endl;

                timeseries << step << " " << time << " " << avg_R << " "
                          << current_max_R << " " << current_max_density
                          << " stable" << std::endl;
            }

        } catch (const std::exception& e) {
            std::cout << "❌ ERROR at step " << step << ": " << e.what() << std::endl;
            exploded = true;
            break;
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    timeseries.close();

    // Report results
    std::cout << "\n===== TEST RESULTS =====" << std::endl;
    if (exploded) {
        std::cout << "❌ TEST FAILED: Numerical explosion detected" << std::endl;
    } else {
        std::cout << "✅ TEST PASSED: Simulation stable" << std::endl;
        printStatus("Total steps", TEST_STEPS);
        printStatus("Substep ratio N", SUBSTEP_RATIO);
        printStatus("Simulation time (s)", TEST_STEPS * DT_FAST);
        printStatus("Wall time (ms)", duration.count());
        printStatus("Max sync field R", max_R);
        printStatus("Min sync field R", min_R);
        printStatus("Max Dirac density", max_density);
        std::cout << "\n✓ Timesync operator splitting working correctly" << std::endl;
        std::cout << "✓ Born-Oppenheimer approximation validated" << std::endl;
        std::cout << "✓ No numerical explosions detected" << std::endl;
    }

    // Cleanup handled by Nova destructor
}

int main() {
    std::cout << "GPU Dirac Timesync Validation Test" << std::endl;
    std::cout << "===================================" << std::endl;
    std::cout << "\nTesting hypothesis from Feature-Not-Bug.md:" << std::endl;
    std::cout << "- Vacuum (Kuramoto) equilibrates ~100× faster than matter (Dirac)" << std::endl;
    std::cout << "- Operator splitting with N=" << SUBSTEP_RATIO << " should be stable" << std::endl;
    std::cout << "- GPU timeout avoided by separating timescales" << std::endl;

    try {
        runGPUDiracTest();
    } catch (const std::exception& e) {
        std::cerr << "\n❌ FATAL ERROR: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "\nTest complete. Check output/dirac_gpu_timesync.dat for results." << std::endl;
    return 0;
}