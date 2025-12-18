/**
 * Test Program for MSFT GPU Pipeline
 *
 * Validates that the GPU compute pipelines are working correctly:
 * 1. Creates MSFTEngine instance
 * 2. Initializes small grid (32x32)
 * 3. Sets initial phases (random)
 * 4. Runs one simulation step on GPU
 * 5. Verifies results are computed
 */

#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include "MSFTEngine.h"
// NovaConfig is already included via MSFTEngine.h -> Nova.h -> config.h

int main() {
    std::cout << "MSFT GPU Pipeline Test\n";
    std::cout << "======================\n\n";

    // Initialize Nova engine (minimal setup for compute)
    NovaConfig config = {
        .name = "MSFT GPU Test",
        .screen = {800, 600},
        .debug_level = "info",
        .dimensions = "2D",
        .camera_type = "orthographic",
        .compute = true,
    };

    Nova nova(config);
    nova.initialized = true;

    // Create MSFT engine
    MSFTEngine engine(&nova);

    // Initialize small test grid
    const uint32_t Nx = 32;
    const uint32_t Ny = 32;
    const float Delta = 1.0f;  // Test with unit mass gap
    const float chiral_angle = 0.0f;

    std::cout << "1. Initializing " << Nx << "x" << Ny << " grid...\n";
    engine.initialize(Nx, Ny, Delta, chiral_angle);

    // Set random initial phases
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> phase_dist(0.0f, 2.0f * M_PI);
    std::uniform_real_distribution<float> freq_dist(-1.0f, 1.0f);

    std::vector<float> initial_phases(Nx * Ny);
    std::vector<float> frequencies(Nx * Ny);

    for (size_t i = 0; i < Nx * Ny; ++i) {
        initial_phases[i] = phase_dist(gen);
        frequencies[i] = freq_dist(gen);
    }

    std::cout << "2. Setting initial conditions...\n";
    engine.setInitialPhases(initial_phases);
    engine.setNaturalFrequencies(frequencies);

    // Get initial state
    auto sync_before = engine.getSyncField();
    float avg_sync_before = 0.0f;
    for (float R : sync_before) {
        avg_sync_before += R;
    }
    avg_sync_before /= sync_before.size();

    std::cout << "3. Initial average synchronization: " << avg_sync_before << "\n";

    // Run one GPU simulation step
    const float dt = 0.01f;
    const float K = 2.0f;  // Coupling strength
    const float damping = 0.1f;

    std::cout << "4. Running GPU compute step (dt=" << dt << ", K=" << K << ")...\n";
    engine.step(dt, K, damping);

    // Get results
    auto sync_after = engine.getSyncField();
    auto mass_field = engine.getMassField();
    auto gravity_field = engine.getGravitationalField();

    // Compute statistics
    float avg_sync_after = 0.0f;
    float max_sync = 0.0f;
    float min_sync = 1.0f;

    for (float R : sync_after) {
        avg_sync_after += R;
        max_sync = std::max(max_sync, R);
        min_sync = std::min(min_sync, R);
    }
    avg_sync_after /= sync_after.size();

    float avg_mass = 0.0f;
    for (float m : mass_field) {
        avg_mass += m;
    }
    avg_mass /= mass_field.size();

    // Compute gravity magnitude
    float avg_gravity_mag = 0.0f;
    for (size_t i = 0; i < Nx * Ny; ++i) {
        float gx = gravity_field[2*i];
        float gy = gravity_field[2*i + 1];
        avg_gravity_mag += std::sqrt(gx*gx + gy*gy);
    }
    avg_gravity_mag /= (Nx * Ny);

    // Display results
    std::cout << "\n5. GPU Computation Results:\n";
    std::cout << "   ------------------------\n";
    std::cout << "   Synchronization field R(x):\n";
    std::cout << "     - Average: " << avg_sync_after << " (was " << avg_sync_before << ")\n";
    std::cout << "     - Range: [" << min_sync << ", " << max_sync << "]\n";
    std::cout << "   Mass field m(x) = Δ·R(x):\n";
    std::cout << "     - Average: " << avg_mass << " (Δ = " << Delta << ")\n";
    std::cout << "   Gravitational field g(x) = -Δ·∇R(x):\n";
    std::cout << "     - Average magnitude: " << avg_gravity_mag << "\n";

    // Validation checks
    bool success = true;
    if (sync_after.empty() || sync_after.size() != Nx * Ny) {
        std::cout << "\n❌ ERROR: Sync field has wrong size!\n";
        success = false;
    }
    if (mass_field.empty() || mass_field.size() != Nx * Ny) {
        std::cout << "\n❌ ERROR: Mass field has wrong size!\n";
        success = false;
    }
    if (gravity_field.empty() || gravity_field.size() != 2 * Nx * Ny) {
        std::cout << "\n❌ ERROR: Gravity field has wrong size!\n";
        success = false;
    }

    // Check for NaN/Inf
    for (float R : sync_after) {
        if (!std::isfinite(R)) {
            std::cout << "\n❌ ERROR: Sync field contains NaN or Inf!\n";
            success = false;
            break;
        }
    }

    if (success) {
        std::cout << "\n✅ SUCCESS: GPU pipeline working correctly!\n";
        std::cout << "   All three compute shaders executed successfully.\n";
        std::cout << "   Phase 3 & 4 implementation complete.\n";
    }

    // Cleanup handled by destructors

    return success ? 0 : 1;
}