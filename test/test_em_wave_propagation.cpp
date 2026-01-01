// test/test_em_wave_propagation.cpp
// Test electromagnetic wave propagation to verify c = 1/√(μ₀ε₀)

#include "SMFTEngine.h"
#include "simulations/SMFTTestRunner.h"
#include "simulations/TestConfig.h"
#include "Nova.h"
#include <iostream>
#include <cmath>
#include <vector>

/**
 * Test G1: Electromagnetic Wave Propagation
 *
 * Goal: Verify c = 1/√(μ₀ε₀) in SMFT electromagnetic fields
 *
 * Method:
 * 1. Initialize Gaussian EM pulse with known wave vector
 * 2. Evolve using Maxwell equations (wave equation)
 * 3. Track pulse position over time
 * 4. Measure phase velocity v = Δx/Δt
 * 5. Compare to theoretical c = 1.0 (natural units)
 *
 * Success criteria: |v - c| / c < 1%
 */

int runManualWaveTest() {
    std::cout << "\n===== Manual EM Wave Propagation Test =====" << std::endl;

    // Initialize Nova
    NovaConfig config{};
    config.name = "EM Wave Test";
    config.screen = {800, 600};
    config.debug_level = "none";
    config.dimensions = "2D";
    config.camera_type = "orbit";
    config.compute = true;

    Nova nova(config);
    if (!nova.initialized) {
        std::cerr << "Failed to initialize Nova" << std::endl;
        return 1;
    }

    // Create SMFTEngine
    SMFTEngine engine(&nova);

    // Initialize grid (256x64 for clear wave propagation)
    const int Nx = 256, Ny = 64;
    const float Delta = 0.0f; // No mass gap for pure EM test
    const float chiral_angle = 0.0f;

    engine.initialize(Nx, Ny, Delta, chiral_angle);

    // Set zero initial phases (no Kuramoto dynamics)
    std::vector<float> theta(Nx * Ny, 0.0f);
    std::vector<float> omega(Nx * Ny, 0.0f);
    engine.setInitialPhases(theta);
    engine.setNaturalFrequencies(omega);

    std::cout << "\nInitial setup:" << std::endl;
    std::cout << "  Grid: " << Nx << " x " << Ny << std::endl;
    std::cout << "  EM coupling: Stückelberg (massless photon)" << std::endl;
    std::cout << std::endl;

    // Initialize Gaussian pulse
    float pulse_center_x = 32.0f;  // Start near left boundary
    float pulse_center_y = 32.0f;
    float pulse_width_x = 8.0f;
    float pulse_width_y = 8.0f;
    float pulse_amplitude = 0.1f;
    float k_x = 10.0f;  // Wave vector for rightward propagation
    float k_y = 0.0f;

    std::cout << "Initializing EM pulse:" << std::endl;
    std::cout << "  Center: (" << pulse_center_x << ", " << pulse_center_y << ")" << std::endl;
    std::cout << "  Width: (" << pulse_width_x << ", " << pulse_width_y << ")" << std::endl;
    std::cout << "  Amplitude: " << pulse_amplitude << std::endl;
    std::cout << "  Wave vector: (" << k_x << ", " << k_y << ")" << std::endl;
    std::cout << "  Expected phase velocity: c = 1.0" << std::endl;
    std::cout << std::endl;

    engine.initializeEMPulse(pulse_center_x, pulse_center_y,
                            pulse_width_x, pulse_width_y,
                            pulse_amplitude, k_x, k_y, "A_y");

    // Measure wave velocity
    const float dt = 0.001f;
    const int measurement_steps = 100;

    std::cout << "Measuring wave velocity over " << measurement_steps
              << " steps (dt = " << dt << ")..." << std::endl;
    std::cout << std::endl;

    float measured_velocity = engine.measureWaveVelocity(dt, measurement_steps);

    // Check result
    float expected_c = 1.0f;
    float relative_error = std::abs(measured_velocity - expected_c) / expected_c;
    bool passed = relative_error < 0.01f;  // 1% tolerance

    std::cout << "\n===== Test Results =====" << std::endl;
    std::cout << "Expected velocity (c): " << expected_c << std::endl;
    std::cout << "Measured velocity: " << measured_velocity << std::endl;
    std::cout << "Relative error: " << (relative_error * 100.0f) << "%" << std::endl;
    std::cout << "Test status: " << (passed ? "PASSED ✓" : "FAILED ✗") << std::endl;

    return passed ? 0 : 1;
}

int runConfigBasedTest() {
    std::cout << "\n===== Config-Based EM Wave Test =====" << std::endl;

    // Load test configuration
    std::string config_path = "config/em_wave_propagation.yaml";
    TestConfig test_config;

    if (!test_config.loadFromYAML(config_path)) {
        std::cerr << "Failed to load config: " << config_path << std::endl;
        return 1;
    }

    // Run test using SMFTTestRunner
    SMFTTestRunner runner(test_config);

    if (!runner.initialize()) {
        std::cerr << "Failed to initialize test runner" << std::endl;
        return 1;
    }

    bool success = runner.run();

    runner.generateReport();

    return success ? 0 : 1;
}

int main(int argc, char* argv[]) {
    std::cout << "=== TEST G1: Electromagnetic Wave Propagation ===" << std::endl;
    std::cout << "Verifying c = 1/√(μ₀ε₀) in natural units" << std::endl;

    // Run manual test (direct API calls)
    int manual_result = runManualWaveTest();

    // Run config-based test (if config exists)
    std::string config_path = "config/em_wave_propagation.yaml";
    std::ifstream config_check(config_path);
    if (config_check.good()) {
        config_check.close();
        std::cout << "\n--- Running config-based test ---" << std::endl;
        int config_result = runConfigBasedTest();

        if (manual_result == 0 && config_result == 0) {
            std::cout << "\n✓✓✓ ALL EM WAVE TESTS PASSED ✓✓✓" << std::endl;
            return 0;
        }
    }

    if (manual_result == 0) {
        std::cout << "\n✓ Manual EM wave test PASSED" << std::endl;
        return 0;
    } else {
        std::cout << "\n✗ EM wave test FAILED" << std::endl;
        return 1;
    }
}