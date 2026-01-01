// test_em_wave_simple.cpp - Standalone EM wave propagation test
#include "SMFTEngine.h"
#include "Nova.h"
#include <iostream>
#include <cmath>

int main() {
    std::cout << "=== EM Wave Propagation Test (G1) ===" << std::endl;
    std::cout << "Goal: Verify c = 1/√(μ₀ε₀) = 1.0 in natural units" << std::endl;
    std::cout << std::endl;

    // Initialize Nova (required for Vulkan)
    NovaConfig config{};
    config.name = "EM Wave Test";
    config.screen = {800, 600};
    config.debug_level = "error";
    config.dimensions = "2D";
    config.camera_type = "orbit";
    config.compute = true;

    Nova nova(config);
    if (!nova.initialized) {
        std::cerr << "Failed to initialize Nova" << std::endl;
        return 1;
    }

    // Initialize SMFT Engine
    SMFTEngine engine(&nova);

    const int Nx = 256, Ny = 64;
    const float Delta = 0.0f;  // No mass gap
    const float chiral_angle = 0.0f;

    engine.initialize(Nx, Ny, Delta, chiral_angle);

    // Zero phases (no Kuramoto dynamics)
    std::vector<float> theta(Nx * Ny, 0.0f);
    std::vector<float> omega(Nx * Ny, 0.0f);
    engine.setInitialPhases(theta);
    engine.setNaturalFrequencies(omega);

    std::cout << "Grid: " << Nx << " x " << Ny << std::endl;
    std::cout << "EM coupling: Stückelberg (massless photon)" << std::endl;
    std::cout << std::endl;

    // Initialize EM pulse
    float center_x = 32.0f;
    float center_y = 32.0f;
    float width_x = 8.0f;
    float width_y = 8.0f;
    float amplitude = 0.1f;
    float k_x = 10.0f;  // Rightward propagation
    float k_y = 0.0f;

    std::cout << "Initializing Gaussian EM pulse:" << std::endl;
    std::cout << "  Position: (" << center_x << ", " << center_y << ")" << std::endl;
    std::cout << "  Width: (" << width_x << ", " << width_y << ")" << std::endl;
    std::cout << "  Amplitude: " << amplitude << std::endl;
    std::cout << "  Wave vector: k = (" << k_x << ", " << k_y << ")" << std::endl;
    std::cout << std::endl;

    engine.initializeEMPulse(center_x, center_y, width_x, width_y,
                            amplitude, k_x, k_y, "A_y");

    // Measure wave velocity
    const float dt = 0.001f;
    const int steps = 100;

    std::cout << "Measuring wave velocity..." << std::endl;
    std::cout << "  Time steps: " << steps << " (dt = " << dt << ")" << std::endl;
    std::cout << "  Total time: " << (steps * dt) << std::endl;
    std::cout << std::endl;

    float velocity = engine.measureWaveVelocity(dt, steps);

    // Validation
    float expected_c = 1.0f;
    float error = std::abs(velocity - expected_c) / expected_c;
    bool passed = (error < 0.01f);  // 1% tolerance

    std::cout << "\n===== Results =====" << std::endl;
    std::cout << "Expected velocity (c): " << expected_c << std::endl;
    std::cout << "Measured velocity: " << velocity << std::endl;
    std::cout << "Relative error: " << (error * 100.0f) << "%" << std::endl;
    std::cout << "Status: " << (passed ? "PASSED ✓" : "FAILED ✗") << std::endl;

    return passed ? 0 : 1;
}
