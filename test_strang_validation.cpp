// test_strang_validation.cpp
// Standalone test to validate Strang splitting energy conservation

#include "ConservativeSolver.h"
#include <iostream>
#include <iomanip>

int main() {
    std::cout << "==========================================" << std::endl;
    std::cout << " Strang Splitting Validation Test" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << std::endl;

    // Configuration
    ConservativeSolver::Config config;
    config.nx = 128;  // Doubled resolution
    config.ny = 128;
    config.nz = 128;
    config.dx = 0.5f;  // Half grid spacing
    config.dt = 0.0005f;  // Smaller timestep (CFL condition)
    config.method = ConservativeSolver::IntegrationMethod::STRANG_SPLITTING;

    ConservativeSolver solver;
    solver.initialize(config);

    // Initialize Gaussian wave packet (smooth initial conditions)
    float x0 = 64.0f;  // Center of 128³ grid
    float y0 = 64.0f;
    float z0 = 64.0f;
    float sigma = 10.0f;  // Doubled to match new grid spacing
    float amplitude = 0.1f;  // Smaller amplitude to reduce nonlinear effects

    std::cout << "Initializing Gaussian wave packet:" << std::endl;
    std::cout << "  Center: (" << x0 << ", " << y0 << ", " << z0 << ")" << std::endl;
    std::cout << "  Sigma: " << sigma << std::endl;
    std::cout << "  Amplitude: " << amplitude << std::endl;
    std::cout << std::endl;

    solver.initializeGaussian(x0, y0, z0, sigma, amplitude);

    float initial_energy = solver.computeTotalEnergy();
    std::cout << "Initial energy: " << initial_energy << std::endl;
    std::cout << std::endl;

    // Run simulation for 20000 steps (same total time: 20000*0.0005 = 10)
    const int num_steps = 20000;
    const int report_interval = 2000;

    std::cout << "Running " << num_steps << " timesteps..." << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << std::endl;

    std::cout << std::setw(10) << "Step"
              << std::setw(20) << "Energy"
              << std::setw(20) << "Drift (%)" << std::endl;
    std::cout << std::string(50, '-') << std::endl;

    for (int step = 1; step <= num_steps; ++step) {
        solver.evolveSineGordon(config.dt);

        if (step % report_interval == 0 || step == num_steps) {
            float current_energy = solver.computeTotalEnergy();
            float drift_percent = std::abs(current_energy - initial_energy) / initial_energy * 100.0f;

            std::cout << std::setw(10) << step
                      << std::setw(20) << current_energy
                      << std::setw(20) << drift_percent << std::endl;
        }
    }

    std::cout << std::endl;

    // Final validation
    std::cout << "==========================================" << std::endl;
    std::cout << " Final Validation" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << std::endl;

    bool energy_passed = solver.validateEnergyConservation(0.0001f);  // 0.01% threshold
    bool reversibility_passed = solver.validateTimeReversibility(1e-4f);

    std::cout << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << " Test Result" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << std::endl;

    if (energy_passed && reversibility_passed) {
        std::cout << "✓ TEST PASSED" << std::endl;
        std::cout << "  - Energy conservation: ✓" << std::endl;
        std::cout << "  - Time reversibility: ✓" << std::endl;
        return 0;
    } else {
        std::cout << "✗ TEST FAILED" << std::endl;
        if (!energy_passed) {
            std::cout << "  - Energy conservation: ✗ FAILED" << std::endl;
        }
        if (!reversibility_passed) {
            std::cout << "  - Time reversibility: ✗ FAILED" << std::endl;
        }
        return 1;
    }
}
