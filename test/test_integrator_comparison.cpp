// test_integrator_comparison.cpp
// Compare Velocity Verlet vs Strang Splitting

#include "ConservativeSolver.h"
#include <iostream>
#include <iomanip>

void runTest(ConservativeSolver::IntegrationMethod method, const char* name) {
    std::cout << "\n=========================================="  << std::endl;
    std::cout << " Testing: " << name << std::endl;
    std::cout << "==========================================" << std::endl;

    ConservativeSolver::Config config;
    config.nx = 64;
    config.ny = 64;
    config.nz = 64;
    config.dx = 1.0f;
    config.dt = 0.005f;
    config.method = method;

    ConservativeSolver solver;
    solver.initialize(config);

    // Initialize Gaussian wave packet
    solver.initializeGaussian(32.0f, 32.0f, 32.0f, 5.0f, 0.1f);

    // Run 2000 steps
    for (int step = 0; step < 2000; ++step) {
        solver.evolveSineGordon(config.dt);
    }

    // Validate
    std::cout << "\nFinal validation:" << std::endl;
    bool energy_ok = solver.validateEnergyConservation(0.0001f);
    bool reversibility_ok = solver.validateTimeReversibility(1e-4f);

    if (energy_ok && reversibility_ok) {
        std::cout << "✓ PASSED" << std::endl;
    } else {
        std::cout << "✗ FAILED" << std::endl;
    }
}

int main() {
    std::cout << "==========================================" << std::endl;
    std::cout << " Integrator Comparison Test" << std::endl;
    std::cout << "==========================================" << std::endl;

    runTest(ConservativeSolver::IntegrationMethod::VELOCITY_VERLET, "Velocity Verlet");
    runTest(ConservativeSolver::IntegrationMethod::STRANG_SPLITTING, "Strang Splitting");
    runTest(ConservativeSolver::IntegrationMethod::RK2_SYMPLECTIC, "RK2 Symplectic");

    return 0;
}
