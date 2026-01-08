// test_grid_resolution.cpp
// Test energy conservation with finer grid resolution

#include "ConservativeSolver.h"
#include <iostream>
#include <iomanip>

void runTest(uint32_t nx, const char* name) {
    std::cout << "\n=========================================="  << std::endl;
    std::cout << " Testing: " << name << " (" << nx << "³ grid)" << std::endl;
    std::cout << "==========================================" << std::endl;

    ConservativeSolver::Config config;
    config.nx = nx;
    config.ny = nx;
    config.nz = nx;
    config.dx = 64.0f / nx;  // Keep domain size constant
    config.dt = 0.005f;
    config.method = ConservativeSolver::IntegrationMethod::STRANG_SPLITTING;
    config.spatial_order = ConservativeSolver::SpatialOrder::SECOND_ORDER;  // Use 2nd-order

    ConservativeSolver solver;
    solver.initialize(config);

    // Initialize Gaussian wave packet (scale sigma with grid)
    float sigma = 5.0f * (64.0f / nx);  // Keep same number of grid points in Gaussian
    solver.initializeGaussian(nx/2.0f * config.dx, nx/2.0f * config.dx, nx/2.0f * config.dx,
                              sigma, 0.1f);

    float E_initial = solver.computeTotalEnergy();

    // Run 2000 steps
    const int num_steps = 2000;
    std::cout << "Running " << num_steps << " steps..." << std::endl;

    for (int step = 0; step < num_steps; ++step) {
        solver.evolveSineGordon(config.dt);
    }

    // Final validation
    float E_final = solver.computeTotalEnergy();
    float final_drift = std::abs(E_final - E_initial) / E_initial * 100.0f;

    std::cout << "  Grid: " << nx << "×" << nx << "×" << nx << std::endl;
    std::cout << "  dx = " << config.dx << std::endl;
    std::cout << "  Initial energy: " << std::fixed << std::setprecision(6) << E_initial << std::endl;
    std::cout << "  Final energy:   " << std::fixed << std::setprecision(6) << E_final << std::endl;
    std::cout << "  Energy drift:   " << std::setprecision(4) << final_drift << "%" << std::endl;

    // Check against GO/NO-GO threshold
    const float threshold = 0.01f;  // 0.01% threshold
    if (final_drift < threshold) {
        std::cout << "  Status: ✓ PASS (drift < " << threshold << "%)" << std::endl;
    } else {
        std::cout << "  Status: ✗ FAIL (drift > " << threshold << "%)" << std::endl;
    }
}

int main() {
    std::cout << "==========================================" << std::endl;
    std::cout << " Grid Resolution Test" << std::endl;
    std::cout << " Objective: <0.01% energy conservation" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << std::endl;
    std::cout << "Theory: Energy error scales as O(dx²)" << std::endl;
    std::cout << "Doubling resolution → 4× better conservation" << std::endl;
    std::cout << std::endl;

    // Test different resolutions
    runTest(64, "Baseline");
    runTest(96, "1.5× resolution");
    runTest(128, "2× resolution");

    std::cout << "\n==========================================" << std::endl;
    std::cout << " Summary" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << "Expected:" << std::endl;
    std::cout << "  64³:  ~0.068% drift (baseline)" << std::endl;
    std::cout << "  96³:  ~0.030% drift (2.25× better)" << std::endl;
    std::cout << "  128³: ~0.017% drift (4× better)" << std::endl;
    std::cout << "==========================================" << std::endl;

    return 0;
}
