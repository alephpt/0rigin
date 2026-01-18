// test_spatial_order_timestep.cpp
// Test 4th-order spatial with different timesteps

#include "ConservativeSolver.h"
#include <iostream>
#include <iomanip>

void runTest(float dt, ConservativeSolver::SpatialOrder order, const char* name) {
    std::cout << "\n=========================================="  << std::endl;
    std::cout << " Testing: " << name << " (dt = " << dt << ")" << std::endl;
    std::cout << "==========================================" << std::endl;

    ConservativeSolver::Config config;
    config.nx = 64;
    config.ny = 64;
    config.nz = 64;
    config.dx = 1.0f;
    config.dt = dt;
    config.method = ConservativeSolver::IntegrationMethod::STRANG_SPLITTING;
    config.spatial_order = order;

    ConservativeSolver solver;
    solver.initialize(config);

    // Initialize Gaussian wave packet
    solver.initializeGaussian(32.0f, 32.0f, 32.0f, 5.0f, 0.1f);

    float E_initial = solver.computeTotalEnergy();

    // Run to same final time (t = 10)
    const int num_steps = static_cast<int>(10.0f / dt);
    std::cout << "Running " << num_steps << " steps (t = 10)..." << std::endl;

    for (int step = 0; step < num_steps; ++step) {
        solver.evolveSineGordon(config.dt);
    }

    // Final validation
    float E_final = solver.computeTotalEnergy();
    float final_drift = std::abs(E_final - E_initial) / E_initial * 100.0f;

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
    std::cout << " 4th-Order Spatial + Timestep Sweep" << std::endl;
    std::cout << " Objective: <0.01% energy conservation" << std::endl;
    std::cout << "==========================================" << std::endl;

    // Baseline: 2nd-order with dt=0.005
    runTest(0.005f, ConservativeSolver::SpatialOrder::SECOND_ORDER,
            "2nd-order (baseline)");

    // 4th-order with various timesteps
    runTest(0.005f, ConservativeSolver::SpatialOrder::FOURTH_ORDER,
            "4th-order");

    runTest(0.002f, ConservativeSolver::SpatialOrder::FOURTH_ORDER,
            "4th-order (smaller dt)");

    runTest(0.001f, ConservativeSolver::SpatialOrder::FOURTH_ORDER,
            "4th-order (very small dt)");

    std::cout << "\n==========================================" << std::endl;
    std::cout << " Summary" << std::endl;
    std::cout << "==========================================" << std::endl;

    return 0;
}
