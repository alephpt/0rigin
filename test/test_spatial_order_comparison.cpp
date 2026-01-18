// test_spatial_order_comparison.cpp
// Compare 2nd-order vs 4th-order spatial discretization for energy conservation

#include "ConservativeSolver.h"
#include <iostream>
#include <iomanip>

void runTest(ConservativeSolver::SpatialOrder order, const char* name) {
    std::cout << "\n=========================================="  << std::endl;
    std::cout << " Testing: " << name << std::endl;
    std::cout << "==========================================" << std::endl;

    ConservativeSolver::Config config;
    config.nx = 64;
    config.ny = 64;
    config.nz = 64;
    config.dx = 1.0f;
    config.dt = 0.005f;
    config.method = ConservativeSolver::IntegrationMethod::STRANG_SPLITTING;
    config.spatial_order = order;

    ConservativeSolver solver;
    solver.initialize(config);

    // Initialize Gaussian wave packet
    solver.initializeGaussian(32.0f, 32.0f, 32.0f, 5.0f, 0.1f);

    float E_initial = solver.computeTotalEnergy();
    std::cout << "Initial energy: " << std::fixed << std::setprecision(6) << E_initial << std::endl;

    // Run 2000 steps
    const int num_steps = 2000;
    std::cout << "Running " << num_steps << " steps..." << std::endl;

    for (int step = 0; step < num_steps; ++step) {
        solver.evolveSineGordon(config.dt);

        // Print energy every 500 steps
        if ((step + 1) % 500 == 0) {
            float E_current = solver.computeTotalEnergy();
            float drift = std::abs(E_current - E_initial) / E_initial * 100.0f;
            std::cout << "  Step " << std::setw(4) << (step + 1)
                      << ": E = " << std::fixed << std::setprecision(6) << E_current
                      << ", drift = " << std::setprecision(4) << drift << "%" << std::endl;
        }
    }

    // Final validation
    std::cout << "\nFinal validation:" << std::endl;
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

    // Time reversibility
    bool reversibility_ok = solver.validateTimeReversibility(1e-4f);

    return;
}

int main() {
    std::cout << "==========================================" << std::endl;
    std::cout << " Spatial Order Comparison Test" << std::endl;
    std::cout << " Objective: <0.01% energy conservation" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << std::endl;
    std::cout << "Configuration:" << std::endl;
    std::cout << "  Grid: 64×64×64" << std::endl;
    std::cout << "  dx = 1.0, dt = 0.005" << std::endl;
    std::cout << "  Integrator: Strang Splitting (T-V-T)" << std::endl;
    std::cout << "  Initial conditions: Gaussian (σ=5, A=0.1)" << std::endl;
    std::cout << "  Duration: 2000 steps (t = 10)" << std::endl;
    std::cout << std::endl;

    // Test 2nd-order baseline
    runTest(ConservativeSolver::SpatialOrder::SECOND_ORDER,
            "2nd-order spatial (6-neighbor Laplacian)");

    std::cout << std::endl;
    std::cout << "==========================================\n" << std::endl;

    // Test 4th-order upgrade
    runTest(ConservativeSolver::SpatialOrder::FOURTH_ORDER,
            "4th-order spatial (12-neighbor Laplacian)");

    std::cout << "\n==========================================" << std::endl;
    std::cout << " Comparison Summary" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << std::endl;
    std::cout << "Expected results:" << std::endl;
    std::cout << "  2nd-order: ~0.068% energy drift (baseline)" << std::endl;
    std::cout << "  4th-order: <0.01% energy drift (target)" << std::endl;
    std::cout << std::endl;
    std::cout << "Improvement factor: ~10-70× expected" << std::endl;
    std::cout << "==========================================" << std::endl;

    return 0;
}
