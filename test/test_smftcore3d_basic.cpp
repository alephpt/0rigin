// test/test_smftcore3d_basic.cpp
#include "SMFTCore3D.h"
#include <iostream>
#include <cassert>
#include <cmath>
#include <chrono>

/**
 * Basic validation test for SMFTCore3D
 *
 * Week 1 Tests:
 *   1. Grid allocation and memory
 *   2. Index mapping (3D <-> 1D)
 *   3. Neighbor access with periodic BCs
 *   4. CPU Kuramoto evolution
 *   5. R-field computation
 */

void test_index_mapping() {
    std::cout << "\n[TEST] Index Mapping Validation" << std::endl;

    SMFTCore3D core;
    SMFTCore3D::Config config;
    config.Nx = 32;
    config.Ny = 32;
    config.Nz = 32;
    core.initialize(config);

    // Test round-trip conversion
    uint32_t test_cases[][3] = {
        {0, 0, 0},      // Origin
        {31, 31, 31},   // Max corner
        {16, 16, 16},   // Center
        {1, 2, 3},      // Arbitrary point
        {31, 0, 0},     // X edge
        {0, 31, 0},     // Y edge
        {0, 0, 31},     // Z edge
    };

    for (const auto& test : test_cases) {
        uint32_t i = test[0], j = test[1], k = test[2];
        uint32_t idx = core.index3D(i, j, k);
        uint32_t i2, j2, k2;
        core.coords3D(idx, i2, j2, k2);

        assert(i == i2 && j == j2 && k == k2);
        std::cout << "  (" << i << "," << j << "," << k << ") -> "
                  << idx << " -> (" << i2 << "," << j2 << "," << k2 << ") ✓" << std::endl;
    }

    std::cout << "[PASS] Index mapping verified" << std::endl;
}

void test_neighbor_access() {
    std::cout << "\n[TEST] Neighbor Access with Periodic BCs" << std::endl;

    SMFTCore3D core;
    SMFTCore3D::Config config;
    config.Nx = 4;  // Small grid for manual verification
    config.Ny = 4;
    config.Nz = 4;
    core.initialize(config);

    // Test corner point (0,0,0) - should wrap around
    auto neighbors = core.getNeighbors(0, 0, 0);

    // Expected neighbors with periodic BC
    assert(neighbors.x_plus == core.index3D(1, 0, 0));
    assert(neighbors.x_minus == core.index3D(3, 0, 0));  // Wraps to 3
    assert(neighbors.y_plus == core.index3D(0, 1, 0));
    assert(neighbors.y_minus == core.index3D(0, 3, 0));  // Wraps to 3
    assert(neighbors.z_plus == core.index3D(0, 0, 1));
    assert(neighbors.z_minus == core.index3D(0, 0, 3));  // Wraps to 3

    std::cout << "  Corner (0,0,0) neighbors verified ✓" << std::endl;

    // Test max corner (3,3,3)
    neighbors = core.getNeighbors(3, 3, 3);
    assert(neighbors.x_plus == core.index3D(0, 3, 3));  // Wraps to 0
    assert(neighbors.x_minus == core.index3D(2, 3, 3));
    assert(neighbors.y_plus == core.index3D(3, 0, 3));  // Wraps to 0
    assert(neighbors.y_minus == core.index3D(3, 2, 3));
    assert(neighbors.z_plus == core.index3D(3, 3, 0));  // Wraps to 0
    assert(neighbors.z_minus == core.index3D(3, 3, 2));

    std::cout << "  Corner (3,3,3) neighbors verified ✓" << std::endl;

    // Test center point (2,2,2)
    neighbors = core.getNeighbors(2, 2, 2);
    assert(neighbors.x_plus == core.index3D(3, 2, 2));
    assert(neighbors.x_minus == core.index3D(1, 2, 2));
    assert(neighbors.y_plus == core.index3D(2, 3, 2));
    assert(neighbors.y_minus == core.index3D(2, 1, 2));
    assert(neighbors.z_plus == core.index3D(2, 2, 3));
    assert(neighbors.z_minus == core.index3D(2, 2, 1));

    std::cout << "  Center (2,2,2) neighbors verified ✓" << std::endl;

    std::cout << "[PASS] Neighbor access with periodic BCs verified" << std::endl;
}

void test_uniform_evolution() {
    std::cout << "\n[TEST] Uniform Phase Evolution" << std::endl;

    SMFTCore3D core;
    SMFTCore3D::Config config;
    config.Nx = 32;
    config.Ny = 32;
    config.Nz = 32;
    config.dt = 0.01f;
    config.coupling_strength = 1.0f;
    core.initialize(config);

    // Initialize uniform phase
    core.initializeUniform(0.0f);

    // Check initial R field (should be 1.0 for uniform)
    float initial_R = core.getAverageR();
    std::cout << "  Initial R = " << initial_R << " (expected ~1.0)" << std::endl;
    assert(std::abs(initial_R - 1.0f) < 1e-5f);

    // Evolve for 100 steps
    auto start = std::chrono::high_resolution_clock::now();

    for (int step = 0; step < 100; ++step) {
        core.evolveKuramotoCPU(config.dt);
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    // Check R field remains high (uniform should stay synchronized)
    float final_R = core.getAverageR();
    std::cout << "  Final R after 100 steps = " << final_R << " (expected ~1.0)" << std::endl;
    std::cout << "  Evolution time: " << duration.count() << " ms" << std::endl;
    assert(final_R > 0.99f);  // Should remain synchronized

    std::cout << "[PASS] Uniform phase evolution verified" << std::endl;
}

void test_random_evolution() {
    std::cout << "\n[TEST] Random Phase Evolution" << std::endl;

    SMFTCore3D core;
    SMFTCore3D::Config config;
    config.Nx = 32;
    config.Ny = 32;
    config.Nz = 32;
    config.dt = 0.01f;
    config.coupling_strength = 2.0f;  // Strong coupling for synchronization
    core.initialize(config);

    // Initialize random phases
    core.initializeRandom(42);

    float initial_R = core.getAverageR();
    std::cout << "  Initial R (random) = " << initial_R << std::endl;

    // Evolve for 500 steps - should see some synchronization
    auto start = std::chrono::high_resolution_clock::now();

    for (int step = 0; step < 500; ++step) {
        core.evolveKuramotoCPU(config.dt);

        if (step % 100 == 99) {
            float R = core.getAverageR();
            std::cout << "  Step " << (step + 1) << ": R = " << R << std::endl;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    float final_R = core.getAverageR();
    std::cout << "  Final R after 500 steps = " << final_R << std::endl;
    std::cout << "  Evolution time: " << duration.count() << " ms" << std::endl;

    // With strong coupling, should see some synchronization increase
    assert(final_R > initial_R);
    std::cout << "  Synchronization increased: "
              << initial_R << " -> " << final_R << " ✓" << std::endl;

    std::cout << "[PASS] Random phase evolution verified" << std::endl;
}

void test_memory_allocation() {
    std::cout << "\n[TEST] Memory Allocation" << std::endl;

    SMFTCore3D core;
    SMFTCore3D::Config config;

    // Test 32³ grid
    config.Nx = 32;
    config.Ny = 32;
    config.Nz = 32;
    core.initialize(config);

    uint32_t expected_points = 32 * 32 * 32;
    assert(core.getTotalPoints() == expected_points);
    std::cout << "  32³ grid: " << expected_points << " points allocated ✓" << std::endl;

    size_t expected_memory = expected_points * sizeof(float);
    std::cout << "  Memory per field: " << expected_memory / 1024 << " KB" << std::endl;
    std::cout << "  Total memory (3 fields): " << (3 * expected_memory) / 1024 << " KB" << std::endl;

    // Verify field sizes
    assert(core.getTheta().size() == expected_points);
    assert(core.getOmega().size() == expected_points);
    assert(core.getRField().size() == expected_points);

    std::cout << "[PASS] Memory allocation verified" << std::endl;
}

void test_boundary_wrapping() {
    std::cout << "\n[TEST] Boundary Wrapping Functions" << std::endl;

    SMFTCore3D core;
    SMFTCore3D::Config config;
    config.Nx = 10;
    config.Ny = 10;
    config.Nz = 10;
    core.initialize(config);

    // Test X wrapping
    assert(core.wrapX(-1) == 9);
    assert(core.wrapX(0) == 0);
    assert(core.wrapX(9) == 9);
    assert(core.wrapX(10) == 0);
    assert(core.wrapX(11) == 1);
    std::cout << "  X wrapping verified ✓" << std::endl;

    // Test Y wrapping
    assert(core.wrapY(-2) == 8);
    assert(core.wrapY(-1) == 9);
    assert(core.wrapY(10) == 0);
    assert(core.wrapY(15) == 5);
    std::cout << "  Y wrapping verified ✓" << std::endl;

    // Test Z wrapping
    assert(core.wrapZ(-3) == 7);
    assert(core.wrapZ(10) == 0);
    assert(core.wrapZ(20) == 0);
    assert(core.wrapZ(25) == 5);
    std::cout << "  Z wrapping verified ✓" << std::endl;

    std::cout << "[PASS] Boundary wrapping verified" << std::endl;
}

int main() {
    std::cout << "=== SMFTCore3D Basic Validation Test ===" << std::endl;
    std::cout << "Week 1 Goals: Core 3D grid functionality" << std::endl;

    try {
        test_memory_allocation();
        test_index_mapping();
        test_boundary_wrapping();
        test_neighbor_access();
        test_uniform_evolution();
        test_random_evolution();

        std::cout << "\n=== ALL TESTS PASSED ===" << std::endl;
        std::cout << "SMFTCore3D Week 1 implementation complete ✓" << std::endl;
        return 0;

    } catch (const std::exception& e) {
        std::cerr << "\n[ERROR] Test failed: " << e.what() << std::endl;
        return 1;
    }
}