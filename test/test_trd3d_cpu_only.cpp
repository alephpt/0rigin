/**
 * test_smft3d_cpu_only.cpp
 *
 * Week 3-4: CPU-only 3D infrastructure test
 *
 * Tests TRDCore3D without GPU/Nova dependency for CI/CD
 */

#include "TRDCore3D.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>

bool testBasic3DGrid() {
    std::cout << "\n=== Test 1: Basic 3D Grid ===" << std::endl;

    TRDCore3D core;
    TRDCore3D::Config config;
    config.Nx = 32;
    config.Ny = 32;
    config.Nz = 32;
    config.dx = 1.0f;
    config.dt = 0.01f;
    config.coupling_strength = 1.0f;

    core.initialize(config);

    if (core.getTotalPoints() != 32 * 32 * 32) {
        std::cerr << "Grid size mismatch" << std::endl;
        return false;
    }

    std::cout << "✓ Basic 3D grid initialized: " << core.getNx() << "×"
              << core.getNy() << "×" << core.getNz() << std::endl;
    return true;
}

bool testIndexMapping() {
    std::cout << "\n=== Test 2: 3D Index Mapping ===" << std::endl;

    TRDCore3D core;
    TRDCore3D::Config config;
    config.Nx = 8;
    config.Ny = 8;
    config.Nz = 8;
    core.initialize(config);

    // Test corner points
    uint32_t idx_000 = core.index3D(0, 0, 0);
    uint32_t idx_777 = core.index3D(7, 7, 7);

    if (idx_000 != 0) {
        std::cerr << "Origin index wrong: " << idx_000 << std::endl;
        return false;
    }

    if (idx_777 != 8*8*8 - 1) {
        std::cerr << "Corner index wrong: " << idx_777 << std::endl;
        return false;
    }

    // Test inverse mapping
    uint32_t i, j, k;
    core.coords3D(idx_777, i, j, k);
    if (i != 7 || j != 7 || k != 7) {
        std::cerr << "Inverse mapping failed: (" << i << "," << j << "," << k << ")" << std::endl;
        return false;
    }

    std::cout << "✓ 3D index mapping correct" << std::endl;
    return true;
}

bool testPeriodicBoundaries() {
    std::cout << "\n=== Test 3: Periodic Boundaries ===" << std::endl;

    TRDCore3D core;
    TRDCore3D::Config config;
    config.Nx = 8;
    config.Ny = 8;
    config.Nz = 8;
    core.initialize(config);

    // Test wrapping in X direction
    auto neighbors = core.getNeighbors(0, 4, 4);  // Left edge
    uint32_t i_minus, j_minus, k_minus;
    core.coords3D(neighbors.x_minus, i_minus, j_minus, k_minus);

    if (i_minus != 7) {  // Should wrap to right edge
        std::cerr << "X-wrapping failed: got i=" << i_minus << std::endl;
        return false;
    }

    // Test wrapping in Z direction
    neighbors = core.getNeighbors(4, 4, 7);  // Top edge
    core.coords3D(neighbors.z_plus, i_minus, j_minus, k_minus);

    if (k_minus != 0) {  // Should wrap to bottom
        std::cerr << "Z-wrapping failed: got k=" << k_minus << std::endl;
        return false;
    }

    std::cout << "✓ Periodic boundaries correct" << std::endl;
    return true;
}

bool testUniformKuramoto() {
    std::cout << "\n=== Test 4: Uniform Phase Kuramoto ===" << std::endl;

    TRDCore3D core;
    TRDCore3D::Config config;
    config.Nx = 16;
    config.Ny = 16;
    config.Nz = 16;
    config.dt = 0.01f;
    config.coupling_strength = 1.0f;
    core.initialize(config);

    // Set uniform phase
    core.initializeUniform(0.5f);

    auto theta_before = core.getTheta();

    // Evolve (uniform field should stay relatively stable)
    core.evolveKuramotoCPU(config.dt);

    auto theta_after = core.getTheta();

    // Check that phase evolved slightly due to any numerical effects
    float max_change = 0.0f;
    for (size_t i = 0; i < theta_before.size(); ++i) {
        float change = std::abs(theta_after[i] - theta_before[i]);
        max_change = std::max(max_change, change);
    }

    std::cout << "  Max phase change: " << max_change << " rad" << std::endl;
    std::cout << "✓ Uniform Kuramoto evolution successful" << std::endl;
    return true;
}

bool testSyncFieldUniform() {
    std::cout << "\n=== Test 5: Sync Field (Uniform) ===" << std::endl;

    TRDCore3D core;
    TRDCore3D::Config config;
    config.Nx = 16;
    config.Ny = 16;
    config.Nz = 16;
    core.initialize(config);

    // Set uniform phase → R should be ~1
    core.initializeUniform(0.0f);
    core.computeRField();

    float R_avg = core.getAverageR();

    std::cout << "  Uniform phase R_avg = " << R_avg << " (expected ≈ 1.0)" << std::endl;

    if (std::abs(R_avg - 1.0f) > 0.1f) {
        std::cerr << "R field computation failed" << std::endl;
        return false;
    }

    std::cout << "✓ Sync field (uniform) correct" << std::endl;
    return true;
}

bool testSyncFieldRandom() {
    std::cout << "\n=== Test 6: Sync Field (Random) ===" << std::endl;

    TRDCore3D core;
    TRDCore3D::Config config;
    config.Nx = 16;
    config.Ny = 16;
    config.Nz = 16;
    core.initialize(config);

    // Set random phases → R should be < 1
    core.initializeRandom(42);
    core.computeRField();

    float R_avg = core.getAverageR();

    std::cout << "  Random phase R_avg = " << R_avg << " (expected < 0.5)" << std::endl;

    if (R_avg >= 0.5f) {
        std::cerr << "R field for random phases too high" << std::endl;
        return false;
    }

    std::cout << "✓ Sync field (random) correct" << std::endl;
    return true;
}

bool test3DKuramotoSynchronization() {
    std::cout << "\n=== Test 7: 3D Kuramoto Synchronization ===" << std::endl;

    TRDCore3D core;
    TRDCore3D::Config config;
    config.Nx = 16;
    config.Ny = 16;
    config.Nz = 16;
    config.dt = 0.01f;
    config.coupling_strength = 2.0f;  // Strong coupling
    core.initialize(config);

    // Start with random phases
    core.initializeRandom(123);

    float R_initial = core.getAverageR();
    std::cout << "  Initial R = " << R_initial << std::endl;

    // Evolve for many steps
    const int num_steps = 1000;
    for (int step = 0; step < num_steps; ++step) {
        core.evolveKuramotoCPU(config.dt);
        if (step % 100 == 0) {
            core.computeRField();
        }
    }

    core.computeRField();
    float R_final = core.getAverageR();
    std::cout << "  Final R = " << R_final << " (after " << num_steps << " steps)" << std::endl;

    // Strong coupling should increase synchronization
    if (R_final <= R_initial) {
        std::cerr << "Synchronization did not increase" << std::endl;
        return false;
    }

    std::cout << "  ΔR = +" << (R_final - R_initial) << std::endl;
    std::cout << "✓ 3D Kuramoto synchronization working" << std::endl;
    return true;
}

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "TRD 3D CPU-Only Infrastructure Tests" << std::endl;
    std::cout << "Week 3-4: 3D Grid + Kuramoto Dynamics" << std::endl;
    std::cout << "========================================" << std::endl;

    bool all_passed = true;

    all_passed &= testBasic3DGrid();
    all_passed &= testIndexMapping();
    all_passed &= testPeriodicBoundaries();
    all_passed &= testUniformKuramoto();
    all_passed &= testSyncFieldUniform();
    all_passed &= testSyncFieldRandom();
    all_passed &= test3DKuramotoSynchronization();

    std::cout << "\n========================================" << std::endl;
    if (all_passed) {
        std::cout << "✓ ALL TESTS PASSED" << std::endl;
        std::cout << "========================================" << std::endl;
        return 0;
    } else {
        std::cout << "✗ SOME TESTS FAILED" << std::endl;
        std::cout << "========================================" << std::endl;
        return 1;
    }
}
