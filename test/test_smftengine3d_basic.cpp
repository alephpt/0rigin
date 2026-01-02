/**
 * test_smftengine3d_basic.cpp
 *
 * Week 3-4: Basic SMFTEngine3D infrastructure test
 *
 * Tests:
 * 1. Initialization of 3D grid
 * 2. Phase field setup
 * 3. CPU Kuramoto evolution
 * 4. R-field computation
 */

#include "SMFTEngine3D.h"
#include "Nova.h"
#include <iostream>
#include <cmath>
#include <vector>

bool testInitialization() {
    std::cout << "\n=== Test 1: SMFTEngine3D Initialization ===" << std::endl;

    // Create Nova instance for GPU access
    NovaConfig config{};
    config.name = "SMFTEngine3D Test";
    config.screen = {800, 600};
    config.debug_level = "none";
    config.dimensions = "3D";
    config.camera_type = "orbit";
    config.compute = true;

    Nova nova(config);
    if (!nova.initialized) {
        std::cerr << "Failed to initialize Nova" << std::endl;
        return false;
    }

    // Create SMFTEngine3D
    SMFTEngine3D engine(&nova);

    // Initialize 32³ grid
    const uint32_t N = 32;
    const float Delta = 1.0f;
    engine.initialize(N, N, N, Delta);

    // Verify dimensions
    if (engine.getNx() != N || engine.getNy() != N || engine.getNz() != N) {
        std::cerr << "Grid dimensions mismatch" << std::endl;
        return false;
    }

    if (engine.getTotalPoints() != N * N * N) {
        std::cerr << "Total points mismatch" << std::endl;
        return false;
    }

    std::cout << "✓ Initialization successful: " << N << "³ grid" << std::endl;
    return true;
}

bool testUniformPhaseEvolution() {
    std::cout << "\n=== Test 2: Uniform Phase Evolution ===" << std::endl;

    NovaConfig config{};
    config.name = "SMFTEngine3D Test";
    config.screen = {800, 600};
    config.debug_level = "none";
    config.dimensions = "3D";
    config.camera_type = "orbit";
    config.compute = true;

    Nova nova(config);
    if (!nova.initialized) return false;

    SMFTEngine3D engine(&nova);
    const uint32_t N = 32;
    engine.initialize(N, N, N, 1.0f);

    // Set uniform phase field
    const uint32_t N_total = N * N * N;
    std::vector<float> theta(N_total, 0.5f);  // θ = 0.5 rad everywhere
    std::vector<float> omega(N_total, 0.0f);  // ω = 0 (no intrinsic frequency)

    engine.setInitialPhases(theta);
    engine.setNaturalFrequencies(omega);

    // Get initial phase field
    auto theta_initial = engine.getPhaseField3D();

    // Verify uniform initialization
    for (size_t i = 0; i < N_total; ++i) {
        if (std::abs(theta_initial[i] - 0.5f) > 1e-5f) {
            std::cerr << "Phase initialization failed at index " << i << std::endl;
            return false;
        }
    }

    // Evolve (uniform field should remain unchanged)
    const float dt = 0.01f;
    const float K = 1.0f;
    const float damping = 0.1f;

    engine.stepKuramoto3D(dt, K, damping);

    auto theta_after = engine.getPhaseField3D();

    // Check that phases evolved (due to damping)
    float max_change = 0.0f;
    for (size_t i = 0; i < N_total; ++i) {
        float change = std::abs(theta_after[i] - theta_initial[i]);
        max_change = std::max(max_change, change);
    }

    std::cout << "  Max phase change: " << max_change << " rad" << std::endl;
    std::cout << "✓ Uniform phase evolution successful" << std::endl;
    return true;
}

bool testSyncFieldComputation() {
    std::cout << "\n=== Test 3: Synchronization Field Computation ===" << std::endl;

    NovaConfig config{};
    config.name = "SMFTEngine3D Test";
    config.screen = {800, 600};
    config.debug_level = "none";
    config.dimensions = "3D";
    config.camera_type = "orbit";
    config.compute = true;

    Nova nova(config);
    if (!nova.initialized) return false;

    SMFTEngine3D engine(&nova);
    const uint32_t N = 32;
    engine.initialize(N, N, N, 1.0f);

    const uint32_t N_total = N * N * N;

    // Test 1: Uniform phase → R = 1
    std::vector<float> theta_uniform(N_total, 0.0f);
    std::vector<float> omega(N_total, 0.0f);

    engine.setInitialPhases(theta_uniform);
    engine.setNaturalFrequencies(omega);
    engine.computeSyncField3D();

    auto R_uniform = engine.getSyncField3D();

    // Check that R ≈ 1 for uniform phase
    float R_avg = 0.0f;
    for (size_t i = 0; i < N_total; ++i) {
        R_avg += R_uniform[i];
    }
    R_avg /= N_total;

    std::cout << "  Uniform phase R_avg = " << R_avg << " (expected ≈ 1.0)" << std::endl;

    if (std::abs(R_avg - 1.0f) > 0.1f) {
        std::cerr << "R field computation failed for uniform phase" << std::endl;
        return false;
    }

    // Test 2: Random phases → R < 1
    std::vector<float> theta_random(N_total);
    for (size_t i = 0; i < N_total; ++i) {
        theta_random[i] = -M_PI + 2.0f * M_PI * (float(i) / N_total);
    }

    engine.setInitialPhases(theta_random);
    engine.computeSyncField3D();

    auto R_random = engine.getSyncField3D();

    R_avg = 0.0f;
    for (size_t i = 0; i < N_total; ++i) {
        R_avg += R_random[i];
    }
    R_avg /= N_total;

    std::cout << "  Random phase R_avg = " << R_avg << " (expected < 1.0)" << std::endl;

    if (R_avg >= 1.0f) {
        std::cerr << "R field computation failed for random phases" << std::endl;
        return false;
    }

    std::cout << "✓ Synchronization field computation successful" << std::endl;
    return true;
}

bool testMassField() {
    std::cout << "\n=== Test 4: Mass Field Computation ===" << std::endl;

    NovaConfig config{};
    config.name = "SMFTEngine3D Test";
    config.screen = {800, 600};
    config.debug_level = "none";
    config.dimensions = "3D";
    config.camera_type = "orbit";
    config.compute = true;

    Nova nova(config);
    if (!nova.initialized) return false;

    SMFTEngine3D engine(&nova);
    const uint32_t N = 16;  // Smaller grid for speed
    const float Delta = 2.5f;  // Vacuum potential
    engine.initialize(N, N, N, Delta);

    const uint32_t N_total = N * N * N;

    // Set uniform phase for R = 1
    std::vector<float> theta(N_total, 0.0f);
    std::vector<float> omega(N_total, 0.0f);

    engine.setInitialPhases(theta);
    engine.setNaturalFrequencies(omega);
    engine.computeSyncField3D();

    // Get mass field m = Δ · R
    auto mass = engine.getMassField3D();

    // Check that m ≈ Delta for uniform phase
    float m_avg = 0.0f;
    for (size_t i = 0; i < N_total; ++i) {
        m_avg += mass[i];
    }
    m_avg /= N_total;

    std::cout << "  Average mass m_avg = " << m_avg << " (expected ≈ " << Delta << ")" << std::endl;

    if (std::abs(m_avg - Delta) > 0.1f * Delta) {
        std::cerr << "Mass field computation failed" << std::endl;
        return false;
    }

    std::cout << "✓ Mass field computation successful" << std::endl;
    return true;
}

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "SMFTEngine3D Basic Infrastructure Tests" << std::endl;
    std::cout << "Week 3-4: GPU Integration + 3D Kuramoto" << std::endl;
    std::cout << "========================================" << std::endl;

    bool all_passed = true;

    all_passed &= testInitialization();
    all_passed &= testUniformPhaseEvolution();
    all_passed &= testSyncFieldComputation();
    all_passed &= testMassField();

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
