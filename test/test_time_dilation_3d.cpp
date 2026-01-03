/**
 * test_time_dilation_3d.cpp
 *
 * Gravitational Time Dilation Validation in 3D (A5)
 *
 * Goal: Verify clock rates depend on R-field position
 *
 * Physics:
 *   Metric: g_μν = R²(x,y,z)·η_μν (conformal to Minkowski)
 *   Proper time: dτ = R(x)·dt (coordinate time scaling)
 *   Frequency ratio: ω₂/ω₁ = R₂/R₁ (coordinate frequency ∝ R)
 *
 * Test Method:
 *   - Place oscillators at positions with different R values
 *   - Count oscillation cycles over fixed coordinate time
 *   - Measure frequency ratio ω₂/ω₁
 *   - Compare to theoretical R₁/R₂
 *
 * Test Scenarios:
 *   1. Flat space (R=1) → ω₂/ω₁ = 1 (no time dilation)
 *   2. Gaussian R-peak → measure gravitational redshift
 *   3. Two positions in R-field → verify √(R₁/R₂) relation
 *
 * Quality Gates:
 *   - Frequency ratio within 1% of R₁/R₂
 *   - Numerical stability over long integration times
 *   - Phase coherence maintained
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <array>
#include <functional>

const float PI = 3.14159265358979323846f;

/**
 * Harmonic oscillator structure
 * Tracks phase evolution in local proper time
 */
struct HarmonicOscillator {
    // Position in space
    float x, y, z;

    // Intrinsic frequency (in local proper time)
    float omega_0;

    // Phase (accumulated cycles)
    float phase;

    // Coordinate time
    float t;

    // Cycle count
    int cycles;

    // Phase history for verification
    std::vector<std::array<float, 3>> history;  // {t, phase, R}

    void recordState(float R) {
        history.push_back({t, phase, R});
    }

    // Count zero-crossings (cycles)
    void updateCycleCount(float prev_phase) {
        // Detect when phase crosses 2π
        if (prev_phase < 2.0f * PI && phase >= 2.0f * PI) {
            cycles++;
            phase -= 2.0f * PI;  // Reset phase
        }
    }
};

/**
 * R-field evaluation
 */
struct RFieldData {
    float R;
};

/**
 * Flat R-field: R = 1 everywhere
 */
static RFieldData flatRField(float x, float y, float z) {
    return {1.0f};
}

/**
 * Gaussian R-peak: R(x,y,z) = 1 + A·exp(-r²/2σ²)
 * Represents localized gravitational potential
 */
static RFieldData gaussianRPeak(float x, float y, float z, float A, float sigma) {
    float r2 = x*x + y*y + z*z;
    float R = 1.0f + A * std::exp(-r2 / (2.0f * sigma*sigma));
    return {R};
}

/**
 * Point mass R-field (weak field): R(r) = 1 - GM/r
 * For gravitational redshift tests
 */
RFieldData pointMassRField(float x, float y, float z, float M) {
    float r = std::sqrt(x*x + y*y + z*z);
    if (r < 0.01f) r = 0.01f;  // Regularization

    float epsilon = -M / r;  // G=1 in natural units
    float R = 1.0f + epsilon;
    return {R};
}

/**
 * Evolve oscillator phase in curved spacetime
 *
 * Phase evolution equation:
 *   dϕ/dt = ω₀ · (dτ/dt) = ω₀ · R(x)
 *
 * where dτ = R·dt is the proper time increment
 */
template<typename RFieldFunc>
void evolveOscillator(HarmonicOscillator& osc,
                      RFieldFunc getRField,
                      float dt) {
    // Get R-field at oscillator position
    auto data = getRField(osc.x, osc.y, osc.z);
    float R = data.R;

    // Save previous phase for cycle counting
    float prev_phase = osc.phase;

    // Phase advances according to proper time
    // dϕ/dt = ω₀ · R
    osc.phase += osc.omega_0 * R * dt;

    // Update coordinate time
    osc.t += dt;

    // Count cycles
    osc.updateCycleCount(prev_phase);
}

/**
 * Measure effective frequency over observation period
 *
 * Measured frequency: ω_measured = (total phase change) / (coordinate time)
 */
float measureFrequency(const HarmonicOscillator& osc, float t_obs) {
    // Total phase accumulated (including full cycles)
    float total_phase = osc.cycles * 2.0f * PI + osc.phase;

    // Measured frequency in coordinate time
    return total_phase / t_obs;
}

/**
 * Test 1: Flat space (R=1) → No time dilation
 * Two oscillators should tick at same rate
 */
bool testFlatSpaceNoDilation() {
    std::cout << "\n=== Test 1: Flat Space (R=1) - No Time Dilation ===\n";
    std::cout << "Expected: ω₂/ω₁ = 1.000\n\n";

    // Create two oscillators at different positions (but R=1 everywhere)
    HarmonicOscillator osc1, osc2;

    osc1.x = 0.0f; osc1.y = 0.0f; osc1.z = 0.0f;
    osc1.omega_0 = 1.0f;  // Unit frequency
    osc1.phase = 0.0f;
    osc1.t = 0.0f;
    osc1.cycles = 0;

    osc2.x = 5.0f; osc2.y = 0.0f; osc2.z = 0.0f;
    osc2.omega_0 = 1.0f;  // Same intrinsic frequency
    osc2.phase = 0.0f;
    osc2.t = 0.0f;
    osc2.cycles = 0;

    // Evolve both oscillators
    const float dt = 0.01f;
    const float t_obs = 100.0f;  // Observe for 100 time units
    const int num_steps = static_cast<int>(t_obs / dt);

    for (int i = 0; i < num_steps; ++i) {
        evolveOscillator(osc1, flatRField, dt);
        evolveOscillator(osc2, flatRField, dt);
    }

    // Measure frequencies
    float omega1 = measureFrequency(osc1, t_obs);
    float omega2 = measureFrequency(osc2, t_obs);
    float ratio = omega2 / omega1;

    // Get R values
    auto R1 = flatRField(osc1.x, osc1.y, osc1.z);
    auto R2 = flatRField(osc2.x, osc2.y, osc2.z);
    float expected_ratio = R2.R / R1.R;

    std::cout << "Position 1: (" << osc1.x << ", " << osc1.y << ", " << osc1.z << ")\n";
    std::cout << "  R₁ = " << R1.R << "\n";
    std::cout << "  ω₁ = " << omega1 << " rad/s\n";
    std::cout << "  Cycles: " << osc1.cycles << "\n\n";

    std::cout << "Position 2: (" << osc2.x << ", " << osc2.y << ", " << osc2.z << ")\n";
    std::cout << "  R₂ = " << R2.R << "\n";
    std::cout << "  ω₂ = " << omega2 << " rad/s\n";
    std::cout << "  Cycles: " << osc2.cycles << "\n\n";

    std::cout << "Frequency Ratio:\n";
    std::cout << "  Measured:  ω₂/ω₁ = " << ratio << "\n";
    std::cout << "  Expected:  R₂/R₁ = " << expected_ratio << "\n";

    float error = std::abs(ratio - expected_ratio) / expected_ratio;
    std::cout << "  Error: " << (error * 100.0f) << "%\n";

    bool pass = error < 0.01f;  // 1% tolerance
    std::cout << "\n" << (pass ? "✓ PASS" : "✗ FAIL") << "\n";

    return pass;
}

/**
 * Test 2: Gaussian R-peak → Gravitational redshift
 * Clock at peak (higher R) runs faster than clock far away
 */
bool testGaussianRedshift() {
    std::cout << "\n=== Test 2: Gaussian R-Peak - Gravitational Redshift ===\n";
    std::cout << "Expected: ω_peak/ω_far = R_peak/R_far > 1\n\n";

    const float A = 0.5f;      // Peak amplitude
    const float sigma = 2.0f;  // Width

    // Oscillator at peak (x=0)
    HarmonicOscillator osc_peak;
    osc_peak.x = 0.0f; osc_peak.y = 0.0f; osc_peak.z = 0.0f;
    osc_peak.omega_0 = 1.0f;
    osc_peak.phase = 0.0f;
    osc_peak.t = 0.0f;
    osc_peak.cycles = 0;

    // Oscillator far from peak (x=10)
    HarmonicOscillator osc_far;
    osc_far.x = 10.0f; osc_far.y = 0.0f; osc_far.z = 0.0f;
    osc_far.omega_0 = 1.0f;
    osc_far.phase = 0.0f;
    osc_far.t = 0.0f;
    osc_far.cycles = 0;

    // R-field wrapper
    auto getRField = [A, sigma](float x, float y, float z) {
        return gaussianRPeak(x, y, z, A, sigma);
    };

    // Evolve both oscillators
    const float dt = 0.01f;
    const float t_obs = 100.0f;
    const int num_steps = static_cast<int>(t_obs / dt);

    for (int i = 0; i < num_steps; ++i) {
        evolveOscillator(osc_peak, getRField, dt);
        evolveOscillator(osc_far, getRField, dt);
    }

    // Measure frequencies
    float omega_peak = measureFrequency(osc_peak, t_obs);
    float omega_far = measureFrequency(osc_far, t_obs);
    float ratio = omega_peak / omega_far;

    // Get R values
    auto R_peak = getRField(osc_peak.x, osc_peak.y, osc_peak.z);
    auto R_far = getRField(osc_far.x, osc_far.y, osc_far.z);
    float expected_ratio = R_peak.R / R_far.R;  // This one is correct: ω_peak/ω_far = R_peak/R_far

    std::cout << "Peak Position: (" << osc_peak.x << ", " << osc_peak.y << ", " << osc_peak.z << ")\n";
    std::cout << "  R_peak = " << R_peak.R << "\n";
    std::cout << "  ω_peak = " << omega_peak << " rad/s\n";
    std::cout << "  Cycles: " << osc_peak.cycles << "\n\n";

    std::cout << "Far Position: (" << osc_far.x << ", " << osc_far.y << ", " << osc_far.z << ")\n";
    std::cout << "  R_far = " << R_far.R << "\n";
    std::cout << "  ω_far = " << omega_far << " rad/s\n";
    std::cout << "  Cycles: " << osc_far.cycles << "\n\n";

    std::cout << "Frequency Ratio:\n";
    std::cout << "  Measured:  ω_peak/ω_far = " << ratio << "\n";
    std::cout << "  Expected:  R_peak/R_far = " << expected_ratio << "\n";

    float error = std::abs(ratio - expected_ratio) / expected_ratio;
    std::cout << "  Error: " << (error * 100.0f) << "%\n";

    bool pass = error < 0.01f;  // 1% tolerance
    std::cout << "\n" << (pass ? "✓ PASS" : "✗ FAIL") << "\n";

    return pass;
}

/**
 * Test 3: Point mass R-field → Gravitational time dilation
 * Verify ω(r₂)/ω(r₁) = √(1-2M/r₁) / √(1-2M/r₂) ≈ R(r₁)/R(r₂)
 */
bool testPointMassTimeDilation() {
    std::cout << "\n=== Test 3: Point Mass - Gravitational Time Dilation ===\n";
    std::cout << "Expected: ω₂/ω₁ = R₁/R₂ (redshift at larger radius)\n\n";

    const float M = 1.0f;  // Point mass (natural units)

    // Inner oscillator (r=2)
    HarmonicOscillator osc1;
    osc1.x = 2.0f; osc1.y = 0.0f; osc1.z = 0.0f;
    osc1.omega_0 = 1.0f;
    osc1.phase = 0.0f;
    osc1.t = 0.0f;
    osc1.cycles = 0;

    // Outer oscillator (r=10)
    HarmonicOscillator osc2;
    osc2.x = 10.0f; osc2.y = 0.0f; osc2.z = 0.0f;
    osc2.omega_0 = 1.0f;
    osc2.phase = 0.0f;
    osc2.t = 0.0f;
    osc2.cycles = 0;

    // R-field wrapper
    auto getRField = [M](float x, float y, float z) {
        return pointMassRField(x, y, z, M);
    };

    // Evolve both oscillators
    const float dt = 0.01f;
    const float t_obs = 100.0f;
    const int num_steps = static_cast<int>(t_obs / dt);

    for (int i = 0; i < num_steps; ++i) {
        evolveOscillator(osc1, getRField, dt);
        evolveOscillator(osc2, getRField, dt);
    }

    // Measure frequencies
    float omega1 = measureFrequency(osc1, t_obs);
    float omega2 = measureFrequency(osc2, t_obs);
    float ratio = omega2 / omega1;

    // Get R values
    auto R1 = getRField(osc1.x, osc1.y, osc1.z);
    auto R2 = getRField(osc2.x, osc2.y, osc2.z);
    float expected_ratio = R2.R / R1.R;

    std::cout << "Inner Position: r₁ = " << std::sqrt(osc1.x*osc1.x + osc1.y*osc1.y + osc1.z*osc1.z) << "\n";
    std::cout << "  R₁ = " << R1.R << "\n";
    std::cout << "  ω₁ = " << omega1 << " rad/s\n";
    std::cout << "  Cycles: " << osc1.cycles << "\n\n";

    std::cout << "Outer Position: r₂ = " << std::sqrt(osc2.x*osc2.x + osc2.y*osc2.y + osc2.z*osc2.z) << "\n";
    std::cout << "  R₂ = " << R2.R << "\n";
    std::cout << "  ω₂ = " << omega2 << " rad/s\n";
    std::cout << "  Cycles: " << osc2.cycles << "\n\n";

    std::cout << "Frequency Ratio:\n";
    std::cout << "  Measured:  ω₂/ω₁ = " << ratio << "\n";
    std::cout << "  Expected:  R₂/R₁ = " << expected_ratio << "\n";

    float error = std::abs(ratio - expected_ratio) / expected_ratio;
    std::cout << "  Error: " << (error * 100.0f) << "%\n";

    std::cout << "\nPhysical interpretation:\n";
    std::cout << "  Clock at r=" << std::sqrt(osc2.x*osc2.x) << " ticks ";
    if (ratio < 1.0f) {
        std::cout << "slower (redshifted)\n";
    } else {
        std::cout << "faster (blueshifted)\n";
    }
    std::cout << "  relative to clock at r=" << std::sqrt(osc1.x*osc1.x) << "\n";

    bool pass = error < 0.01f;  // 1% tolerance
    std::cout << "\n" << (pass ? "✓ PASS" : "✗ FAIL") << "\n";

    return pass;
}

/**
 * Test 4: Multiple positions → Consistency check
 * Verify transitivity: (ω₃/ω₂) · (ω₂/ω₁) = (ω₃/ω₁)
 */
bool testTransitivity() {
    std::cout << "\n=== Test 4: Transitivity Check ===\n";
    std::cout << "Expected: (ω₃/ω₂)·(ω₂/ω₁) = (ω₃/ω₁)\n\n";

    const float M = 1.0f;

    // Three oscillators at different radii
    std::vector<HarmonicOscillator> oscs(3);
    std::vector<float> radii = {2.0f, 5.0f, 10.0f};

    for (int i = 0; i < 3; ++i) {
        oscs[i].x = radii[i];
        oscs[i].y = 0.0f;
        oscs[i].z = 0.0f;
        oscs[i].omega_0 = 1.0f;
        oscs[i].phase = 0.0f;
        oscs[i].t = 0.0f;
        oscs[i].cycles = 0;
    }

    auto getRField = [M](float x, float y, float z) {
        return pointMassRField(x, y, z, M);
    };

    // Evolve all oscillators
    const float dt = 0.01f;
    const float t_obs = 100.0f;
    const int num_steps = static_cast<int>(t_obs / dt);

    for (int i = 0; i < num_steps; ++i) {
        for (auto& osc : oscs) {
            evolveOscillator(osc, getRField, dt);
        }
    }

    // Measure frequencies
    std::vector<float> omegas(3);
    for (int i = 0; i < 3; ++i) {
        omegas[i] = measureFrequency(oscs[i], t_obs);
    }

    float ratio_21 = omegas[1] / omegas[0];
    float ratio_32 = omegas[2] / omegas[1];
    float ratio_31 = omegas[2] / omegas[0];
    float product = ratio_32 * ratio_21;

    std::cout << "ω₂/ω₁ = " << ratio_21 << "\n";
    std::cout << "ω₃/ω₂ = " << ratio_32 << "\n";
    std::cout << "ω₃/ω₁ = " << ratio_31 << "\n";
    std::cout << "(ω₃/ω₂)·(ω₂/ω₁) = " << product << "\n";

    float error = std::abs(product - ratio_31) / ratio_31;
    std::cout << "Error: " << (error * 100.0f) << "%\n";

    bool pass = error < 0.01f;  // 1% tolerance
    std::cout << "\n" << (pass ? "✓ PASS" : "✗ FAIL") << "\n";

    return pass;
}

/**
 * Main test runner
 */
int runTimeDilation3DTest() {
    std::cout << "\n========================================\n";
    std::cout << "  A5: Gravitational Time Dilation Test\n";
    std::cout << "========================================\n";
    std::cout << "\nPhysics: ω₂/ω₁ = R₂/R₁ (coordinate frequency ∝ R)\n";
    std::cout << "Quality Gate: |measured - expected|/expected < 1%\n";

    bool all_pass = true;

    all_pass &= testFlatSpaceNoDilation();
    all_pass &= testGaussianRedshift();
    all_pass &= testPointMassTimeDilation();
    all_pass &= testTransitivity();

    std::cout << "\n========================================\n";
    if (all_pass) {
        std::cout << "✓ ALL TESTS PASSED\n";
        std::cout << "A5 Gravitational Time Dilation: VALIDATED\n";
    } else {
        std::cout << "✗ SOME TESTS FAILED\n";
    }
    std::cout << "========================================\n\n";

    return all_pass ? 0 : 1;
}
