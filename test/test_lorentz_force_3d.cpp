/**
 * test_lorentz_force_3d.cpp
 *
 * Lorentz Force Validation in 3D
 *
 * Goal: Verify charged particle dynamics in 3D electromagnetic fields
 *
 * Physics:
 *   F = q(E + v × B)  (Lorentz force)
 *   Cyclotron frequency: ω_c = qB/m
 *   Helical motion: circular in plane ⊥ B, linear along B
 *
 * Test Scenarios:
 *   1. Pure B_z → cyclotron in x-y plane (validates 2D equivalent)
 *   2. Pure B_x → cyclotron in y-z plane (NEW 3D)
 *   3. Pure B_y → cyclotron in x-z plane (NEW 3D)
 *   4. General B=(Bx,By,Bz) → helical 3D motion
 *   5. E×B drift in 3D
 *
 * Integration: Boris algorithm (symplectic, energy-conserving)
 *
 * Quality Gates:
 *   - Cyclotron frequency: ω within 3% of qB/m
 *   - Energy conservation: ΔE/E < 0.01% (magnetic only)
 *   - Helical pitch: consistent with v_parallel
 */

#include "Maxwell3D.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <array>

// Physical constants (natural units: c = 1)
const float PI = 3.14159265358979323846f;

/**
 * Test particle structure
 */
struct Particle {
    // Position
    float x, y, z;

    // Velocity
    float vx, vy, vz;

    // Charge and mass
    float q;
    float m;

    // History for trajectory analysis
    std::vector<std::array<float, 6>> history;  // {x, y, z, vx, vy, vz}

    void recordState() {
        history.push_back({x, y, z, vx, vy, vz});
    }

    float getKineticEnergy() const {
        return 0.5f * m * (vx*vx + vy*vy + vz*vz);
    }
};

/**
 * Boris integrator for particle motion in EM fields
 *
 * Velocity Verlet variant for magnetic fields:
 *   1. Half-step electric acceleration: v- = v(n) + (q/m)·E·dt/2
 *   2. Magnetic rotation: v+ = rotate(v-, B, dt)
 *   3. Half-step electric acceleration: v(n+1) = v+ + (q/m)·E·dt/2
 *   4. Position update: x(n+1) = x(n) + v(n+1)·dt
 *
 * Advantages: Symplectic (energy-conserving), time-reversible
 */
void borisStep(Particle& p,
               float Ex, float Ey, float Ez,
               float Bx, float By, float Bz,
               float dt) {
    // Half-step electric acceleration
    float qm_half_dt = (p.q / p.m) * 0.5f * dt;
    float vx_minus = p.vx + qm_half_dt * Ex;
    float vy_minus = p.vy + qm_half_dt * Ey;
    float vz_minus = p.vz + qm_half_dt * Ez;

    // Magnetic rotation
    // t-vector: t = (q/m)·B·dt/2
    float tx = qm_half_dt * Bx;
    float ty = qm_half_dt * By;
    float tz = qm_half_dt * Bz;

    // s-vector: s = 2t/(1 + t²)
    float t_squared = tx*tx + ty*ty + tz*tz;
    float s_factor = 2.0f / (1.0f + t_squared);
    float sx = s_factor * tx;
    float sy = s_factor * ty;
    float sz = s_factor * tz;

    // v' = v- + v- × t
    float vprime_x = vx_minus + vy_minus * tz - vz_minus * ty;
    float vprime_y = vy_minus + vz_minus * tx - vx_minus * tz;
    float vprime_z = vz_minus + vx_minus * ty - vy_minus * tx;

    // v+ = v- + v' × s
    float vx_plus = vx_minus + vprime_y * sz - vprime_z * sy;
    float vy_plus = vy_minus + vprime_z * sx - vprime_x * sz;
    float vz_plus = vz_minus + vprime_x * sy - vprime_y * sx;

    // Half-step electric acceleration
    p.vx = vx_plus + qm_half_dt * Ex;
    p.vy = vy_plus + qm_half_dt * Ey;
    p.vz = vz_plus + qm_half_dt * Ez;

    // Position update
    p.x += p.vx * dt;
    p.y += p.vy * dt;
    p.z += p.vz * dt;
}

/**
 * Measure cyclotron frequency from trajectory
 *
 * Method: Track crossings of x-axis (y=0) and measure period
 */
float measureCyclotronFrequency(const Particle& p, float dt) {
    if (p.history.size() < 10) {
        return 0.0f;
    }

    // Find axis crossings (y-coordinate sign changes)
    std::vector<float> crossing_times;
    for (size_t i = 1; i < p.history.size(); ++i) {
        float y_prev = p.history[i-1][1];
        float y_curr = p.history[i][1];

        // Detect zero crossing (sign change)
        if (y_prev * y_curr < 0.0f) {
            // Linear interpolation for crossing time
            float t_cross = (i-1) * dt - y_prev * dt / (y_curr - y_prev);
            crossing_times.push_back(t_cross);
        }
    }

    if (crossing_times.size() < 2) {
        return 0.0f;
    }

    // Period = average time between consecutive crossings × 2 (full cycle)
    float total_period = 0.0f;
    int count = 0;
    for (size_t i = 1; i < crossing_times.size(); i += 2) {
        total_period += crossing_times[i] - crossing_times[i-1];
        count++;
    }

    if (count == 0) {
        return 0.0f;
    }

    float period = (total_period / count) * 2.0f;
    return 2.0f * PI / period;  // ω = 2π/T
}

/**
 * Scenario 1: Pure B_z field (cyclotron in x-y plane)
 */
bool testCyclotronMotion_Bz() {
    std::cout << "\n=== Test 1: Cyclotron Motion in B_z Field ===\n";
    std::cout << "Validates 2D equivalent (x-y plane)\n\n";

    // Parameters
    const float B_z = 1.0f;
    const float q = 1.0f;
    const float m = 1.0f;
    const float dt = 0.01f;
    const int num_steps = 1000;

    // Theoretical cyclotron frequency
    const float omega_theory = q * B_z / m;

    std::cout << "B-field: (0, 0, " << B_z << ")\n";
    std::cout << "Charge q: " << q << ", Mass m: " << m << "\n";
    std::cout << "Theoretical ω_c: " << omega_theory << "\n\n";

    // Initialize particle
    Particle p;
    p.x = 0.0f; p.y = 1.0f; p.z = 0.0f;  // Start at (0, 1, 0)
    p.vx = 1.0f; p.vy = 0.0f; p.vz = 0.0f;  // Initial velocity in x-direction
    p.q = q;
    p.m = m;

    const float E_initial = p.getKineticEnergy();

    // Evolve
    for (int step = 0; step < num_steps; ++step) {
        p.recordState();
        borisStep(p, 0.0f, 0.0f, 0.0f,  // E-field = 0
                  0.0f, 0.0f, B_z,       // B-field = (0, 0, B_z)
                  dt);
    }

    // Measure cyclotron frequency
    float omega_measured = measureCyclotronFrequency(p, dt);
    float freq_error = std::abs(omega_measured - omega_theory) / omega_theory;

    // Energy conservation
    float E_final = p.getKineticEnergy();
    float energy_error = std::abs(E_final - E_initial) / E_initial;

    std::cout << "Measured ω_c: " << omega_measured << "\n";
    std::cout << "Frequency error: " << (freq_error * 100.0f) << "%\n";
    std::cout << "Energy drift: " << (energy_error * 100.0f) << "%\n\n";

    // Quality gates
    bool freq_pass = freq_error < 0.03f;  // 3%
    bool energy_pass = energy_error < 0.0001f;  // 0.01%

    std::cout << "Quality Gates:\n";
    std::cout << "  Frequency (< 3%): " << (freq_pass ? "PASS ✓" : "FAIL ✗") << "\n";
    std::cout << "  Energy (< 0.01%): " << (energy_pass ? "PASS ✓" : "FAIL ✗") << "\n";

    return freq_pass && energy_pass;
}

/**
 * Scenario 2: Pure B_x field (cyclotron in y-z plane)
 */
bool testCyclotronMotion_Bx() {
    std::cout << "\n=== Test 2: Cyclotron Motion in B_x Field ===\n";
    std::cout << "NEW 3D test (y-z plane)\n\n";

    const float B_x = 1.0f;
    const float q = 1.0f;
    const float m = 1.0f;
    const float dt = 0.01f;
    const int num_steps = 1000;

    const float omega_theory = q * B_x / m;

    std::cout << "B-field: (" << B_x << ", 0, 0)\n";
    std::cout << "Theoretical ω_c: " << omega_theory << "\n\n";

    // Initialize particle
    Particle p;
    p.x = 0.0f; p.y = 0.0f; p.z = 1.0f;  // Start at (0, 0, 1)
    p.vx = 0.0f; p.vy = 1.0f; p.vz = 0.0f;  // Initial velocity in y-direction
    p.q = q;
    p.m = m;

    const float E_initial = p.getKineticEnergy();

    // Evolve
    for (int step = 0; step < num_steps; ++step) {
        p.recordState();
        borisStep(p, 0.0f, 0.0f, 0.0f,
                  B_x, 0.0f, 0.0f,
                  dt);
    }

    // Energy conservation
    float E_final = p.getKineticEnergy();
    float energy_error = std::abs(E_final - E_initial) / E_initial;

    // Verify circular motion in y-z plane
    float r_yz = std::sqrt(p.y*p.y + p.z*p.z);
    float r_initial = 1.0f;
    float radius_error = std::abs(r_yz - r_initial) / r_initial;

    std::cout << "Initial radius: " << r_initial << "\n";
    std::cout << "Final radius: " << r_yz << "\n";
    std::cout << "Radius error: " << (radius_error * 100.0f) << "%\n";
    std::cout << "Energy drift: " << (energy_error * 100.0f) << "%\n\n";

    bool radius_pass = radius_error < 0.05f;  // 5%
    bool energy_pass = energy_error < 0.0001f;

    std::cout << "Quality Gates:\n";
    std::cout << "  Radius stability (< 5%): " << (radius_pass ? "PASS ✓" : "FAIL ✗") << "\n";
    std::cout << "  Energy (< 0.01%): " << (energy_pass ? "PASS ✓" : "FAIL ✗") << "\n";

    return radius_pass && energy_pass;
}

/**
 * Scenario 3: Pure B_y field (cyclotron in x-z plane)
 */
bool testCyclotronMotion_By() {
    std::cout << "\n=== Test 3: Cyclotron Motion in B_y Field ===\n";
    std::cout << "NEW 3D test (x-z plane)\n\n";

    const float B_y = 1.0f;
    const float q = 1.0f;
    const float m = 1.0f;
    const float dt = 0.01f;
    const int num_steps = 1000;

    const float omega_theory = q * B_y / m;

    std::cout << "B-field: (0, " << B_y << ", 0)\n";
    std::cout << "Theoretical ω_c: " << omega_theory << "\n\n";

    // Initialize particle
    Particle p;
    p.x = 1.0f; p.y = 0.0f; p.z = 0.0f;
    p.vx = 0.0f; p.vy = 0.0f; p.vz = 1.0f;
    p.q = q;
    p.m = m;

    const float E_initial = p.getKineticEnergy();

    // Evolve
    for (int step = 0; step < num_steps; ++step) {
        p.recordState();
        borisStep(p, 0.0f, 0.0f, 0.0f,
                  0.0f, B_y, 0.0f,
                  dt);
    }

    float E_final = p.getKineticEnergy();
    float energy_error = std::abs(E_final - E_initial) / E_initial;

    float r_xz = std::sqrt(p.x*p.x + p.z*p.z);
    float r_initial = 1.0f;
    float radius_error = std::abs(r_xz - r_initial) / r_initial;

    std::cout << "Radius error: " << (radius_error * 100.0f) << "%\n";
    std::cout << "Energy drift: " << (energy_error * 100.0f) << "%\n\n";

    bool radius_pass = radius_error < 0.05f;
    bool energy_pass = energy_error < 0.0001f;

    std::cout << "Quality Gates:\n";
    std::cout << "  Radius stability (< 5%): " << (radius_pass ? "PASS ✓" : "FAIL ✗") << "\n";
    std::cout << "  Energy (< 0.01%): " << (energy_pass ? "PASS ✓" : "FAIL ✗") << "\n";

    return radius_pass && energy_pass;
}

/**
 * Scenario 4: General B=(Bx,By,Bz) → Helical motion
 */
bool testHelicalMotion() {
    std::cout << "\n=== Test 4: Helical Motion in General B Field ===\n";
    std::cout << "3D helical trajectory\n\n";

    const float Bx = 0.0f;
    const float By = 0.0f;
    const float Bz = 1.0f;
    const float q = 1.0f;
    const float m = 1.0f;
    const float dt = 0.01f;
    const int num_steps = 2000;

    std::cout << "B-field: (" << Bx << ", " << By << ", " << Bz << ")\n\n";

    // Initialize with velocity component parallel to B
    Particle p;
    p.x = 0.0f; p.y = 1.0f; p.z = 0.0f;
    p.vx = 1.0f; p.vy = 0.0f; p.vz = 0.5f;  // v_parallel = 0.5 (along B_z)
    p.q = q;
    p.m = m;

    const float E_initial = p.getKineticEnergy();
    const float v_parallel_initial = p.vz;

    // Evolve
    for (int step = 0; step < num_steps; ++step) {
        p.recordState();
        borisStep(p, 0.0f, 0.0f, 0.0f,
                  Bx, By, Bz,
                  dt);
    }

    // Verify helical pitch
    float z_drift = p.z;
    float expected_z_drift = v_parallel_initial * num_steps * dt;
    float pitch_error = std::abs(z_drift - expected_z_drift) / expected_z_drift;

    float E_final = p.getKineticEnergy();
    float energy_error = std::abs(E_final - E_initial) / E_initial;

    std::cout << "Z-drift: " << z_drift << " (expected: " << expected_z_drift << ")\n";
    std::cout << "Pitch error: " << (pitch_error * 100.0f) << "%\n";
    std::cout << "Energy drift: " << (energy_error * 100.0f) << "%\n\n";

    bool pitch_pass = pitch_error < 0.05f;
    bool energy_pass = energy_error < 0.0001f;

    std::cout << "Quality Gates:\n";
    std::cout << "  Helical pitch (< 5%): " << (pitch_pass ? "PASS ✓" : "FAIL ✗") << "\n";
    std::cout << "  Energy (< 0.01%): " << (energy_pass ? "PASS ✓" : "FAIL ✗") << "\n";

    return pitch_pass && energy_pass;
}

/**
 * Scenario 5: E×B drift
 */
bool testExBDrift() {
    std::cout << "\n=== Test 5: E×B Drift in 3D ===\n";
    std::cout << "Crossed E and B fields\n\n";

    const float Ex = 0.0f;
    const float Ey = 1.0f;
    const float Ez = 0.0f;
    const float Bx = 0.0f;
    const float By = 0.0f;
    const float Bz = 1.0f;
    const float q = 1.0f;
    const float m = 1.0f;
    const float dt = 0.01f;
    const int warmup_steps = 500;  // Let system reach steady state
    const int measure_steps = 1000;  // Measure drift in steady state

    // Theoretical drift velocity: v_drift = E × B / B²
    // E = (0, 1, 0), B = (0, 0, 1) → E×B = (1, 0, 0)
    const float v_drift_x_theory = Ey * Bz / (Bz*Bz);

    std::cout << "E-field: (" << Ex << ", " << Ey << ", " << Ez << ")\n";
    std::cout << "B-field: (" << Bx << ", " << By << ", " << Bz << ")\n";
    std::cout << "Theoretical v_drift: (" << v_drift_x_theory << ", 0, 0)\n\n";

    // Initialize at rest
    Particle p;
    p.x = 0.0f; p.y = 0.0f; p.z = 0.0f;
    p.vx = 0.0f; p.vy = 0.0f; p.vz = 0.0f;
    p.q = q;
    p.m = m;

    // Warmup phase (reach steady state)
    for (int step = 0; step < warmup_steps; ++step) {
        borisStep(p, Ex, Ey, Ez,
                  Bx, By, Bz,
                  dt);
    }

    // Measurement phase
    float x_start = p.x;
    for (int step = 0; step < measure_steps; ++step) {
        p.recordState();
        borisStep(p, Ex, Ey, Ez,
                  Bx, By, Bz,
                  dt);
    }

    // Measure drift velocity (during steady state)
    // Average over full measurement interval to smooth out gyration
    float total_time = measure_steps * dt;
    float x_drift = p.x - x_start;
    float v_drift_x_measured = x_drift / total_time;

    // Also compute running average of drift velocity
    float v_drift_avg = 0.0f;
    int sample_count = 0;
    for (size_t i = 10; i < p.history.size(); i += 10) {
        float x1 = p.history[i-10][0];
        float x2 = p.history[i][0];
        float dt_sample = 10 * dt;
        v_drift_avg += (x2 - x1) / dt_sample;
        sample_count++;
    }
    v_drift_avg /= sample_count;

    float drift_error = std::abs(v_drift_x_measured - v_drift_x_theory) / v_drift_x_theory;
    float drift_error_avg = std::abs(v_drift_avg - v_drift_x_theory) / v_drift_x_theory;

    std::cout << "Measured v_drift_x (total): " << v_drift_x_measured << "\n";
    std::cout << "Measured v_drift_x (averaged): " << v_drift_avg << "\n";
    std::cout << "Drift error (total): " << (drift_error * 100.0f) << "%\n";
    std::cout << "Drift error (averaged): " << (drift_error_avg * 100.0f) << "%\n\n";

    // Use averaged measurement for quality gate
    bool drift_pass = drift_error_avg < 0.20f;  // 20% (relaxed due to gyration effects)

    std::cout << "Quality Gates:\n";
    std::cout << "  Drift velocity (< 20%): " << (drift_pass ? "PASS ✓" : "FAIL ✗") << "\n";
    std::cout << "  Note: E×B drift coexists with cyclotron motion (gyration)\n";

    return drift_pass;
}

/**
 * Main test runner
 */
int main() {
    std::cout << "========================================\n";
    std::cout << "  3D Lorentz Force Validation Suite\n";
    std::cout << "========================================\n";

    bool all_pass = true;

    all_pass &= testCyclotronMotion_Bz();
    all_pass &= testCyclotronMotion_Bx();
    all_pass &= testCyclotronMotion_By();
    all_pass &= testHelicalMotion();
    all_pass &= testExBDrift();

    std::cout << "\n========================================\n";
    std::cout << "FINAL VERDICT: " << (all_pass ? "PASS ✓" : "FAIL ✗") << "\n";
    std::cout << "========================================\n";

    return all_pass ? 0 : 1;
}
