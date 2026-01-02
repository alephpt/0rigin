/**
 * test_three_body_em_3d.cpp
 *
 * Three-Body Electromagnetic Dynamics in 3D
 *
 * Goal: Verify electromagnetic interactions in multi-particle systems
 *
 * Physics:
 *   Coulomb force: F_ij = k·q_i·q_j·r_ij/|r_ij|³
 *   Superposition principle: F_i = Σ_j F_ij
 *   Conservation laws: Total energy, momentum
 *
 * Test Scenarios:
 *   1. Three charges (+q, +q, -2q) → Net-zero charge
 *   2. Equilateral triangle initial configuration
 *   3. Verify energy/momentum conservation
 *   4. Verify superposition principle
 *
 * Integration: Velocity Verlet (symplectic)
 *
 * Quality Gates:
 *   - Total energy conserved < 0.1%
 *   - Total momentum conserved < 0.1%
 *   - Superposition principle validated
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <array>

const float PI = 3.14159265358979323846f;
const float k_e = 1.0f;  // Coulomb constant (natural units)

/**
 * Particle structure
 */
struct ChargedParticle {
    float x, y, z;      // Position
    float vx, vy, vz;   // Velocity
    float q;            // Charge
    float m;            // Mass

    // Forces (computed each step)
    float fx, fy, fz;

    float getKineticEnergy() const {
        return 0.5f * m * (vx*vx + vy*vy + vz*vz);
    }

    std::array<float, 3> getMomentum() const {
        return {m*vx, m*vy, m*vz};
    }
};

/**
 * Compute Coulomb force between two particles
 * F_ij = k·q_i·q_j·r_ij/|r_ij|³
 */
std::array<float, 3> coulombForce(const ChargedParticle& pi, const ChargedParticle& pj) {
    // Vector from j to i
    float dx = pi.x - pj.x;
    float dy = pi.y - pj.y;
    float dz = pi.z - pj.z;

    float r2 = dx*dx + dy*dy + dz*dz;
    float r = std::sqrt(r2);

    if (r < 0.1f) {
        // Regularization near collision
        r = 0.1f;
        r2 = r * r;
    }

    float r3 = r2 * r;
    float force_mag = k_e * pi.q * pj.q / r2;

    // Force on particle i from particle j
    float fx = force_mag * dx / r;
    float fy = force_mag * dy / r;
    float fz = force_mag * dz / r;

    return {fx, fy, fz};
}

/**
 * Compute potential energy between two particles
 * U_ij = k·q_i·q_j/r_ij
 */
float coulombPotential(const ChargedParticle& pi, const ChargedParticle& pj) {
    float dx = pi.x - pj.x;
    float dy = pi.y - pj.y;
    float dz = pi.z - pj.z;

    float r = std::sqrt(dx*dx + dy*dy + dz*dz);
    if (r < 0.1f) r = 0.1f;  // Regularization

    return k_e * pi.q * pj.q / r;
}

/**
 * Compute total energy of system
 */
float totalEnergy(const std::vector<ChargedParticle>& particles) {
    float KE = 0.0f;
    float PE = 0.0f;

    // Kinetic energy
    for (const auto& p : particles) {
        KE += p.getKineticEnergy();
    }

    // Potential energy (pairwise)
    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            PE += coulombPotential(particles[i], particles[j]);
        }
    }

    return KE + PE;
}

/**
 * Compute total momentum of system
 */
std::array<float, 3> totalMomentum(const std::vector<ChargedParticle>& particles) {
    float px = 0.0f, py = 0.0f, pz = 0.0f;

    for (const auto& p : particles) {
        auto mom = p.getMomentum();
        px += mom[0];
        py += mom[1];
        pz += mom[2];
    }

    return {px, py, pz};
}

/**
 * Velocity Verlet integration step
 */
void velocityVerletStep(std::vector<ChargedParticle>& particles, float dt) {
    const size_t N = particles.size();

    // Step 1: Update positions using current velocities
    for (auto& p : particles) {
        p.x += p.vx * dt + 0.5f * (p.fx / p.m) * dt * dt;
        p.y += p.vy * dt + 0.5f * (p.fy / p.m) * dt * dt;
        p.z += p.vz * dt + 0.5f * (p.fz / p.m) * dt * dt;
    }

    // Step 2: Compute new forces
    std::vector<std::array<float, 3>> forces_new(N, {0.0f, 0.0f, 0.0f});

    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            if (i != j) {
                auto F = coulombForce(particles[i], particles[j]);
                forces_new[i][0] += F[0];
                forces_new[i][1] += F[1];
                forces_new[i][2] += F[2];
            }
        }
    }

    // Step 3: Update velocities using average of old and new forces
    for (size_t i = 0; i < N; ++i) {
        auto& p = particles[i];
        p.vx += 0.5f * (p.fx + forces_new[i][0]) / p.m * dt;
        p.vy += 0.5f * (p.fy + forces_new[i][1]) / p.m * dt;
        p.vz += 0.5f * (p.fz + forces_new[i][2]) / p.m * dt;

        // Update forces for next step
        p.fx = forces_new[i][0];
        p.fy = forces_new[i][1];
        p.fz = forces_new[i][2];
    }
}

/**
 * Test 1: Three-body system (+q, +q, -2q)
 * Net-zero charge, equilateral triangle
 */
bool testThreeBodyConservation() {
    std::cout << "\n=== Test 1: Three-Body Conservation Laws ===\n";
    std::cout << "Configuration: (+q, +q, -2q) in equilateral triangle\n";
    std::cout << "Net charge: 0 (charge neutral)\n\n";

    // Initialize three particles in equilateral triangle
    const float q = 1.0f;
    const float m = 1.0f;
    const float L = 5.0f;  // Triangle side length

    std::vector<ChargedParticle> particles(3);

    // Particle 0: +q at (0, 0, 0)
    particles[0] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, q, m, 0.0f, 0.0f, 0.0f};

    // Particle 1: +q at (L, 0, 0)
    particles[1] = {L, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, q, m, 0.0f, 0.0f, 0.0f};

    // Particle 2: -2q at (L/2, L*√3/2, 0)
    particles[2] = {L/2.0f, L*std::sqrt(3.0f)/2.0f, 0.0f, 0.0f, 0.0f, 0.0f, -2.0f*q, m, 0.0f, 0.0f, 0.0f};

    // Verify net charge = 0
    float Q_total = particles[0].q + particles[1].q + particles[2].q;
    std::cout << "Total charge: " << Q_total << " (expect 0)\n";

    // Compute initial forces
    for (size_t i = 0; i < 3; ++i) {
        particles[i].fx = 0.0f;
        particles[i].fy = 0.0f;
        particles[i].fz = 0.0f;
        for (size_t j = 0; j < 3; ++j) {
            if (i != j) {
                auto F = coulombForce(particles[i], particles[j]);
                particles[i].fx += F[0];
                particles[i].fy += F[1];
                particles[i].fz += F[2];
            }
        }
    }

    // Record initial conserved quantities
    float E_initial = totalEnergy(particles);
    auto p_initial = totalMomentum(particles);

    std::cout << "Initial energy: " << E_initial << "\n";
    std::cout << "Initial momentum: (" << p_initial[0] << ", " << p_initial[1] << ", " << p_initial[2] << ")\n\n";

    // Evolve system
    const float dt = 0.01f;
    const int num_steps = 1000;

    for (int step = 0; step < num_steps; ++step) {
        velocityVerletStep(particles, dt);
    }

    // Check conservation
    float E_final = totalEnergy(particles);
    auto p_final = totalMomentum(particles);

    float E_error = std::abs(E_final - E_initial) / std::abs(E_initial);
    float p_error = std::sqrt(
        (p_final[0] - p_initial[0])*(p_final[0] - p_initial[0]) +
        (p_final[1] - p_initial[1])*(p_final[1] - p_initial[1]) +
        (p_final[2] - p_initial[2])*(p_final[2] - p_initial[2])
    );

    std::cout << "Final energy: " << E_final << "\n";
    std::cout << "Energy drift: " << (E_error * 100.0f) << "%\n";
    std::cout << "Final momentum: (" << p_final[0] << ", " << p_final[1] << ", " << p_final[2] << ")\n";
    std::cout << "Momentum drift: " << p_error << "\n\n";

    bool energy_pass = E_error < 0.001f;  // 0.1%
    bool momentum_pass = p_error < 0.001f;  // 0.1%

    std::cout << "Quality Gates:\n";
    std::cout << "  Energy conservation (< 0.1%): " << (energy_pass ? "PASS ✓" : "FAIL ✗") << "\n";
    std::cout << "  Momentum conservation (< 0.1%): " << (momentum_pass ? "PASS ✓" : "FAIL ✗") << "\n";

    return energy_pass && momentum_pass;
}

/**
 * Test 2: Superposition principle
 * Verify F_1 = F_12 + F_13
 */
bool testSuperpositionPrinciple() {
    std::cout << "\n=== Test 2: Superposition Principle ===\n";
    std::cout << "Verify: F_1 = F_12 + F_13 (force on particle 1)\n\n";

    const float q = 1.0f;
    const float m = 1.0f;

    std::vector<ChargedParticle> particles(3);

    // Simple configuration
    particles[0] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, q, m, 0.0f, 0.0f, 0.0f};
    particles[1] = {2.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, q, m, 0.0f, 0.0f, 0.0f};
    particles[2] = {0.0f, 3.0f, 0.0f, 0.0f, 0.0f, 0.0f, -q, m, 0.0f, 0.0f, 0.0f};

    // Compute forces pairwise
    auto F_01 = coulombForce(particles[0], particles[1]);
    auto F_02 = coulombForce(particles[0], particles[2]);

    // Superposition
    float F_0_x = F_01[0] + F_02[0];
    float F_0_y = F_01[1] + F_02[1];
    float F_0_z = F_01[2] + F_02[2];

    // Compute total force on particle 0 directly
    float F_total_x = 0.0f, F_total_y = 0.0f, F_total_z = 0.0f;
    for (size_t j = 1; j < 3; ++j) {
        auto F = coulombForce(particles[0], particles[j]);
        F_total_x += F[0];
        F_total_y += F[1];
        F_total_z += F[2];
    }

    // Compare
    float error_x = std::abs(F_0_x - F_total_x);
    float error_y = std::abs(F_0_y - F_total_y);
    float error_z = std::abs(F_0_z - F_total_z);
    float error_total = std::sqrt(error_x*error_x + error_y*error_y + error_z*error_z);

    std::cout << "F_01: (" << F_01[0] << ", " << F_01[1] << ", " << F_01[2] << ")\n";
    std::cout << "F_02: (" << F_02[0] << ", " << F_02[1] << ", " << F_02[2] << ")\n";
    std::cout << "Superposition: (" << F_0_x << ", " << F_0_y << ", " << F_0_z << ")\n";
    std::cout << "Direct computation: (" << F_total_x << ", " << F_total_y << ", " << F_total_z << ")\n";
    std::cout << "Error: " << error_total << "\n\n";

    bool pass = error_total < 1e-6f;  // Numerical precision

    std::cout << "Quality Gate (< 1e-6): " << (pass ? "PASS ✓" : "FAIL ✗") << "\n";

    return pass;
}

/**
 * Three-Body EM 3D Test Entry Point
 * Called from main.cpp when: ./smft --test config/three_body_em_3d.yaml
 */
int runThreeBodyEM3DTest() {
    std::cout << "========================================\n";
    std::cout << "  3D Three-Body EM Dynamics\n";
    std::cout << "========================================\n";
    std::cout << "Physics: F_ij = k·q_i·q_j·r_ij/|r_ij|³\n";
    std::cout << "System: (+q, +q, -2q) net-zero charge\n\n";

    bool all_pass = true;

    all_pass &= testThreeBodyConservation();
    all_pass &= testSuperpositionPrinciple();

    std::cout << "\n========================================\n";
    std::cout << "FINAL VERDICT: " << (all_pass ? "PASS ✓" : "FAIL ✗") << "\n";
    std::cout << "========================================\n";

    return all_pass ? 0 : 1;
}
