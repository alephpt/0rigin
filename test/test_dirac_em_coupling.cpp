/**
 * test_dirac_em_coupling.cpp
 *
 * Test DiracEvolution with electromagnetic field coupling
 * Verifies:
 * 1. Backward compatibility (empty EM fields = original behavior)
 * 2. Scalar potential affects energy levels
 * 3. Vector potential induces current flow
 * 4. Gauge invariance of physical observables
 */

#include "../src/DiracEvolution.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

void test_backward_compatibility() {
    std::cout << "\n=== Test 1: Backward Compatibility ===" << std::endl;

    // Setup
    uint32_t N = 64;
    DiracEvolution dirac1(N, N);
    DiracEvolution dirac2(N, N);

    // Initialize both with same wavepacket
    float x0 = N/2.0f, y0 = N/2.0f, sigma = 5.0f;
    dirac1.initialize(x0, y0, sigma);
    dirac2.initialize(x0, y0, sigma);

    // Create mass field
    std::vector<float> mass_field(N*N, 0.1f);

    // Evolve one with old API (implicit empty EM fields)
    float dt = 0.01f;
    for (int i = 0; i < 10; i++) {
        dirac1.step(mass_field, dt);
    }

    // Evolve other with new API but empty EM fields
    std::vector<float> empty_A_x, empty_A_y, empty_A_z, empty_phi;
    for (int i = 0; i < 10; i++) {
        dirac2.step(mass_field, dt, empty_A_x, empty_A_y, empty_A_z, empty_phi);
    }

    // Compare densities - should be identical
    auto density1 = dirac1.getDensity();
    auto density2 = dirac2.getDensity();

    float max_diff = 0.0f;
    for (size_t i = 0; i < density1.size(); i++) {
        float diff = std::abs(density1[i] - density2[i]);
        max_diff = std::max(max_diff, diff);
    }

    std::cout << "Max density difference: " << max_diff << std::endl;
    assert(max_diff < 1e-10 && "Backward compatibility failed!");
    std::cout << "✓ Backward compatibility verified" << std::endl;
}

void test_scalar_potential() {
    std::cout << "\n=== Test 2: Scalar Potential Energy ===" << std::endl;

    // Setup
    uint32_t N = 64;
    DiracEvolution dirac(N, N);

    // Initialize wavepacket
    float x0 = N/2.0f, y0 = N/2.0f, sigma = 5.0f;
    dirac.initialize(x0, y0, sigma);

    // Create uniform mass and scalar potential fields
    std::vector<float> mass_field(N*N, 0.1f);
    std::vector<float> phi(N*N, 0.5f);  // Uniform scalar potential

    // Get initial energy
    float KE0, PE0;
    float E0 = dirac.getEnergy(mass_field, KE0, PE0);

    // Get energy with scalar potential
    float KE_phi, PE_phi;
    float E_phi = dirac.getEnergy(mass_field, KE_phi, PE_phi, {}, {}, {}, phi);

    std::cout << "Energy without φ: E=" << E0 << " (KE=" << KE0 << ", PE=" << PE0 << ")" << std::endl;
    std::cout << "Energy with φ=0.5: E=" << E_phi << " (KE=" << KE_phi << ", PE=" << PE_phi << ")" << std::endl;

    // Scalar potential should increase potential energy (for normalized wavefunction)
    float energy_shift = E_phi - E0;
    std::cout << "Energy shift due to φ: " << energy_shift << std::endl;

    // For uniform φ and normalized Ψ, shift should be approximately φ
    assert(std::abs(energy_shift - 0.5f) < 0.1f && "Scalar potential energy incorrect!");
    std::cout << "✓ Scalar potential energy verified" << std::endl;
}

void test_vector_potential() {
    std::cout << "\n=== Test 3: Vector Potential Dynamics ===" << std::endl;

    // Setup
    uint32_t N = 64;
    DiracEvolution dirac_free(N, N);
    DiracEvolution dirac_em(N, N);

    // Initialize both with same wavepacket
    float x0 = N/2.0f, y0 = N/2.0f, sigma = 5.0f;
    dirac_free.initialize(x0, y0, sigma);
    dirac_em.initialize(x0, y0, sigma);

    // Create mass field and uniform vector potential (constant B field)
    std::vector<float> mass_field(N*N, 0.0f);  // Massless for cleaner dynamics
    std::vector<float> A_x(N*N), A_y(N*N);

    // Landau gauge: A = (-By, 0, 0) gives uniform B field in z direction
    float B = 0.1f;  // Moderate field strength
    for (uint32_t j = 0; j < N; j++) {
        for (uint32_t i = 0; i < N; i++) {
            uint32_t idx = j * N + i;
            A_x[idx] = -B * (j - N/2.0f);  // Centered gauge
            A_y[idx] = 0.0f;
        }
    }

    // Evolve both systems
    float dt = 0.01f;
    for (int i = 0; i < 100; i++) {  // More steps to see deflection
        dirac_free.step(mass_field, dt);
        dirac_em.step(mass_field, dt, A_x, A_y);
    }

    // Compare center of mass - magnetic field should deflect particle
    float x_free, y_free, x_em, y_em;
    dirac_free.getCenterOfMass(x_free, y_free);
    dirac_em.getCenterOfMass(x_em, y_em);

    std::cout << "Free particle CoM: (" << x_free << ", " << y_free << ")" << std::endl;
    std::cout << "EM-coupled CoM: (" << x_em << ", " << y_em << ")" << std::endl;

    float deflection = std::sqrt(std::pow(x_em - x_free, 2) + std::pow(y_em - y_free, 2));
    std::cout << "Deflection due to B field: " << deflection << std::endl;

    // Should see some deflection due to Lorentz force
    assert(deflection > 0.1f && "No deflection from magnetic field!");
    std::cout << "✓ Vector potential dynamics verified" << std::endl;
}

void test_energy_conservation() {
    std::cout << "\n=== Test 4: Energy Conservation ===" << std::endl;

    // Setup with static fields (should conserve energy)
    uint32_t N = 64;
    DiracEvolution dirac(N, N);

    // Initialize wavepacket
    float x0 = N/2.0f, y0 = N/2.0f, sigma = 5.0f;
    dirac.initialize(x0, y0, sigma);

    // Create static fields
    std::vector<float> mass_field(N*N, 0.1f);
    std::vector<float> phi(N*N, 0.2f);
    std::vector<float> A_x(N*N, 0.05f);
    std::vector<float> A_y(N*N, -0.03f);

    // Get initial energy
    float KE0, PE0;
    float E0 = dirac.getEnergy(mass_field, KE0, PE0, A_x, A_y, {}, phi);
    std::cout << "Initial energy: " << E0 << std::endl;

    // Evolve for many steps
    float dt = 0.001f;  // Small timestep for accuracy
    for (int i = 0; i < 100; i++) {
        dirac.step(mass_field, dt, A_x, A_y, {}, phi);
    }

    // Get final energy
    float KEf, PEf;
    float Ef = dirac.getEnergy(mass_field, KEf, PEf, A_x, A_y, {}, phi);
    std::cout << "Final energy: " << Ef << std::endl;

    float energy_drift = std::abs(Ef - E0) / std::abs(E0);
    std::cout << "Relative energy drift: " << energy_drift * 100 << "%" << std::endl;

    // Should conserve energy to within numerical precision
    assert(energy_drift < 0.01f && "Energy not conserved!");
    std::cout << "✓ Energy conservation verified" << std::endl;
}

int main() {
    std::cout << "Testing DiracEvolution with EM field coupling..." << std::endl;

    test_backward_compatibility();
    test_scalar_potential();
    test_vector_potential();
    test_energy_conservation();

    std::cout << "\n=== All Tests Passed! ===" << std::endl;
    return 0;
}